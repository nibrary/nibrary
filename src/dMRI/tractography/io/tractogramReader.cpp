#include <cstdint>
#include <cstring>
#include <iterator>
#include <type_traits>
#include <cctype>
#include <atomic>

#include "tractogramReader.h"
#include "tractogramField.h"
#include "base/dataTypeHandler.h"

using namespace NIBR;

TractogramReader::TractogramReader(const std::string& _fileName, std::size_t _preloadCount)
        : fileName(_fileName),
          preloadCount(_preloadCount)
{
	std::memset(imgDims, 0, sizeof(imgDims));
	std::memset(pixDims, 0, sizeof(pixDims));
	std::memset(voxOrdr, 0, sizeof(voxOrdr));
	std::memset(ijk2xyz, 0, sizeof(ijk2xyz));
	std::memset(xyz2ijk, 0, sizeof(xyz2ijk));
	ijk2xyz[3][3] = 1.0f; // Ensure valid affine matrix base
	xyz2ijk[3][3] = 1.0f;
	if (!initReader()) {
		cleanupMemory();
		throw std::runtime_error("Failed to initialize TractogramReader for file: " + _fileName);
	}
	disp(MSG_DEBUG, "TractogramReader initialized successfully for %s", fileName.c_str());
}

TractogramReader::~TractogramReader() {
	cleanupMemory();
	disp(MSG_DEBUG, "TractogramReader destroyed for %s", fileName.c_str());
}

void TractogramReader::cleanupMemory() {
	
	std::lock_guard<std::mutex> lock(reader_mutex);

	disp(MSG_DEBUG, "Cleaning up memory for %s", fileName.c_str());
	
	if (file != NULL) {
		fclose(file);
		file = NULL;
		disp(MSG_DEBUG, "Closed file handle for %s", fileName.c_str());
	}

	for (auto bf : batchFile) {
		if (bf != NULL) {
			fclose(bf);
			bf = NULL;
			disp(MSG_DEBUG, "Closed file batch file handle");
		}
	}
	batchFile.clear();

	if (len != NULL)			delete[] len;
	if (streamlinePos != NULL) 	delete[] streamlinePos;

	fileName.erase();
	fileDescription.erase();

	resetCache();

	numberOfStreamlines = 0;
	numberOfPoints 	 	= 0;
}


bool NIBR::TractogramReader::initReader() {

	std::lock_guard<std::mutex> lock(reader_mutex);

	file = fopen(fileName.c_str(), "rb");

	if (file == NULL) {
		disp(MSG_ERROR, "Failed to open tractogram file: %s", fileName.c_str());
		return false;
	}

	const std::size_t strLength = 256;
	char dummy[strLength];

	std::string extension = getFileExtension(fileName);

	disp(MSG_DEBUG,"File extension: %s", extension.c_str());

	// --- TCK Format ---
	if (extension == "tck") {

		fileFormat 		= TCK;
		fileDescription = SGNTR();

		std::string tmps;
		long pos = 0;
		do {
			fgets(dummy, strLength, file);
			tmps = std::string(dummy);

			std::size_t column = tmps.find_last_of(":");

			if (column != std::string::npos) {
				std::string key = tmps.substr(0, column);

				if (key == "count") {	
					numberOfStreamlines = std::atoi((tmps.substr(column + 1)).c_str());
				} 
				if (key == "file") {
					pos = std::atoi((tmps.substr(column + 3)).c_str());
				} 
			}

		} while (tmps != "END\n");

		numberOfPoints = 0;
		
		if (numberOfStreamlines > 0) {

			len = new uint32_t[numberOfStreamlines];
			streamlinePos = new long[numberOfStreamlines];
			streamlinePos[0] = pos;
			len[0] = 0;
			fseek(file, pos, SEEK_SET);

			float  tmp[3] = { 0,0,0 };
			std::size_t ind = 0;

			while (!feof(file)) {
				std::fread(tmp, sizeof(float), 3, file);
				len[ind]++;
				
				if (std::isnan(tmp[0])) {
					len[ind]--;
					numberOfPoints += len[ind];
					ind++;
					if (ind == numberOfStreamlines)
						break;
					else {
						streamlinePos[ind] = ftell(file);
						len[ind] = 0;
					}
				}

				// if (NIBR::VERBOSE()>VERBOSE_INFO) {
                //     std::cout << "                                                                                                                 " << '\r' << std::flush; // clean terminal
                //     std::cout << std::fixed << std::setprecision (2)  << preamble << message << ": " << float(ind)/float(numberOfStreamlines)*100 << "%" <<"\033[0m" << '\r' << std::flush;
                // }

			}

		} else {
			numberOfStreamlines = 0;
			len 				= NULL;
			streamlinePos 		= NULL;
		}

	}

	if (extension == "trk") {
		
		fileFormat 		= TRK;
		fileDescription = SGNTR();

		trkFileStruct trkFileHeader;

		std::fread(&trkFileHeader, sizeof(trkFileStruct), 1, file);

		for(int i=0;i<4;++i)
			for(int j=0;j<4;++j)
				ijk2xyz[i][j] = trkFileHeader.vox_to_ras[i][j];

		for(int i=0;i<3;++i) {
			imgDims[i] = trkFileHeader.dim[i];
			pixDims[i] = trkFileHeader.voxel_size[i];
		}
		std::strcpy(voxOrdr,"LAS");

		bool allZero = true;

		for(int i=0;i<3;++i)
			for(int j=0;j<3;++j)
				if (ijk2xyz[i][j] != 0)
					allZero = false;

		if ( allZero || (trkFileHeader.version != 2) ) {
			disp(MSG_FATAL,"Trk file version must be 2. Old trk files are not accepted.");
			return false;
		}

		inverseAffine(ijk2xyz, xyz2ijk);

		n_scalars_trk 		= trkFileHeader.n_scalars;
		n_properties_trk	= trkFileHeader.n_properties;

		numberOfPoints      = 0;
		numberOfStreamlines	= trkFileHeader.n_count;

		if (numberOfStreamlines > 0) {

			len 				= new uint32_t[numberOfStreamlines];
			streamlinePos 		= new long[numberOfStreamlines];
			
			int tmp;
			std::fread(&tmp, sizeof(int), 1, file);

			len[0] 				= tmp;
			streamlinePos[0] 	= ftell(file);
			numberOfPoints     += tmp;

			for (std::size_t i = 1; i < numberOfStreamlines; i++) {
				long offset = sizeof(float) * (len[i-1]*(3+n_scalars_trk)+n_properties_trk);
				std::fseek(file, offset, SEEK_CUR);
				std::fread(&tmp, sizeof(int), 1, file);

				len[i] 			 = tmp;
				streamlinePos[i] = ftell(file);
				numberOfPoints  += tmp;

				// if (NIBR::VERBOSE()>VERBOSE_INFO) {
				// 	std::cout << "                                                                                                                 " << '\r' << std::flush; // clean terminal
				// 	std::cout << std::fixed << std::setprecision (2)  << preamble << message << ": " << float(i)/float(numberOfStreamlines)*100 << "%" <<"\033[0m" << '\r' << std::flush;
				// }

			}

		} else {
			numberOfStreamlines = 0;
			len 				= NULL;
			streamlinePos 		= NULL;
		}

	}

	if (extension == "vtk") {

		std::string vtkFormat;

		fgets(dummy, strLength, file);                        // vtk version
		
		if (std::string(dummy).find("3")==std::string::npos) {
			disp(MSG_WARN,"Vtk version is not supported: %s", dummy);
		}

		disp(MSG_DEBUG,"Vtk version is: %s", dummy);

		fgets(dummy, strLength, file);                        // file description
		fileDescription = std::string(dummy);
		fileDescription = fileDescription.substr(0, fileDescription.size() - 1);
		std::fscanf(file, "%s ", dummy);                     // ascii or binary
		vtkFormat = std::string(dummy);
		fgets(dummy, strLength, file);                        // always DATASET POLYDATA, we skip to check this for now
		std::fscanf(file, "%*s %zu %*s ", &numberOfPoints);  // number of points and datatype, we assume datatype is float and skip checking this
		long posData = ftell(file);

		disp(MSG_DEBUG,"numberOfPoints: %d", numberOfPoints);

		if ((vtkFormat == "ascii") || (vtkFormat == "ASCII")) fileFormat = VTK_ASCII;
		else if ((vtkFormat == "binary") || (vtkFormat == "BINARY")) fileFormat = VTK_BINARY;
		else return false;

		if (numberOfPoints > 0) {

			// Skip points
			if (fileFormat == VTK_BINARY) {
				std::fseek(file, sizeof(float) * numberOfPoints * 3, SEEK_CUR);
				int tmp = std::fgetc(file); if (tmp != '\n') std::ungetc(tmp, file); // Make sure to go end of the line
			}

			if (fileFormat == VTK_ASCII) {
				for (std::size_t i = 0; i < numberOfPoints; i++)
					std::fscanf(file, "%*f %*f %*f ");
			}

			// Get number of streamlines, lengths and streamlinePos
			auto checkLine = (std::fscanf(file, "LINES %zu %*d ", &numberOfStreamlines) == 1);

			if ((numberOfStreamlines < 1) || !checkLine) {

				disp(MSG_DEBUG,"No lines in file");

				numberOfStreamlines = 0;
				len 				= NULL;
				streamlinePos 		= NULL;

			} else {

				disp(MSG_DEBUG,"numberOfStreamlines: %d", numberOfStreamlines);

				len 			 = new uint32_t[numberOfStreamlines];
				streamlinePos 	 = new long[numberOfStreamlines];
				streamlinePos[0] = posData;

				if (fileFormat == VTK_BINARY) {

					int tmp;
					std::fread(&tmp, sizeof(int), 1, file);
					swapByteOrder(tmp);
					std::fseek(file, tmp * sizeof(int), SEEK_CUR);
					len[0] = tmp;

					for (std::size_t i = 1; i < numberOfStreamlines; i++) {
						streamlinePos[i] = streamlinePos[i - 1] + long(sizeof(float) * len[i - 1] * 3);
						std::fread(&tmp, sizeof(int), 1, file);
						swapByteOrder(tmp);
						std::fseek(file, tmp * sizeof(int), SEEK_CUR);
						len[i] = tmp;
					}
				}

				if (fileFormat == VTK_ASCII) {
					for (std::size_t i = 0; i < numberOfStreamlines; i++) {
						std::fscanf(file, "%u ", &len[i]);
						for (uint32_t j = 0; j < (len[i] - 1); j++)
							std::fscanf(file, "%*d ");
						std::fscanf(file, "%*d ");
						// disp(MSG_DEBUG,"len: %d", len[i]);
					}

					fseek(file, posData, SEEK_SET);

					for (std::size_t i = 0; i < (numberOfStreamlines - 1); i++) {
						for (std::size_t j = 0; j < len[i]; j++)
							std::fscanf(file, "%*f %*f %*f ");
						streamlinePos[i + 1] = ftell(file);

					}
				}

			}

		}
		else {
			numberOfStreamlines = 0;
			len 				= NULL;
			streamlinePos 		= NULL;
		}

	}

	if (numberOfStreamlines > 0) {
		std::fseek(file, streamlinePos[numberOfStreamlines-1], 		SEEK_SET);
		std::fseek(file, sizeof(float)*len[numberOfStreamlines-1]*3,SEEK_CUR);
	}
	endPosOfStreamlines = ftell(file);

	for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {

		FILE* bf = fopen(fileName.c_str(), "rb");

		if (bf == NULL) {
			disp(MSG_ERROR, "Failed to open tractogram file for parallel reading.");
			return false;
		}

		batchFile.push_back(bf);

	}


	return true;
}

float** NIBR::TractogramReader::readStreamline(std::size_t n) {

	std::vector<Point> points_vec = readStreamlinePoints(n);

	if (points_vec.empty()) {
		disp(MSG_WARN, "readStreamline: Failed to read points for streamline %zu.", n);
		return NULL;
	}

	float** points = NULL;

	try {
		std::size_t numPoints = points_vec.size();
		points = new float*[numPoints];

		for (std::size_t n = 0; n < numPoints; ++n) {
			points[n] = new float[3];
			points[n][0] = points_vec[n].x;
			points[n][1] = points_vec[n].y;
			points[n][2] = points_vec[n].z;
		}
	} catch (const std::bad_alloc& e) {
		disp(MSG_FATAL, "Memory allocation failed in readStreamline for streamline %zu: %s", n, e.what());
		// Accept some memory leak in the case of error
		delete[] points;
		return NULL;
	}

	return points;

}

std::vector<std::vector<float>> NIBR::TractogramReader::readStreamlineVector(std::size_t n) {

	std::vector<Point> points_vec = readStreamlinePoints(n);

	std::vector<std::vector<float>> v;

	if (points_vec.empty()) {
		disp(MSG_WARN, "readStreamlineVector: Failed to read points for streamline %zu.", n);
		return v;
	}

	// Reserve space for efficiency
	try {
		v.reserve(points_vec.size());
	} catch (const std::bad_alloc& e) {
		disp(MSG_ERROR, "Failed to reserve memory for %u points for streamline %zu: %s", points_vec.size(), n, e.what());
		return std::vector<std::vector<float>>();
	}

	for (const auto& p : points_vec) {
		v.emplace_back(std::vector<float>{p.x, p.y, p.z});
	}

	return v;

}


std::vector<Point> TractogramReader::readStreamlinePoints(std::size_t n) {

	std::vector<Point> points;

	std::lock_guard<std::mutex> lock(reader_mutex);

	if (preloadCount > 0) {

		std::size_t requestedBatchStart = (n / preloadCount) * preloadCount;

		// Check if the required batch is currently loaded
		if (requestedBatchStart != preloadedStartIndex || preloadedStreamlines.empty()) {

			missedCacheCounter++;

			if (missedCacheCounter < MT::MAXNUMBEROFTHREADS()) {
				disp(MSG_DEBUG, "Cache load failed for streamline %zu. Falling back to direct read.", n);
				if (!readStreamlinePointsFromFile(file, n, points)) {
					disp(MSG_ERROR, "Direct read fallback also failed for streamline %zu.", n);
					return {};
				}
				return points;
			} else {
				disp(MSG_DEBUG, "Too many missed caches. Attempting to load batch starting at %zu.", n, requestedBatchStart);
				if (!loadBatchContaining(n)) {
					disp(MSG_ERROR, "Failed to load required cache batch for streamline %zu.", n);
					return {};
				}
				disp(MSG_DEBUG, "Cache loaded successfully for batch starting at %zu.", requestedBatchStart);
			}
		}

		// Cache hit
		std::size_t indexInCache = n - preloadedStartIndex;
		if (indexInCache < preloadedStreamlines.size()) {
			disp(MSG_DEBUG, "Cache hit for streamline %zu (index %zu in cache).", n, indexInCache);
			// **Return a COPY from the cache** for thread safety after lock release
			 points = preloadedStreamlines[indexInCache];
		} else {
			// This indicates a logic error if loadBatchContaining_Locked succeeded
			disp(MSG_FATAL, "Logic error: Streamline %zu not found in cache (size %zu) after loading batch starting at %zu.",
				 n, preloadedStreamlines.size(), preloadedStartIndex);
			// Clear cache to indicate error state?
			resetCache();
			return {}; // Return empty on logic error
		}

	} else {
		// No Preloading
		if (!readStreamlinePointsFromFile(file, n, points)) {
			disp(MSG_ERROR, "Read failed for streamline %zu.", n);
			return {};
		}
	}

	return points;
}



std::vector<std::vector<std::vector<float>>> NIBR::TractogramReader::read() {

	std::vector<std::vector<std::vector<float>>> out;

	if (numberOfStreamlines==0) 
		return out;

	for (std::size_t n = 0; n < numberOfStreamlines; n++) {
		out.emplace_back(readStreamlineVector(n));
	}

	return out;

}

bool NIBR::TractogramReader::readStreamlinePointsFromFile(FILE* bf, std::size_t n, std::vector<Point>& points) {

	points.clear();

	uint32_t lineLen = len[n];

	try {
		points.reserve(lineLen);
	} catch (const std::bad_alloc& e) {
		disp(MSG_ERROR, "Failed to reserve memory for %u points for streamline %zu: %s", lineLen, n, e.what());
		return false;
	}

	if (fseek(bf, streamlinePos[n], SEEK_SET) != 0) {
		disp(MSG_ERROR, "Failed to seek to streamline %zu position %ld (errno %d).", n, streamlinePos[n], errno);
		return false;
	}

	float tmp_p[3];

	try {
		switch (fileFormat) {
			case TCK: {
				for (uint32_t i = 0; i < lineLen; ++i) {
					if (fread(tmp_p, sizeof(float), 3, bf) != 3) {
						 disp(MSG_ERROR, "TCK: Read error reading point %u/%u for streamline %zu at offset %ld.", i + 1, lineLen, n, ftell(bf));
						 return false;
					}
					// TODO: Handle TCK endianness if dataType indicated non-native
					points.emplace_back(Point{tmp_p[0], tmp_p[1], tmp_p[2]});
				}
				break;
			}
			case TRK: {
				// Buffer to read point (coords + scalars)
				std::size_t floats_per_point = 3 + n_scalars_trk;
				float* point_buffer = new float[floats_per_point];

				for (uint32_t i = 0; i < lineLen; ++i) {
					if (fread(point_buffer, sizeof(float), floats_per_point, bf) != floats_per_point) {
						disp(MSG_ERROR, "TRK: Read error reading data for point %u/%u for streamline %zu at offset %ld.", i + 1, lineLen, n, ftell(bf));
						delete[] point_buffer;
						return false;
					}
					point_buffer[0] -= 0.5f;
					point_buffer[1] -= 0.5f;
					point_buffer[2] -= 0.5f;
					
					// Applies transform only to the coordinate part (first 3 floats)
					applyTransform(point_buffer,ijk2xyz);
					points.emplace_back(Point{point_buffer[0], point_buffer[1], point_buffer[2]});
				}

				delete[] point_buffer;

				break;
			}
			case VTK_BINARY: {
				for (uint32_t i = 0; i < lineLen; ++i) {
					if (fread(tmp_p, sizeof(float), 3, bf) != 3) {
						disp(MSG_ERROR, "VTK Binary: Read error reading point %u/%u for streamline %zu at offset %ld.", i + 1, lineLen, n, ftell(bf));
						return false;
					}
					// Swap bytes after reading if needed (assuming VTK standard is Big Endian)
					NIBR::swapByteOrder(tmp_p[0]);
					NIBR::swapByteOrder(tmp_p[1]);
					NIBR::swapByteOrder(tmp_p[2]);
					points.emplace_back(Point{tmp_p[0], tmp_p[1], tmp_p[2]});
				}
				break;
			}
			case VTK_ASCII: {
				// Seeking should have placed us at the start of point data for this line
				for (uint32_t i = 0; i < lineLen; ++i) {
					if (fscanf(bf, "%f %f %f", &tmp_p[0], &tmp_p[1], &tmp_p[2]) != 3) {
						disp(MSG_ERROR, "VTK ASCII: Read error scanning point %u/%u for streamline %zu at offset %ld.", i + 1, lineLen, n, ftell(bf));
						// Check for EOF vs other errors
						if (feof(bf)) disp(MSG_ERROR, " (EOF reached unexpectedly)");
						return false;
					}
					points.emplace_back(Point{tmp_p[0], tmp_p[1], tmp_p[2]});
					// Consume potential trailing whitespace/newline robustly
					int c;
					while ((c = fgetc(bf)) != EOF && isspace(c)) {}
					if (c != EOF) ungetc(c, bf);
				}
				break;
			}
			default:
				// Should not happen if initReader worked correctly
				disp(MSG_FATAL, "readStreamlinePointsFromFile: Invalid file format encountered.");
				return false;
		}
	} catch (const std::exception& e) {
		// Catch potential exceptions from vector operations (e.g., bad_alloc in emplace_back?)
		disp(MSG_ERROR, "Exception during file read for streamline %zu: %s", n, e.what());
		return false;
	}

	// Final check: Ensure the correct number of points were added
	if (points.size() != lineLen) {
		 disp(MSG_ERROR, "Internal error: Read %zu points for streamline %zu, but expected %u.", points.size(), n, lineLen);
		 return false; // Data inconsistency
	}

	return true; // Success

}

void TractogramReader::resetCache() {
	preloadedStreamlines.clear();
	preloadedStreamlines.shrink_to_fit();
	preloadedStartIndex = std::numeric_limits<std::size_t>::max();
	disp(MSG_DEBUG, "Preload cache reset.");
}

bool TractogramReader::loadBatchContaining(std::size_t streamlineIndex) {
	// Assumes reader_mutex is held by caller

	// Should not be called if preloading is disabled, but check anyway
	if (preloadCount == 0 || file == NULL || len == NULL || streamlinePos == NULL) {
		disp(MSG_ERROR, "loadBatchContaining called unexpectedly (preload disabled or reader not ready).");
		return false;
	}

	// Calculate the start index of the batch containing streamlineIndex
	std::size_t batchStart = (streamlineIndex / preloadCount) * preloadCount;

	// Check if the requested batch is already loaded (redundant check, caller usually does this, but safe)
	if (batchStart == preloadedStartIndex && !preloadedStreamlines.empty()) {
		disp(MSG_DEBUG, "Batch starting at %zu already loaded.", batchStart);
		return true; // Already loaded
	}

	disp(MSG_DEBUG, "Loading cache batch starting at index %zu (for requested index %zu)...", batchStart, streamlineIndex);

	// Determine how many streamlines to actually load in this batch
	std::size_t countToLoad = std::min(preloadCount, numberOfStreamlines - batchStart);

	if (countToLoad == 0) {
		// This can happen if streamlineIndex >= numberOfStreamlines
		resetCache(); // Clear cache if requesting beyond the end
		disp(MSG_WARN, "Attempted to load batch starting at or beyond end of file (%zu >= %zu)", batchStart, numberOfStreamlines);
		return false; // Nothing to load
	}

	// --- Allocate space in the cache ---
	std::vector<std::vector<Point>> newBatchData; // Load into temporary first
	try {
		newBatchData.resize(countToLoad);
	} catch (const std::bad_alloc& e) {
		disp(MSG_FATAL, "Failed to allocate memory for temporary preload batch: %s", e.what());
		resetCache(); // Reset main cache state
		return false; // Indicate failure
	}

	// --- Read the streamlines for this batch into the temporary vector ---
	/*
	bool success = true;

	std::atomic<int64_t> mt_failed_streamline_idx(-1);

	auto batchRead = [&] (MT::TASK task) {

		if (mt_failed_streamline_idx.load(std::memory_order_relaxed) > -1) {
            return;
        }

		std::size_t currentIndex 	= batchStart + task.no;
        int fileHandleIndex 		= task.threadId % batchFile.size();
        std::vector<Point> tmp;


		if (!readStreamlinePointsFromFile(batchFile[fileHandleIndex], currentIndex, tmp)) {
			int64_t expected = -1;
			mt_failed_streamline_idx.compare_exchange_strong(expected, (int64_t)currentIndex, std::memory_order_relaxed);
			return;
		}

		MT::PROC_MX().lock();
		newBatchData[task.no] = std::move(tmp);
		MT::PROC_MX().unlock();
	};
	MT::MTRUN(countToLoad,batchRead);

	if (mt_failed_streamline_idx.load() > -1) {
		disp(MSG_ERROR, "Failed to read streamline %d into preload batch.", mt_failed_streamline_idx.load());
		success = false;
	}
	*/


	bool success = true;
	for (std::size_t i = 0; i < countToLoad; ++i) {
		std::size_t currentIndex = batchStart + i;
		// Use the locked file reading function
		if (!readStreamlinePointsFromFile(file, currentIndex, newBatchData[i])) {
			disp(MSG_ERROR, "Failed to read streamline %zu into preload batch.", currentIndex);
			success = false;
		}
	}

	// --- If successful, move the temporary data to the main cache ---
	if (success) {
		 try {
			 // Move assignment is efficient
			 preloadedStreamlines = std::move(newBatchData);
			 preloadedStartIndex = batchStart; // Update the index of the loaded batch
			 disp(MSG_DEBUG, "Successfully preloaded %zu streamlines starting at index %zu.", preloadedStreamlines.size(), batchStart);
		 } catch (const std::bad_alloc& e) {
			 // Should be unlikely with move, but handle defensively
			 disp(MSG_FATAL, "Failed to move preloaded data to cache: %s", e.what());
			 resetCache(); // Reset main cache state
			 success = false;
		 }
	} else {
		disp(MSG_ERROR, "Failed to load batch starting at index %zu.", batchStart);
		// Ensure main cache is reset if loading failed
		resetCache();
	}

	missedCacheCounter = 0;

	return success;
}


std::vector<NIBR::TractogramField> NIBR::TractogramReader::findTractogramFields() {

	std::vector<NIBR::TractogramField> fieldList;

    if (fileFormat == VTK_ASCII) {
        disp(MSG_ERROR,"Can't read fields from ASCII files");
        return fieldList;
    }

    auto input = file;
    const std::size_t strLength = 256;
	char dummy[strLength];

    std::fseek(input, endPosOfStreamlines, SEEK_SET);
    int tmp = std::fgetc(input);
    if (tmp != '\n') std::ungetc(tmp,input);
    std::fgets(dummy,strLength,input); // Skip the line about line number
    std::fseek(input, sizeof(int)*(numberOfStreamlines+numberOfPoints),SEEK_CUR);
    tmp = std::fgetc(input);
    if (tmp != '\n') std::ungetc(tmp,input);

    char* name = new char[128];
    char* type = new char[128];
    int   dimension;
    bool  cellDataFound  = false;
    bool  pointDataFound = false;
    int   tmpi;

	auto toUpperCase = [&] (const std::string& input) -> std::string {
		std::string result = input;
		for (char &c : result) {
			c = std::toupper(static_cast<unsigned char>(c));
		}
		return result;
	};
    
    while(feof(input) == 0) {
        
        std::fgets(dummy,strLength,input);
        
        if (std::string(dummy).find("CELL_DATA")!=std::string::npos) {
            disp(MSG_DEBUG,"Cell data found");
            cellDataFound   = true;
            pointDataFound  = false;
            continue;
        }
    
        if (std::string(dummy).find("POINT_DATA")!=std::string::npos)  {
            disp(MSG_DEBUG,"Point data found");
            cellDataFound  = false;
            pointDataFound = true;
            continue;
        }

        if (std::string(dummy).find("SCALARS")!=std::string::npos) {
            std::sscanf(dummy,"SCALARS %s %s %d\n", name, type, &dimension);
            std::fgets(dummy,strLength,input);

            disp(MSG_DEBUG,"Found field %s", name);

            if (cellDataFound==true) {
                TractogramField f = {STREAMLINE_OWNER,std::string(name),NIBR::getTypeId(toUpperCase(std::string(type))),dimension,NULL};
                fieldList.push_back(f);
                if (std::string(type)=="float") std::fseek(input,sizeof(float)*numberOfStreamlines*dimension,SEEK_CUR);
                if (std::string(type)=="int")   std::fseek(input,sizeof(int)*numberOfStreamlines*dimension,SEEK_CUR);
            }
            if (pointDataFound==true) {
                TractogramField f = {POINT_OWNER,std::string(name),NIBR::getTypeId(toUpperCase(std::string(type))),dimension,NULL};
                fieldList.push_back(f);
                if (std::string(type)=="float") std::fseek(input,sizeof(float)*numberOfPoints*dimension,SEEK_CUR);
                if (std::string(type)=="int")   std::fseek(input,sizeof(int)*numberOfPoints*dimension,SEEK_CUR);
            }
            tmpi = std::fgetc(input); if (tmpi != '\n') std::ungetc(tmpi,input); // Make sure to go end of the line
        }

    }

    delete[] name;
    delete[] type;

    disp(MSG_DEBUG,"Found %d fields", fieldList.size());
    
    return fieldList;

}


void NIBR::TractogramReader::printInfo() {

	disp(MSG_INFO,"Tractogram info");
	std::cout << "\033[32m";

	std::cout << "File name: " << fileName << std::endl;

	std::cout << "Format:                   ";
	if (fileFormat==NIBR::TCK)              std::cout << "tck"          << std::endl << std::flush;
	else if (fileFormat==NIBR::TRK)         std::cout << "trk"          << std::endl << std::flush;
	else if (fileFormat==NIBR::VTK_ASCII)   std::cout << "vtk (ascii)"  << std::endl << std::flush;
	else if (fileFormat==NIBR::VTK_BINARY)  std::cout << "vtk (binary)" << std::endl << std::flush;
	else {
		std::cout << "unknown"      << std::endl << std::flush;
		return;
	}

	std::cout << "Description:              "  << fileDescription      << std::endl << std::flush;
	std::cout << "Streamline count:         "  << numberOfStreamlines  << std::endl << std::flush;
	std::cout << "Number of points:         "  << numberOfPoints       << std::endl << std::flush;

	std::size_t total = 0;
	for (std::size_t i=0; i<numberOfStreamlines; i++) {
		total += len[i];
	}
	std::cout << "Number of points (check): "  << total << std::endl << std::flush;

	if (fileFormat==NIBR::TRK) {

		std::cout << "ijk2xyz: " << std::endl << std::flush;
		
		std::cout << "   " << ijk2xyz[0][0] << " " << ijk2xyz[0][1] << " " << ijk2xyz[0][2] << " " << ijk2xyz[0][3] << std::endl << std::flush;
		std::cout << "   " << ijk2xyz[1][0] << " " << ijk2xyz[1][1] << " " << ijk2xyz[1][2] << " " << ijk2xyz[1][3] << std::endl << std::flush;
		std::cout << "   " << ijk2xyz[2][0] << " " << ijk2xyz[2][1] << " " << ijk2xyz[2][2] << " " << ijk2xyz[2][3] << std::endl << std::flush;
		std::cout << "   " << ijk2xyz[3][0] << " " << ijk2xyz[3][1] << " " << ijk2xyz[3][2] << " " << ijk2xyz[3][3] << std::endl << std::flush;

	}

	auto fields = findTractogramFields();

	std::cout << "Field count: " << fields.size() << std::endl << std::flush;

	int i = 1;
	for (auto f : fields) {
		std::cout << "   " << i++ << ". " << f.name.c_str() << ": " << ((f.owner == POINT_OWNER) ? "Point" : "Streamline") << " field with dim " << f.dimension << std::endl << std::flush;
	}

	std::cout << "\033[0m";

	return;


}