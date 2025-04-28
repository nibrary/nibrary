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




float** NIBR::TractogramReader::readStreamlineRawFromFile(FILE* bf, std::size_t n) {

	float** points;
	points = new float* [len[n]];

	fseek(bf, streamlinePos[n], SEEK_SET);

	if (fileFormat == TCK) {
		float tmp;
		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, bf);
				points[i][j] = tmp;
			}
		}
	}


	if (fileFormat == TRK) {
		float tmp, p_tmp[3];
		
		for (uint32_t i = 0; i < len[n]; i++) {

			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, bf);
				p_tmp[j] = tmp - 0.5f;
			}
			applyTransform(points[i],p_tmp,ijk2xyz);
			
			for (int j = 0; j < n_scalars_trk; j++)
				std::fread(&tmp, sizeof(float), 1, bf);
		}
	}
	

	if (fileFormat == VTK_BINARY) {
		float tmp;
		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, bf);
				swapByteOrder(tmp);
				points[i][j] = tmp;
			}
		}
	}

	if (fileFormat == VTK_ASCII) {

		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];
			std::fscanf(bf, "%f %f %f ", &points[i][0], &points[i][1], &points[i][2]);
		}

	}

	return points;

}




bool TractogramReader::loadBatchContaining(std::size_t streamlineIndex) {

    // Calculate the start index of the batch
    std::size_t batchStart = (streamlineIndex / preloadCount) * preloadCount;

    // Check if already loaded
    if (batchStart == preloadedStartIndex && !preloadedStreamlines.empty()) {
        return true;
    }

    disp(MSG_DETAIL, "Loading cache batch starting at index %zu (for requested index %zu)...", batchStart, streamlineIndex);

    // Determine actual number to load
    std::size_t countToLoad = std::min(preloadCount, numberOfStreamlines - batchStart);
    if (countToLoad == 0) {
        resetCache(); // Clear previous cache if requesting beyond end
        return false; // Index out of bounds
    }

    // --- Read the streamlines for this batch ---
    // Create a temporary vector to hold the newly read float** pointers
    std::vector<float**> newBatchData(countToLoad, nullptr); // Initialize with nullptrs
    bool success = true;

    // --- Fallback: Sequential Loading (if parallel is complex/not implemented) ---
    for (std::size_t i = 0; i < countToLoad; ++i) {
        std::size_t currentIndex = batchStart + i;
        // Use the primary file handle for sequential loading
        float** streamline_ptr = readStreamlineRawFromFile(file, currentIndex);

        // Check if read failed for a non-empty streamline
        if (!streamline_ptr && len[currentIndex] > 0) {
            disp(MSG_ERROR, "Failed to read streamline %zu into preload batch.", currentIndex);
            success = false;
            break; // Stop loading this batch
        }
        newBatchData[i] = streamline_ptr; // Store pointer (nullptr if empty)
    }
    // --- END Sequential Loading ---


    // --- If successful, replace the old cache with the new batch ---
    if (success) {
         try {
             // First, clear the old cache (deletes old float** data)
             resetCache();
             // Now, move the newly read data into the main cache
             preloadedStreamlines = std::move(newBatchData);
             preloadedStartIndex = batchStart; // Update the index of the loaded batch
             disp(MSG_DEBUG, "Successfully preloaded %zu streamlines starting at index %zu.", preloadedStreamlines.size(), batchStart);
         } catch (const std::bad_alloc& e) {
             // Should be unlikely with move, but handle defensively
             disp(MSG_FATAL, "Failed to move preloaded data to cache: %s", e.what());
             // Cleanup the data we allocated in newBatchData before it was moved
             for(size_t i=0; i<newBatchData.size(); ++i) deleteStreamlineRaw(newBatchData[i], (batchStart+i < numberOfStreamlines) ? len[batchStart+i] : 0);
             resetCache(); // Ensure main cache is reset
             success = false;
         }
    } else {
        disp(MSG_ERROR, "Failed to load batch starting at index %zu.", batchStart);
        // Cleanup the partially/failed loaded data in newBatchData
        for(size_t i=0; i<newBatchData.size(); ++i) deleteStreamlineRaw(newBatchData[i], (batchStart+i < numberOfStreamlines) ? len[batchStart+i] : 0);
        // Ensure main cache is still in a clean state (should have been reset if old data existed)
        resetCache();
    }

    return success;
}

// Safely delete memory allocated for a float** streamline
void TractogramReader::deleteStreamlineRaw(float** points, uint32_t numPoints) {
    if (!points) return;
    for (uint32_t i = 0; i < numPoints; ++i) {
        delete[] points[i]; // Delete each float[3] array
    }
    delete[] points; // Delete the array of float* pointers
}

void TractogramReader::resetCache() {
    // Assumes reader_mutex is held
    // disp(MSG_DEBUG, "Resetting preload cache. Deleting %zu cached streamlines.", preloadedStreamlines.size());
    // if (len && !preloadedStreamlines.empty() && preloadedStartIndex < numberOfStreamlines) {
    //     std::size_t countInCache = preloadedStreamlines.size();
    //     for (std::size_t i = 0; i < countInCache; ++i) {
    //         std::size_t streamlineIdx = preloadedStartIndex + i;
    //         // Check index validity before accessing len
    //         if (streamlineIdx < numberOfStreamlines) {
    //              deleteStreamlineRaw(preloadedStreamlines[i], len[streamlineIdx]);
    //         } else {
    //              // Should not happen, but handle defensively
    //              deleteStreamlineRaw(preloadedStreamlines[i], 0); // Pass 0 length if index is bad
    //         }
    //     }
    // } else {
    //     // Fallback if len is not available or cache is empty/invalid index
    //     for (float** pts : preloadedStreamlines) {
    //          // We don't know the length, which is risky.
    //          // Assume deleteStreamlineRaw handles nullptr points[i] gracefully if needed.
    //          // This path indicates a potential logic error elsewhere.
    //          if (pts) delete[] pts; // Delete only the outer array if inner is unknown
    //     }
    // }

	// Caller clears memory
    preloadedStreamlines.clear();
    preloadedStreamlines.shrink_to_fit(); // Optional: release memory
    preloadedStartIndex = std::numeric_limits<std::size_t>::max(); // Reset start index
    disp(MSG_DEBUG, "Preload cache reset complete.");
}




float** TractogramReader::readStreamline(std::size_t n) {
    // Acquire lock for the entire operation
    std::lock_guard<std::mutex> lock(reader_mutex);

    if (n >= numberOfStreamlines) {
        disp(MSG_WARN,"Requested streamline index %zu is out of bounds (0-%zu).", n, numberOfStreamlines > 0 ? numberOfStreamlines-1 : 0);
        return nullptr;
    }

    // --- Preloading Logic (inside lock) ---
    if (preloadCount > 0) {
        std::size_t requestedBatchStart = (n / preloadCount) * preloadCount;

		// Check if the required batch is currently loaded
		if (requestedBatchStart != preloadedStartIndex || preloadedStreamlines.empty()) {

			missedCacheCounter++;

			if (missedCacheCounter < MT::MAXNUMBEROFTHREADS()) {
				disp(MSG_DEBUG, "Cache load failed for streamline %zu. Falling back to direct read.", n);
				float** out = readStreamlineRawFromFile(file,n);
				if (!out && len[n]>0) {
					disp(MSG_ERROR, "Direct read fallback also failed for streamline %zu.", n);
					return NULL;
				}
				return out;
			} else {
				disp(MSG_DEBUG, "Too many missed caches. Attempting to load batch starting at %zu.", n, requestedBatchStart);
				if (!loadBatchContaining(n)) {
					disp(MSG_ERROR, "Failed to load required cache batch for streamline %zu.", n);
					return NULL;
				}
				disp(MSG_DEBUG, "Cache loaded successfully for batch starting at %zu.", requestedBatchStart);
			}
		}

        // Cache hit (or cache miss just loaded the batch)
        std::size_t indexInCache = n - preloadedStartIndex;
        if (indexInCache < preloadedStreamlines.size()) {
             disp(MSG_DEBUG, "Cache hit for streamline %zu.", n);
            // Return pointer directly from cache. CAUTION: Cache owns this memory.
            return preloadedStreamlines[indexInCache];
        } else {
            // This indicates a logic error if loadBatchContaining_Locked succeeded
            disp(MSG_FATAL,"Logic error: Streamline %zu not found in cache after loading batch starting %zu.", n, preloadedStartIndex);
            return nullptr; // Return null on error
        }

    } else {
        float** out = readStreamlineRawFromFile(file,n);
		if (!out && len[n]>0) {
			disp(MSG_ERROR, "Direct read fallback also failed for streamline %zu.", n);
			return NULL;
		}
		return out;
    }
}


std::vector<Point> TractogramReader::readStreamlinePoints(std::size_t n) {
    std::vector<Point> points_vec;
    
    // Call readStreamline which handles locking and caching
    float** points_raw = readStreamline(n); 

    // Need to re-acquire lock briefly to get length safely
    uint32_t lineLen = 0;
    { // Short lock scope
        std::lock_guard<std::mutex> lock(reader_mutex);
        if (len == nullptr || n >= numberOfStreamlines) {
             // Error already logged by readStreamline or index is invalid
             return points_vec; // Return empty
        }
        lineLen = len[n];
    } // Lock released


    if (points_raw != nullptr && lineLen > 0) {
        try {
            points_vec.reserve(lineLen);
            for (uint32_t i = 0; i < lineLen; ++i) {
                if (points_raw[i] != nullptr) { // Check inner pointer
                    points_vec.push_back({points_raw[i][0], points_raw[i][1], points_raw[i][2]});
                } else {
                    // Handle potential null inner pointer if allocation failed partially
                    disp(MSG_ERROR, "Null inner pointer encountered during conversion for streamline %zu, point %u.", n, i);
                    points_vec.clear(); // Invalidate result
                    return points_vec;
                }
            }
        } catch (const std::bad_alloc& e) {
            disp(MSG_ERROR, "Memory allocation failed during Point conversion for streamline %zu: %s", n, e.what());
            points_vec.clear(); // Return empty on error
        }
    } else if (lineLen > 0 && points_raw == nullptr) {
         // Log error if readStreamline failed for a non-empty streamline
         disp(MSG_WARN, "readStreamlinePoints: Failed to get raw data for streamline %zu.", n);
    }
    // If lineLen is 0, points_raw might be null, return empty vector which is correct.

    return points_vec;
}


// Reads streamline 'n' and converts to std::vector<std::vector<float>>.
std::vector<std::vector<float>> TractogramReader::readStreamlineVector(std::size_t n) {
	std::vector<std::vector<float>> points_vec;
   
   // Call readStreamline which handles locking and caching
   float** points_raw = readStreamline(n); 

   // Need to re-acquire lock briefly to get length safely
   uint32_t lineLen = 0;
   { // Short lock scope
	   std::lock_guard<std::mutex> lock(reader_mutex);
	   if (len == nullptr || n >= numberOfStreamlines) {
			return points_vec; // Return empty
	   }
	   lineLen = len[n];
   } // Lock released


   if (points_raw != nullptr && lineLen > 0) {
	   try {
		   points_vec.reserve(lineLen);
		   for (uint32_t i = 0; i < lineLen; ++i) {
				if (points_raw[i] != nullptr) {
				   points_vec.push_back({points_raw[i][0], points_raw[i][1], points_raw[i][2]});
				} else {
				   disp(MSG_ERROR, "Null inner pointer encountered during vector<float> conversion for streamline %zu, point %u.", n, i);
				   points_vec.clear();
				   return points_vec;
				}
		   }
	   } catch (const std::bad_alloc& e) {
		   disp(MSG_ERROR, "Memory allocation failed during vector<float> conversion for streamline %zu: %s", n, e.what());
		   points_vec.clear();
	   }
   } else if (lineLen > 0 && points_raw == nullptr) {
		disp(MSG_WARN, "readStreamlineVector: Failed to get raw data for streamline %zu.", n);
   }

   return points_vec;
}

// Reads all streamlines into a vector of vectors of float vectors.
std::vector<std::vector<std::vector<float>>> TractogramReader::read() {
    std::vector<std::vector<std::vector<float>>> allStreamlines;
     if (numberOfStreamlines == 0) return allStreamlines;

    disp(MSG_INFO, "Reading all %zu streamlines into vector<vector<float>>...", numberOfStreamlines);
    try {
        allStreamlines.reserve(numberOfStreamlines);
        for (std::size_t n = 0; n < numberOfStreamlines; ++n) {
            allStreamlines.push_back(readStreamlineVector(n)); // Uses caching internally
        }
     } catch (const std::bad_alloc& e) {
         disp(MSG_FATAL, "Memory allocation failed while reading all streamlines (Vectors): %s", e.what());
         allStreamlines.clear();
     }
    disp(MSG_INFO, "Finished reading all streamlines.");
    return allStreamlines;
}
