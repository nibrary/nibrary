#include "dMRI/tractography/io/tractogramReader.h"
#include "tractogramField.h"
#include <cstdint>
#include <cstring>
#include <iterator>
#include "expat/expat.h"

using namespace NIBR;
			
NIBR::TractogramReader::TractogramReader() {
	fileName 			= "";
	fileDescription 	= "";
	file 				= NULL;
	usePreload 			= false;
	streamlineBuffer 	= NULL;
	fileBuffer 			= NULL;

}

NIBR::TractogramReader::TractogramReader(std::string _fileName) {
	fileName 			= "";
	fileDescription 	= "";
	file 				= NULL;
	usePreload 			= false;
	streamlineBuffer 	= NULL;
	fileBuffer 			= NULL;
	initReader(_fileName);
}

NIBR::TractogramReader::TractogramReader(std::string _fileName, bool _usePreload) {
	fileName 			= "";
	fileDescription		= "";
	file 				= NULL;
	usePreload 			= _usePreload;
	streamlineBuffer 	= NULL;
	fileBuffer 			= NULL;
	initReader(_fileName,usePreload);
}

NIBR::TractogramReader::TractogramReader(const TractogramReader& obj) {
	this->copyFrom(obj);
}

NIBR::TractogramReader::~TractogramReader() {

	if (finishedLoading && streamlineBuffer != NULL) {
		for (std::size_t n = 0; n < numberOfStreamlines; n++) {
			if (fileBuffer == NULL) {
				for (std::size_t l = 0; l < len[n]; l++) {
					delete[] streamlineBuffer->at(n)[l];
				}
			}
			delete[] streamlineBuffer->at(n);
		}
	}

	if (fileBuffer != NULL) delete[] fileBuffer;
	fileBuffer = NULL;

	fileName.erase();
	fileDescription.erase();
	
	if (file != NULL)			fclose(file);
	if (len != NULL)			delete[] len;
	if (streamlinePos != NULL) 	delete[] streamlinePos;

	
	file 						= NULL;
	len 						= NULL;
	streamlinePos 				= NULL;

	if (streamlineBuffer != NULL) {
		streamlineBuffer->clear();
		delete streamlineBuffer;
	}
	streamlineBuffer 	= NULL;

}



void NIBR::TractogramReader::copyFrom(const TractogramReader& obj) {
	fileName 					= obj.fileName;
	fileDescription 			= obj.fileDescription;
	fileFormat	 				= obj.fileFormat;
	numberOfPoints 				= obj.numberOfPoints;
	numberOfStreamlines 		= obj.numberOfStreamlines;
	len 						= obj.len;
	streamlinePos 				= obj.streamlinePos;

	n_scalars_trk 		= obj.n_scalars_trk;
	n_properties_trk 	= obj.n_properties_trk;
	for(int i=0;i<4;++i) {
		for(int j=0;j<4;++j){
			ijk2xyz[i][j] = obj.ijk2xyz[i][j];
			xyz2ijk[i][j] = obj.xyz2ijk[i][j];
		}
	}

	file = fopen(fileName.c_str(), "rb");

	streamlineBuffer 	= obj.streamlineBuffer;
	usePreload 			= obj.usePreload;
	finishedLoading 	= obj.finishedLoading;
	fileBuffer 			= obj.fileBuffer;
}

void NIBR::TractogramReader::destroyCopy() {
	fileName.erase();
	fileDescription.erase();
	if (file != NULL) fclose(file);
	file 		  		= NULL;
	len 		  		= NULL;
	streamlinePos 		= NULL;
	streamlineBuffer 	= NULL;
	usePreload 			= false;
	finishedLoading 	= false;
	fileBuffer			= NULL;
}

bool NIBR::TractogramReader::initReader(std::string _fileName,bool _usePreload) {

	usePreload = _usePreload;

	if (usePreload) {

		streamlineBuffer = new std::vector<float**>();

		if (!initReader(_fileName)) 
			return false;

		return readToMemory();

	} else {
		return initReader(_fileName);
	}
	
	return true;
}

bool NIBR::TractogramReader::initReader(std::string _fileName) {

	fileName 	= _fileName;
	file 		= fopen(fileName.c_str(), "rb");

	if (file == NULL)
		return false;

	const std::size_t strLength = 256;
	char dummy[strLength];

	std::string extension = getFileExtension(fileName);

	disp(MSG_DEBUG,"File extension: %s", extension.c_str());


	if (extension == "vtp") {
        
        fileFormat      = VTP;
        fileDescription = SGNTR();

        if (!parseVTPHeaderWithExpat()) { // <-- CALL EXPAT VERSION
            fclose(file); file = NULL; return false;
        }

        // --- Check endianness and set swap flag ---
        vtp_needs_swap_ = (is_system_little_endian() != vtp_is_little_endian_);
        disp(MSG_DEBUG, "VTP Byte Order: %s, Needs Swap: %s", 
            vtp_is_little_endian_ ? "LittleEndian" : "BigEndian", 
            vtp_needs_swap_ ? "Yes" : "No");
        // ----------------------------------------

        disp(MSG_DEBUG, "VTP Parsed: Points=%zu, Lines=%zu, P_off=%ld, C_off=%ld, O_off=%ld", 
             numberOfPoints, numberOfStreamlines, points_xml_offset_, connectivity_xml_offset_, offsets_xml_offset_);

        // Read Offsets and Connectivity into memory
        if (!readVTPAppendedData(offsets_xml_offset_, vtp_offsets_)) {
            disp(MSG_FATAL, "Failed to read VTP offsets data.");
            fclose(file); file = NULL; return false;
        }
        if (!readVTPAppendedData(connectivity_xml_offset_, vtp_connectivity_)) {
            disp(MSG_FATAL, "Failed to read VTP connectivity data.");
            fclose(file); file = NULL; return false;
        }
        
        // ... (Calculate 'len' array as before) ...
        len = new uint32_t[numberOfStreamlines];
        len[0] = vtp_offsets_[0];
        for (size_t i = 1; i < numberOfStreamlines; ++i) {
            len[i] = vtp_offsets_[i] - vtp_offsets_[i - 1];
        }

        // Calculate the absolute start position for reading points
        uint64_t points_size_bytes = 0;
        fseek(file, appended_data_start_pos_ + points_xml_offset_, SEEK_SET);
        fread(&points_size_bytes, sizeof(uint64_t), 1, file);
        if (vtp_needs_swap_) swapByteOrder(points_size_bytes); // Swap size too!
        vtp_points_file_offset_ = ftell(file);

    }

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

	return true;
}

float** NIBR::TractogramReader::readStreamline(std::size_t n) {

	if (finishedLoading) {
		return streamlineBuffer->at(n);
	}

	float** points;
	points = new float* [len[n]];

	fseek(file, streamlinePos[n], SEEK_SET);

	if (fileFormat == VTP) {

        int64_t start_conn_idx = (n == 0) ? 0 : vtp_offsets_[n-1];

        for (uint32_t i = 0; i < len[n]; i++) {
            points[i] = new float[3];
            int64_t point_idx = vtp_connectivity_[start_conn_idx + i];
            if (!readVTPPoints(point_idx, points[i])) {
                // Handle error - clean up and return nullptr
                disp(MSG_ERROR, "Failed to read VTP point %lld for streamline %zu", (long long)point_idx, n);
                for(uint32_t j=0; j<=i; ++j) delete[] points[j];
                delete[] points;
                return nullptr;
            }
        }

    }

	if (fileFormat == TCK) {
		float tmp;
		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, file);
				points[i][j] = tmp;
			}
		}
	}

	if (fileFormat == TRK) {
		float tmp, p_tmp[3];
		
		for (uint32_t i = 0; i < len[n]; i++) {

			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, file);
				p_tmp[j] = tmp - 0.5f;
			}
			applyTransform(points[i],p_tmp,ijk2xyz);
			
			for (int j = 0; j < n_scalars_trk; j++)
				std::fread(&tmp, sizeof(float), 1, file);
		}
	}
	

	if (fileFormat == VTK_BINARY) {
		float tmp;
		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];

			for (int j = 0; j < 3; j++) {
				std::fread(&tmp, sizeof(float), 1, file);
				swapByteOrder(tmp);
				points[i][j] = tmp;
			}
		}
	}

	if (fileFormat == VTK_ASCII) {

		for (uint32_t i = 0; i < len[n]; i++) {
			points[i] = new float[3];
			std::fscanf(file, "%f %f %f ", &points[i][0], &points[i][1], &points[i][2]);
		}

	}

	return points;

}

void NIBR::TractogramReader::deleteStreamline(float** streamline, std::size_t n) {
	if (finishedLoading == false) {
		for (uint32_t i=0; i<len[n]; i++)
			delete[] streamline[i];
		delete[] streamline;
	}
}

std::vector<std::vector<float>> NIBR::TractogramReader::readStreamlineVector(std::size_t n) {

	std::vector<std::vector<float>> out;
	out.reserve(len[n]);

	float** tmp = readStreamline(n);

	for (uint32_t i = 0; i < len[n]; i++) {
		out.push_back({tmp[i][0],tmp[i][1],tmp[i][2]});
	}

	if (!finishedLoading) {
		for (std::size_t l = 0; l < len[n]; l++) {
			delete[] tmp[l];
		}
		delete[] tmp;
	}

	return out;

}

std::vector<Point> NIBR::TractogramReader::readStreamlinePoints(std::size_t n) {

	std::vector<Point> out;
	out.reserve(len[n]);

	float** tmp = readStreamline(n);

	for (uint32_t i = 0; i < len[n]; i++) {
		out.push_back({tmp[i][0],tmp[i][1],tmp[i][2]});
	}

	if (!finishedLoading) {
		for (std::size_t l = 0; l < len[n]; l++) {
			delete[] tmp[l];
		}
		delete[] tmp;
	}

	return out;

}

bool NIBR::TractogramReader::readToMemory() {

	streamlineBuffer->resize(numberOfStreamlines);

	finishedLoading = false;

	if (numberOfStreamlines==0) {
		finishedLoading = true;
		return true;
	}

	if (fileFormat == VTK_BINARY) {

		fileBuffer = new float[numberOfPoints * 3];

		fseek(file, streamlinePos[0], SEEK_SET);

		size_t floatsRead = std::fread(fileBuffer, sizeof(float), numberOfPoints * 3, file);

		if (floatsRead != (numberOfPoints * 3)) {
            NIBR::disp(MSG_ERROR, "Failed to read the entire tractogram file into memory. Expected %ld, got %zu floats.", (numberOfPoints * 3), floatsRead);
			delete[] fileBuffer;
            finishedLoading = false;
            return false;
        }

		for (size_t i = 0; i < (numberOfPoints * 3); ++i) {
			swapByteOrder(fileBuffer[i]);
		}

		float s = 100.0f / float(numberOfStreamlines);
		float* currentPtr = fileBuffer;

		for (size_t n = 0; n < numberOfStreamlines; n++) {
			streamlineBuffer->at(n) = new float*[len[n]];
			for (uint32_t l = 0; l < len[n]; l++) {
				streamlineBuffer->at(n)[l] = currentPtr;
				currentPtr += 3;
			}
			if ((n>0) && (NIBR::VERBOSE()>=VERBOSE_INFO)) std::cout << "\033[A\r\033[K" << std::flush;
			NIBR::disp(MSG_INFO,"Preloading streamlines: %.2f%%", (n+1)*s);
		}
		if (NIBR::VERBOSE()>=VERBOSE_INFO) std::cout << "\033[A\r\033[K" << std::flush;
		NIBR::disp(MSG_INFO,"Preloading streamlines: 100%%");
		
	} else {

		float s = 100.0f / float(numberOfStreamlines);

		for (size_t n = 0; n < numberOfStreamlines; n++) {
			streamlineBuffer->at(n) = readStreamline(n);
			if ((n>0) && (NIBR::VERBOSE()>=VERBOSE_INFO)) std::cout << "\033[A\r\033[K" << std::flush;
			NIBR::disp(MSG_INFO,"Preloading streamlines: %.2f%%", (n+1)*s);
		}
		if (NIBR::VERBOSE()>=VERBOSE_INFO) std::cout << "\033[A\r\033[K" << std::flush;
		NIBR::disp(MSG_INFO,"Preloading streamlines: 100%%");

	}

	finishedLoading = true;
	
	return true;

}

bool NIBR::TractogramReader::readNextBatch(size_t batchSize, std::vector<std::vector<std::vector<float>>>& batch_out) {

    if (usePreload) {
        disp(MSG_WARN, "readNextBatch called with preloading enabled. Consider disabling preload for memory efficiency.");
    }

    if (file == NULL && !usePreload) {
        disp(MSG_ERROR, "Reader not initialized or file closed.");
        return false;
    }

    if (currentStreamlineIdx >= numberOfStreamlines) {
        return false;
    }

    batch_out.clear();
    batch_out.reserve(batchSize);

    size_t streamlinesRead = 0;
    while (streamlinesRead < batchSize && currentStreamlineIdx < numberOfStreamlines) {
        
        batch_out.push_back(readStreamlineVector(currentStreamlineIdx));
        
        currentStreamlineIdx++;
        streamlinesRead++;
    }

    return (streamlinesRead > 0);
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

	std::vector<NIBR::TractogramField> fields = findTractogramFields(*this);

	std::cout << "Field count: " << fields.size() << std::endl << std::flush;

	int i = 1;
	for (auto f : fields) {
		std::cout << "   " << i++ << ". " << f.name.c_str() << ": " << ((f.owner == POINT_OWNER) ? "Point" : "Streamline") << " field with dim " << f.dimension << std::endl << std::flush;
	}

	std::cout << "\033[0m";

	return;


}

// VTP reader functions

// Structure to hold state and results during Expat parsing
struct VTPExpatUserData {
    NIBR::TractogramReader* reader; 
    std::string current_DataArray_Name;
    bool        parsing_ok;
    bool        found_appended_tag;
    XML_Parser  parser; // <-- Holds the parser instance
};

// Helper to process attributes - NOW TAKES VTPExpatUserData*
static void process_attributes(VTPExpatUserData* data, const XML_Char *name, const XML_Char **atts) {
    
    NIBR::TractogramReader* reader = data->reader; // Get reader from data
    XML_Parser parser = data->parser;             // Get parser from data

    long long* num_points_ptr   = reinterpret_cast<long long*>(&reader->numberOfPoints);
    long long* num_lines_ptr    = reinterpret_cast<long long*>(&reader->numberOfStreamlines);
    long* points_off_ptr   = &reader->points_xml_offset_;
    long* conn_off_ptr     = &reader->connectivity_xml_offset_;
    long* offs_off_ptr     = &reader->offsets_xml_offset_;
    bool* is_little_ptr    = &reader->vtp_is_little_endian_;
    std::string* current_name_ptr = &data->current_DataArray_Name;

    if (strcmp(name, "VTKFile") == 0) {
        for (int i = 0; atts[i]; i += 2) {
            if (strcmp(atts[i], "byte_order") == 0) {
                *is_little_ptr = (strcmp(atts[i+1], "BigEndian") != 0);
            }
        }
    } else if (strcmp(name, "Piece") == 0) {
        for (int i = 0; atts[i]; i += 2) {
            if (strcmp(atts[i], "NumberOfPoints") == 0) *num_points_ptr = std::stoll(atts[i+1]);
            if (strcmp(atts[i], "NumberOfLines") == 0)  *num_lines_ptr  = std::stoll(atts[i+1]);
        }
    } else if (strcmp(name, "DataArray") == 0) {
        std::string format = "";
        long offset = -1;
        std::string arr_name = "";

        for (int i = 0; atts[i]; i += 2) {
            if (strcmp(atts[i], "Name") == 0)   arr_name = atts[i+1];
            if (strcmp(atts[i], "format") == 0) format = atts[i+1];
            if (strcmp(atts[i], "offset") == 0) offset = std::stoll(atts[i+1]);
        }

        if (format != "appended") {
            NIBR::disp(MSG_FATAL, "VTP Reader currently only supports format=\"appended\". Found %s.", format.c_str());
            data->parsing_ok = false;
            XML_StopParser(parser, XML_FALSE); // <-- USE LOCAL PARSER
            return;
        }

        *current_name_ptr = arr_name;

        if (arr_name == "Points")       *points_off_ptr = offset;
        if (arr_name == "connectivity") *conn_off_ptr   = offset;
        if (arr_name == "offsets")      *offs_off_ptr   = offset;
        
    } else if (strcmp(name, "AppendedData") == 0) {
        data->found_appended_tag = true;
        XML_StopParser(parser, XML_TRUE); // <-- USE LOCAL PARSER
    }
}

// Expat Start Element Handler - NOW PASSES 'data' TO HELPER
static void XMLCALL startElementHandler(void *userData, const XML_Char *name, const XML_Char **atts) {
    VTPExpatUserData* data = reinterpret_cast<VTPExpatUserData*>(userData);
    
    if (!data->parsing_ok || data->found_appended_tag) return;

    data->reader->vtp_current_tag_ = name;
    process_attributes(data, name, atts); // <-- Pass 'data'
}

// Expat End Element Handler
static void XMLCALL endElementHandler(void *userData, const XML_Char *name) {
    VTPExpatUserData* data = reinterpret_cast<VTPExpatUserData*>(userData);
    data->reader->vtp_current_tag_ = "";
}

// Expat Character Data Handler
static void XMLCALL characterDataHandler(void *userData, const XML_Char *s, int len) {
    // No action needed
}

// Basic VTP XML header parser.
bool NIBR::TractogramReader::parseVTPHeaderWithExpat() {

    XML_Parser parser = XML_ParserCreate(NULL);
    if (!parser) {
        disp(MSG_FATAL, "Couldn't allocate memory for Expat parser.");
        return false;
    }

    VTPExpatUserData userData;
    userData.reader = this;
    userData.parsing_ok = true;
    userData.found_appended_tag = false;
    userData.parser = parser; // <-- Store parser in userData

    XML_SetUserData(parser, &userData);
    XML_SetElementHandler(parser, startElementHandler, endElementHandler);
    XML_SetCharacterDataHandler(parser, characterDataHandler);

    fseek(file, 0, SEEK_SET);

    const int buffer_size = 4096;
    char buffer[buffer_size];
    bool done = false;

    disp(MSG_DEBUG, "Starting Expat VTP header parsing...");

    do {
        int bytes_read = fread(buffer, 1, buffer_size, file);
        done = bytes_read < buffer_size;

        if (XML_Parse(parser, buffer, bytes_read, done) == XML_STATUS_ERROR) {
            if (XML_GetErrorCode(parser) != XML_ERROR_ABORTED) {
                disp(MSG_FATAL, "Expat XML parse error: %s at line %lu",
                    XML_ErrorString(XML_GetErrorCode(parser)),
                    XML_GetCurrentLineNumber(parser));
                userData.parsing_ok = false;
            }
            break; 
        }

    } while (!done && userData.parsing_ok && !userData.found_appended_tag);

    long expat_index   = XML_GetCurrentByteIndex(parser);

    XML_ParserFree(parser); // Free the parser

    if (!userData.parsing_ok) return false;
    if (!userData.found_appended_tag) {
        disp(MSG_FATAL, "VTP parsing finished without finding <AppendedData> tag.");
        return false;
    }

    // Find the '_' marker starting from where Expat stopped.
    fseek(file, expat_index, SEEK_SET);
    char c;
    while ((c = fgetc(file)) != EOF) {
        if (c == '_') {
            appended_data_start_pos_ = ftell(file);
            disp(MSG_DEBUG, "Found AppendedData marker '_' at position %ld", appended_data_start_pos_);
            return true;
        }
    }

    disp(MSG_FATAL, "Could not find '_' marker after <AppendedData> tag.");
    return false;
}


bool NIBR::TractogramReader::readVTPPoints(int64_t index, float* point) {
    if (vtp_points_file_offset_ <= 0) return false; // Use member var
    
    long pos = vtp_points_file_offset_ + (index * 3 * sizeof(float)); // Use member var
    fseek(file, pos, SEEK_SET);
    
    if (fread(point, sizeof(float), 3, file) != 3) return false;

    if (vtp_needs_swap_) { // Use member var
        swapByteOrder(point[0]);
        swapByteOrder(point[1]);
        swapByteOrder(point[2]);
    }
    return true;
}