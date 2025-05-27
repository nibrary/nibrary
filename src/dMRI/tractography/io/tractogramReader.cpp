#include "dMRI/tractography/io/tractogramReader.h"
#include "tractogramField.h"
#include <cstdint>
#include <cstring>
#include <iterator>

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
