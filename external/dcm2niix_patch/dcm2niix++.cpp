#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#define strtok_portable(str, delim, context) strtok_s(str, delim, context)

bool isPathFolder(const std::string& dirPath) {
    DWORD attributes = GetFileAttributesA(dirPath.c_str());
    return (attributes != INVALID_FILE_ATTRIBUTES && (attributes & FILE_ATTRIBUTE_DIRECTORY));
}

std::string getFirstFilenameInDirectory(const std::string& dirPath) {
    std::string searchPath = dirPath + "\\*";
    WIN32_FIND_DATAA findFileData;
    HANDLE hFind = FindFirstFileA(searchPath.c_str(), &findFileData);

    if (hFind == INVALID_HANDLE_VALUE) {
        return "";
    }

    std::string firstFilename;
    do {
        if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            firstFilename = findFileData.cFileName;
            break;
        }
    } while (FindNextFileA(hFind, &findFileData) != 0);

    FindClose(hFind);
    return firstFilename;
}

#else
#include <dirent.h>
#define strtok_portable(str, delim, context) strtok_r(str, delim, context)

bool isPathFolder(const std::string& dirPath) {
    DIR* dir = opendir(dirPath.c_str());
    if (dir != nullptr) {
        closedir(dir);
        return true;
    }
    return false;
}

std::string getFirstFilenameInDirectory(const std::string& dirPath) {
    DIR* dir = opendir(dirPath.c_str());
    if (dir == nullptr) {
        return "";
    }

    std::string firstFilename;
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_REG) {
            firstFilename = entry->d_name;
            break;
        }
    }

    closedir(dir);
    return firstFilename;
}

#endif

std::string dirname(const std::string& path) {
    if (path.empty()) return ".";

    std::size_t lastSlashPos = path.find_last_of("/\\");

    if (lastSlashPos != std::string::npos) {
        return path.substr(0, lastSlashPos + 1);
    }

	// Return current path
	#ifdef _WIN32
	return ".\\";
	#endif

	return "./";
    
}

std::string basename(const std::string& path) {
    if (path.empty()) return "";
    
    std::size_t lastSlashPos = path.find_last_of("/\\");

    if (lastSlashPos == std::string::npos) {
        return path; // No slash found, the whole path is the basename
    }
    
	return path.substr(lastSlashPos + 1);
}

#include "nii_dicom.h"
#include "nii_dicom_batch.h"

// #include "../dcm2niix/console/nii_dicom.h"
// #include "../dcm2niix/console/nii_dicom_batch.h"

#include "dcm2niix++.h"

dcm2niix::dcm2niix()
{
	TDCMopts* tmp = new TDCMopts();
	setDefaultOpts(tmp,NULL);
	opts = static_cast<void*>(tmp);
	updateFullPath();
}

dcm2niix::~dcm2niix()
{
	clear();
}

void dcm2niix::updateFullPath()
{

	std::string modifiedFolderPath = folderPath;
    if (!modifiedFolderPath.empty() && (modifiedFolderPath.back() == '/' || modifiedFolderPath.back() == '\\')) {
        modifiedFolderPath.pop_back();
    }

    if (!modifiedFolderPath.empty()) {
       #ifdef _WIN32
       fullPath = modifiedFolderPath + "\\" + baseFileName;
       #else
       fullPath = modifiedFolderPath + "/" + baseFileName;
       #endif 
    } else if (modifiedFolderPath.empty()  && !baseFileName.empty()) {
	   #ifdef _WIN32
       fullPath = ".\\" + baseFileName;
       #else
       fullPath = "./" + baseFileName;
       #endif
    } else {
       fullPath = ".";
	}

	if (fullPath_charp != NULL) {
		delete[] fullPath_charp;
		fullPath_charp = NULL;
	}
	
	fullPath_charp 	= new char[fullPath.length() + 1];
	strcpy(fullPath_charp, fullPath.c_str());
}

bool dcm2niix::setInputPath(std::string inputPath)
{
	isFolder = isPathFolder(inputPath);
	
	if (isFolder) {
		folderPath 	 = inputPath;
		baseFileName = getFirstFilenameInDirectory(folderPath);
	} else {
		folderPath 	 = dirname(inputPath);
		baseFileName = basename(inputPath);
	}
	
	updateFullPath();
	
	return isDICOMfile(fullPath_charp);
}

bool dcm2niix::toNii()
{
	if (fullPath == ".") return false;

	std::string opts_str = "v=0";
	this->setOpts(opts_str.c_str());

	// std::cout << "readDICOM: " << fullPath_charp << std::endl << std::flush;
	struct TDICOMdata tdicomData = readDICOM(fullPath_charp);
	// std::cout << "readDICOM...Done" << std::endl << std::flush;

	TDCMopts& tdcmOpts = *(TDCMopts*)(opts);

	double seriesNo = (double)tdicomData.seriesUidCrc;
	if (tdcmOpts.isIgnoreSeriesInstanceUID)
		seriesNo = (double)tdicomData.seriesNum;

	// set TDCMopts to convert just one series
	tdcmOpts.seriesNumber[0] = seriesNo;
	tdcmOpts.numSeries 	 	 = 1;
	tdcmOpts.isCreateBIDS    = false;
	tdcmOpts.isOnlyBIDS      = false;

	int out;

	if (isFolder) {
		out = nii_loadDirCore(tdcmOpts.indir, &tdcmOpts);
	} else {
		out = singleDICOM(&tdcmOpts, fullPath_charp);
	}

	return ( out == EXIT_SUCCESS); 
}

// return nifti header saved in MRIFSSTRUCT
nifti_1_header dcm2niix::getNiiHeader() {
	return nii_getMrifsStruct()->hdr0;
}

// return image data saved in MRIFSSTRUCT
const unsigned char* dcm2niix::getMRIimg() {
	return nii_getMrifsStruct()->imgM;
}



std::vector<float> dcm2niix::getbvals()
{
	auto gradTable = nii_getMrifsStruct()->tdti;
	int  dirCount  = nii_getMrifsStruct()->numDti;
	
	std::vector<float> out;
	out.reserve(dirCount);

	for (int n = 0; n < dirCount; n++) {
		out.push_back(gradTable[n].V[0]);
	}

	return out;
}

std::vector<std::array<float, 3>> dcm2niix::getbvecs()
{
	auto gradTable = nii_getMrifsStruct()->tdti;
	int  dirCount  = nii_getMrifsStruct()->numDti;
	
	std::vector<std::array<float, 3>> out;
	out.reserve(dirCount);

	for (int n = 0; n < dirCount; n++) {
		std::array<float,3> tmp;

		tmp[0] = gradTable[n].V[1];
		tmp[1] = gradTable[n].V[2];
		tmp[2] = gradTable[n].V[3];

		out.push_back(tmp);
	}

	return out;
}


void dcm2niix::clear() {
	// std::cout << "clearing: " << std::endl << std::flush;
	if (opts != NULL) {
		delete static_cast<TDCMopts*>(opts);
		opts = NULL;
	}
	// std::cout << "cleared opts " << std::endl << std::flush;

	baseFileName.clear();
	folderPath.clear();
	fullPath.clear();
	// std::cout << "cleared strings " << std::endl << std::flush;
	
	if (fullPath_charp != NULL) {
		delete[] fullPath_charp;
		fullPath_charp = NULL;
	}

	// std::cout << "cleared fullPath_charp " << std::endl << std::flush;

	nii_clrMrifsStructVector();
	// std::cout << "cleared nii_clrMrifsStructVector " << std::endl << std::flush;

	nii_clrMrifsStruct();
	// std::cout << "reseted nii_clrMrifsStruct " << std::endl << std::flush;
	

	// std::cout << "Done" << std::endl << std::flush;
}


void dcm2niix::setOpts(const char *dcm2niixopts) {

	// std::cout << "Setting options" << std::endl << std::flush;

	TDCMopts& tdcmOpts = *static_cast<TDCMopts*>(opts);

	// std::cout << "strcpy(tdcmOpts.indir, folderPath.c_str())" << std::endl << std::flush;
	// std::cout << "size: " << sizeof(tdcmOpts.indir) << std::endl << std::flush;
	// std::cout << "folderPath: " << folderPath << std::endl << std::flush;
	// std::cout << "folderPath_char: " << folderPath.c_str() << std::endl << std::flush;
	// std::cout << "tdcmOpts.indir: " << tdcmOpts.indir << std::endl << std::flush;

	strcpy(&tdcmOpts.indir[0], folderPath.c_str());

	// dcmunpack actually uses seriesDescription, set FName = `printf %04d.$descr $series`
	// change it from "%4s.%p" to "%4s.%d"
	// std::cout << "strcpy(tdcmOpts.filename)" << std::endl << std::flush;
	strcpy(&tdcmOpts.filename[0], "%4s.%d");

	if (dcm2niixopts != NULL) {

		// std::cout << "Setting dcm2niixopts options" << std::endl << std::flush;

		char *restOpts = (char *)malloc(strlen(dcm2niixopts) + 1);
		memset(restOpts, 0, strlen(dcm2niixopts) + 1);
		memcpy(restOpts, dcm2niixopts, strlen(dcm2niixopts));

		char *nextOpt = strtok_portable((char *)dcm2niixopts, ",", &restOpts);
		while (nextOpt != NULL) {
			char *k = nextOpt;
			char *v = strchr(nextOpt, '=');
			if (v != NULL)
				*v = '\0';
			v++; // move past '='

			// skip leading white spaces
			while (*k == ' ')
				k++;

			if (strcmp(k, "b") == 0) {
				if (*v == 'n' || *v == 'N' || *v == '0')
					tdcmOpts.isCreateBIDS = false;
				else if (*v == 'i' || *v == 'I') {
					tdcmOpts.isCreateBIDS = false;
					tdcmOpts.isOnlyBIDS = true;
				} else if (*v == 'y' || *v == 'Y')
					tdcmOpts.isCreateBIDS = true;
			} else if (strcmp(k, "ba") == 0)
				tdcmOpts.isAnonymizeBIDS = (*v == 'n' || *v == 'N') ? false : true;
			else if (strcmp(k, "f") == 0)
				strcpy(&tdcmOpts.filename[0], v);
			else if (strcmp(k, "i") == 0)
				tdcmOpts.isIgnoreDerivedAnd2D = (*v == 'y' || *v == 'Y') ? true : false;
			else if (strcmp(k, "m") == 0) {
				if (*v == 'n' || *v == 'N' || *v == '0')
					tdcmOpts.isForceStackSameSeries = 0;
				else if (*v == 'y' || *v == 'Y' || *v == '1')
					tdcmOpts.isForceStackSameSeries = 1;
				else if (*v == '2')
					tdcmOpts.isForceStackSameSeries = 2;
				else if (*v == 'o' || *v == 'O')
					tdcmOpts.isForceStackDCE = false;
			} else if (strcmp(k, "v") == 0) {
				if (*v == 'n' || *v == 'N' || *v == '0')
					tdcmOpts.isVerbose = 0;
				else if (*v == 'h' || *v == 'H' || *v == '2')
					tdcmOpts.isVerbose = 2;
				else
					tdcmOpts.isVerbose = 1;
			} else if (strcmp(k, "o") == 0)
				strcpy(&tdcmOpts.outdir[0], v);
			else if (strcmp(k, "t") == 0)
				tdcmOpts.isCreateText = (*v == 'y' || *v == 'Y') ? true : false;
			else if (strcmp(k, "p") == 0) {
				if (*v == 'n' || *v == 'N' || *v == '0')
					tdcmOpts.isPhilipsFloatNotDisplayScaling = false;
			} else if (strcmp(k, "x") == 0) {
				if (*v == 'y' || *v == 'Y')
					tdcmOpts.isCrop = true;
				else if (*v == 'i' || *v == 'I') {
					tdcmOpts.isRotate3DAcq = false;
					tdcmOpts.isCrop = false;
				}
			} else {
				printf("[WARN] dcm2niix option %s=%s skipped\n", k, v);
			}

			nextOpt = strtok_portable(NULL, ",", &restOpts);
		}

		// std::cout << "Setting dcm2niixopts options... Done" << std::endl << std::flush;
	}

	// std::cout << "Setting options... Done" << std::endl << std::flush;
}



