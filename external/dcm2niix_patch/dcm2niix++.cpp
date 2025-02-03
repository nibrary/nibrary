#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <iostream>
#include <string.h>
#include <windows.h>

char* dirname(char* path) {
    if (path == nullptr || strlen(path) == 0) {
        return path; // Or return "." for an empty path
    }

    char* lastSlash = strrchr(path, '\\'); // Find the last backslash

    if (lastSlash != nullptr) {
        *lastSlash = '\0'; // Terminate the string at the last backslash
        return path;
    } else {
        return "."; // No backslash found, return current directory "."
    }
}

#define strtok_portable(str, delim, context) strtok_s(str, delim, context)

#else
#include <libgen.h> // For dirname on non-Windows systems
#define strtok_portable(str, delim, context) strtok_r(str, delim, context)
#endif


#include "nii_dicom.h"
#include "nii_dicom_batch.h"

// #include "../dcm2niix/console/nii_dicom.h"
// #include "../dcm2niix/console/nii_dicom_batch.h"

#include "dcm2niix++.h"

dcm2niix::dcm2niix()
{
	TDCMopts* tmp = new TDCMopts();
	void* nibr_tdcmopts = (void*)(tmp);
	fileName = "";
}

dcm2niix::~dcm2niix()
{
	clear();
}

bool dcm2niix::setFileName(std::string inputFileName)
{
	fileName = inputFileName;
	return isDICOMfile(fileName.c_str());
}

bool dcm2niix::toNii()
{
	if (fileName == "") return false;

	std::string opts_str = "v=0";
	this->setOpts(fileName.c_str(),opts_str.c_str());

	struct TDICOMdata tdicomData = readDICOM((char *)(fileName.c_str()));

	TDCMopts& tdcmOpts = *(TDCMopts*)(opts);

	double seriesNo = (double)tdicomData.seriesUidCrc;
	if (tdcmOpts.isIgnoreSeriesInstanceUID)
		seriesNo = (double)tdicomData.seriesNum;

	// set TDCMopts to convert just one series
	tdcmOpts.seriesNumber[0] = seriesNo;
	tdcmOpts.numSeries = 1;

	return (nii_loadDirCore(tdcmOpts.indir, &tdcmOpts) == EXIT_SUCCESS);
}

// return nifti header saved in MRIFSSTRUCT
nifti_1_header* dcm2niix::getNiiHeader() {
	MRIFSSTRUCT* mrifsStruct = nii_getMrifsStruct();
	return &mrifsStruct->hdr0;
}

// return image data saved in MRIFSSTRUCT
const unsigned char *dcm2niix::getMRIimg() {
	MRIFSSTRUCT* mrifsStruct = nii_getMrifsStruct();
	return mrifsStruct->imgM;
}

void dcm2niix::clear() {
	if (opts != NULL) {
		delete (TDCMopts*)(opts);
		opts = NULL;
	}
	fileName = "";
	nii_clrMrifsStruct();
	nii_clrMrifsStructVector();
}


void dcm2niix::setOpts(const char *dcmindir, const char *dcm2niixopts) {

	TDCMopts& tdcmOpts = *(TDCMopts*)(opts);

	setDefaultOpts(&tdcmOpts, NULL);

	if (dcmindir != NULL)
		strcpy(tdcmOpts.indir, dcmindir);

	// dcmunpack actually uses seriesDescription, set FName = `printf %04d.$descr $series`
	// change it from "%4s.%p" to "%4s.%d"
	strcpy(tdcmOpts.filename, "%4s.%d");

	if (dcm2niixopts != NULL) {

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
				strcpy(tdcmOpts.filename, v);
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
				strcpy(tdcmOpts.outdir, v);
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


	}

}



