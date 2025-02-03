#ifndef DCM2NIIXPP_H
#define DCM2NIIXPP_H

#ifdef USING_DCM2NIIXFSWRAPPER
#include "nifti1.h"
#else
#include "niftilib/nifti2/nifti1.h"
#endif

#include <string>

class dcm2niix {

  	public:

		dcm2niix();		
		~dcm2niix();

		bool setFileName(std::string inputFileName);
		bool toNii();

		nifti_1_header* 	 getNiiHeader();
		const unsigned char* getMRIimg();

		void clear();

	private:

		std::string fileName;
		void* opts;

		void setOpts(const char *dcmindir, const char *dcm2niixopts = NULL);

};

#endif
