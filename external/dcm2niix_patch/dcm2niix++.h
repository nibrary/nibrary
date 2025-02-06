#ifndef DCM2NIIXPP_H
#define DCM2NIIXPP_H

#ifdef USING_DCM2NIIXFSWRAPPER
#include "nifti1.h"
#else
#include "niftilib/nifti2/nifti1.h"
#endif

#include <string>
#include <vector>
#include <array>

class dcm2niix {

  	public:

		dcm2niix();		
		~dcm2niix();

		bool setInputPath(std::string path);
		bool toNii();

		nifti_1_header 	 	    			getNiiHeader();
		const unsigned char*    			getMRIimg();

		std::vector<float>      			getbvals();		// not tested!!
		std::vector<std::array<float, 3>>   getbvecs(); 	// not tested!!

		void clear();

	private:

		bool        isFolder{false};
		std::string baseFileName{""};
		std::string folderPath{"."};
		std::string fullPath{"."};
		char* 		fullPath_charp{NULL};
		void* 		opts{NULL};	

		void 		updateFullPath();
		void 		setOpts(const char *dcm2niixopts = NULL);

};

#endif
