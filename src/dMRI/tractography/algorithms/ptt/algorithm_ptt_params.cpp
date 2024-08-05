#include "../../tracker/tracker.h"
#include "algorithm_ptt_params.h"
#include "math/disc.h"

using namespace NIBR;


bool Params_PTT::update() {

	if (pttIsReady)
		return true;

	if ((img_FOD==NULL) || ((img_FOD!=NULL) &&  ( (img_FOD->filePath != img_FOD_path) || (img_FOD->getSphereFileName() != fod_sphere_path) ))) {

		if (img_FOD!=NULL) delete img_FOD;

		if (fod_sphere_path=="")
			img_FOD 			= new FOD_Image(img_FOD_path);
		else
			img_FOD 			= new FOD_Image(img_FOD_path,fod_sphere_path);

		orderOfDirections = convertOrderOfDirections(orderOfDirectionsTextInput);

		img_FOD->setOrderOfDirections (orderOfDirections);
		img_FOD->setdiscretizationFlag(fodDiscretization);

		if (fodIsSym) 
			img_FOD->fodIsSym(); 
		else 
			img_FOD->fodIsAsym();

		if(!img_FOD->read()) {
			disp(MSG_ERROR, "Can't read FOD image");
			return false;		
		}

		smallestPixDim = img_FOD->smallestPixDim;

		// if (img_FOD->imgDims[0]>1)
		//	return false;
		
	}

	setDefaults();

	pttIsReady = true;
	
	return true;

}

void Params_PTT::clear() {

	pttIsReady = false;

	// FOD internal
	if (img_FOD!=NULL) {
		delete img_FOD;
		img_FOD = NULL;
	}

	fodSphere.clear();
	fodIsSpheresliced           = false;
	fodIsSym					= true;
	smallestPixDim              = NAN;
	orderOfDirections           = XYZ;

	// FOD options in Trekker
	img_FOD_path                = "";
	fod_sphere_path         	= "";
	fodDiscretization           = true;	
	orderOfDirectionsTextInput  = "";
	
	// Tracking options
	stepSize                    = NAN;
	minRadiusOfCurvature        = NAN;
	minDataSupport              = NAN;
	dataSupportExponent         = NAN;
	weakLinkThresh              = NAN;

	// Sampling options
	maxEstInterval              = -1;
	initMaxEstTrials            = -1;
	propMaxEstTrials            = -1;
	triesPerRejectionSampling   = -1;
	useBestAtInit               = false;
	useLegacySampling           = false;
	samplingQuality             = -1;

	// Probe options
	probeCount                  = NAN;
    probeRadius                 = NAN;
    probeLength                 = NAN;
	probeQuality                = NAN;

	// Output options
	saveFrame                   = false;
	outputStep                  = NAN;

	// Derived parameters
    maxCurvature                = NAN;
	checkWeakLinks              = NAN;
	modMinDataSupport           = NAN;
}



void Params_PTT::setDefaults() {

	// Handle stepSize and outputStep
	if (std::isnan(stepSize) 			 || (stepSize<=0)) 				stepSize 			 = smallestPixDim*DEFAULT_PTT_STEPSIZE_IN_PIXELDIM;
	if (std::isnan(outputStep) 			 || (outputStep<=0)) 			outputStep 			 = smallestPixDim*DEFAULT_PTT_OUTPUTSTEPSIZE_IN_PIXELDIM;

	// Handle minRadiusOfCurvature
	if (std::isnan(minRadiusOfCurvature) || (minRadiusOfCurvature<=0))	minRadiusOfCurvature = smallestPixDim*DEFAULT_PTT_MINRADIUSOFCURVATURE_IN_PIXELDIM;

	// Handle minDataSupport
	if (std::isnan(minDataSupport) 		 || (minDataSupport<0.0)) 		minDataSupport 	   	 = DEFAULT_PTT_MINDATASUPPORT;

	// Handle dataSupportExponent
	if (std::isnan(dataSupportExponent)  || (dataSupportExponent<0.0))  dataSupportExponent  = DEFAULT_PTT_DATASUPPORTEXPONENT;

	// Handle weak link checking
	if (std::isnan(weakLinkThresh) 		 || (weakLinkThresh<=0)) 		weakLinkThresh 		 = DEFAULT_PTT_WEAKLINKTHRESH;

	// Handle maxEstInterval, initMaxEstTrials, propMaxEstTrials, triesPerRejectionSampling and samplingQuality
	if (maxEstInterval<=0.0)			maxEstInterval 			  = DEFAULT_PTT_MAXESTINTERVAL;
	if (initMaxEstTrials<=0.0)			initMaxEstTrials 		  = DEFAULT_PTT_INITMAXESTTRIALS;
	if (propMaxEstTrials<=0.0)			propMaxEstTrials 		  = DEFAULT_PTT_PROPMAXESTTRIALS;
	if (triesPerRejectionSampling<=0.0)	triesPerRejectionSampling = DEFAULT_PTT_TRIESPERREJECTIONSAMPLING;
	if (samplingQuality<=0.0)			samplingQuality 		  = DEFAULT_PTT_SAMPLINGQUALITY;

	// Handle probeCount and probeRadius
	if (probeCount<1)   probeCount = NAN;
	if (probeRadius<0)  probeRadius= NAN;

	if (isnan(probeCount) && isnan(probeRadius)) {
		probeCount = DEFAULT_PTT_PROBECOUNT;
		probeRadius= DEFAULT_PTT_PROBERADIUS_IN_PIXELDIM*smallestPixDim;
	} else if (isnan(probeCount) && probeRadius!=NAN) {
		probeCount = DEFAULT_PTT_PROBECOUNT_WHEN_THEREIS_PROBERADIUS;
	} else if (!isnan(probeCount) && isnan(probeRadius)) {
		probeRadius= DEFAULT_PTT_PROBERADIUS_IN_PIXELDIM*smallestPixDim;
	}

	if (probeCount==1) {
		probeRadius=0;
	} else {
		if (probeRadius==0) {
			probeCount=1;
		} else {
			if (probeRadius>minRadiusOfCurvature) {
				probeRadius = minRadiusOfCurvature;
			}
		}
	}

	// Handle probeLength and probeQuality
	if (std::isnan(probeLength)  || (probeLength<=0.0))							probeLength  = smallestPixDim*DEFAULT_PTT_PROBELENGTH_IN_PIXELDIM;
	if (std::isnan(probeQuality) || (probeQuality<=1.0) || (probeQuality>100))	probeQuality = DEFAULT_PTT_PROBEQUALITY;


	// Derived parameters
	maxCurvature = 1/minRadiusOfCurvature;
	if (maxCurvature<1e-4) maxCurvature = 1e-4;

	checkWeakLinks = (weakLinkThresh>0) ? true : false;
	modMinDataSupport = std::pow(minDataSupport,dataSupportExponent);

	// Prep CDF domain
	auto fill_cdfk1k2 = [&](auto& verts) {
		cdfk1k2.reserve(cdfVertCnt);
		for (int n = 0; n < cdfVertCnt; n++) {
			float k1 = verts[n][0] * maxCurvature;
			float k2 = verts[n][1] * maxCurvature;
			cdfk1k2.push_back(std::make_pair(k1,k2));
		}
	};

	auto fill_cdfFace = [&](auto& faces) {
		cdfFace.reserve(cdfFaceCnt);
		for (int n = 0; n < cdfFaceCnt; n++) {
			cdfFace.push_back({faces[n][0],faces[n][1],faces[n][2]});
		}
	};



	switch (samplingQuality) {
		case 1: {cdfVertCnt = DISC_1_VERT_CNT; cdfFaceCnt = DISC_1_FACE_CNT; fill_cdfk1k2(DISC_1_VERT); fill_cdfFace(DISC_1_FACE); break;}
		case 2: {cdfVertCnt = DISC_2_VERT_CNT; cdfFaceCnt = DISC_2_FACE_CNT; fill_cdfk1k2(DISC_2_VERT); fill_cdfFace(DISC_2_FACE); break;}
		case 3: {cdfVertCnt = DISC_3_VERT_CNT; cdfFaceCnt = DISC_3_FACE_CNT; fill_cdfk1k2(DISC_3_VERT); fill_cdfFace(DISC_3_FACE); break;}
		case 4: {cdfVertCnt = DISC_4_VERT_CNT; cdfFaceCnt = DISC_4_FACE_CNT; fill_cdfk1k2(DISC_4_VERT); fill_cdfFace(DISC_4_FACE); break;}
		case 5: {cdfVertCnt = DISC_5_VERT_CNT; cdfFaceCnt = DISC_5_FACE_CNT; fill_cdfk1k2(DISC_5_VERT); fill_cdfFace(DISC_5_FACE); break;}
		case 6: {cdfVertCnt = DISC_6_VERT_CNT; cdfFaceCnt = DISC_6_FACE_CNT; fill_cdfk1k2(DISC_6_VERT); fill_cdfFace(DISC_6_FACE); break;}
		case 7: {cdfVertCnt = DISC_7_VERT_CNT; cdfFaceCnt = DISC_7_FACE_CNT; fill_cdfk1k2(DISC_7_VERT); fill_cdfFace(DISC_7_FACE); break;}
		default: break;
	}

}




void Params_PTT::print() {

	if (NIBR::VERBOSE()<VERBOSE_INFO) {
        return;
    }

	disp(MSG_INFO,"PTT OPTIONS");

	std::cout << "\033[32m";
    
    std::cout << "algorithm            : parallel transport tracker (ptt)"  << std::endl;
	std::cout << "fod                  : "  << img_FOD_path << std::endl;
	if (fodDiscretization)
		std::cout << "fodDiscretization    : ON "  << std::endl;
	else
		std::cout << "fodDiscretization    : OFF " << std::endl;

	if (fodIsSpheresliced) {
		std::cout << "fod spherical domain : "  << fod_sphere_path << std::endl;
		if (fodIsSym)
			std::cout << "fodIsSym             : ON "  << std::endl;
		else
			std::cout << "fodIsSym             : OFF " << std::endl;
	}

	std::cout << "stepSize             : "  << to_string_with_precision(stepSize,4) << std::endl;
	std::cout << "writeStepSize        : "  << to_string_with_precision(outputStep,4) << std::endl;
	std::cout << "minRadiusOfCurvature : "  << to_string_with_precision(minRadiusOfCurvature,4) << std::endl;

	std::cout << "minDataSupport       : "  << to_string_with_precision(minDataSupport,4) << std::endl;
	std::cout << "dataSupportExponent  : "  << to_string_with_precision(dataSupportExponent,4) << std::endl;

    std::cout << "ignoreWeakLinks      : "  << to_string_with_precision(weakLinkThresh,4) << std::endl;
    std::cout << "maxEstInterval       : "  << maxEstInterval 				<< std::endl;
	std::cout << "maxSamplingPerStep   : "  << triesPerRejectionSampling 	<< std::endl;
	std::cout << "initMaxEstTrials     : "  << initMaxEstTrials 			<< std::endl;
	std::cout << "propMaxEstTrials     : "  << propMaxEstTrials 			<< std::endl;

	if (useBestAtInit)
		std::cout << "useBestAtInit        : ON "  << std::endl;
	else
		std::cout << "useBestAtInit        : OFF " << std::endl;

	if (useLegacySampling)
		std::cout << "useLegacySampling    : ON "  << std::endl;
	else
		std::cout << "useLegacySampling    : OFF " << std::endl;

	std::cout << "samplingQuality      : "  << samplingQuality 			    << std::endl;

	std::cout << "probeLength          : "  << to_string_with_precision(probeLength,4) << std::endl;
	std::cout << "probeRadius          : "  << to_string_with_precision(probeRadius,4) << std::endl;
	std::cout << "probeCount           : "  << int(probeCount)	       		<< std::endl;
	std::cout << "probeQuality         : "  << int(probeQuality)	   		<< std::endl;

	if (NIBR::VERBOSE()==VERBOSE_DEBUG) {
		img_FOD->printInfo();
	}

	std::cout << "\033[0m";

}
