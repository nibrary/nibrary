#include "../../tracker/tracker.h"
#include "algorithm_ptt_params.h"
#include "math/disc.h"
#include "base/vectorOperations.h"


using namespace NIBR;

bool Params_PTT::update() {

	if (pttIsReady)
		return true;

	if ((img_FOD==NULL) || ((img_FOD!=NULL) &&  ( (img_FOD->filePath != img_FOD_path) || (img_FOD->getSphereFileName() != fod_sphere_path) ))) {

		if (img_FOD!=NULL) {
			delete img_FOD;
			img_FOD = NULL;
		}

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

	if ( (img_param_mask==NULL) || ((img_param_mask!=NULL) && (img_param_mask->filePath != img_param_mask_path)) ) {
		
		delete_img(img_param_mask);

		if (img_param_mask_path != "") {
			
			disp(MSG_INFO, "Reading parameter mask image");
			img_param_mask = new Image<bool>(img_param_mask_path);

			if(!img_param_mask->read()) {
				disp(MSG_ERROR, "Can't read parameter mask image");
				return false;		
			}

			if (img_param_mask->numberOfDimensions != 3) {
				disp(MSG_ERROR, "Parameter mask has to be a 3D image");
				return false;		
			}

			img_param_mask->setInterpolationMethod(NEAREST);
		}

	}

	if(!update_param_img(outputStep_img,				outputStep_img_path,				"outputStep", 					toRefresh, img_param_mask)) return false;
	if(!update_param_img(stepSize_img,					stepSize_img_path,					"stepSize", 					toRefresh, img_param_mask)) return false;

	if(!update_param_img(minRadiusOfCurvature_img,		minRadiusOfCurvature_img_path,		"minRadiusOfCurvature", 		toRefresh, img_param_mask)) return false;	
	if(!update_param_img(minDataSupport_img,			minDataSupport_img_path,			"minDataSupport", 				toRefresh, img_param_mask)) return false;
	if(!update_param_img(dataSupportExponent_img,		dataSupportExponent_img_path,		"dataSupportExponent", 			toRefresh, img_param_mask)) return false;

	if(!update_param_img(maxEstInterval_img,			maxEstInterval_img_path,			"maxEstInterval", 				toRefresh, img_param_mask)) return false;
	if(!update_param_img(initMaxEstTrials_img,			initMaxEstTrials_img_path,			"initMaxEstTrials", 			toRefresh, img_param_mask)) return false;
	if(!update_param_img(propMaxEstTrials_img,			propMaxEstTrials_img_path,			"propMaxEstTrials", 			toRefresh, img_param_mask)) return false;
	if(!update_param_img(triesPerRejectionSampling_img,	triesPerRejectionSampling_img_path,	"triesPerRejectionSampling", 	toRefresh, img_param_mask)) return false;

	if(!update_param_img(probeLength_img,				probeLength_img_path,				"probeLength", 					toRefresh, img_param_mask)) return false;
	if(!update_param_img(probeQuality_img,				probeQuality_img_path,				"probeQuality", 				toRefresh, img_param_mask)) return false;
	if(!update_param_img(probeCount_img,				probeCount_img_path,				"probeCount", 					toRefresh, img_param_mask)) return false;
	if(!update_param_img(probeRadius_img,				probeRadius_img_path,				"probeRadius", 					toRefresh, img_param_mask)) return false;

	setDefaults();

	pttIsReady = true;

	return true;

}

void Params_PTT::clear() {

	pttIsReady = false;

	// FOD internal
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
	weakLinkThresh              = NAN;

	outputStep_global      		= NAN;
	stepSize_global             = NAN;

	minRadiusOfCurvature_global = NAN;
	minDataSupport_global       = NAN;
	dataSupportExponent_global  = NAN;
	

	// Sampling options
	useBestAtInit               = false;
	useLegacySampling           = false;
	samplingQuality             = -1;

	maxEstInterval_global              = -1;
	initMaxEstTrials_global            = -1;
	propMaxEstTrials_global            = -1;
	triesPerRejectionSampling_global   = -1;

	// Probe options
	probeCount_global      		= NAN;
    probeRadius_global       	= NAN;
    probeLength_global       	= NAN;
	probeQuality_global       	= NAN;	

	// Derived parameters
	checkWeakLinks              = NAN;

    maxCurvature_global         = NAN;
	modMinDataSupport_global    = NAN;

	if (img_FOD!=NULL) {
		disp(MSG_DETAIL,"Cleaning memory for %s", img_FOD->filePath.c_str());
		delete img_FOD;
		img_FOD = NULL;
	}

	delete_img(img_param_mask);

	delete_img(outputStep_img);
	delete_img(stepSize_img);

	delete_img(minRadiusOfCurvature_img);
	delete_img(minDataSupport_img);
	delete_img(dataSupportExponent_img);

	delete_img(maxEstInterval_img);
	delete_img(initMaxEstTrials_img);
	delete_img(propMaxEstTrials_img);
	delete_img(triesPerRejectionSampling_img);

	delete_img(probeCount_img);
	delete_img(probeRadius_img);
	delete_img(probeLength_img);
	delete_img(probeQuality_img);

	delete_img(img_peaks);
	delete_img(img_bias);
	

	disp(MSG_DEBUG,"PTT parameters cleared.");

}



void Params_PTT::setDefaults() {

	// Handle outputStep, stepSize, and minRadiusOfCurvature
	if (std::isnan(outputStep_global)     		|| (outputStep_global<=0)			) 	outputStep_global 	 		= smallestPixDim*DEFAULT_PTT_OUTPUTSTEPSIZE_IN_PIXELDIM;
	if (std::isnan(stepSize_global) 	  		|| (stepSize_global<=0)   			) 	stepSize_global 	 		= smallestPixDim*DEFAULT_PTT_STEPSIZE_IN_PIXELDIM;
	
	if (std::isnan(minRadiusOfCurvature_global) || (minRadiusOfCurvature_global<=0)	)	minRadiusOfCurvature_global = smallestPixDim*DEFAULT_PTT_MINRADIUSOFCURVATURE_IN_PIXELDIM;	
	maxCurvature_global = 1.0f / minRadiusOfCurvature_global;
	if (maxCurvature_global < 1e-4) maxCurvature_global = 1e-4;

	// Handle minDataSupport, dataSupportExponent, and weakLinkThresh
	if (std::isnan(minDataSupport_global) 		 || (minDataSupport_global<0.0)		) 	minDataSupport_global 		= DEFAULT_PTT_MINDATASUPPORT;
	if (std::isnan(dataSupportExponent_global)   || (dataSupportExponent_global<0.0))   dataSupportExponent_global  = DEFAULT_PTT_DATASUPPORTEXPONENT;
	if (std::isnan(weakLinkThresh) 		 || (weakLinkThresh<=0)			)	weakLinkThresh 		 = DEFAULT_PTT_WEAKLINKTHRESH;

	modMinDataSupport_global = std::pow(minDataSupport_global,dataSupportExponent_global);	
	checkWeakLinks 			 = (weakLinkThresh>0) ? true : false;

	// Handle maxEstInterval, initMaxEstTrials, propMaxEstTrials, triesPerRejectionSampling, and samplingQuality
	if (maxEstInterval_global			 <= 0.0)	maxEstInterval_global 			  = DEFAULT_PTT_MAXESTINTERVAL;
	if (initMaxEstTrials_global		  	 <= 0.0)	initMaxEstTrials_global 		  = DEFAULT_PTT_INITMAXESTTRIALS;
	if (propMaxEstTrials_global		 	 <= 0.0)	propMaxEstTrials_global 		  = DEFAULT_PTT_PROPMAXESTTRIALS;
	if (triesPerRejectionSampling_global <= 0.0)	triesPerRejectionSampling_global  = DEFAULT_PTT_TRIESPERREJECTIONSAMPLING;

	// Handle probeCount, probeRadius, probeLength, and probeQuality
	if (std::isnan(probeLength_global)  || (probeLength_global<=0.0))								probeLength_global  = smallestPixDim*DEFAULT_PTT_PROBELENGTH_IN_PIXELDIM;
	if (std::isnan(probeQuality_global) || (probeQuality_global<=1.0) || (probeQuality_global>100))	probeQuality_global = DEFAULT_PTT_PROBEQUALITY;
	
	if (probeCount_global<1)   probeCount_global = NAN;
	if (probeRadius_global<0)  probeRadius_global= NAN;

	if (isnan(probeCount_global) && isnan(probeRadius_global)) {
		probeCount_global  = DEFAULT_PTT_PROBECOUNT;
		probeRadius_global = DEFAULT_PTT_PROBERADIUS_IN_PIXELDIM*smallestPixDim;
	} else if (isnan(probeCount_global) && probeRadius_global!=NAN) {
		probeCount_global = DEFAULT_PTT_PROBECOUNT_WHEN_THEREIS_PROBERADIUS;
	} else if (!isnan(probeCount_global) && isnan(probeRadius_global)) {
		probeRadius_global= DEFAULT_PTT_PROBERADIUS_IN_PIXELDIM*smallestPixDim;
	}

	if (probeCount_global==1.0f) {
		probeRadius_global=0.0f;
	} else {
		if (probeRadius_global==0.0f) {
			probeCount_global=1.0f;
		} else {
			if (probeRadius_global>minRadiusOfCurvature_global) {
				probeRadius_global = minRadiusOfCurvature_global;
			}
		}
	}

	if ((probeRadius_img == NULL) && (probeRadius_global == 0.0f)) toRefresh.erase("probeCount"); // Don't refresh probeCount
	if ((probeCount_img  == NULL) && (probeCount_global  == 1.0f)) toRefresh.erase("probeRadius"); // Don't refresh probeRadius

	// Handle samplingQuality
	if (samplingQuality<=0.0) samplingQuality = DEFAULT_PTT_SAMPLINGQUALITY;

	// Prep CDF domain

	cdfCurvatures.clear();

	float maxMaxCurvature = maxCurvature_global;

	if (minRadiusOfCurvature_img != NULL) {

		const auto& img = minRadiusOfCurvature_img;

		float minMinRadiusOfCurvature = minRadiusOfCurvature_global;

		for (int n = 0; n < img->numel; n++) {
			if ((img->data[n] > 0) && (img->data[n] < minMinRadiusOfCurvature))
				minMinRadiusOfCurvature = img->data[n];
		}

		if (minMinRadiusOfCurvature < 1e-4) minMinRadiusOfCurvature = 0.001;

		maxMaxCurvature = 1.0f / minMinRadiusOfCurvature;
		if (maxMaxCurvature < 1e-4) maxMaxCurvature = 1e-4;

		cdfCurvatures = linspace(0.0001f, maxMaxCurvature, 1000);

	} else {
		cdfCurvatures.push_back(maxMaxCurvature);
	}


	auto fill_cdfk1k2 = [&](auto& verts) {

		cdfk1k2_global.clear();
		cdfk1k2_global.reserve(cdfCurvatures.size());
		
		for (auto curv : cdfCurvatures) {

			std::vector<std::pair<float,float>> tmp;

			for (int n = 0; n < cdfVertCnt; n++) {
				float k1 = verts[n][0] * curv;
				float k2 = verts[n][1] * curv;
				tmp.push_back(std::make_pair(k1,k2));
			}

			cdfk1k2_global.push_back(tmp);

		}

	};

	auto fill_cdfFace = [&](auto& faces) {
		cdfFace_global.clear();
		cdfFace_global.reserve(cdfFaceCnt);
		for (int n = 0; n < cdfFaceCnt; n++) {
			cdfFace_global.push_back({faces[n][0],faces[n][1],faces[n][2]});
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

	auto dispParam = [&](std::string pstr, Image<float>* img, float pval, int precision) {

		std::cout << pstr << to_string_with_precision(pval,precision);

		if (img != NULL)
			std::cout << "," << img->filePath << std::endl;
		else
			std::cout << std::endl;

	};

	auto dispParamInt = [&](std::string pstr, Image<int>* img, float pval) {

		std::cout << pstr << pval;

		if (img != NULL)
			std::cout << "," << img->filePath << std::endl;
		else
			std::cout << std::endl;

	};

	dispParam("stepSize             : ", stepSize_img, 	 			stepSize_global,   			 4);
	dispParam("writeStepSize        : ", outputStep_img, 			outputStep_global, 			 4);
	
	dispParam("minRadiusOfCurvature : ", minRadiusOfCurvature_img, 	minRadiusOfCurvature_global, 4);
	dispParam("minDataSupport       : ", minDataSupport_img, 		minDataSupport_global, 		 4);
	dispParam("dataSupportExponent  : ", dataSupportExponent_img, 	dataSupportExponent_global,  4);

    std::cout << "ignoreWeakLinks      : "  << to_string_with_precision(weakLinkThresh,4) << std::endl;

	dispParamInt("maxEstInterval       : ", maxEstInterval_img, 	 		maxEstInterval_global);
	dispParamInt("maxSamplingPerStep   : ", triesPerRejectionSampling_img, 	triesPerRejectionSampling_global);
	dispParamInt("initMaxEstTrials     : ", initMaxEstTrials_img, 	 		initMaxEstTrials_global);
	dispParamInt("propMaxEstTrials     : ", propMaxEstTrials_img, 	 		propMaxEstTrials_global);

	if (useBestAtInit)
		std::cout << "useBestAtInit        : ON "  << std::endl;
	else
		std::cout << "useBestAtInit        : OFF " << std::endl;

	if (useLegacySampling)
		std::cout << "useLegacySampling    : ON "  << std::endl;
	else
		std::cout << "useLegacySampling    : OFF " << std::endl;

	std::cout << "samplingQuality      : "  << samplingQuality 			    << std::endl;

	dispParam("probeLength          : ", probeLength_img, 		probeLength_global, 	4);
	dispParam("probeRadius          : ", probeRadius_img, 		probeRadius_global, 	4);
	dispParam("probeCount           : ", probeCount_img, 		probeCount_global, 		0);
	dispParam("probeQuality         : ", probeQuality_img, 		probeQuality_global, 	0);

	if (NIBR::VERBOSE()==VERBOSE_DEBUG) {
		img_FOD->printInfo();
	}

	std::cout << "\033[0m";

}
