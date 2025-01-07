#include "../../tracker/tracker.h"
#include "math/disc.h"
#include "algorithm_ptt.h"
#include "algorithm_ptt_params.h"
#include "ptf.h"

using namespace NIBR;

TrackWith_PTT::TrackWith_PTT() {
	curve 				= new PTF(this);
	initial_curve 		= new PTF(this);

	cdfVertVal.resize(TRACKER::params_ptt.cdfVertCnt);
	cdf.resize(TRACKER::params_ptt.cdfFaceCnt);

	// Initialize parameter to global ones
	outputStep 					= TRACKER::params_ptt.outputStep_global;
	stepSize   					= TRACKER::params_ptt.stepSize_global;

	minRadiusOfCurvature 		= TRACKER::params_ptt.minRadiusOfCurvature_global;
	minDataSupport 		 		= TRACKER::params_ptt.minDataSupport_global;
	dataSupportExponent  		= TRACKER::params_ptt.dataSupportExponent_global;

	maxEstInterval 				= TRACKER::params_ptt.maxEstInterval_global;
	initMaxEstTrials 			= TRACKER::params_ptt.initMaxEstTrials_global;
	propMaxEstTrials 			= TRACKER::params_ptt.propMaxEstTrials_global;
	triesPerRejectionSampling 	= TRACKER::params_ptt.triesPerRejectionSampling_global;

	probeCount   				= TRACKER::params_ptt.probeCount_global;
	probeRadius  				= TRACKER::params_ptt.probeRadius_global;
	probeLength  				= TRACKER::params_ptt.probeLength_global;
	probeQuality 				= TRACKER::params_ptt.probeQuality_global;

	maxCurvature 	  			= TRACKER::params_ptt.maxCurvature_global;
	modMinDataSupport 			= TRACKER::params_ptt.modMinDataSupport_global;

	cdfFace 					= &TRACKER::params_ptt.cdfFace_global;

	setCDF();

	parametersAreReady 			= false;
	
}

TrackWith_PTT::~TrackWith_PTT() {
	if (initial_curve!=NULL) 	delete 	 initial_curve;
	if (curve!=NULL) 			delete   curve;
}

void TrackWith_PTT::setSeed() {
	parametersAreReady = false;
	fetchParams(thread->seed_coordinates);
	curve->setPosition(thread->seed_coordinates);
	curve->refreshParams();
}

void TrackWith_PTT::setCDF() {
	auto maxCurvIdx = findFirstGreaterIndex(TRACKER::params_ptt.cdfCurvatures,maxCurvature);
	cdfk1k2 = &TRACKER::params_ptt.cdfk1k2_global[maxCurvIdx];
}

void TrackWith_PTT::reset() {
	if (initial_curve!=NULL) 	delete 	 initial_curve;
	if (curve!=NULL) 			delete   curve;
	curve 			= new PTF(this);
	initial_curve 	= new PTF(this);

	// Reset parameters
	outputStep 					= TRACKER::params_ptt.outputStep_global;
	stepSize   					= TRACKER::params_ptt.stepSize_global;

	minRadiusOfCurvature 		= TRACKER::params_ptt.minRadiusOfCurvature_global;
	minDataSupport 		 		= TRACKER::params_ptt.minDataSupport_global;
	dataSupportExponent  		= TRACKER::params_ptt.dataSupportExponent_global;

	maxEstInterval 				= TRACKER::params_ptt.maxEstInterval_global;
	initMaxEstTrials 			= TRACKER::params_ptt.initMaxEstTrials_global;
	propMaxEstTrials 			= TRACKER::params_ptt.propMaxEstTrials_global;
	triesPerRejectionSampling 	= TRACKER::params_ptt.triesPerRejectionSampling_global;

	probeCount   				= TRACKER::params_ptt.probeCount_global;
	probeRadius  				= TRACKER::params_ptt.probeRadius_global;
	probeLength  				= TRACKER::params_ptt.probeLength_global;
	probeQuality 				= TRACKER::params_ptt.probeQuality_global;

	maxCurvature 	  			= TRACKER::params_ptt.maxCurvature_global;
	modMinDataSupport 			= TRACKER::params_ptt.modMinDataSupport_global;

	setCDF();
	
	parametersAreReady 			= false;
}

void TrackWith_PTT::append() {

    thread->walker->streamline->push_back({curve->p[0],curve->p[1],curve->p[2]});
	// disp(MSG_DEBUG,"thread->walker->streamline[%d]=[%.2f,%.2f,%.2f]",thread->walker->streamline->size()-1,curve->p[0],curve->p[1],curve->p[2]);

	/*
	if (TRACKER::params_ptt.saveFrame) {
		tangent.push_back({curve->F[0][0],curve->F[0][1],curve->F[0][2]});
		k1axis.push_back ({curve->F[1][0],curve->F[1][1],curve->F[1][2]});
		k2axis.push_back ({curve->F[2][0],curve->F[2][1],curve->F[2][2]});
		k1.push_back(curve->k1);
		k2.push_back(curve->k2);
		likelihood.push_back(curve->likelihood);
	}
	*/

}

float TrackWith_PTT::writeStepSize() {
	return outputStep;
}

int TrackWith_PTT::appendInterval() {
	return std::max( int(outputStep / stepSize),1);
}

Propagation_Decision TrackWith_PTT::flip() {

	curve->getFlippedCopy(initial_curve);
	parametersAreReady = false;
	fetchParams(curve->p);
	curve->refreshParams();

	std::reverse(thread->walker->streamline->begin(),thread->walker->streamline->end());
	
	/*
	if (TRACKER::params_ptt.saveFrame) {
		std::reverse(tangent.begin(),tangent.end());
		std::reverse(k1axis.begin(),k1axis.end());
		std::reverse(k2axis.begin(),k2axis.end());
		std::reverse(k1.begin(),k1.end());
		std::reverse(k2.begin(),k2.end());
		std::reverse(likelihood.begin(),likelihood.end());
	}
	*/

	// Already find the next curve! 
	// Because next call for propagate needs to know the curve in order to walk.
	
	// Estimate posterior
    posteriorMax = curve->getInitPosteriorMax();

	// disp(MSG_DEBUG,"Flipped curve:");
	// curve->print();

	// Rejection sample
	return rejectionSample();

}

void TrackWith_PTT::fetchParams(float *p) {

	// If no parameter mask is provided, we will use the global parameters, which were initially set
	if (TRACKER::params_ptt.img_param_mask == NULL) 
		return;

	// Refresh parameters if needed
	if ((*TRACKER::params_ptt.img_param_mask)(p)>0) {

		for (auto s : TRACKER::params_ptt.toRefresh) {

			if (s == "outputStep") {
				outputStep = (*TRACKER::params_ptt.outputStep_img)(p);
			} 
			else if (s == "stepSize") {
				stepSize   = (*TRACKER::params_ptt.stepSize_img)(p);
			}
			else if (s == "minRadiusOfCurvature") {
				minRadiusOfCurvature = (*TRACKER::params_ptt.minRadiusOfCurvature_img)(p);
				maxCurvature = 1.0f / minRadiusOfCurvature;
				if (maxCurvature < 1e-4) maxCurvature = 1e-4;
				setCDF();
			}
			else if (s == "minDataSupport") {
				minDataSupport = (*TRACKER::params_ptt.minDataSupport_img)(p);
			}
			else if (s == "dataSupportExponent") {
				dataSupportExponent = (*TRACKER::params_ptt.dataSupportExponent_img)(p);
			}
			else if (s == "maxEstInterval") {
				maxEstInterval = (*TRACKER::params_ptt.maxEstInterval_img)(p);
			}
			else if (s == "initMaxEstTrials") {
				initMaxEstTrials = (*TRACKER::params_ptt.initMaxEstTrials_img)(p);
			}
			else if (s == "propMaxEstTrials") {
				propMaxEstTrials = (*TRACKER::params_ptt.propMaxEstTrials_img)(p);
			}
			else if (s == "triesPerRejectionSampling") {
				triesPerRejectionSampling = (*TRACKER::params_ptt.triesPerRejectionSampling_img)(p);
			}
			else if (s == "probeLength") {
				probeLength  = (*TRACKER::params_ptt.probeLength_img)(p);
			}
			else if (s == "probeQuality") {
				probeQuality = (*TRACKER::params_ptt.probeQuality_img)(p);
			}
			else if (s == "probeCount") {
				probeCount 	 = (*TRACKER::params_ptt.probeCount_img)(p);
			}
			else if (s == "probeRadius") {
				probeRadius  = (*TRACKER::params_ptt.probeRadius_img)(p);
			}

		}

		if ((TRACKER::params_ptt.toRefresh.find("minDataSupport") != TRACKER::params_ptt.toRefresh.end()) || (TRACKER::params_ptt.toRefresh.find("dataSupportExponent") != TRACKER::params_ptt.toRefresh.end())) {
			modMinDataSupport = std::pow(minDataSupport,dataSupportExponent);
		}

		if ((TRACKER::params_ptt.toRefresh.find("probeCount") != TRACKER::params_ptt.toRefresh.end()) || (TRACKER::params_ptt.toRefresh.find("probeRadius") != TRACKER::params_ptt.toRefresh.end())) {

			if (probeCount==1.0f) {
				probeRadius=0.0f;
			} else {
				if (probeRadius==0.0f) {
					probeCount=1.0f;
				} else {
					if (probeRadius > minRadiusOfCurvature) {
						probeRadius = minRadiusOfCurvature;
					}
				}
			}
		}

		parametersAreReady = false;

	} else {

		if (parametersAreReady) return;

		for (auto s : TRACKER::params_ptt.toRefresh) {

			if (s == "outputStep") {
				outputStep = TRACKER::params_ptt.outputStep_global;
			} 
			else if (s == "stepSize") {
				stepSize   = TRACKER::params_ptt.stepSize_global;
			}
			else if (s == "minRadiusOfCurvature") {
				minRadiusOfCurvature = TRACKER::params_ptt.minRadiusOfCurvature_global;
				maxCurvature 		 = TRACKER::params_ptt.maxCurvature_global;
				setCDF();
			}
			else if (s == "minDataSupport") {
				minDataSupport = TRACKER::params_ptt.minDataSupport_global;
			}
			else if (s == "dataSupportExponent") {
				dataSupportExponent = TRACKER::params_ptt.dataSupportExponent_global;
			}
			else if (s == "maxEstInterval") {
				maxEstInterval = TRACKER::params_ptt.maxEstInterval_global;
			}
			else if (s == "initMaxEstTrials") {
				initMaxEstTrials = TRACKER::params_ptt.initMaxEstTrials_global;
			}
			else if (s == "propMaxEstTrials") {
				propMaxEstTrials = TRACKER::params_ptt.propMaxEstTrials_global;
			}
			else if (s == "triesPerRejectionSampling") {
				triesPerRejectionSampling = TRACKER::params_ptt.triesPerRejectionSampling_global;
			}
			else if (s == "probeLength") {
				probeLength = TRACKER::params_ptt.probeLength_global;
			}
			else if (s == "probeQuality") {
				probeQuality = TRACKER::params_ptt.probeQuality_global;
			}
			else if (s == "probeCount") {
				probeCount = TRACKER::params_ptt.probeCount_global;
			}
			else if (s == "probeRadius") {
				probeRadius = TRACKER::params_ptt.probeRadius_global;
			}

		}

		modMinDataSupport = TRACKER::params_ptt.modMinDataSupport_global; // This is faster to update like this than to check whether minDataSupport or dataSupportExponent needed refreshing

		parametersAreReady = true;

	}

}
