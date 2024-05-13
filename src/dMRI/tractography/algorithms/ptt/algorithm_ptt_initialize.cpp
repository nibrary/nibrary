#include "../../tracker/tracker.h"
#include "algorithm_ptt.h"

using namespace NIBR;

Propagation_Decision TrackWith_PTT::initialize() {

	// Sample initial curve by rejection sampling
	int   tries  	  = 0;
	bool  failed 	  = false;
	float dataSupport = 0;

    
	// Initial max estimate
	posteriorMax 	= 0;

	for (tries=0; tries < TRACKER::params_ptt.initMaxEstTrials; tries++) {
		curve->getInitCandidate(thread->seed_init_direction);
		if (curve->likelihood > posteriorMax) {
			posteriorMax = curve->likelihood;
			initial_curve->copy(curve);			// Saved for useBestAtInit
		}
	}

	posteriorMax = std::pow(posteriorMax*DEFAULT_PTT_MAXPOSTESTCOMPENS,TRACKER::params_ptt.dataSupportExponent);
	// disp(MSG_DEBUG,"posteriorMax: %f", posteriorMax);

	if (TRACKER::params_ptt.useBestAtInit) {

		dataSupport = std::pow(initial_curve->likelihood,TRACKER::params_ptt.dataSupportExponent);

		// Skip rejection sampling for initialization
		if (dataSupport < TRACKER::params_ptt.modMinDataSupport) {
			failed = true;
		} else {
			curve->copy(initial_curve);
			curve->saveLastVal();
		}

	} else {

		// Do rejection sampling for initialization
		for (tries=0; tries<TRACKER::params_ptt.triesPerRejectionSampling; tries++) {

			curve->getInitCandidate(thread->seed_init_direction);

			dataSupport = std::pow(curve->likelihood,TRACKER::params_ptt.dataSupportExponent);

			if (dataSupport > posteriorMax) {
				failed = true;
				break;
			} else if ((TRACKER::params_ptt.modMinDataSupport<=dataSupport) && (curve->doRandomThings.uniform_01()*posteriorMax <= dataSupport )) { // Equal helps to sample extrema
                // This candidate is now selected and it will be propagated
				initial_curve->copy(curve);
				curve->saveLastVal();
				break;
			}

		}
		
        if (tries==TRACKER::params_ptt.triesPerRejectionSampling)
            return PROP_STOP;

	}
	
	stepCounter = 0;
	initial_curve->setInitPosteriorMax(posteriorMax);

    if (failed) {
		// disp(MSG_DEBUG,"Initialization failed");
        return PROP_FAIL;
	} else {
		// disp(MSG_DEBUG,"Initialization successful. posteriorMax: %f, dataSupport: %f", posteriorMax, dataSupport);
		// disp(MSG_DEBUG,"Initial curve:");
		// initial_curve->print();
    	return PROP_CONTINUE;
	}

}

