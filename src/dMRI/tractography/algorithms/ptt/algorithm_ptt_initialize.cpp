#include "../../tracker/tracker.h"
#include "algorithm_ptt.h"
#include "ptf.h"

using namespace NIBR;

Propagation_Decision TrackWith_PTT::initialize() {

	// Sample initial curve by rejection sampling
	int   tries  	  = 0;
	bool  failed 	  = false;
	float dataSupport = 0;

    
	// Initial max estimate
	posteriorMax 	= 0;

	for (tries=0; tries < initMaxEstTrials; tries++) {
		curve->getInitCandidate(trackingThread->seed_init_direction);
		if (curve->likelihood > posteriorMax) {
			posteriorMax = curve->likelihood;
			initial_curve->copy(curve);			// Saved for useBestAtInit
		}
	}

	posteriorMax = std::pow(posteriorMax*DEFAULT_PTT_MAXPOSTESTCOMPENS,dataSupportExponent);
	// disp(MSG_DEBUG,"posteriorMax: %f", posteriorMax);

	if (TRACKER::params_ptt.useBestAtInit) {

		dataSupport = std::pow(initial_curve->likelihood,dataSupportExponent);

		// Skip rejection sampling for initialization
		if (dataSupport < modMinDataSupport) {
			failed = true;
		} else {
			// curve is selected, and next it will be propagated-first value can be reset here
			curve->resetFirstVal();
		}

	} else {

		// Do rejection sampling for initialization
		for (tries=0; tries<triesPerRejectionSampling; tries++) {

			curve->getInitCandidate(trackingThread->seed_init_direction);

			dataSupport = std::pow(curve->likelihood,dataSupportExponent);

			if (dataSupport > posteriorMax) {
				failed = true;
				break;
			} else if ((modMinDataSupport<=dataSupport) && (curve->doRandomThings.uniform_01()*posteriorMax <= dataSupport )) { // Equal helps to sample extrema
                // This candidate is now selected and it will be propagated
				initial_curve->copy(curve);
				// curve is selected, and next it will be propagated-first value can be reset here
				curve->resetFirstVal();
				break;
			}

		}
		
        if (tries==triesPerRejectionSampling)
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

