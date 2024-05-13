#include "../../tracker/tracker.h"
#include "algorithm_ptt.h"

using namespace NIBR;

void TrackWith_PTT::estimatePosteriorMax() {

	posteriorMax = 0;

	for (int i=0; i<TRACKER::params_ptt.propMaxEstTrials; i++) {
		curve->getCandidate();
		if (curve->likelihood > posteriorMax)
			posteriorMax = curve->likelihood;
	}

	posteriorMax = std::pow(posteriorMax*DEFAULT_PTT_MAXPOSTESTCOMPENS,TRACKER::params_ptt.dataSupportExponent);
	// disp(MSG_DEBUG,"posteriorMax: %f", posteriorMax);

}

Propagation_Decision TrackWith_PTT::rejectionSample() {

	for (int tries=0; tries<TRACKER::params_ptt.triesPerRejectionSampling; tries++) {

		curve->getCandidate();

		float dataSupport = std::pow(curve->likelihood,TRACKER::params_ptt.dataSupportExponent);

		if (dataSupport > posteriorMax) {
			// disp(MSG_DEBUG,"curve->likelihood (failed): %f", curve->likelihood);
			return PROP_FAIL;
		} else if ((TRACKER::params_ptt.modMinDataSupport<=dataSupport) && (curve->doRandomThings.uniform_01()*posteriorMax <= dataSupport )) { // Equal helps to sample extrema
            // This candidate is now selected and it will be propagated
            curve->saveLastVal();
			// disp(MSG_DEBUG,"posteriorMax: %f, dataSupport: %f", posteriorMax, dataSupport);
			return PROP_CONTINUE;
		}

	}
	
	// if (tries==TRACKER::params_ptt.triesPerRejectionSampling)
	// disp(MSG_DEBUG,"Prop stop");
	return PROP_STOP;

}


// New sampling strategy picks a random k1-k2 from the cumulative distribution function
Propagation_Decision TrackWith_PTT::sampleFromCDF() {

	auto& p = TRACKER::params_ptt;

	// Calculate data support for each vertex of the k1-k2 disc 
	for (int n = 0; n < p.cdfVertCnt; n++) {
		cdfVertVal[n] = std::pow(curve->calcDataSupport(p.cdfk1k2[n].first,p.cdfk1k2[n].second),TRACKER::params_ptt.dataSupportExponent);
	}

	// Calculate the data support for each face (sum of face vertices)
	float cumSum = 0.0f;
	for (int n = 0; n < p.cdfFaceCnt; n++) {

		auto& f = p.cdfFace[n];

		if ((cdfVertVal[f[0]] >= p.modMinDataSupport) && (cdfVertVal[f[1]] >= p.modMinDataSupport) && (cdfVertVal[f[2]] >= p.modMinDataSupport) ) {

			float faceDataSupport = cdfVertVal[f[0]] + cdfVertVal[f[1]] + cdfVertVal[f[2]];

			cumSum += faceDataSupport;

			cdf[n] = cumSum;			

		} else {

			cdf[n] = FLT_MIN;

		}

	}

	if (cumSum == 0.0f) {
		return PROP_STOP;
	}

	// Sample face
	for (int tries=0; tries<TRACKER::params_ptt.triesPerRejectionSampling; tries++) {

		float randVal = curve->doRandomThings.uniform_01()*cumSum;
		
		int n = 0;
		for (n = 0; n < p.cdfFaceCnt; n++) {
			if (cdf[n] >= randVal)
				break;
		}

		auto& f = p.cdfFace[n];

		float r1 = curve->doRandomThings.uniform_01();
		float r2 = curve->doRandomThings.uniform_01();

		if (r1 + r2 > 1) {
			r1 = 1 - r1;
			r2 = 1 - r2;
		}

		float k1t = p.cdfk1k2[f[0]].first  + r1 * (p.cdfk1k2[f[1]].first -p.cdfk1k2[f[0]].first ) + r2 * (p.cdfk1k2[f[2]].first -p.cdfk1k2[f[0]].first );
		float k2t = p.cdfk1k2[f[0]].second + r1 * (p.cdfk1k2[f[1]].second-p.cdfk1k2[f[0]].second) + r2 * (p.cdfk1k2[f[2]].second-p.cdfk1k2[f[0]].second);

		float dataSupport = std::pow(curve->calcDataSupport(k1t,k2t),TRACKER::params_ptt.dataSupportExponent);

		if (TRACKER::params_ptt.modMinDataSupport<=dataSupport) {
            // This candidate is now selected and it will be propagated
            curve->saveLastVal();
			// disp(MSG_DEBUG,"posteriorMax: %f, dataSupport: %f", posteriorMax, dataSupport);
			return PROP_CONTINUE;
		}

	}

	return PROP_STOP;

	

}


Propagation_Decision TrackWith_PTT::propagate() {
    
    // Take a step forward
	// disp(MSG_DEBUG,"Before walk");
	// curve->print();
	curve->walk();
	// disp(MSG_DEBUG,"After walk");
	// curve->print();
	stepCounter++;
    
	if (TRACKER::params_ptt.useLegacySampling) {

		// Estimate posterior
		if (stepCounter%TRACKER::params_ptt.maxEstInterval==0) {
			estimatePosteriorMax();
		}

		// Rejection sample
		return rejectionSample();

	} else {

		// Sample from CDF
		return sampleFromCDF();

	}

}
