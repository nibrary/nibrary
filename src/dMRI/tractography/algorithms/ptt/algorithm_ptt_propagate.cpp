#include "../../tracker/tracker.h"
#include "algorithm_ptt.h"
#include "ptf.h"

using namespace NIBR;

void TrackWith_PTT::estimatePosteriorMax() {

	posteriorMax = 0;

	for (int i=0; i<propMaxEstTrials; i++) {
		curve->getCandidate();
		if (curve->likelihood > posteriorMax)
			posteriorMax = curve->likelihood;
	}

	posteriorMax = std::pow(posteriorMax*DEFAULT_PTT_MAXPOSTESTCOMPENS,dataSupportExponent);
	// disp(MSG_DEBUG,"posteriorMax: %f", posteriorMax);

}

Propagation_Decision TrackWith_PTT::rejectionSample() {

	for (int tries=0; tries<triesPerRejectionSampling; tries++) {

		curve->getCandidate();

		float dataSupport = std::pow(curve->likelihood,dataSupportExponent);

		if (dataSupport > posteriorMax) {
			// disp(MSG_DEBUG,"curve->likelihood (failed): %f", curve->likelihood);
			return PROP_FAIL;
		} else if ((modMinDataSupport<=dataSupport) && (curve->doRandomThings.uniform_01()*posteriorMax <= dataSupport )) { // Equal helps to sample extrema
            // This candidate is now selected and it will be propagated
            curve->resetFirstVal();
			// disp(MSG_DEBUG,"posteriorMax: %f, dataSupport: %f", posteriorMax, dataSupport);
			return PROP_CONTINUE;
		}

	}
	
	// if (tries==triesPerRejectionSampling)
	// disp(MSG_DEBUG,"Prop stop");
	return PROP_STOP;

}


// New sampling strategy picks a random k1-k2 from the cumulative distribution function
Propagation_Decision TrackWith_PTT::sampleFromCDF() {

	auto& p = TRACKER::params_ptt;

	// Calculate data support for each vertex of the k1-k2 disc 
	for (int n = 0; n < p.cdfVertCnt; n++) {
		cdfVertVal[n] = std::pow(curve->calcDataSupport(cdfk1k2->at(n).first,cdfk1k2->at(n).second),dataSupportExponent);
	}

	// Calculate the data support for each face (sum of face vertices)
	float cumSum = 0.0f;
	for (int n = 0; n < p.cdfFaceCnt; n++) {

		auto& f = cdfFace->at(n);

		if ((cdfVertVal[f[0]] >= modMinDataSupport) && (cdfVertVal[f[1]] >= modMinDataSupport) && (cdfVertVal[f[2]] >= modMinDataSupport) ) {

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
	for (int tries=0; tries<triesPerRejectionSampling; tries++) {

		float randVal = curve->doRandomThings.uniform_01()*cumSum;
		
		int n = 0;
		for (n = 0; n < p.cdfFaceCnt; n++) {
			if (cdf[n] >= randVal)
				break;
		}

		auto& f = cdfFace->at(n);

		float r1 = curve->doRandomThings.uniform_01();
		float r2 = curve->doRandomThings.uniform_01();

		if (r1 + r2 > 1) {
			r1 = 1 - r1;
			r2 = 1 - r2;
		}

		float k1t = cdfk1k2->at(f[0]).first  + r1 * (cdfk1k2->at(f[1]).first -cdfk1k2->at(f[0]).first ) + r2 * (cdfk1k2->at(f[2]).first -cdfk1k2->at(f[0]).first );
		float k2t = cdfk1k2->at(f[0]).second + r1 * (cdfk1k2->at(f[1]).second-cdfk1k2->at(f[0]).second) + r2 * (cdfk1k2->at(f[2]).second-cdfk1k2->at(f[0]).second);

		float dataSupport = std::pow(curve->calcDataSupport(k1t,k2t),dataSupportExponent);

		if (modMinDataSupport<=dataSupport) {
            // This candidate is now selected and it will be propagated
            curve->resetFirstVal();
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
	fetchParams(curve->p);
	curve->refreshParams();

	// disp(MSG_DEBUG,"After walk");
	// curve->print();
	stepCounter++;
    
	if (TRACKER::params_ptt.useLegacySampling) {

		// Estimate posterior
		if (stepCounter%maxEstInterval==0) {
			estimatePosteriorMax();
		}

		// Rejection sample
		return rejectionSample();

	} else {

		// Sample from CDF
		return sampleFromCDF();

	}

}
