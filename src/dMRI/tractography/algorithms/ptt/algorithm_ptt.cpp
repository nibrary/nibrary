#include "../../tracker/tracker.h"
#include "algorithm_ptt.h"
#include "algorithm_ptt_params.h"
#include "ptf.h"

using namespace NIBR;

TrackWith_PTT::TrackWith_PTT() {	
	curve 			= new PTF();
	initial_curve 	= new PTF();

	cdfVertVal.resize(TRACKER::params_ptt.cdfVertCnt);
	cdf.resize(TRACKER::params_ptt.cdfFaceCnt);
	
}

TrackWith_PTT::~TrackWith_PTT() {
	if (initial_curve!=NULL) 	delete 	 initial_curve;
	if (curve!=NULL) 			delete   curve;
}

void TrackWith_PTT::setSeed() {
	curve->setPosition(thread->seed_coordinates);
}

void TrackWith_PTT::reset() {
	if (initial_curve!=NULL) 	delete 	 initial_curve;
	if (curve!=NULL) 			delete   curve;
	curve 			= new PTF();
	initial_curve 	= new PTF();
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
	return TRACKER::params_ptt.outputStep;
}

int TrackWith_PTT::appendInterval() {
	return std::max( int(TRACKER::params_ptt.outputStep / TRACKER::params_ptt.stepSize),1);
}

Propagation_Decision TrackWith_PTT::flip() {

	curve->getFlippedCopy(initial_curve);

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
