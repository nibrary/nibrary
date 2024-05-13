#pragma once

#include "../tractographyAlgorithm.h"
#include "ptf.h"

namespace NIBR {

class TrackWith_PTT : public TractographyAlgorithm {

public:

	TrackWith_PTT();
	~TrackWith_PTT();

	virtual Propagation_Decision initialize();
	virtual Propagation_Decision propagate();
	virtual Propagation_Decision flip();

	virtual void  setSeed();
	virtual void  reset();
	virtual void  append();
	virtual float writeStepSize();
	virtual int   appendInterval();

private:

	PTF		*curve 		   = NULL;
	PTF 	*initial_curve = NULL;

	void 	 				estimatePosteriorMax();
	Propagation_Decision 	rejectionSample();
	Propagation_Decision 	sampleFromCDF();	

	float 	 posteriorMax;
	int      stepCounter;

	std::vector<float> cdfVertVal; // Data support for each vertex of the k1-k2 disc
	std::vector<float> cdf;        // The cumulative distribution function (evaluated on faces)

	// std::vector<Point> tangent;
	// std::vector<Point> k1axis;
	// std::vector<Point> k2axis;

	// std::vector<float> k1;
	// std::vector<float> k2;
	// std::vector<float> likelihood;
    
};

}