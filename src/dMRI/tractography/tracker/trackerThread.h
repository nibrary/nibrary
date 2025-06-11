#pragma once

#include "base/nibr.h"
#include "../pathway/pathway.h"
#include "../io/tractogramWriter.h"

namespace NIBR {

class TractographyAlgorithm;

class TrackingThread {
public:

	TrackingThread(int _seed_no);
	~TrackingThread();

	int                     id;
	int 					threadId;

	TractographyAlgorithm  *method;
	Walker             	   *walker;
	Streamline 				streamline;
	
	float*             		seed_coordinates;
	float*             		seed_init_direction;

	void                    init();
	void                    reset();
	void                    clear();
	bool 		            track(TractogramWriter* writer = NULL); // returns true if tracking was successful. It no writer is provided, then saves in internal tractogram.

};

}
