#include "../algorithms/ptt/algorithm_ptt.h"
#include "../utility/resampleStreamline.h"

#include "trackerThread.h"
#include "tracker.h"

#include <cstdint>

using namespace NIBR;

TrackingThread::TrackingThread(int _id) 
{
	id 		 = _id;
	threadId = _id % TRACKER::nThreads;
	init();
}

TrackingThread::~TrackingThread() 
{
	clear();
}

void TrackingThread::reset() 
{
	clear();
	init();
}

void TrackingThread::init()
{
	switch (TRACKER::algorithm) {
	case PTT:
		method	= new TrackWith_PTT();
		break;
	case LOCAL_PROBABILISTIC:
		// method	= new TrackWith_Local_Probabilistic(this);
		break;
	default :
		break;
	}

	method->setThread(this);
	walker = NULL;

	seed_coordinates = new float[3];

	if ( (TRACKER::seed.getMode()==SEED_LIST_WITH_DIRECTIONS) 			||
		 (TRACKER::seed.getMode()==SEED_SURFACE_MASK_WITH_DIRECTIONS) 	||
		 (TRACKER::seed.getMode()==SEED_SURFACE_RS_WITH_DIRECTIONS) ) {
		seed_init_direction = new float[3];
	} else {
		seed_init_direction = NULL;
	}
}

void TrackingThread::clear() 
{
	if (method				!= NULL)	{delete   method; 			   method = NULL;}
	if (walker				!= NULL)	{delete   walker;			   walker = NULL;}
	if (seed_coordinates    != NULL)	{delete[] seed_coordinates;	   seed_coordinates 	= NULL;}
	if (seed_init_direction != NULL)	{delete[] seed_init_direction; seed_init_direction 	= NULL;}
	streamline.clear();
}


bool hitTimeLimit(Walker* walker) 
{
	if (!TRACKER::isWithinTimeLimits()) {
		walker->action 						= DISCARD;
		walker->discardingReason 			= REACHED_TIME_LIMIT;
		disp(MSG_DEBUG,"Reached run time limit");
		return true;
	}

	return false;
}


bool TrackingThread::track()
{

	// track starts with reseting the object, i.e., everything is cleared.
	disp(MSG_DETAIL, "Starting tracker: %d", id);
	reset();

	int trialNo = 0;

	walker = TRACKER::pw.createWalker(&streamline,trialNo); // walker = NULL after reset(), so this is OK.

	disp(MSG_DETAIL, "Getting seed");

	switch(TRACKER::seed.getSeed(seed_coordinates,seed_init_direction,threadId))	{
		case SEED_OK:
			disp(MSG_DEBUG, "SEED_OK");
			disp(MSG_DEBUG, "Tracker no: %d", id);
			disp(MSG_DETAIL, "Seed coordinate : [%.2f,%.2f,%.2f]", seed_coordinates[0], 	seed_coordinates[1], 	seed_coordinates[2]   );
			if (seed_init_direction!=NULL) {
				disp(MSG_DEBUG, "seed direction   : [%.2f,%.2f,%.2f]", seed_init_direction[0], 	seed_init_direction[1], seed_init_direction[2]);
			}
			break;

		case SEED_REACHED_MAX_COUNT:
			disp(MSG_DEBUG, "SEED_REACHED_MAX_COUNT");
			return false;

		case SEED_ERROR:
			disp(MSG_DEBUG, "SEED_ERROR");
			return false;
	}

	while ((walker->action != KEEP) && (trialNo<TRACKER::seed.trials) && (!TRACKER::countIsReached) && (!TRACKER::runtimeLimitReached) && (!TRACKER::idletimeLimitReached)) {

		method->reset();
		method->setSeed();

        Propagation_Decision algDecision = PROP_STOP;
		int propCounter = 0;
				
		// This is a good place to check the runTimeLimit
        if (!hitTimeLimit(walker)) {
			disp(MSG_DEBUG, "Initializing tracker, trial: %d", trialNo);
			algDecision = method->initialize();
        }
		
		// Algorithm could not initialize at the this point
		if (algDecision != PROP_CONTINUE) {

			disp(MSG_DEBUG, "Tracker failed to initialize, trial: %d", trialNo);

			walker->action 		  = FAIL;
			walker->failingReason = FAILED_BY_THE_ALGORITHM_AT_INITIALIZATION;

		// Algorithm initialized correctly, now check the pathway rules and iterate
		} else {

			disp(MSG_DEBUG, "Tracker is initialized, trial: %d", trialNo);

			walker->seedInd = 0;
			walker->begInd  = 0;
			walker->endInd  = 0;
			method->append(); // Append the seed point in the algorithm
			// TRACKER::pw.printWalker(walker);

			if (TRACKER::pw.directionality==ONE_SIDED) {

				walker->terminationReasonSideA = SEED_POINT;

				if (TRACKER::pw.checkSeed(walker,side_B)!=DISCARD) {
					// This will set to discard if seed point is inside a discard_if_inside region
					TRACKER::pw.checkEndsInside(&walker->streamline->at(walker->seedInd).x,walker);
				} else {
					disp(MSG_DEBUG, "checkSeed failed, discarding trial: %d", trialNo);
				}

			} else {
				
				// Checks the initialization and sets the seedInd of the walker in step
				// step can return CONTINUE or DISCARD
				if (TRACKER::pw.checkSeed(walker) != DISCARD) {
					disp(MSG_DEBUG, "checkSeed success, seed point is appended, trial: %d", trialNo);
				} else {
					disp(MSG_DEBUG, "checkSeed failed, discarding trial: %d", trialNo);
				}

			}

			bool appended = false;

			auto appendAndCheckPw = [&]()->void {
				if (!appended) {
					method->append();
					TRACKER::pw.checkWalker(walker);
					walker->endInd = streamline.size()-1;
					appended = true;
				}
			};

			auto algProp = [&]()->void {

				// First propagate the algorithm and then move the point
				// Note: algDecision is for the next propagation!
				// i.e. if initialization is successful, that means the method can propagate to the next point.
				// Therefore, we immediately move to the next point. THEN check whether the algorithm can move to the next point.
				// So, if algDecision below is PROP_FAIL, that means we can propagate but this is the last one.
				// If algDecision is PROP_STOP, that means we can propagate but we need to stop there.
				// We are doing it this way because after the algDecision, we still need to check the pw rules.
				algDecision = method->propagate();
				propCounter++;
				appended = false;

				// TRACKER::pw.printWalker(walker);

				if ((propCounter % method->appendInterval()) == 0) {
					appendAndCheckPw();
				}

				// PROP_FAIL case
				if (algDecision == PROP_FAIL) {
					walker->action 		  = FAIL;
					walker->failingReason = FAILED_BY_THE_ALGORITHM;
				}


			};

			// Track first side
			while ( (walker->action == CONTINUE) && (algDecision == PROP_CONTINUE) && !hitTimeLimit(walker) ) {
				algProp();
				// disp(MSG_DEBUG, "propCounter: %d",propCounter);
			}

			if ( (walker->action != DISCARD) && (walker->action != FAIL) ) {
				appendAndCheckPw();
			}

			if (TRACKER::pw.directionality==ONE_SIDED) {

				// ONE_SIDED tracking
				TRACKER::pw.tieSecondEnd(walker); // walker->action is either DISCARD or KEEP

				if (TRACKER::pw.skipSeedROI) {
					// Changes the beginning position of the streamline.
					// Sets walker->action to DISCARD if streamline is TOO_SHORT
					TRACKER::pw.skipSeed(walker,false); 
				}

			} else {		

				// TWO_SIDED tracking

				TRACKER::pw.tieFirstEnd(walker); // walker->action is either DISCARD or STOP

				if ((walker->action != DISCARD) && (walker->action != FAIL)) {

					disp(MSG_DEBUG, "Starting tracking of the second end");

					// Flip and get the algorithm decision for whether a propagation can be made towards the other end or not
					algDecision = method->flip();

					// Flipping was not possible, e.g. because posteriorMax was exceeded during flip when using ptt algorithm
					if (algDecision==PROP_FAIL) {

						walker->action 		  = FAIL;
						walker->failingReason = FAILED_BY_THE_ALGORITHM;
						disp(MSG_DEBUG, "Flipping was not possible");

					} else {

						// seedInd is now at the end
						walker->seedInd = walker->streamline->size()-1;
						TRACKER::pw.flipSide(walker);

						// Track other side
						while ( (walker->action == CONTINUE) && (algDecision == PROP_CONTINUE) && !hitTimeLimit(walker) ) {
							algProp();
							// disp(MSG_DEBUG, "propCounter: %d",propCounter);
						}
						
						if ( (walker->action != DISCARD) && (walker->action != FAIL) ) {
							appendAndCheckPw();
						}

						TRACKER::pw.tieSecondEnd(walker); // walker->action is either DISCARD or KEEP

						// TWO_SIDED tracking is now finished.
						// If walker->action is DISCARD, we will check the rules also with the switch sides.
						// if (walker->action == DISCARD) {

						// 	// We will not check if discardingReason is TOO_SHORT, TOO_LONG or REACHED_TIME_LIMIT
						// 	switch (walker->discardingReason) {
						// 		case TOO_SHORT:
						// 		case TOO_LONG:
						// 		case REACHED_TIME_LIMIT: break;
						// 		default: {
						// 			TRACKER::pw.softReset(walker);
						// 			if (TRACKER::pw.checkSeed(walker)!=DISCARD) {
						// 				TRACKER::pw.walkSecondEnd(walker);
						// 				TRACKER::pw.tieFirstEnd(walker);
						// 				if (walker->action!=DISCARD) {
						// 					TRACKER::pw.flipSide(walker);
						// 					TRACKER::pw.walkFirstEnd(walker);
						// 					TRACKER::pw.tieSecondEnd(walker);
						// 				}
						// 			}
						// 			break;
						// 		}
						// 	}
						// }

					}
				}

				

			}

		}

		if ( (walker->action != DISCARD) && (walker->action != FAIL)) walker->action = KEEP;
        
        // At this point tracking of the streamline is complete and one of the following decisions is made:
        // - DISCARD
        // - FAIL
		// - KEEP

		if (walker->action == DISCARD) {

			switch (walker->discardingReason) {
			case TOO_SHORT:
				disp(MSG_DEBUG, "TOO_SHORT");
				TRACKER::log_discard_TOO_SHORT.fetch_add(1);
				break;
			case TOO_LONG:
				disp(MSG_DEBUG, "TOO_LONG");
				TRACKER::log_discard_TOO_LONG.fetch_add(1);
				break;
			case DISCARD_REGION_REACHED:
				disp(MSG_DEBUG, "DISCARD_REGION_REACHED");
				TRACKER::log_discard_DISCARD_ROI_REACHED.fetch_add(1);
				break;
			case REQUIRED_ROI_NOT_MET:
				disp(MSG_DEBUG, "REQUIRED_ROI_NOT_MET");
				TRACKER::log_discard_REQUIRED_ROI_NOT_MET.fetch_add(1);
				break;
			case REQUIRED_ORDER_NOT_MET:
				disp(MSG_DEBUG, "REQUIRED_ROI_ORDER_NOT_MET");
				TRACKER::log_discard_REQUIRED_ROI_ORDER_NOT_MET.fetch_add(1);
				break;
			case CANT_MEET_STOP_CONDITION:
				disp(MSG_DEBUG, "CANT_MEET_STOP_CONDITION");
				TRACKER::log_discard_CANT_MEET_STOP_CONDITION.fetch_add(1);
				break;
			case ENDED_INSIDE_DISCARD_ROI:
				disp(MSG_DEBUG, "ENDED_INSIDE_DISCARD_ROI");
				TRACKER::log_discard_ENDED_INSIDE_DISCARD_ROI.fetch_add(1);
				break;
			case REACHED_TIME_LIMIT:
				disp(MSG_DEBUG, "REACHED_TIME_LIMIT");
				TRACKER::log_discard_REACHED_TIME_LIMIT.fetch_add(1);
				break;
			default:
				break;
			}

		}

		if (walker->action == FAIL) {

			switch (walker->failingReason) {
			case FAILED_BY_THE_ALGORITHM_AT_INITIALIZATION:
				disp(MSG_DEBUG, "FAILED_BY_THE_ALGORITHM_AT_INITIALIZATION");
				TRACKER::log_failed_BY_THE_ALGORITHM_AT_INITIALIZATION.fetch_add(1);
				break;
			case FAILED_BY_THE_ALGORITHM:
				disp(MSG_DEBUG, "FAILED_BY_THE_ALGORITHM");
				TRACKER::log_failed_BY_THE_ALGORITHM.fetch_add(1);
				break;
			default:
				break;
			}

		}

		if (walker->action != KEEP) {
			delete walker;
			streamline.clear();
			walker = TRACKER::pw.createWalker(&streamline,trialNo+1);
		}
		
		trialNo++;

	}


	if (walker->action == KEEP) {

		std::vector<std::vector<float>> outTrk;
		for (auto p : streamline) {
			outTrk.push_back({p.x,p.y,p.z});
		}

		{
			std::lock_guard<std::mutex> lock(NIBR::MT::PROC_MX());

			TRACKER::countIsReached = (TRACKER::tractogram.size()==size_t(TRACKER::seed.getCount()));

			if (!TRACKER::countIsReached) {

				TRACKER::tractogram.push_back(outTrk);

				TRACKER::lastTime = std::chrono::steady_clock::now();
			
				if      ((walker->terminationReasonSideA==MAX_LENGTH_REACHED)      || (walker->terminationReasonSideB==MAX_LENGTH_REACHED)     )
					TRACKER::log_success_REACHED_MAXLENGTH_LIMIT.fetch_add(1);
				else if ((walker->terminationReasonSideA==MIN_DATASUPPORT_REACHED) || (walker->terminationReasonSideB==MIN_DATASUPPORT_REACHED))
					TRACKER::log_success_REACHED_MINDATASUPPORT_LIMIT.fetch_add(1);
				else
					TRACKER::log_success_SATISFIED_PATHWAY_RULES.fetch_add(1);

			}

			TRACKER::countIsReached = (TRACKER::tractogram.size()==size_t(TRACKER::seed.getCount()));

		}

		disp(MSG_DETAIL, "Tracked %d in %d trials.", id, trialNo);

		return true;

	}

	disp(MSG_DETAIL, "Failed tracking %d in %d trials.", id, trialNo);
	return false;	

}
