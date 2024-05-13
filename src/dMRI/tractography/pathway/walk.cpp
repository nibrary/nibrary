#include "pathway.h"
#include "pathwayRule.h"

using namespace NIBR;

// For offline mode - Walk the walker
//
// The following lists all the possible scenerios that are handled:
//   1. with seed, one_sided,   either ->     OK, inOrder allowed and can stop   (internally isSided is set to true. seed is set to side_A, either rules are converted to side_B.)
//   2. with seed, one_sided, side_A/B -> NOT OK
//   3. with seed, two_sided,   either ->     OK, inOrder allowed but can't stop (isSided = false)
//   4. with seed, two_sided, side_A/B ->     OK, inOrder allowed and can stop   (isSided = true )
//   5.   no seed, one_sided           -> NOT OK
//   6.   no seed, two_sided,   either ->     OK, inOrder allowed but can't stop (internally isSided is set to true. Ends are set as A/B, depending on whichever propagation direction satisfies the pathway rules.)
//   7.   no seed, two_sided, side_A/B -> NOT OK
//
void NIBR::Pathway::walk(NIBR::Walker* walker) {

    // disp(MSG_DEBUG, "Processing: %d",walker->ind);

    // First check length thresholds if possible
    if ( (canStop()==false) && ( (minLength>0) || (maxLength<FLT_MAX) ) ) {

        walker->begInd = 0;
        walker->endInd = walker->streamline->size()-1;

        float length = getStreamlineLength(walker);

        if (length<minLength) {
            walker->trackedLength    = length;
            walker->action           = DISCARD;
            walker->discardingReason = TOO_SHORT;
            return;
        }

        if (length>maxLength) {
            walker->trackedLength    = length;
            walker->action           = DISCARD;
            walker->discardingReason = TOO_LONG;
            return;
        }

    }

    // Tracking has no seed so it is always considered as two_sided.
    // Seed point is assigned as one end of the streamline, and it is assigned as side_A. 
    // The other side is side_B. i.e. w->either option is NOT possible.
    if (hasOneSeed()==false) {

        // If there is no seed, there can't be stop_at_exit/entry so we always preserve the whole streamline.
        // Streamline can either be kept or discarded
        walker->begInd = 0;
        walker->endInd = walker->streamline->size()-1;

        // To speed things up, quickly discard for discard_if_ends_inside, discard_if_enters and !req_end_inside
        if (checkEndsInside(&walker->streamline->front().x,walker)) return;
        if (checkEndsInside(&walker->streamline->back().x, walker)) return;

        // Starting from the seed, walk towards the end of the streamline and get endInd
        walker->seedInd = 0;
        walker->terminationReasonSideA = SEED_POINT;        // Side_A is terminated at the seed
        if (checkSeed(walker,side_B) != DISCARD) {          // Side_B is tracked

            for (size_t i = 1; i < walker->streamline->size(); i++) {
                // disp(MSG_DEBUG,"Checking: %d -> %d",i-1,i);
                if (checkWalker(walker,i-1,i) != CONTINUE) {
                    // disp(MSG_DEBUG,"Break at: %d",i);
                    break;
                }
            }

            tieSecondEnd(walker);
            
            if (walker->action == KEEP)
                return;

        }
        
        // If walker is discarded because of required order was not met, then we walk also from the other end
        if (walker->discardingReason == REQUIRED_ORDER_NOT_MET) {
            softReset(walker);

            walker->seedInd = walker->streamline->size()-1;
            walker->terminationReasonSideA = SEED_POINT;    // Side_A is terminated at the seed
            if (checkSeed(walker,side_B)==DISCARD) return;  // Side_B is tracked

            for (size_t i = walker->streamline->size()-1; i > 0; i--) {
                // disp(MSG_DEBUG,"Checking: %d -> %d",i,i-1);
                if (checkWalker(walker,i,i-1) != CONTINUE) {
                    // disp(MSG_DEBUG,"Break at: %d",i);
                    break;
                }
            }

            tieSecondEnd(walker);
        }
        
        return;

    } else {

        // In the case where seed is provided

        // First get a valid seed point on the given seed region.
        // This functions add this seed point in the streamline if needed, which is mostly the case.
        // We will try with as many seeds as seedTrials, so we loop with do while

        int trials = 0;

        do {

            disp(MSG_DEBUG,"\n\nTRIAL: %d", trials);
            softReset(walker);         

            if (!getSeedInd(walker)) {
                walker->action           = DISCARD;
                walker->discardingReason = SEED_NOT_FOUND;
                return;
            }

            // disp(MSG_DEBUG,"Seed index: %d", walker->seedInd);
            // disp(MSG_DEBUG,"Streamline size: %d", walker->streamline->size());
            // disp(MSG_DEBUG,"Seed pos: [%.2f,%.2f,%.2f]", walker->streamline->at(walker->seedInd).x,walker->streamline->at(walker->seedInd).y,walker->streamline->at(walker->seedInd).z);
            
            // Tracking has seed and it is one_sided.
            // In the output:
            //     begInd shows the index of the seed point
            //     endInd shows the other end
            // Therefore begInd can be bigger than endInd
            // Seed point is one end of the streamline, and it is assigned as side_A immediately. The other side is side_B.
            // inOrder is allowed and can stop. Both begInd and endInd can ben floating, since the seed point can be in the middle of a segment
            if (directionality==ONE_SIDED) {

                // Starting from the seed, walk towards the end of the streamline and get endInd
                // endInd can be floating if termination is due to a stop rule
                walker->terminationReasonSideA = SEED_POINT;    // Side_A is terminated at the seed
                if (checkSeed(walker,side_B)==DISCARD) continue;

                // If seed point is inside a discard_if_inside region, immediately discard.
                // Notice that we are doing this here, because discard_if_inside is B_pulled so needs to be checked for walker side_B. The check is therefore done for the seed point.
                if (checkEndsInside(&walker->streamline->at(walker->seedInd).x,walker)) 
                    continue;

                // We first try to walk towards the first end, i.e. seed -> walker.endInd
                // disp(MSG_DEBUG, "Walking the first end");
                walkFirstEnd(walker);

                // If the seed point is inserted we need to shift the begInd accordingly
                // begInd is not an integer if the seed is found in the middle of a segment
                if (walker->seedInserted) {
                    float prevLen = dist(walker->streamline->at(walker->seedInd),walker->streamline->at(walker->seedInd-1));
                    float nextLen = dist(walker->streamline->at(walker->seedInd),walker->streamline->at(walker->seedInd+1));
                    float corrLen = prevLen / (prevLen + nextLen);
                    walker->begInd = walker->seedInd - 1 + corrLen;
                } else {
                    walker->begInd = walker->seedInd;
                }
                tieSecondEnd(walker);
                // disp(MSG_DEBUG, "begInd: %.2f", walker->begInd);
                // disp(MSG_DEBUG, "endInd: %.2f", walker->endInd);

                if (walker->action == KEEP) {
                    if (skipSeedROI && skipSeed(walker,false)) {
                        return; // KEEP case
                    }
                    return; // KEEP case
                }


                // The seed -> walker.endInd did not meet the needs so we try seed -> walker.begInd other way around.
                // disp(MSG_DEBUG, "Walking the second end");
                softReset(walker);
                
                // Starting from the seed, walk towards the beginning of the streamline and get begInd
                // begInd can be floating if termination is due to a stop rule
                walker->terminationReasonSideA = SEED_POINT;      // Side_A is terminated at the seed
                if (checkSeed(walker,side_B)==DISCARD) continue;
                walkSecondEnd(walker);
                
                // If the seed point is inserted we need to shift the endInd accordingly
                // endInd is not an integer if the seed is found in the middle of a segment
                if (walker->seedInserted) {
                    float prevLen = dist(walker->streamline->at(walker->seedInd),walker->streamline->at(walker->seedInd-1));
                    float nextLen = dist(walker->streamline->at(walker->seedInd),walker->streamline->at(walker->seedInd+1));
                    float corrLen = prevLen / (prevLen + nextLen);
                    walker->endInd = walker->seedInd - 1 + corrLen;
                } else {
                    walker->endInd = walker->seedInd;
                }

                // We will switch values of begInd and endInd here so the seed point is at the begInd
                std::swap(walker->endInd,walker->begInd);

                tieSecondEnd(walker);

                if (walker->action == KEEP) {
                    if (skipSeedROI) {
                        if (skipSeedROI && skipSeed(walker,true)) {
                            return; // KEEP case
                        }
                        return; // KEEP case
                    }

                }

                return; // KEEP case
                
            } 


            if (checkSeed(walker)!=DISCARD) {
                disp(MSG_DEBUG,"Getting endInd");
                walkFirstEnd(walker);
                if (walker->side==side_A)       disp(MSG_DEBUG,"side: A, endInd: %.2f", walker->endInd);
                else if (walker->side==side_B)  disp(MSG_DEBUG,"side: B, endInd: %.2f", walker->endInd);
                else                            disp(MSG_DEBUG,"side: either, endInd: %.2f", walker->endInd);
                disp(MSG_DEBUG,"Tying first end");
                tieFirstEnd(walker);
                if (walker->action!=DISCARD) {
                    disp(MSG_DEBUG,"Flipping");
                    flipSide(walker);
                    disp(MSG_DEBUG,"Getting begInd");
                    walkSecondEnd(walker);
                    if (walker->side==side_A)       disp(MSG_DEBUG,"side: A, begInd: %.2f", walker->begInd);
                    else if (walker->side==side_B)  disp(MSG_DEBUG,"side: B, begInd: %.2f", walker->begInd);
                    else                            disp(MSG_DEBUG,"side: either, begInd: %.2f", walker->begInd);
                    disp(MSG_DEBUG,"Tying second end");
                    tieSecondEnd(walker);
                }
            }

            if (isnan(walker->begInd)) disp(MSG_FATAL,"NAN begInd, id: %d",walker->ind);
            if (isnan(walker->endInd)) disp(MSG_FATAL,"NAN endInd, id: %d",walker->ind);

            disp(MSG_DEBUG,"[begInd,endInd]: %.2f-%.2f",walker->begInd,walker->endInd);


            if (walker->action == KEEP)
                return;

            disp(MSG_DEBUG,"Trying the other way around");
            softReset(walker);
            if (checkSeed(walker)!=DISCARD) {
                disp(MSG_DEBUG,"Getting endInd");
                walkSecondEnd(walker);
                if (walker->side==side_A)       disp(MSG_DEBUG,"side: A, endInd: %.2f", walker->endInd);
                else if (walker->side==side_B)  disp(MSG_DEBUG,"side: B, endInd: %.2f", walker->endInd);
                else                            disp(MSG_DEBUG,"side: either, endInd: %.2f", walker->endInd);
                disp(MSG_DEBUG,"Tying first end");
                tieFirstEnd(walker);
                if (walker->action!=DISCARD) {
                    disp(MSG_DEBUG,"Flipping");
                    flipSide(walker);
                    disp(MSG_DEBUG,"Getting begInd");
                    walkFirstEnd(walker);
                    if (walker->side==side_A)       disp(MSG_DEBUG,"side: A, begInd: %.2f", walker->begInd);
                    else if (walker->side==side_B)  disp(MSG_DEBUG,"side: B, begInd: %.2f", walker->begInd);
                    else                            disp(MSG_DEBUG,"side: either, begInd: %.2f", walker->begInd);
                    disp(MSG_DEBUG,"Tying second end");
                    tieSecondEnd(walker);
                }
            }

            if (isnan(walker->begInd)) disp(MSG_FATAL,"NAN begInd, id: %d",walker->ind);
            if (isnan(walker->endInd)) disp(MSG_FATAL,"NAN endInd, id: %d",walker->ind);

            disp(MSG_DEBUG,"[begInd,endInd]: %.2f-%.2f",walker->begInd,walker->endInd);

            if (walker->action == KEEP)
                return;
            
            trials++;
        } 
        while (trials < seedTrials);

    }
}