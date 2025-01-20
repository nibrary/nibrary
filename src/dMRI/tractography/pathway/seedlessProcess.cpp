#include "pathway.h"
#include "pathwayRule.h"
#include "../utility/streamline_operators.h"

using namespace NIBR;

// With seedlessProcess we handle cases with:
// - no seed
// - no one_sided
// - no stop rules
// - no stopAtMax
// - no req_exit
// - no discard_if_exits
//
// With these we allow:
// - minLength
// - maxLength
// - in_order
// - req_end_inside             (has to be defined with A/B)
// - discard_if_ends_inside     (has to be defined with A/B)
// - req_entry                  (cannot have A/B)
// - discard_if_enters          (cannot have A/B)

void NIBR::Pathway::seedlessProcess(NIBR::Walker* w) {

    // disp(MSG_DEBUG, "Applying seedless filter");
    // disp(MSG_DEBUG, "Processing: %d",w->ind);

    // First check length thresholds if possible
    if ( (minLength>0) || (maxLength<FLT_MAX) ) {

        float length = getStreamlineLength(*(w->streamline));

        if (length<minLength) {
            w->trackedLength    = length;
            w->action           = DISCARD;
            w->discardingReason = TOO_SHORT;
            return;
        }

        if (length>maxLength) {
            w->trackedLength    = length;
            w->action           = DISCARD;
            w->discardingReason = TOO_LONG;
            return;
        }
        
    }
	

    // First check end point rules
    float firstPoint[3];
    firstPoint[0]  = w->streamline->front().x;
    firstPoint[1]  = w->streamline->front().y;
    firstPoint[2]  = w->streamline->front().z;

    float lastPoint[3];
    lastPoint[0]  = w->streamline->back().x;
    lastPoint[1]  = w->streamline->back().y;
    lastPoint[2]  = w->streamline->back().z;

    bool discard = false;
    bool keep    = true;

    auto checkEnds = [&](float* A_end, float* B_end)->void {

        discard = false;
        keep    = true;

        for (size_t n=0; n<prules.size(); n++) {

            bool inside_A_end = isPointInsideRule(A_end,n);
            bool inside_B_end = isPointInsideRule(B_end,n);
        
            if ((inside_A_end && (prules[n].side == side_A)) ||
                (inside_B_end && (prules[n].side == side_B)) ) 
            {
                if (prules[n].type == discard_if_ends_inside)   {discard = true; return;}
            }

            if ((!inside_A_end && (prules[n].side == side_A)) ||
                (!inside_B_end && (prules[n].side == side_B)) ) 
            {
                if (prules[n].type == req_end_inside)           {keep = false; return;}
            }

        }

    };

    // Check whether setting first point is in side_A and last point is in side_B works
    checkEnds(firstPoint,lastPoint); 

    // We try the other way around, with last point on side_B and first point on side_B
    if (discard || !keep) {       
        checkEnds(lastPoint,firstPoint); 
    }

    if (discard) {
        w->action           = DISCARD;
        w->discardingReason = ENDED_INSIDE_DISCARD_ROI;
        return;
    }

    if (!keep) {
        w->action           = DISCARD;
        w->discardingReason = REQUIRED_ROI_NOT_MET;
        return;
    }

    // Check ends for discard_if_enters
    for (size_t n=0; n<prules.size(); n++) {

        bool isInside = isPointInsideRule(firstPoint,n) || isPointInsideRule(lastPoint,n);

        if ((prules[n].type == discard_if_enters) && isInside) { 
            w->action 			= DISCARD;
            w->discardingReason = DISCARD_REGION_REACHED;
            disp(MSG_DEBUG,"Rule %d. DISCARD - DISCARD_REGION_REACHED");
            return;
        }

    }


    // Check end for req_entry
    int reqOrder = 0;

    auto checkReqEntry = [&](int n)->NIBR::WalkerAction {

        if (satisfy_requirements_in_order == IN_ORDER) {
            if (order_of_prules[reqOrder] != n) {
                w->action           = DISCARD;
                w->discardingReason = REQUIRED_ORDER_NOT_MET;
                disp(MSG_DEBUG,"Rule %d. DISCARD - REQUIRED_ORDER_NOT_MET", n);
                return DISCARD;
            } else {
                reqOrder++;
            }
        }

        w->entry_status[n]  = entered;
        w->isDone[n]        = true;

        return CONTINUE;

    };

    auto checkEndForReqEntry = [&](float* p)->NIBR::WalkerAction {
        
        reqOrder = 0; // Let's initialize this here to start following the reqOrder

        for (size_t n=0; n<prules.size(); n++) {

            bool isInside = isPointInsideRule(p,n);

            if ((prules[n].type == req_entry) && isInside) { 

                if (checkReqEntry(n) == DISCARD) 
                    return DISCARD;

            }

        }

        return CONTINUE;

    };

    
    auto checkSegment = [&](int b, int e)->NIBR::WalkerAction {

        bool isSegmentReady = false;

        for (size_t n = 0; n < prules.size(); n++) {

            if (w->isDone[n]) 
                continue;

            disp(MSG_DEBUG,"Rule %d. Segment %d - %d. Checking...", n, b, e);

            // Prepare segment
            if (isSegmentReady == false) {
                
                w->segment.beg   = &(w->streamline->at(b).x);
                w->segment.end   = &(w->streamline->at(e).x);

                vec3sub(w->segment.dir,w->segment.end,w->segment.beg);
                w->segment.len   = norm(w->segment.dir);
                w->segStopLength = w->segment.len;
                normalize(w->segment.dir);

                isSegmentReady = true;
            }

            setEntryStatus(w, n);   // Returns false only when stop rules can't be met, which is not a case we need to consider here

            if (w->entry_status[n] != entered) // notEnteredYet, notExitedYet, exited states are not relevant in this case
			    continue;

            if (prules[n].type == discard_if_enters) {
                w->action 			= DISCARD;
                w->discardingReason = DISCARD_REGION_REACHED;
                disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - DISCARD_REGION_REACHED", n, b, e);
                return DISCARD;
            }

            if (prules[n].type == req_entry) {
                if (checkReqEntry(n) == DISCARD) 
                    return DISCARD;
            }

        }

        return CONTINUE;

    };


    // Walk one way
    if (checkEndForReqEntry(firstPoint) == CONTINUE) {

        for (size_t i = 1; i < w->streamline->size(); i++) {
            if (checkSegment(i-1,i) != CONTINUE) 
                break;
        }

    }

    
    // Walk the other way if DISCARD due to REQUIRED_ORDER_NOT_MET
    if ((w->action == DISCARD) && (w->discardingReason == REQUIRED_ORDER_NOT_MET)) {

        for (size_t n = 0; n < prules.size(); n++) {
            w->entry_status[n] = notEnteredYet;
            w->isDone[n]       = false;
        }

        if (checkEndForReqEntry(lastPoint) == CONTINUE) {
            for (size_t i = w->streamline->size()-1; i > 0; i--) {
                if (checkSegment(i,i-1) != CONTINUE) 
                    break;
            }
        }

    }


    // Check all req_entry rules were met
    for (size_t n = 0; n < prules.size(); n++) {
        if ((prules[n].type == req_entry) && (!w->isDone[n])) {
            w->action 			= DISCARD;
            w->discardingReason = REQUIRED_ROI_NOT_MET;
        }
    }


    if (w->action == DISCARD) 
        return;
    
    // Streamline will not be cropped anywhere
    w->terminationReasonSideA = MIN_DATASUPPORT_REACHED;
    w->terminationReasonSideB = MIN_DATASUPPORT_REACHED;

    w->begInd = 0;
    w->endInd = w->streamline->size()-1;

    w->action = KEEP;

    return;
    
}