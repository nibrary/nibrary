#include "pathway.h"

using namespace NIBR;

// When the second end is tied, the action is set to either DISCARD or KEEP
NIBR::Walker *NIBR::Pathway::tieSecondEnd(NIBR::Walker *w)
{

    // printWalker(w);

    // Check if too short
	if (w->trackedLength < minLength) {
		w->action = DISCARD;
		w->discardingReason = TOO_SHORT;
		return w;
	}

    if (tieFirstEnd(w)->action == DISCARD) return w;

    auto sideKeeper    = w->side;
    auto actionKeeper  = w->action;
    auto endKeeper     = w->segment.end;

    // If the current side is A, and there is no termination reason set for side B
    // this means, side B was terminated due to reaching low data support
    // then set the side to side_B, and tie the rules
    // note that in this case beginning of streamline is end of side_B
    if ((sideKeeper == side_A) && (w->terminationReasonSideB == TERMINATIONREASON_NOTSET) ) {

        w->side = side_B;        
        w->terminationReasonSideB = MIN_DATASUPPORT_REACHED;

        // Check if side_B ends inside a discard region
        w->segment.end  = &(w->streamline->at(0).x);
        if (tieDiscardRules(w)->action == DISCARD) return w;

        // Check if side_B satisfies required rules
        if (tieRequireRules(w)->action == DISCARD) return w;

    }

    // If the current side is B, and there is no termination reason set for side A
    // this means, side A was terminated due to reaching low data support
    // then set the side to side_A, and tie the rules
    // note that in this case beginning of streamline is end of side_A
    if ((sideKeeper == side_B) && (w->terminationReasonSideA == TERMINATIONREASON_NOTSET) ) {

        w->side = side_A;        
        w->terminationReasonSideA = MIN_DATASUPPORT_REACHED;

        // Check if side_A ends inside a discard region
        w->segment.end  = &(w->streamline->at(0).x);
        if (tieDiscardRules(w)->action == DISCARD) return w;

        // Check if side_A satisfies required rules
        if (tieRequireRules(w)->action == DISCARD) return w;

    }

    w->action       = actionKeeper;
    w->segment.end  = endKeeper;

    // If the current side is "either", 
    // then we first set the side to side_A and check whether that satisfies everything
    // if not then we set the side to side_B and check whether that satisfies everything
    // if these options don't work then we discard

    if (sideKeeper == either) {

        bool satisfiesAll = true;

        w->side                   = side_A;
        w->terminationReasonSideA = MIN_DATASUPPORT_REACHED;
        w->segment.end            = &(w->streamline->at(0).x);
        if (satisfiesAll && (tieDiscardRules(w)->action == DISCARD) ) satisfiesAll = false;
        if (satisfiesAll && (tieRequireRules(w)->action == DISCARD) ) satisfiesAll = false;

        if (satisfiesAll) {
            w->side                   = side_B;
            w->terminationReasonSideB = MIN_DATASUPPORT_REACHED;
            w->segment.end            = endKeeper;
            if (satisfiesAll && (tieDiscardRules(w)->action == DISCARD) ) satisfiesAll = false;
            if (satisfiesAll && (tieRequireRules(w)->action == DISCARD) ) satisfiesAll = false;
        }

        // Try the other way around
        if (!satisfiesAll) {

            satisfiesAll = true;

            w->side                   = side_A;
            w->terminationReasonSideA = MIN_DATASUPPORT_REACHED;
            w->segment.end            = endKeeper;
            if (satisfiesAll && (tieDiscardRules(w)->action == DISCARD) ) satisfiesAll = false;
            if (satisfiesAll && (tieRequireRules(w)->action == DISCARD) ) satisfiesAll = false;

            if (satisfiesAll) {
                w->side                   = side_B;
                w->terminationReasonSideB = MIN_DATASUPPORT_REACHED;
                w->segment.end            = &(w->streamline->at(0).x);
                if (satisfiesAll && (tieDiscardRules(w)->action == DISCARD) ) satisfiesAll = false;
                if (satisfiesAll && (tieRequireRules(w)->action == DISCARD) ) satisfiesAll = false;
            }

        }

        // This means w->action was set to DISCARD
        if (!satisfiesAll) return w;

    }
    

	// Final verification of pathway rules
	std::vector<bool>::iterator lbit = w->isDone.begin();
	for (std::vector<PathwayRule>::iterator it = prules.begin(); it != prules.end(); ++it)
	{
		if ( ((it->type == req_entry) || (it->type == req_exit) || (it->type == req_end_inside)) && ((*lbit) == false))
		{
			w->action = DISCARD;
			w->discardingReason = REQUIRED_ROI_NOT_MET;
			return w;
		}
		lbit++;
	}

	w->action      = KEEP;
	
	return w;
}