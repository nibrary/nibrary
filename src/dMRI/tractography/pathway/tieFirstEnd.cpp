#include "pathway.h"

using namespace NIBR;

// When the first end is tied, the action is set to either DISCARD or STOP
// Tying the first end does not change w->side
// i.e. if w->side = either, it will check as either, and it will remain that way
NIBR::Walker *NIBR::Pathway::tieFirstEnd(NIBR::Walker *w)
{

	if (w->action == DISCARD) {
        disp(MSG_DEBUG,"Discarding streamline");
		return w;
    }

	// If the action was not already set to STOP until here, this means that fiber tracking stoped due to reaching min data support
    // And if w->side has not been set until here, leave the termination reason as TERMINATIONREASON_NOTSET,
    // the side and termination reason will be second while tying the second end.
	if (w->action != STOP)
	{
        disp(MSG_DEBUG,"First end reached min data support");

		if (w->side == side_A)
			w->terminationReasonSideA = MIN_DATASUPPORT_REACHED;
		else if (w->side == side_B)
			w->terminationReasonSideB = MIN_DATASUPPORT_REACHED;

		w->action = STOP;
		w->segCrosLength = 1.0;
	}

	if (tieDiscardRules(w)->action == DISCARD) return w;	
    if (tieRequireRules(w)->action == DISCARD) return w;

    disp(MSG_DEBUG,"First end is tied successfully");

	return w;
}