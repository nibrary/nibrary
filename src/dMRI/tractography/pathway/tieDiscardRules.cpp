#include "pathway.h"

using namespace NIBR;

// Discard if it ends inside a discard region
NIBR::Walker *NIBR::Pathway::tieDiscardRules(NIBR::Walker *w)
{

	for (int n = 0; n < ruleCnt; n++)
	{
		if (prules[n].type == discard_if_ends_inside)
		{

            disp(MSG_DEBUG,"Checking if first end ends within rule %d",n);

            if ((prules[n].side == either) || (w->side == either) ||
				     (((w->side == side_A) && (prules[n].side == side_A)) ||
				   	  ((w->side == side_B) && (prules[n].side == side_B))))
				{

                    // if (isPointInsideRule(w->segment.end,n) || isPointAtEdgeOfRule(w->segment.end,n,EPS4)) {
                    
                    if (isPointInsideRule(w->segment.end,n)) {
                        disp(MSG_DEBUG,"End point is inside the discard region");
                        w->action = DISCARD;
                        w->discardingReason = ENDED_INSIDE_DISCARD_ROI;
                        return w;
                    } else {
                        disp(MSG_DEBUG,"End point is outside the discard region");
                    }

                }
		}
	}

    return w;

}