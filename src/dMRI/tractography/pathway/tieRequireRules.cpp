#include "pathway.h"

using namespace NIBR;

NIBR::Walker *NIBR::Pathway::tieRequireRules(NIBR::Walker *w)
{
    // Check if pathway rules are OK on this side
	if (w->side == side_A)
	{
        disp(MSG_DEBUG,"Checking side A");
		std::vector<bool>::iterator lbit = w->isDone.begin();
		for (std::vector<PathwayRule>::iterator it = prules.begin(); it != prules.end(); ++it)
		{
			if ((it->side == side_A) && ((it->type == req_entry) || (it->type == req_exit)) && ((*lbit) == false))
			{
                disp(MSG_DEBUG,"Does not meet rule %d for side A", it->uniqueId);
				w->action = DISCARD;
				w->discardingReason = REQUIRED_ROI_NOT_MET;
				return w;
			}

			if ((it->side == side_A) && (it->type == req_end_inside))
			{

				if (!isPointInsideRule(w->segment.end,it->uniqueId) && !isPointAtEdgeOfRule(w->segment.end,it->uniqueId,EPS2))
				// if (!isPointInsideRule(w->segment.end,it->uniqueId))
				{
					disp(MSG_DEBUG,"End point is not inside the required region");
					w->action = DISCARD;
					w->discardingReason = REQUIRED_ROI_NOT_MET;
					(*lbit) = false;
					return w;
				} else {
					disp(MSG_DEBUG,"End point is inside the required region");
					(*lbit) = true;
				}
			}

			lbit++;
		}
        disp(MSG_DEBUG,"Side A is OK");

	}
	
	else if (w->side == side_B)
	{
		
        disp(MSG_DEBUG,"Checking side B");
		std::vector<bool>::iterator lbit = w->isDone.begin();
		for (std::vector<PathwayRule>::iterator it = prules.begin(); it != prules.end(); ++it)
		{
			if ((it->side == side_B) && ((it->type == req_entry) || (it->type == req_exit)) && ((*lbit) == false))
			{
                disp(MSG_DEBUG,"Does not meet rule %d for side B", it->uniqueId);
				w->action = DISCARD;
				w->discardingReason = REQUIRED_ROI_NOT_MET;
				return w;
			}

			if ((it->side == side_B) && (it->type == req_end_inside))
			{
				if (!isPointInsideRule(w->segment.end,it->uniqueId) && !isPointAtEdgeOfRule(w->segment.end,it->uniqueId,EPS2))
				// if (!isPointInsideRule(w->segment.end,it->uniqueId))
				{
					disp(MSG_DEBUG,"End point is not inside the required region");
					w->action = DISCARD;
					w->discardingReason = REQUIRED_ROI_NOT_MET;
					(*lbit) = false;
					return w;
				} else {
					disp(MSG_DEBUG,"End point is inside the required region");
					(*lbit) = true;
				}
			}

			lbit++;
		}
        disp(MSG_DEBUG,"Side B is OK");
	
	}
	
	else {

		for (int n = 0; n < ruleCnt; n++)
		{
			if ((prules[n].type == req_end_inside) && (isPointInsideRule(w->segment.end,n) || isPointAtEdgeOfRule(w->segment.end,n,EPS2)))
			// if ((prules[n].type == req_end_inside) && isPointInsideRule(w->segment.end,n))
			{

				if (prules[n].side == either) {

					w->isDone[n] = true;

				} else {

					if (w->side == either)
						w->side = prules[n].side;

					if (w->side == side_A) {
						
						if (satisfy_requirements_in_order == IN_ORDER) {
							if (order_of_side_A_prules[w->sideAorder] != n) {
								w->action           = DISCARD;
								w->discardingReason = REQUIRED_ORDER_NOT_MET;
								return w;
							} else
								w->sideAorder++;	
						}

						w->isDone[n] = true;
					}

					if (w->side == side_B) {
						
						if (satisfy_requirements_in_order == IN_ORDER) {
							if (order_of_side_B_prules[w->sideBorder] != n) {
								w->action           = DISCARD;
								w->discardingReason = REQUIRED_ORDER_NOT_MET;
								return w;
							} else
								w->sideBorder++;	
						}

						w->isDone[n] = true;
					}
				}

			}
		}

	}
	
	for (int n = 0; n < ruleCnt; n++)
	{
		if (prules[n].type == req_end_inside)
		{

            disp(MSG_DEBUG,"Checking if first end ends within rule %d",n);

            if ((prules[n].side == either) ||        (w->side == either)  ||
				     (((w->side == side_A) && (prules[n].side == side_A)) ||
				   	  ((w->side == side_B) && (prules[n].side == side_B))))
				{

                    if (!isPointInsideRule(w->segment.end,n) && !isPointAtEdgeOfRule(w->segment.end,n,EPS4)) {
					// if (!isPointInsideRule(w->segment.end,n)) {
                        disp(MSG_DEBUG,"End point is not inside the required region");
						w->action = DISCARD;
						w->discardingReason = REQUIRED_ROI_NOT_MET;
                        return w;
                    } else {
                        disp(MSG_DEBUG,"End point is inside the required region");
						w->isDone[n] = true;
                    }

                }
		}
	}

	return w;

}