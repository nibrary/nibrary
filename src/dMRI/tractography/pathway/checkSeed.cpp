#include "pathway.h"

using namespace NIBR;

// This function returns either DISCARD or CONTINUE, and w->side is set to side_A/_B at the seed immediately.
// This function handles the following:
//    1. no seed, one_sided -> YES
//    2. no seed, two_sided -> YES
//    3.    seed, one_sided -> YES
//    4.    seed, two_sided -> NO. This case is handled by checkSeed(Walker *w)
NIBR::WalkerAction NIBR::Pathway::checkSeed(NIBR::Walker *w, NIBR::Tracking_Side side)
{
	// Prepare segment
	w->segment.beg    = &(w->streamline->at(w->seedInd).x);
    w->segment.end    = &(w->streamline->at(w->seedInd).x);
	w->segment.len    = 0.0;
	w->segment.dir[0] = 1.0;
	w->segment.dir[1] = 0.0;
	w->segment.dir[2] = 0.0;
	w->segCrosLength  = 1.0;

	w->side = side;
	// disp(MSG_DEBUG,"Walker side: %d ", w->side);

	for (int n = 0; n < ruleCnt; n++) {

		// disp(MSG_DEBUG,"checkSeed for prule: %d, side: %d",n,prules[n].side);

		// If there is seed, that is not allowed on edges, and if this is for the seed rule
		// unless the seed is surf_src when surfIs2D
		// then discard
		if (   hasOneSeed() && noEdgeSeeds && (n == theOneSeed)					    &&
		    !((prules[n].src == surf_src) && surfIs2D[n])    &&
			   isPointAtEdgeOfRule(&(w->streamline->at(w->seedInd).x),n,EPS3)) {

			w->action = DISCARD;
			w->discardingReason = IMPROPER_SEED;
			// disp(MSG_DEBUG,"Found edge seed");
			return DISCARD;

		}

		if ((prules[n].type) == seed) continue;

		bool inside = isPointInsideRule(&(w->streamline->at(w->seedInd).x),n);

		// If the seed is not already inside, if there is seed, that is not allowed on edges, and if this is not for the seed rule
		// if the seed point is close to the surf_src when surfIs2D
		// then consider inside
		if ( !inside && hasOneSeed() && noEdgeSeeds && (n != theOneSeed)			&&
		     ((prules[n].src == surf_src) && surfIs2D[n])   	&&
			   isPointAtEdgeOfRule(&(w->streamline->at(w->seedInd).x),n,EPS8)) {
			inside = true;
		}

		if (!inside) continue;

		w->entry_status[n] = entered;
		// disp(MSG_DEBUG,"Entered: %d, side: %d",n,prules[n].side);

		switch (prules[n].type) {

			case discard_seed: 
			{
				w->action = DISCARD;
				w->discardingReason = DISCARD_SEED;
				return DISCARD;	
			}
			
			case req_entry: 
			{
				if (prules[n].side == either) {

					w->isDone[n] = true;
					// disp(MSG_DEBUG,"Rule %d is done", n);

				} else if (prules[n].side == side) {

					if (satisfy_requirements_in_order == IN_ORDER) {
						
						if (w->side == side_A) {
							if (order_of_side_A_prules[w->sideAorder] != n) {
								w->action           = DISCARD;
								w->discardingReason = REQUIRED_ORDER_NOT_MET;
								return DISCARD;
							} else
								w->sideAorder++;
						}

						if (w->side == side_B) {
							if (order_of_side_B_prules[w->sideBorder] != n) {
								w->action           = DISCARD;
								w->discardingReason = REQUIRED_ORDER_NOT_MET;
								return DISCARD;
							} else
								w->sideBorder++;
						}

					}

					w->isDone[n] = true;
					// disp(MSG_DEBUG,"Rule %d is done", n);
				}

				break;
			}

			case req_exit: 
			{
				break;
			}

			// We will only handle the case that the next step is also inside seed. 
			// This is a reasonable assumption, and otherwise, this case will be very hard to handle.		
			case stop_at_entry:
			case stop_before_entry:
			case stop_after_entry:
			{
				w->action = DISCARD;
				w->discardingReason = IMPROPER_SEED;
				return DISCARD;
			}

			case stop_at_exit:
			case stop_before_exit:
			case stop_after_exit:
			{
				break;
			}

			// We will only handle the case that the next step is also inside seed. This is a reasonable assumption, and otherwise, this case will be very hard to handle.
			case discard_if_enters: 
			{
				w->action = DISCARD;
				w->discardingReason = DISCARD_REGION_REACHED;
				return DISCARD;
			}

			case discard_if_exits: 
			{
				break;
			}

			// This case is checked at the end
			case req_end_inside:
			case discard_if_ends_inside: 
			{
				break;
			}

			default:
			{
				break;
			}

		}

	}

	w->action = CONTINUE;
	return CONTINUE;
}



// This function returns either DISCARD or CONTINUE, and w->side can be set to either, side_A or side_B depending on the situation.
// This function handles the following:
//    1. no seed, one_sided -> NO.  This case is handled by checkSeed(Walker *w, Tracking_Side side)
//    2. no seed, two_sided -> NO.  This case is handled by checkSeed(Walker *w, Tracking_Side side)
//    3.    seed, one_sided -> NO.  This case is handled by checkSeed(Walker *w, Tracking_Side side)
//    4.    seed, two_sided -> YES. This case is handled by checkSeed(Walker *w)
//
// The following assignments can be made by the function.
//    1. If rules have side "either", then w->side also always have "either" side.
//    2. If rules have A/B side, then w->side can either remain as "either" or the side of the rule.

NIBR::WalkerAction NIBR::Pathway::checkSeed(NIBR::Walker *w)
{
	// Prepare segment
	w->segment.beg    = &(w->streamline->at(w->seedInd).x);
    w->segment.end    = &(w->streamline->at(w->seedInd).x);
	w->segment.len    = 0.0;
	w->segment.dir[0] = 1.0;
	w->segment.dir[1] = 0.0;
	w->segment.dir[2] = 0.0;
	w->segCrosLength  = 1.0;

	disp(MSG_DEBUG,"Walker side: %d", w->side);

	for (int n = 0; n < ruleCnt; n++) {

		disp(MSG_DEBUG,"checkSeed for prule: %d, side: %d",n,prules[n].side);

		// If there is seed, that is not allowed on edges, and if this is for the seed rule
		// unless the seed is surf_src when surfIs2D
		// then discard
		if (   hasOneSeed() && noEdgeSeeds && (n == theOneSeed)						&&
		    !((prules[n].src == surf_src) && surfIs2D[n])   	&&
			   isPointAtEdgeOfRule(&(w->streamline->at(w->seedInd).x),n,EPS3)) {

			w->action = DISCARD;
			w->discardingReason = IMPROPER_SEED;
			disp(MSG_DEBUG,"Found edge seed");
			return DISCARD;

		}

		if ((prules[n].type) == seed) continue;

		bool inside = isPointInsideRule(&(w->streamline->at(w->seedInd).x),n);

		// If the seed is not already inside, if there is seed, that is not allowed on edges, and if this is not for the seed rule
		// if the seed point is close to the surf_src when surfIs2D
		// then consider inside
		if ( !inside && hasOneSeed() && noEdgeSeeds && (n != theOneSeed)			&&
		     ((prules[n].src == surf_src) && surfIs2D[n])  	&&
			   isPointAtEdgeOfRule(&(w->streamline->at(w->seedInd).x),n,EPS8)) {
			inside = true;
		}

		if (!inside) continue;

		disp(MSG_DEBUG,"INSIDE");

		switch (prules[n].type) {

			case discard_seed: 
			{
				w->action = DISCARD;
				w->discardingReason = DISCARD_SEED;
				return DISCARD;	
			}
		
			case req_entry: 
			{

				if (prules[n].side == either) {
					w->entry_status[n] = entered;
					w->isDone[n]       = true;
					disp(MSG_DEBUG,"Entered side: %d", w->side);
					break;
				}
				
				if (w->side == either) { // Notice that w->side might be set to A/B by a previous rule.
					w->side = prules[n].side;
				}
				
				if (prules[n].side == w->side) {

					w->entry_status[n] = entered;
					disp(MSG_DEBUG,"Entered side: %d", w->side);

					if (satisfy_requirements_in_order == IN_ORDER) {
						
						if (w->side == side_A) {
							if (order_of_side_A_prules[w->sideAorder] != n) {
								w->action           = DISCARD;
								w->discardingReason = REQUIRED_ORDER_NOT_MET;
								return DISCARD;
							} else
								w->sideAorder++;
						}

						if (w->side == side_B) {
							if (order_of_side_B_prules[w->sideBorder] != n) {
								w->action           = DISCARD;
								w->discardingReason = REQUIRED_ORDER_NOT_MET;
								return DISCARD;
							} else
								w->sideBorder++;
						}

					}

					w->isDone[n] = true;
					break;
				}

				break;
			}

			// For the following cases, we only set the entry status and/or the side if needed
			case req_exit:
			case discard_if_exits:
			case stop_before_exit:
			case stop_after_exit:
			case stop_at_exit:
			{

				if (prules[n].side == either) {
					w->entry_status[n] = entered;
					disp(MSG_DEBUG,"Entered side: %d", w->side);
					break;
				}
				
				if (w->side == either) {
					w->side = prules[n].side;
				}
				
				if (prules[n].side == w->side) {
					w->entry_status[n] = entered;
					disp(MSG_DEBUG,"Entered side: %d", w->side);
					break;
				}

				break;
			}

			// We will only handle the case that the next step is also inside seed. 
			// This is a reasonable assumption, and otherwise, this case will be very hard to handle.
			case stop_at_entry:
			case stop_before_entry:
			case stop_after_entry:
			{
				w->action = DISCARD;
				w->discardingReason = IMPROPER_SEED;
				disp(MSG_DEBUG,"Stopped at seed. Discarding with IMPROPER_SEED.");
				return DISCARD;
			}

			// In these cases we will discard if needed. 
			case discard_if_enters:
			{
				w->action = DISCARD;
				w->discardingReason = DISCARD_REGION_REACHED;
				return DISCARD;	
			}

			// This case is checked at the end
			case req_end_inside:
			case discard_if_ends_inside: 
			{
				break;
			}

			default:
			{
				break;
			}

		}

	}

	w->action = CONTINUE;
	return CONTINUE;
}