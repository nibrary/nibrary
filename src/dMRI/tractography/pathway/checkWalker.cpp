#include "pathway.h"

using namespace NIBR;

NIBR::WalkerAction NIBR::Pathway::checkWalker(NIBR::Walker *w) {
	disp(MSG_DEBUG, "Checking between streamline indices %d - %d",int(w->streamline->size())-2, int(w->streamline->size())-1);
	return checkWalker(w,w->streamline->size()-2,w->streamline->size()-1);
}


// Return the action for the walker and the reason for that decision
// The returned action is either CONTINUE, DISCARD or STOP
NIBR::WalkerAction NIBR::Pathway::checkWalker(NIBR::Walker *w, int b, int e)
{

	// Prepare segment
	w->segment.beg   = &(w->streamline->at(b).x);
    w->segment.end   = &(w->streamline->at(e).x);
	vec3sub(w->segment.dir,w->segment.end,w->segment.beg);
    w->segment.len   = norm(w->segment.dir);
    normalize(w->segment.dir);
	w->segCrosLength = 1.0;

	for (int n = 0; n < ruleCnt; n++) {

		if (prules[n].type==seed) continue;

		disp(MSG_DEBUG,"Rule %d. Segment %d - %d. Checking...", ruleCnt, b, e);

		// No need to continue if the rule is done already
		if (w->isDone[n])
			continue;

		// disp(MSG_DEBUG,"Rule: %d, check 1", ruleCnt);
		// disp(MSG_DEBUG,"  w->side: %d, rule->side: %d", w->side, prules[n].side);

		// Continue for the following cases:
		//    1. prules.side = either
		//    2. w->side     = either
		//    3. prules.side = side_A and w->size_A too
		//    4. prules.side = side_B and w->size_B too
		if ((prules[n].side!=either) && (w->side != either) && (prules[n].side != w->side))
			continue;

		// disp(MSG_DEBUG,"Rule: %d, check 2", ruleCnt);

		// Check and set entry status. Discard if stop conditions can't be satisfied.
		// disp(MSG_DEBUG,"Checking %d - %d", b,e);
		if(!setEntryStatus(w, n)) {
			w->action 			= DISCARD;
            w->discardingReason = CANT_MEET_STOP_CONDITION;
			disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - CANT_MEET_STOP_CONDITION", ruleCnt, b, e);
			return DISCARD;
		}

		// disp(MSG_DEBUG,"Rule: %d, check 3", ruleCnt);

		// There is nothing to do if the there is no entry or exit
		if ( (w->entry_status[n] == notEnteredYet) || (w->entry_status[n] == notExitedYet))
			continue;

		// disp(MSG_DEBUG,"   Entered or exited: %d - %d", b,e);
		disp(MSG_DEBUG,"Rule %d. Segment %d - %d. Entered or exited", ruleCnt, b, e);

		

		// First handle the two discard cases:
		//    1. discard_if_enters 
		//    2. discard_if_exits
		// discard_if_ends_inside is handled at the end of the streamline.
		//
		// Therefore, no further check is necessary. If prules[n].side != w->side, code never comes until here.
		if (   (prules[n].type == discard_if_enters) || 
		     ( (prules[n].type == discard_if_exits ) && (w->entry_status[n]==exited) ) ) {

				w->action 			= DISCARD;
				w->discardingReason = DISCARD_REGION_REACHED;
				disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - DISCARD_REGION_REACHED", ruleCnt, b, e);
				return DISCARD;

		}

		// disp(MSG_DEBUG,"Rule: %d, check 4", ruleCnt);

		// Then handle the four remaining cases one by one:
		// 	  1. req_entry
		//    2. req_exit
		//    3. stop_at_entry, stop_before_entry, stop_after_entry
		//    4. stop_at_exit, stop_before_exit, stop_after_exit

		switch (prules[n].type) {

		// Both prules[n]->side and w->side can be "either"
		case req_entry:
		case req_exit:
		{

			// disp(MSG_DEBUG,"Rule: %d, req_entry, req_exit", ruleCnt);
			// if (prules[n].side == either) disp(MSG_DEBUG,"side: either");
			// if (prules[n].side == side_A) disp(MSG_DEBUG,"side: A");
			// if (prules[n].side == side_B) disp(MSG_DEBUG,"side: B");


			if ( ((prules[n].type==req_entry) && (w->entry_status[n]==entered)) ||
				 ((prules[n].type==req_exit)  && (w->entry_status[n]==exited )) ) {

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
								disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - REQUIRED_ORDER_NOT_MET", ruleCnt, b, e);
								return DISCARD;
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
								disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - REQUIRED_ORDER_NOT_MET", ruleCnt, b, e);
								return DISCARD;
							} else
								w->sideBorder++;	
						}

						w->isDone[n] = true;
					}
				}

			}


			break;
		}

		case stop_before_entry:
		case stop_at_entry:
		case stop_after_entry:
		case stop_before_exit:
		case stop_at_exit:
		case stop_after_exit:
		{
			if ((((prules[n].type==stop_before_entry) || (prules[n].type==stop_at_entry) || (prules[n].type==stop_after_entry)) && (w->entry_status[n]==entered)) ||
				(((prules[n].type==stop_before_exit ) || (prules[n].type==stop_at_exit ) || (prules[n].type==stop_after_exit )) && (w->entry_status[n]==exited )))
			{
				if (w->side == either)
					w->side = prules[n].side;

				if (w->side == side_A) {
					w->terminationReasonSideA = STOP_ROI_REACHED;
					w->action = STOP;
					w->isDone[n] = true;
				} else if (w->side == side_B) {
					w->terminationReasonSideB = STOP_ROI_REACHED;
					w->action = STOP;            
					w->isDone[n] = true;
				}

			}

		}

		default: { break; }
		
		}

	}


	if (w->action == STOP) {
		disp(MSG_DEBUG,"Shortened segment by %.12f%% and %.12f mm", w->segCrosLength, w->segment.len*(1.0f-w->segCrosLength));
		w->segment.len   *= w->segCrosLength;
		w->segCrosLength  = 1.0f;
		// Don't modifying segment end unless tracking
		if (isTracking) {
			w->segment.end[0] = w->segment.beg[0] + w->segment.dir[0] * w->segment.len;
			w->segment.end[1] = w->segment.beg[1] + w->segment.dir[1] * w->segment.len;
			w->segment.end[2] = w->segment.beg[2] + w->segment.dir[2] * w->segment.len;
		}
	}
	
	w->trackedLength += w->segment.len;


	// Check atMaxLength case
	// Instead of the precise value of maxLength, we will use maxLength-EPS4
	// This makes sure that if one runs discard atMaxLength at a later run with the same maxLength value, 
	// this streamline will not be discarded in some cases due to floating point errors
	if ( w->trackedLength > (maxLength-EPS4) ) {

		if (atMaxLength == ATMAXLENGTH_DISCARD) {
			w->action 			= DISCARD;
			w->discardingReason = TOO_LONG;
			disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - TOO_LONG", ruleCnt, b, e);
			return DISCARD;
		}

		// Handle ATMAXLENGTH_STOP
		else {			

			if (w->side == side_A)
				w->terminationReasonSideA = MAX_LENGTH_REACHED;
			else
				w->terminationReasonSideB = MAX_LENGTH_REACHED;

			w->action = STOP;

			float shorten = w->trackedLength - (maxLength-EPS4);

			// This could happen but it is not possible to shorten more than the segment length.
			if (shorten > w->segment.len) { 
				w->action 			= DISCARD;
				w->discardingReason = CANT_MEET_STOP_CONDITION;
				disp(MSG_DEBUG,"Rule %d. Segment %d - %d. DISCARD - CANT_MEET_STOP_CONDITION (Segment too short)", ruleCnt, b, e);
				return DISCARD;
			}

			w->segCrosLength   = 1.0f - shorten / w->segment.len;

			// Don't modifying segment end unless tracking
			if (isTracking) {
				disp(MSG_DEBUG,"Shortened segment by %.12f%% and %.12f mm", w->segCrosLength, shorten);
				w->segment.len    -= shorten;
				w->segCrosLength   = 1.0f;
				w->segment.end[0]  = w->segment.beg[0] + w->segment.dir[0]*w->segment.len;
				w->segment.end[1]  = w->segment.beg[1] + w->segment.dir[1]*w->segment.len;
				w->segment.end[2]  = w->segment.beg[2] + w->segment.dir[2]*w->segment.len;
				w->trackedLength  -= shorten;
			}

		}

	}


	if (w->action == STOP) {
		disp(MSG_DEBUG,"Rule %d. Segment %d - %d. STOP", ruleCnt, b, e);
		return STOP;
	}
	
	w->action = CONTINUE;
	disp(MSG_DEBUG,"Rule %d. Segment %d - %d. CONTINUE", ruleCnt, b, e);
	return CONTINUE;
}
