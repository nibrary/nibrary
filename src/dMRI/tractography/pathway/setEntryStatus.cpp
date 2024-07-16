#include "pathway.h"

// Check and set entry status. Discard if stop conditions can't be satisfied.
bool NIBR::Pathway::setEntryStatus(NIBR::Walker* w, int ruleNo) {

    // if (VERBOSE()==VERBOSE_DEBUG) {
    //     switch (w->entry_status[ruleNo]) {
    //     case notEnteredYet: disp(MSG_DEBUG,"Before setEntryStatus (%d): notEnteredYet",ruleNo);   break;
    //     case entered:       disp(MSG_DEBUG,"Before setEntryStatus (%d): entered",ruleNo);         break;
    //     case exited:        disp(MSG_DEBUG,"Before setEntryStatus (%d): exited",ruleNo);          break;
    //     case notExitedYet:  disp(MSG_DEBUG,"Before setEntryStatus (%d): notExitedYet",ruleNo);    break;
    //     }
    // }

    switch (w->entry_status[ruleNo]) {
    case notEnteredYet: w->entry_status[ruleNo] = isSegmentEntering(w, ruleNo) ? entered : notEnteredYet;   break;
    case entered:       w->entry_status[ruleNo] = isSegmentExiting (w, ruleNo) ? exited  : notExitedYet;    break;
    case exited:        w->entry_status[ruleNo] = isSegmentEntering(w, ruleNo) ? entered : exited;          break;
    case notExitedYet:  w->entry_status[ruleNo] = isSegmentExiting (w, ruleNo) ? exited  : notExitedYet;    break;
    }

    // if (VERBOSE()==VERBOSE_DEBUG) {
    //     switch (w->entry_status[ruleNo]) {
    //     case notEnteredYet: disp(MSG_DEBUG,"After setEntryStatus (%d): notEnteredYet",ruleNo);   break;
    //     case entered:       disp(MSG_DEBUG,"After setEntryStatus (%d): entered",ruleNo);         break;
    //     case exited:        disp(MSG_DEBUG,"After setEntryStatus (%d): exited",ruleNo);          break;
    //     case notExitedYet:  disp(MSG_DEBUG,"After setEntryStatus (%d): notExitedYet",ruleNo);    break;
    //     }
    // }

    if (VERBOSE()==VERBOSE_DEBUG) {
        switch (w->entry_status[ruleNo]) {
        case notEnteredYet: disp(MSG_DEBUG,"   notEnteredYet");   break;
        case entered:       disp(MSG_DEBUG,"   entered");         break;
        case exited:        disp(MSG_DEBUG,"   exited");          break;
        case notExitedYet:  disp(MSG_DEBUG,"   notExitedYet");    break;
        }
    }

    // Slightly move back before the stop rule
    if ( ((w->entry_status[ruleNo] == exited ) && (prules[ruleNo].type==stop_before_exit)) ||
         ((w->entry_status[ruleNo] == entered) && (prules[ruleNo].type==stop_before_entry)) ) {

        switch (srcType[ruleNo]) {

            // Nothing to do in these cases - these rules can't accept stop
            case undef_src:
            case res_pnt_src: {
                return true;
            }

            // Surface case
            // Move 1 micrometer before the stop rule
            case sph_src: 
            case surf_src: {

                // disp(MSG_DEBUG,"  Stopping streamline");
                w->segCrosLength -= EPS3/w->segment.len;

                if (w->segCrosLength <= 0) {
                    disp(MSG_DEBUG,"  Can't stop streamline: segment can't be shorter");
                    return false;
                }

                break;

            }
            
            // Image case
            // Move one downsampleFactor before the stop rule
            case img_mask_src:
            case img_label_src:
            case img_pvf_src: {

                // disp(MSG_DEBUG,"  Stopping streamline");
                float downsampleFactor = w->segment.len * maxSegSizeScaler[ruleNo];

                if (downsampleFactor > 1) {
                    float s = w->segment.len / float(std::ceil(downsampleFactor));
                    w->segCrosLength -= s/w->segment.len;
                } else {
                    w->segCrosLength  = 0;
                }

                break;
          
            }

        }

        return true;

    }


    // Slightly move forward after the stop rule
    if ( ((w->entry_status[ruleNo] == exited ) && (prules[ruleNo].type==stop_after_exit)) ||
         ((w->entry_status[ruleNo] == entered) && (prules[ruleNo].type==stop_after_entry)) ) {

        switch (srcType[ruleNo]) {

            // Nothing to do in these cases
            case undef_src:
            case res_pnt_src: {
                return true;
            }

            // Surface case
            // Move 1 micrometer after the stop rule
            case sph_src: 
            case surf_src: {

                // disp(MSG_DEBUG,"  Stopping streamline");
                w->segCrosLength += EPS3/w->segment.len; 
                
                if (w->segCrosLength>1) {
                    disp(MSG_DEBUG,"  Can't stop streamline: segment can't be longer");
                    return false;
                }

                break;

            }
            
            // Image case
            // Nothing needed since the end point is very likely to be not exactly on the border but already beyond the stop rule
            // i.e. when the last point is checked, it will always show to be beyond the stopping rule - it will not appear 50% inside/outside
            case img_mask_src:
            case img_label_src:
            case img_pvf_src: {
                break;
            }

        }

        return true;

    }

    return true;

}