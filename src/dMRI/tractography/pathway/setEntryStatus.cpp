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

    auto confirmIsInside = [&]()->bool {
        float p[3];
        p[0] = w->segment.beg[0] + w->segment.dir[0] * w->segCrosLength * w->segment.len;
        p[1] = w->segment.beg[1] + w->segment.dir[1] * w->segCrosLength * w->segment.len;
        p[2] = w->segment.beg[2] + w->segment.dir[2] * w->segCrosLength * w->segment.len;
        return isPointInsideRule(&p[0],ruleNo);
    };

    auto confirmIsOutside = [&]()->bool {
        float p[3];
        p[0] = w->segment.beg[0] + w->segment.dir[0] * w->segCrosLength * w->segment.len;
        p[1] = w->segment.beg[1] + w->segment.dir[1] * w->segCrosLength * w->segment.len;
        p[2] = w->segment.beg[2] + w->segment.dir[2] * w->segCrosLength * w->segment.len;
        return !isPointInsideRule(&p[0],ruleNo);
    };

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
            // Move 1 micrometer before the stop rule when segCrosLength is precise
            case sph_src: 
            case surf_src: 
            case img_mask_src:
            case img_label_src: {
                disp(MSG_DEBUG,"  Stopping before exit/entry (by %.12f)", w->segCrosLength);
                w->segCrosLength -= EPS3 / w->segment.len;
                break;
            }

            // Image pvf case
            // Move one downsampleFactor before the stop rule
            case img_pvf_src: {
                float downsampleFactor = w->segment.len * maxSegSizeScaler[ruleNo];

                if (downsampleFactor > 1) {
                    float s = w->segment.len / float(std::ceil(downsampleFactor));
                    w->segCrosLength -= s / w->segment.len;
                }

                break;
          
            }

        }

        if (w->segCrosLength <= 0.0) {
            disp(MSG_DEBUG,"  Can't stop streamline: segment can't be shorter");
            // wait("...");
            return false;
        }

        if ((prules[ruleNo].type==stop_before_exit) && (!confirmIsInside())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not inside");
            // wait("...");
            return false;
        }

        if ((prules[ruleNo].type==stop_before_entry) && (!confirmIsOutside())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not outside");
            // wait("...");
            return false;
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
            // Move 1 micrometer after the stop rule when segCrosLength is precise
            case sph_src: 
            case surf_src: 
            case img_mask_src:
            case img_label_src: {
                disp(MSG_DEBUG,"  Stopping after exit/entry (by %.12f)", w->segCrosLength);
                w->segCrosLength += EPS3 / w->segment.len;
                break;
            }

            // Image pvf case
            // Nothing needed since the end point is very likely to be not exactly on the border but already beyond the stop rule
            // i.e. when the last point is checked, it will always show to be beyond the stopping rule - it will not appear 50% inside/outside
            case img_pvf_src: {
                break;
            }

        }

        if (w->segCrosLength>1.0f) {
            disp(MSG_DEBUG,"  Can't stop streamline: segment can't be longer");
            // wait("...");
            return false;
        }

        if ((prules[ruleNo].type==stop_after_exit) && (!confirmIsOutside())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not outside");
            // wait("...");
            return false;
        }

        if ((prules[ruleNo].type==stop_after_entry) && (!confirmIsInside())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not inside");
            // wait("...");
            return false;
        }

        return true;

    }

    return true;

}