#include "pathway.h"

// Check and set entry status. Discard if stop conditions can't be satisfied.
bool NIBR::Pathway::setEntryStatus(NIBR::Walker* w, int ruleNo) {

    bool isCrossing = false;
    float crossLen  = 0.0f;

    switch (w->entry_status[ruleNo]) {
        case notEnteredYet: {std::tie(isCrossing,crossLen) = isSegmentEntering(w->segment, ruleNo); w->entry_status[ruleNo] = isCrossing ? entered : notEnteredYet;   break;}
        case entered:       {std::tie(isCrossing,crossLen) =  isSegmentExiting(w->segment, ruleNo); w->entry_status[ruleNo] = isCrossing ? exited  : notExitedYet;    break;}
        case exited:        {std::tie(isCrossing,crossLen) = isSegmentEntering(w->segment, ruleNo); w->entry_status[ruleNo] = isCrossing ? entered : exited;          break;}
        case notExitedYet:  {std::tie(isCrossing,crossLen) =  isSegmentExiting(w->segment, ruleNo); w->entry_status[ruleNo] = isCrossing ? exited  : notExitedYet;    break;}
    }

    if (VERBOSE()==VERBOSE_DEBUG) {
        switch (w->entry_status[ruleNo]) {
        case notEnteredYet: disp(MSG_DEBUG,"   notEnteredYet");   break;
        case entered:       disp(MSG_DEBUG,"   entered");         break;
        case exited:        disp(MSG_DEBUG,"   exited");          break;
        case notExitedYet:  disp(MSG_DEBUG,"   notExitedYet");    break;
        }
    }


    // Stop exactly at exit/entry
    if ( ((w->entry_status[ruleNo] == exited ) && (prules[ruleNo].type==stop_at_exit)) ||
         ((w->entry_status[ruleNo] == entered) && (prules[ruleNo].type==stop_at_entry)) ) {

            disp(MSG_DEBUG,"  Stopping at exit/entry (by %.12f)", crossLen);
            w->segStopLength = crossLen;
            return true;

    }


    auto exitingOrEnteringWrapper = [&](bool testType)->bool {

        LineSegment checkSeg;
        checkSeg.beg    = new float[3]; 
        checkSeg.end    = new float[3];
        checkSeg.len    = w->segStopLength;

        for (int i = 0; i < 3; i++) {
            checkSeg.beg[i] = w->segment.beg[i];
            checkSeg.dir[i] = w->segment.dir[i];
            checkSeg.end[i] = checkSeg.beg[i] + checkSeg.dir[i] * checkSeg.len;    
        }

        auto [isCrossing, atLength] = testType ? isSegmentExiting(checkSeg, ruleNo) : isSegmentEntering(checkSeg, ruleNo);
        
        return isCrossing;
    };

    auto isExiting  = [&]()->bool { return exitingOrEnteringWrapper(true);  };
    auto isEntering = [&]()->bool { return exitingOrEnteringWrapper(false); };

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
            case surf_src: 
            case img_mask_src:
            case img_label_src: {
                disp(MSG_DEBUG,"  Stopping before exit/entry (by %.12f)", crossLen);
                w->segStopLength = crossLen - EPS3;
                break;
            }

            // Image pvf case
            // Move one downsampleFactor before the stop rule
            case img_pvf_src: {
                disp(MSG_DEBUG,"  Stopping before exit/entry (by %.12f)", crossLen);
                w->segStopLength = crossLen - miniSegment[ruleNo];
                break;
          
            }

        }

        if (w->segStopLength <= 0.0f) w->segStopLength = EPS3;

        if ((prules[ruleNo].type==stop_before_exit) && (isExiting())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not inside");
            // wait("...");
            return false;
        }

        if ((prules[ruleNo].type==stop_before_entry) && (isEntering())) {
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
            // Move 1 micrometer after the stop rule
            case sph_src: 
            case surf_src: 
            case img_mask_src:
            case img_label_src: {
                disp(MSG_DEBUG,"  Stopping after exit/entry (by %.12f)", crossLen);
                w->segStopLength = crossLen + EPS3;
                break;
            }

            // Image pvf case
            // Nothing needed since the end point is very likely to be not exactly on the border but already beyond the stop rule
            // i.e. when the last point is checked, it will always show to be beyond the stopping rule - it will not appear 50% inside/outside
            case img_pvf_src: {
                disp(MSG_DEBUG,"  Stopping after exit/entry (by %.12f)", crossLen);
                w->segStopLength = crossLen + miniSegment[ruleNo];
                break;
            }

        }

        if (w->segStopLength >= w->segment.len) w->segStopLength = w->segment.len - EPS3;

        if ((prules[ruleNo].type==stop_after_exit) && (!isExiting())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not outside");
            // wait("...");
            return false;
        }

        if ((prules[ruleNo].type==stop_after_entry) && (!isEntering())) {
            disp(MSG_DEBUG,"  Can't stop streamline: shorter segment is not inside");
            // wait("...");
            return false;
        }

        return true;

    }

    return true;

}