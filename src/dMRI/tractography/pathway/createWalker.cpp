#include "pathway.h"

using namespace NIBR;

// Create a walker
NIBR::Walker* NIBR::Pathway::createWalker(std::vector<Point>* streamline) {
    
    NIBR::Walker* walker            = new NIBR::Walker();
    walker->streamline              = streamline;
    walker->ind                     = -1;
    walker->side                    = either;
    walker->sideAorder              = 0;
    walker->sideBorder              = 0;
    walker->action                  = CONTINUE;
    walker->seedInd                 = -1;
    walker->seedRange               = std::vector<float>();
    walker->seedInserted            = false;
    walker->begInd                  = 0;
    walker->endInd                  = 0;
    walker->segment.beg             = NULL;
    walker->segment.end             = NULL;
    walker->segment.len             = 0;
    walker->segment.dir[0]          = 1;
    walker->segment.dir[1]          = 0;
    walker->segment.dir[2]          = 0;
    walker->segCrosLength           = 1;
    walker->trackedLength           = 0;
    
    walker->terminationReasonSideA  = TERMINATIONREASON_NOTSET;
    walker->terminationReasonSideB  = TERMINATIONREASON_NOTSET;
    walker->discardingReason        = DISCARDINGREASON_NOTSET;
    walker->failingReason           = FAILREASON_NOTSET;
    walker->successReason           = SUCCESSREASON_NOTSET;
    
    for (int n=0; n<ruleCnt; n++) {
        walker->entry_status.push_back(notEnteredYet);
        walker->isDone.push_back(false);
    }
    
    return walker;
}

NIBR::Walker* NIBR::Pathway::createWalker(std::vector<Point>* streamline, int ind) {
    NIBR::Walker* walker = createWalker(streamline);
    walker->ind          = ind;   
    return walker;
}

void NIBR::Pathway::softReset(NIBR::Walker* walker) {

    walker->side                    = either;
    walker->sideAorder              = 0;
    walker->sideBorder              = 0;
    walker->action                  = CONTINUE;
    walker->trackedLength           = 0;
    
    walker->terminationReasonSideA  = TERMINATIONREASON_NOTSET;
    walker->terminationReasonSideB  = TERMINATIONREASON_NOTSET;
    walker->discardingReason        = DISCARDINGREASON_NOTSET;
    walker->failingReason           = FAILREASON_NOTSET;
    walker->successReason           = SUCCESSREASON_NOTSET;
    
    for (int n=0; n<ruleCnt; n++) {
        walker->entry_status[n] = notEnteredYet;
        walker->isDone[n]       = false;
    }

}

