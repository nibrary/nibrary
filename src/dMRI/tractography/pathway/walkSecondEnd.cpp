#include "pathway.h"
#include "pathwayRule.h"

using namespace NIBR;

// Starting from the seed, walk towards the beginning of the streamline and get begInd
void NIBR::Pathway::walkSecondEnd(NIBR::Walker* walker) 
{

    int n = walker->seedInd;

    if (walker->action == STOP) {
        walker->begInd = n;
        return;
    }

    // If seedInd is the first index, then only prepare the segment
    if (n != 0) {
        for (n=walker->seedInd; n>0; n--) {
            if (checkWalker(walker,n,n-1) != CONTINUE)
                break;
        }
    } else {
        walker->segment.beg   = walker->streamline->at(n).data();
        walker->segment.end   = walker->streamline->at(n).data();
        walker->segment.len   = 0;
        walker->segment.dir[0]= 0;
        walker->segment.dir[1]= 0;
        walker->segment.dir[2]= 0;
    }

    // Set the begInd if STOP is reached
    if (walker->action == STOP) {

        double segStopFrac = walker->segStopLength / walker->segment.len;

        walker->begInd = n - segStopFrac;

        // If the seed point is inserted we need to shift the begInd accordingly
        if ((walker->seedInserted) && (n==walker->seedInd)) {

            float prevLen = dist(walker->streamline->at(n),walker->streamline->at(n-1));
            float nextLen = dist(walker->streamline->at(n),walker->streamline->at(n+1));
            float crsLen  = prevLen * segStopFrac;
            float corrLen = (prevLen - crsLen) / (prevLen + nextLen);

            walker->begInd = n - 1.0 + corrLen;

        }

    } else {
        walker->begInd = 0;
    }

}