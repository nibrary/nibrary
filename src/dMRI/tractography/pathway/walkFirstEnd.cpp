#include "pathway.h"
#include "pathwayRule.h"

using namespace NIBR;

// Starting from the seed, walk towards the end of the streamline and get endInd
void NIBR::Pathway::walkFirstEnd(NIBR::Walker* walker) {

    int n = walker->seedInd;

    if (walker->action == STOP) {
        walker->endInd = n;
        return;
    }

    // If seedInd is the last index, then only prepare the segment
    if (n != int(walker->streamline->size()-1)) {

        // Starting from the seed, walk towards the end of the streamline
        for (n=walker->seedInd; n<int(walker->streamline->size()-1); n++) {
            if (checkWalker(walker,n,n+1) != CONTINUE)
                break;
        }

    } else {
        walker->segment.beg   = &(walker->streamline->at(n).x);
        walker->segment.end   = &(walker->streamline->at(n).x);
        walker->segment.len   = 0;
        walker->segment.dir[0]= 0;
        walker->segment.dir[1]= 0;
        walker->segment.dir[2]= 0;
    }

    // Set the endInd if STOP is reached
    if (walker->action == STOP) {

        walker->endInd = n + walker->segCrosLength;
        
        // If the seed point is inserted we need to shift the endInd accordingly
        if (walker->seedInserted) {

            if (n == walker->seedInd) {

                float prevLen  = dist(walker->streamline->at(n),walker->streamline->at(n-1));
                float nextLen  = dist(walker->streamline->at(n),walker->streamline->at(n+1));
                float crsLen   = nextLen * walker->segCrosLength;
                float corrLen  = (prevLen + crsLen) / (prevLen + nextLen);

                walker->endInd = n - 1.0f + corrLen;

            } else {
                walker->endInd = n - 1.0f + walker->segCrosLength;
            }

        }

    } else {
        walker->endInd = (walker->seedInserted) ? walker->streamline->size()-2 : walker->streamline->size()-1;
    }

    // disp(MSG_DEBUG,"n: %d, walker->segCrosLength]: %.2f", n, walker->segCrosLength);
        
}