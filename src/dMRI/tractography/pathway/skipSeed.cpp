#include "pathway.h"

using namespace NIBR;

bool NIBR::Pathway::skipSeed(NIBR::Walker *walker, bool reverseDir)
{
    if (walker->action == DISCARD) return false;

    if (!reverseDir) {

        // Tracked length needs to be recomputed
        walker->trackedLength = 0;

        // Descending search
        float ind    = std::ceil(walker->endInd);
        float res    = ind - walker->endInd;

        if (res>0) {

            // Last segment is not integer. This needs to be handled separately.
            walker->segment.beg = &(walker->streamline->at(ind).x);
            walker->segment.end = &(walker->streamline->at(ind-1).x);
            prepSegment(walker);

            // Already the last segment is inside the seed region
            if (isSegmentEntering(walker,theOneSeed)) {

                walker->trackedLength = walker->segCrosLength*walker->segment.len;

                if (walker->segCrosLength > res) {
                    // Part of the last segment is inside the seed region, we will cut that part.
                    walker->begInd        = ind - walker->segCrosLength;

                    if (walker->trackedLength < minLength) {
                        walker->action           = DISCARD;
                        walker->discardingReason = TOO_SHORT;
                        return false;
                    }

                    // The streamline is 1 segment long. But it is longer than the minLength threshold. So we still keep it.

                    // In the case of tracking, where only reverseDir = false, we modify the streamline too that directly changes the output.
                    // Note that segment is running in descending direction, so segment.end is in the seed, and segment.beg is outside.
                    if (isTracking) {
                        walker->segment.end[0] = walker->segment.beg[0] + walker->segment.dir[0]*((walker->segCrosLength)*walker->segment.len);
                        walker->segment.end[1] = walker->segment.beg[1] + walker->segment.dir[1]*((walker->segCrosLength)*walker->segment.len);
                        walker->segment.end[2] = walker->segment.beg[2] + walker->segment.dir[2]*((walker->segCrosLength)*walker->segment.len);
                        walker->streamline->erase(walker->streamline->begin(), walker->streamline->begin() + ind - 1);
                    }

                    return true; // KEEP case

                } else {
                    // The whole last segment is inside the seed region
                    walker->action           = DISCARD;
                    walker->discardingReason = SEED_NOT_FOUND;
                    return false;
                }

            }

            ind = ind - 1;

        }

        // We can now check the remaining streamline for the seed region
        for (float desInd=ind; desInd>walker->begInd; desInd--) {

            walker->segment.beg = &(walker->streamline->at(desInd).x);
            walker->segment.end = &(walker->streamline->at(desInd-1).x);
            prepSegment(walker);

            // disp(MSG_DEBUG,"ind: %.2f, beg: %.2f, end: %.2f", desInd, walker->segment.beg[0], walker->segment.end[0]);
            if (isSegmentEntering(walker,theOneSeed)) {

                // disp(MSG_DEBUG,"  shift: %.2f", walker->segCrosLength);

                walker->begInd         = desInd - walker->segCrosLength;
                walker->trackedLength += (walker->segCrosLength)*walker->segment.len;

                if (walker->seedInserted)
                    walker->begInd -= 1;
                
                if (walker->trackedLength < minLength) {
                    walker->action           = DISCARD;
                    walker->discardingReason = TOO_SHORT;
                    return false;
                }

                // In the case of tracking, where only reverseDir = false, we modify the streamline too that directly changes the output.
                // Note that segment is running in descending direction, so segment.end is in the seed, and segment.beg is outside.
                if (isTracking) {
                    walker->segment.end[0] = walker->segment.beg[0] + walker->segment.dir[0]*((walker->segCrosLength)*walker->segment.len);
                    walker->segment.end[1] = walker->segment.beg[1] + walker->segment.dir[1]*((walker->segCrosLength)*walker->segment.len);
                    walker->segment.end[2] = walker->segment.beg[2] + walker->segment.dir[2]*((walker->segCrosLength)*walker->segment.len);
                    walker->streamline->erase(walker->streamline->begin(), walker->streamline->begin() + desInd - 1);
                }

                return true;  // KEEP case
            } else {
                walker->trackedLength += walker->segment.len;
            }

        }

    } else {

        // Tracked length needs to be recomputed
        walker->trackedLength = 0;

        // Ascending search
        float ind    = std::floor(walker->begInd);
        float res    = walker->begInd - ind;

        if (res>0) {

            // Last segment is not integer. This needs to be handled separately.
            walker->segment.beg = &(walker->streamline->at(ind).x);
            walker->segment.end = &(walker->streamline->at(ind+1).x);
            prepSegment(walker);

            // Already the last segment is inside the seed region
            if (isSegmentEntering(walker,theOneSeed)) {

                walker->trackedLength = walker->segCrosLength*walker->segment.len;

                if (walker->segCrosLength > res) {
                    // Part of the last segment is inside the seed region, we will cut that part.
                    walker->begInd        = ind + walker->segCrosLength;

                    if (walker->seedInserted)
                        walker->begInd -= 1;

                    if (walker->trackedLength < minLength) {
                        walker->action           = DISCARD;
                        walker->discardingReason = TOO_SHORT;
                        return false;
                    }

                    // The streamline is 1 segment long. But it is longer than the minLength threshold. So we still keep it.
                    return true; // KEEP case

                } else {
                    // The whole last segment is inside the seed region
                    walker->action           = DISCARD;
                    walker->discardingReason = SEED_NOT_FOUND;
                    return false;
                }

            }

            ind = ind + 1;

        }

        // We can now check the remaining streamline for the seed region
        for (float ascInd=ind; ascInd<walker->endInd; ascInd--) {

            walker->segment.beg = &(walker->streamline->at(ascInd).x);
            walker->segment.end = &(walker->streamline->at(ascInd+1).x);
            prepSegment(walker);
            
            if (isSegmentEntering(walker,theOneSeed)) {
                walker->endInd         = ascInd + walker->segCrosLength;
                walker->trackedLength += (walker->segCrosLength)*walker->segment.len;
                
                if (walker->trackedLength < minLength) {
                    walker->action           = DISCARD;
                    walker->discardingReason = TOO_SHORT;
                    return false;
                }

                return true; // KEEP case
            } else {
                walker->trackedLength += walker->segment.len;
            }
        }

    }

    walker->action           = DISCARD;
    walker->discardingReason = SEED_NOT_FOUND;
    return false;

}