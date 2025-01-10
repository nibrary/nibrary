#include "pathway.h"
#include "pathwayRule.h"
#include <cmath>
#include <ostream>
#include <tuple>

using namespace NIBR;

// Computes the segment.len and segment.dir for the current walker segment
void NIBR::Pathway::prepSegment(NIBR::Walker* walker) 
{
    vec3sub(walker->segment.dir,walker->segment.end,walker->segment.beg);
    walker->segment.len = norm(walker->segment.dir);
    normalize(walker->segment.dir);
}

void NIBR::Pathway::flipSide(NIBR::Walker* walker) 
{

    for (int n=0; n<ruleCnt; n++) {
        walker->entry_status[n] = notEnteredYet;
    }

	if (walker->side == side_A)
		walker->side = side_B;
	else if (walker->side == side_B)
		walker->side = side_A;

    checkSeed(walker,walker->side);

}

// To speed things up, quickly discard for discard_if_ends_inside, discard_if_enters and !req_end_inside
bool NIBR::Pathway::checkEndsInside(float* p, NIBR::Walker* walker) {

    for (int n=0; n<ruleCnt; n++) {
        
        if ( (prules[n].side == either) || (walker->side == either) || (walker->side == prules[n].side) ) {

            if ( (prules[n].type == discard_if_ends_inside) && (isPointInsideRule(p,n)) ) {
                walker->action           = DISCARD;
                walker->discardingReason = ENDED_INSIDE_DISCARD_ROI;
                return true;
            }

            if ( (prules[n].type == discard_if_enters)     && (isPointInsideRule(p,n)) ) {
                walker->action           = DISCARD;
                walker->discardingReason = DISCARD_REGION_REACHED;
                return true;
            }

            if ( (prules[n].type == req_end_inside) && (!isPointInsideRule(p,n)) ) {
                walker->action           = DISCARD;
                walker->discardingReason = REQUIRED_ROI_NOT_MET;
                return true;
            }
            
        }

    }
    return false;
}
