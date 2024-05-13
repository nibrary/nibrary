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

float NIBR::Pathway::getStreamlineLength(Walker* walker) 
{

    if(walker->begInd > walker->endInd) std::swap(walker->begInd,walker->endInd);

    int   intBeg   = (int)walker->begInd;
    float fractBeg = walker->begInd-intBeg;

    int   intEnd   = (int)walker->endInd;
    float fractEnd = walker->endInd-intEnd;

    float  length = 0;
    float* beg;
    float* end;

    // Handle first segment fraction if exists
    if (fractBeg!=0) {
        beg    = new float[3];
        beg[0] = walker->streamline->at(intBeg).x * (1 - fractBeg) + walker->streamline->at(intBeg+1).x * (fractBeg);
        beg[1] = walker->streamline->at(intBeg).y * (1 - fractBeg) + walker->streamline->at(intBeg+1).y * (fractBeg);
        beg[2] = walker->streamline->at(intBeg).z * (1 - fractBeg) + walker->streamline->at(intBeg+1).z * (fractBeg);

        end    = new float[3];
        end[0] = walker->streamline->at(intBeg+1).x;
        end[1] = walker->streamline->at(intBeg+1).y;
        end[2] = walker->streamline->at(intBeg+1).z;

        length += dist(beg,end);
        delete[] beg;
        delete[] end;
    }

    // Handle full segments
    for (int j=intBeg; j<intEnd; j++) {
        beg     = &walker->streamline->at(j).x;
        end     = &walker->streamline->at(j+1).x;
        length += dist(beg,end);
    }

    // Handle last segment fraction if exists
    if (fractEnd!=0) {

        beg    = new float[3];
        beg[0] = walker->streamline->at(intEnd).x * (1 - fractEnd) + walker->streamline->at(intEnd+1).x * (fractEnd);
        beg[1] = walker->streamline->at(intEnd).y * (1 - fractEnd) + walker->streamline->at(intEnd+1).y * (fractEnd);
        beg[2] = walker->streamline->at(intEnd).z * (1 - fractEnd) + walker->streamline->at(intEnd+1).z * (fractEnd);

        end    = new float[3];
        end[0] = walker->streamline->at(intEnd+1).x;
        end[1] = walker->streamline->at(intEnd+1).y;
        end[2] = walker->streamline->at(intEnd+1).z;

        length += dist(beg,end);
        delete[] beg;
        delete[] end;
    }

    return length;

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
