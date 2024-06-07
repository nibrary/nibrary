#include "surface.h"
#include "findSegmentTriangleIntersection.h"
#include <cfloat>
#include <tuple>

// Returns <isBegInside,isEndInside,intersectingFaceInd,distance>
// distance=NAN and intersectionFaceInd= INT_MIN, if there is no intersection
std::tuple<bool,bool,int,float> NIBR::Surface::intersect(LineSegment* seg) 
{  

    // Compute ijk of segment.beg and segment.end
    float p0[3];
    maskAndBoundary.to_ijk(seg->beg,p0);
    int A[3] = {int(std::round(p0[0])), int(std::round(p0[1])), int(std::round(p0[2]))};
    // disp(MSG_DEBUG,"A: [%d,%d,%d]",A[0],A[1],A[2]);

    float p1[3];
    maskAndBoundary.to_ijk(seg->end,p1);
    int B[3] = {int(std::round(p1[0])), int(std::round(p1[1])), int(std::round(p1[2]))};
    // disp(MSG_DEBUG,"B: [%d,%d,%d]",B[0],B[1],B[2]);

    float dist        = NAN;
    float minDist     = FLT_MAX;
    int   intFaceInd  = INT_MIN;

    // Returns <doesIntersect,towardsOutside,faceInd,dist>
    auto checkFaceDistSign=[&](int i, int j, int k)->std::tuple<bool,bool,int,float> {

        dist    = NAN;
        minDist = FLT_MAX;

        bool doesIntersect  = false;
        bool towardsOutside = false;

        // disp(MSG_DEBUG,"Checking faces");
        for (auto faceInd : grid[i][j][k]) {
            
            // disp(MSG_DEBUG,"faceInd: %d", faceInd);

            findSegmentTriangleIntersection(this, faceInd, seg->beg, seg->dir, seg->len, NULL, &dist);

            // disp(MSG_DEBUG,"faceInd: %d, dist: %f", faceInd, dist);

            if ( (!isnan(dist)) && (dist>=0) && (dist<minDist) ) { // equality is at the mesh surface

                minDist       = dist;
                intFaceInd    = faceInd;
                doesIntersect = true;
                
                if (dot(&seg->dir[0],this->normalsOfFaces[faceInd])>0)
                    towardsOutside = true;
                else
                    towardsOutside = false;
            }

        }
        // disp(MSG_DEBUG,"Face check done");

        dist = (minDist==FLT_MAX) ? NAN : minDist;

        // if (doesIntersect) {
        //     if (towardsOutside)
        //         disp(MSG_DEBUG,"Outward intersection at faceInd: %d, dist: %.2f",intFaceInd,dist);
        //     else
        //         disp(MSG_DEBUG,"Inward intersection at faceInd: %d, dist: %.2f",intFaceInd,dist);
        // }

        return std::make_tuple(doesIntersect,towardsOutside,intFaceInd,dist);
    };

    // Check segment intersection. Operates on voxel A.
    // Returns <isBegInside,isEndInside,intFaceInd,dist>
    auto isInside=[&]()->std::tuple<bool,bool,int,float> {

        float val = maskAndBoundary(A[0],A[1],A[2]);

        if (val==OUTSIDE) {

            // Voxel is OUTSIDE
            disp(MSG_DEBUG, "Outside.");
            return std::make_tuple(false,false,INT_MIN,NAN);

        } else if (val==INSIDE) {

            // Voxel is INSIDE
            disp(MSG_DEBUG, "Inside.");
            return std::make_tuple(true,true,INT_MIN,NAN);

        } else {

            // Voxel is on the BOUNDARY
            disp(MSG_DEBUG, "Boundary.");
            auto [doesIntersect,towardsOutside,intFaceInd,dist] = checkFaceDistSign(A[0],A[1],A[2]);

            // There is intersection
            if (doesIntersect) {
                disp(MSG_DEBUG, "There is intersection.");
                if (this->isClosed()) {
                    // Segment goes outside the region
                    if (towardsOutside) {
                        disp(MSG_DEBUG, "in -> out");
                        return std::make_tuple(true,false,intFaceInd,dist);
                    } else {
                        disp(MSG_DEBUG, "out -> in");
                        return std::make_tuple(false,true,intFaceInd,dist);
                    }
                    
                } else {
                    return std::make_tuple(false,false,intFaceInd,dist);
                }
            }

            // No intersection
            disp(MSG_DEBUG, "There is no intersection.");
            bool endIsInside = this->isPointInside(seg->end);
            return std::make_tuple(endIsInside,endIsInside,INT_MAX,NAN);

        }

    };

    

    // Segment does not leave the voxel
    if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
        return isInside();
    }



    // Segment leaves the voxel
    
    // Find length and direction of segment in grid space
    double pt0[3], pt1[3];
    for (int m=0;m<3;m++) {
        pt0[m] = p0[m];
        pt1[m] = p1[m];
    }

    double dir[3], length, t;
    vec3sub(dir,pt1,pt0);
    length = norm(dir);
    vec3scale(dir,1.0/length);

    std::tuple<bool,bool,int,float> interCheck;

    while (length>0.0) {

        interCheck = isInside();

        if (!isnan(std::get<3>(interCheck))) {
            return interCheck;
        }

        if (!rayTraceVoxel(A,pt0,dir,t)) 
            t = 0.0; // so t is not NAN

        if (t>length) break;

        t += EPS4; // so that p0 enters the next voxel

        for (int m=0;m<3;m++) {
            pt0[m] += t*dir[m];
            A[m]    = std::round(pt0[m]);
        }

        length -= t;

    }

    // End of segment reached
    return interCheck;

}
