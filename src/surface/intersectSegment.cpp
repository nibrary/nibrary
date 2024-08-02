#include "surface.h"
#include "findSegmentTriangleIntersection.h"
#include <cfloat>
#include <tuple>

// Checks if a segment intersects the surface or not. 
// <isSegBegInside,isSegEndInside,distFromSegBegToMesh,intersectingFaceIndex,intersectionIsInsideToOutside>
// distFromBegToMesh is NAN if segment is not intersecting the mesh.
std::tuple<bool,bool,double,int,bool> NIBR::Surface::intersectSegment(LineSegment* seg) 
{  

    // Compute ijk of segment.beg and segment.end
    double p0[3];
    maskAndBoundary.to_ijk(seg->beg,p0);
    int A[3] = {int(std::round(p0[0])), int(std::round(p0[1])), int(std::round(p0[2]))};
    // disp(MSG_DEBUG,"A: [%d,%d,%d]",A[0],A[1],A[2]);

    double p1[3];
    maskAndBoundary.to_ijk(seg->end,p1);
    int B[3] = {int(std::round(p1[0])), int(std::round(p1[1])), int(std::round(p1[2]))};
    // disp(MSG_DEBUG,"B: [%d,%d,%d]",B[0],B[1],B[2]);

    double  dist        = NAN;
    double  minDist     = DBL_MAX;
    int     intFaceInd  = INT_MIN;

    std::set<int> facesDone;

    // Returns <doesIntersect,towardsOutside,faceInd,dist>
    auto checkFaceDistSign=[&](int i, int j, int k)->std::tuple<double,int,bool> {

        dist    = NAN;
        minDist = DBL_MAX;

        bool doesIntersect  = false;
        bool towardsOutside = false;

        // disp(MSG_DEBUG,"Checking faces");
        for (auto faceInd : grid[i][j][k]) {

            if (facesDone.find(faceInd) != facesDone.end()) 
                continue;
            
            facesDone.insert(faceInd);

            // disp(MSG_DEBUG,"faceInd: %d", faceInd);
            
            findSegmentTriangleIntersection(this, faceInd, seg->beg, seg->end, NULL, &dist);

            // disp(MSG_DEBUG,"faceInd: %d, dist: %f", faceInd, dist);

            if ( (!isnan(dist)) && (dist<minDist) ) {
                
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

        dist = (minDist==DBL_MAX) ? NAN : minDist;

        if (doesIntersect) {
            return std::make_tuple(dist,intFaceInd,towardsOutside);
        } else {
            return std::make_tuple(NAN,INT_MIN,false);
        }

    };

    // Check segment intersection. Operates on voxel A.
    // Returns <isBegInside,isEndInside,dist,intFaceInd,towardsOutside>
    bool computedBegAndEnd = false;
    bool begIsInside       = false;
    bool endIsInside       = false;

    auto isInside=[&]()->std::tuple<bool,bool,double,int,bool> {

        int8_t val = maskAndBoundary(A[0],A[1],A[2]);

        if (val==OUTSIDE) {

            // Voxel is OUTSIDE
            disp(MSG_DEBUG, "Outside.");
            return std::make_tuple(false,false,NAN,INT_MIN,false);

        } else if (val==INSIDE) {

            // Voxel is INSIDE
            disp(MSG_DEBUG, "Inside.");
            return std::make_tuple(true,true,NAN,INT_MIN,false);

        } else {

            // Voxel is on the BOUNDARY
            disp(MSG_DEBUG, "Boundary.");
            if (!computedBegAndEnd) {
                begIsInside = this->isPointInside(seg->beg);
                endIsInside = this->isPointInside(seg->end);
                computedBegAndEnd = true;
            }

            auto [dist,faceId,towardOutside] = checkFaceDistSign(A[0],A[1],A[2]); 
            return std::make_tuple(begIsInside,endIsInside,dist,faceId,towardOutside);

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

    std::tuple<bool,bool,double,int,bool> interCheck(false,false,NAN,INT_MIN,false);

    int voxCnt = 0;
    while (length>0.0) {

        disp(MSG_DEBUG, "Doing vox: %d", voxCnt++);

        interCheck = isInside();

        if (!isnan(std::get<2>(interCheck))) {
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
