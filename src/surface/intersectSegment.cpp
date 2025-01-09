#include "surface.h"
#include "findSegmentTriangleIntersection.h"
#include <cfloat>
#include <tuple>

// Checks if a segment intersects the surface or not. 
// <isSegBegInside,isSegEndInside,distFromSegBegToMesh,intersectingFaceIndex,intersectionIsInsideToOutside,boundaryTransitionDist>
// distFromBegToMesh      ≠ NAN if the segment is intersecting the mesh. Then intersectingFaceIndex and intersectionIsInsideToOutside are valid.
// boundaryTransitionDist ≠ NAN if the segment is transitioning through the boundary without intersection.
// i.e. a segment can transition from the boundary to outside or the other way around, without intersecting the mesh
// this can happen for example when the segment beg is within SURFTHICKNESS and the end is outside on the opposite of the normal
std::tuple<bool,bool,double,int,bool,double> NIBR::Surface::intersectSegment(LineSegment* seg) 
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

    
    // Returns <isBegInside,isEndInside,dist,intFaceInd,towardsOutside,boundaryTransitionDist>
    bool    begIsInside             = false;
    bool    endIsInside             = false;
    double  dist                    = NAN;
    int     intFaceInd              = INT_MIN;
    bool    towardsOutside          = false;
    double  boundaryTransitionDist  = NAN;

    bool    computedBegAndEnd       = false;
    std::set<int> facesDone;
    double  minDist = DBL_MAX;

    // Check segment intersection. Operates on voxel A.
    // Returns <doesIntersect,towardsOutside,faceInd,dist>
    auto checkFaceDistSign=[&](int i, int j, int k)->std::tuple<double,int,bool> {

        dist    = NAN;
        minDist = DBL_MAX;

        bool doesIntersect  = false;

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
                
                if (dot(seg->dir,this->normalsOfFaces[faceInd])>0)
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

    auto isInside=[&]()->std::tuple<bool,bool,double,int,bool,double> {

        int8_t val = maskAndBoundary(A[0],A[1],A[2]);

        if (val==OUTSIDE) {

            // Voxel is OUTSIDE
            // disp(MSG_DEBUG, "Outside.");
            begIsInside = endIsInside = false;
            return std::make_tuple(begIsInside,endIsInside,NAN,INT_MIN,false,NAN);

        } else if (val==INSIDE) {

            // Voxel is INSIDE
            // disp(MSG_DEBUG, "Inside.");
            begIsInside = endIsInside = true;
            return std::make_tuple(begIsInside,endIsInside,NAN,INT_MIN,false,NAN);

        } else {

            // Voxel is on the BOUNDARY
            // disp(MSG_DEBUG, "Boundary.");
            if (!computedBegAndEnd) {
                begIsInside = this->isPointInside(seg->beg);
                endIsInside = this->isPointInside(seg->end);
                computedBegAndEnd = true;
            }

            auto [dist,faceId,towardOutside] = checkFaceDistSign(A[0],A[1],A[2]); 
            return std::make_tuple(begIsInside,endIsInside,dist,faceId,towardOutside,NAN);

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

    // int voxCnt = 0;
    while (length>0.0) {

        // disp(MSG_DEBUG, "Doing vox: %d", voxCnt++);

        isInside();

        if (!isnan(dist)) break;

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

    // Handle boundary transition cases in the case of 2D-interpreted surfaces
    auto exitedBoundary = [&](double x)->bool {
        float p[3];
        p[0] = seg->beg[0] + dir[0] * x;
        p[1] = seg->beg[1] + dir[1] * x;
        p[2] = seg->beg[2] + dir[2] * x;
        return !isPointInside(&p[0]);
    };

    if (interpretAs2D) {

        if (!isnan(dist)) {
            // Segment went from outside to outside crossing the mesh towardsOutside
            // This means it must have first entered somewhere along the way before intersecting the mesh
            if (!begIsInside && !endIsInside && towardsOutside) {
                bool foundPointOutOfBoundary = false;
                for (double d = double(dist)-EPS7; d > 0; d -= EPS7) {
                    if (exitedBoundary(d)) {
                        boundaryTransitionDist  = d+EPS7; // is inside the boundary
                        foundPointOutOfBoundary = true;
                        break;
                    }
                }
                if (!foundPointOutOfBoundary) {
                    disp(MSG_ERROR,"Expected point along the boundary was not found");
                }
                // disp(MSG_DEBUG,"boundaryTransitionDist towards outside: %.12f", boundaryTransitionDist);
            }
        } else {

            // Segment went from boundary to outside without crossing the mesh
            if (begIsInside && !endIsInside) {
                bool foundPointOutOfBoundary = false;
                for (double d = EPS7; d <= length; d += EPS7) {
                    if (exitedBoundary(d)) {
                        boundaryTransitionDist = d; // is outside the boundary
                        foundPointOutOfBoundary = true;
                        break;
                    }
                }
                if (!foundPointOutOfBoundary) {
                    disp(MSG_ERROR,"Expected point along the boundary was not found");
                }
                // disp(MSG_DEBUG,"boundaryTransitionDist towards outside: %.12f", boundaryTransitionDist);
            }

            // Segment went from outside to the boundary without crossing the mesh
            if (!begIsInside && endIsInside) {
                bool foundPointOutOfBoundary = false;
                for (double d = double(length)-EPS7; d > 0; d -= EPS7) { // endIsInside = true, so we go back towards to the beg until p is outside. We then add EPS7, so the output is inside.
                    if (exitedBoundary(d)) {
                        boundaryTransitionDist = d + EPS7; // is inside the boundary
                        foundPointOutOfBoundary = true;
                        break;
                    }
                }
                if (!foundPointOutOfBoundary) {
                    disp(MSG_ERROR,"Expected point along the boundary was not found");
                }
                // disp(MSG_DEBUG,"boundaryTransitionDist towards inside: %.12f", boundaryTransitionDist);
            }

            // If the segment leaves the boundary find the transition point
            if (begIsInside && endIsInside) {
                for (double d = EPS7; d <= SURFTHICKNESS; d += EPS7) {
                    if (exitedBoundary(d)) {
                        boundaryTransitionDist = d; // is outside the boundary
                        break;
                    }
                }
            }

            // If the segment enters the boundary without crossing
            // we will assume that segment never enters the boundary
            // if (!begIsInside && !endIsInside) { ... }

        }



    }

    /*
    if (begIsInside) {disp(MSG_DEBUG,"begIsInside: YES");} else {disp(MSG_DEBUG,"begIsInside: NO");} 
    if (endIsInside) {disp(MSG_DEBUG,"endIsInside: YES");} else {disp(MSG_DEBUG,"endIsInside: NO");} 
    disp(MSG_DEBUG,"dist:       %.12f", dist);
    disp(MSG_DEBUG,"intFaceInd: %d",    intFaceInd);
    if (towardsOutside) {disp(MSG_DEBUG,"towardsOutside: YES");} else {disp(MSG_DEBUG,"towardsOutside: NO");}
    if (boundaryTransitionDist) {disp(MSG_DEBUG,"boundaryTransitionDist: YES");} else {disp(MSG_DEBUG,"boundaryTransitionDist: NO");}
    */

    return std::make_tuple(begIsInside,endIsInside,dist,intFaceInd,towardsOutside,boundaryTransitionDist);

}
