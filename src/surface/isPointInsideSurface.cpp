#include "surface.h"
#include "math/sphere.h"
#include "surface2imageMapper.h"
#include <cmath>

// isPointInside_basedOnWindingNumber(p) uses libigl for computing winding number and determining whether p is inside the mesh or not.
// isPointInside(p) does an initial check using the mask image - only if necessary (p is on the boundary) then the winding number is computed.
// When testing large numbers of points, using the mask and grid noticable increases the computation speed.

void NIBR::Surface::enablePointCheck(float gridRes) {

    disp(MSG_DEBUG,"enablePointCheck()");

    if (enabledPointCheck) return;

    prepIglAABBTree();
    calcNormalsOfFaces();

    pointCheckGridRes = gridRes;
    mapSurface2Image(this,&maskAndBoundary,pointCheckGridRes,NULL,&grid,MASK_WITH_BOUNDARY);
    maskAndBoundary.setInterpolationMethod(NEAREST);

    if (compClosedAndOpen[0].nv > 0) compClosedAndOpen[0].prepIglAABBTree();
    if (compClosedAndOpen[1].nv > 0) compClosedAndOpen[1].prepIglAABBTree();

    enabledPointCheck = true;
    disp(MSG_DEBUG,"Done enablePointCheck()");

}

bool NIBR::Surface::isPointInside_basedOnWindingNumber(float* p) {

    Eigen::MatrixXd point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    return (igl::fast_winding_number(fwn_bvh,2,point) > 0.5);
}

// OUTSIDE  -> dist = DBL_MIN
// INSIDE   -> dist = DBL_MAX
// BOUNDARY -> dist = dist to surface mesh
bool NIBR::Surface::isPointInside(float* p, double& dist, int& faceInd, float* closestPoint) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return false;
    }

    float ijk[3];
    maskAndBoundary.to_ijk(p,ijk);
    int64_t A[3] = {int64_t(std::round(ijk[0])), int64_t(std::round(ijk[1])), int64_t(std::round(ijk[2]))};

    int8_t vox = (maskAndBoundary)(A[0],A[1],A[2]);

    auto isOnBoundary = [&]()->bool {
        double v[3];
        dist = squaredDistToPoint(p, faceInd, closestPoint);
        vec3sub(&v[0],vertices[faces[faceInd][0]],p);
        dist = (dot(normalsOfFaces[faceInd],&v[0]) > 0.0) ? std::sqrt(dist) : -std::sqrt(dist);

        // disp(MSG_DEBUG, "dist: %.12f",dist);

        if (dist < 0.0)             return false;
        if (dist <= SURFTHICKNESS)  return true;
        
        dist = -dist;
        return false;
        
    };

    if (vox == OUTSIDE) {dist = DBL_MIN; return false;}

    if (interpretAs2D == false) {

        if (vox == INSIDE)  {dist = DBL_MAX; return true;}

        // The point is close to the boundary

        // If the surface has closed components
        if (compClosedAndOpen[0].nv > 0) {
            if (this->compClosedAndOpen[0].isPointInside_basedOnWindingNumber(p)) {
                dist = DBL_MAX;
                return true;
            }
        }

        // If the surface has open components,
        // the point is inside if it is within boundary
        if (compClosedAndOpen[1].nv > 0) {
            return isOnBoundary();
        }

        dist = DBL_MIN;
        return false;

    }

    return isOnBoundary();

}

// Returns true if a point is inside a surface mesh
bool NIBR::Surface::isPointInside(float* p) {
    double dist;
    int faceInd;
    float closestPoint[3];
    return isPointInside(p, dist, faceInd, &closestPoint[0]);
}

// Returns the location of the point:
// OUTSIDE   0
// INSIDE    1
// BOUNDARY  2
short NIBR::Surface::whereIsPoint(float* p) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return false;
    }

    float ijk[3];
    maskAndBoundary.to_ijk(p,ijk);
    int64_t A[3] = {int64_t(std::round(ijk[0])), int64_t(std::round(ijk[1])), int64_t(std::round(ijk[2]))};

    int8_t vox = (maskAndBoundary)(A[0],A[1],A[2]);

    auto isOnBoundary = [&]()->bool {
        int faceInd;
        double v[3];
        double dist = squaredDistToPoint(p, faceInd);
        vec3sub(&v[0],vertices[faces[faceInd][0]],p);
        dist = (dot(normalsOfFaces[faceInd],&v[0]) > 0.0) ? std::sqrt(dist) : -std::sqrt(dist);
        return ((dist <= SURFTHICKNESS) && (dist > 0.0));
    };

    if (interpretAs2D == false) {

        if (vox == OUTSIDE) return OUTSIDE;
        if (vox == INSIDE)  return INSIDE;

        // Check if the point is on the boundary
        if (isOnBoundary()) return BOUNDARY;

        // If the surface has closed components
        if (compClosedAndOpen[0].nv > 0) {
            return compClosedAndOpen[0].isPointInside_basedOnWindingNumber(p);
        }

        return OUTSIDE;
    }

    return (isOnBoundary()) ? BOUNDARY : OUTSIDE;

}
