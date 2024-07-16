#include "surface.h"
#include "math/sphere.h"
#include "findSegmentTriangleIntersection.h"
#include "surface2imageMapper.h"
#include <cmath>

// isPointInside_basedOnWindingNumber(p) uses libigl for computing winding number and determining whether p is inside the mesh or not.
// isPointInside(p) does an initial check using the mask image - only if necessary (p is on the boundary) then the winding number is computed.
// When testing large numbers of points, using the mask and grid noticable increases the computation speed.

void NIBR::Surface::enablePointCheck(float gridRes, bool interpretAs2D) {

    disp(MSG_DEBUG,"enablePointCheck()");

    if (enabledPointCheck) return;

    pointCheckGridRes = gridRes;
    mapSurface2Image(this,&maskAndBoundary,pointCheckGridRes,NULL,&grid,MASK_WITH_BOUNDARY);
    maskAndBoundary.setInterpolationMethod(NEAREST);

    if (interpretAs2D) {
        for (int n = 0; n < maskAndBoundary.numel; n++) {
            if (maskAndBoundary.data[n] == INSIDE) 
                maskAndBoundary.data[n] = OUTSIDE;
        }
    }

    enabledPointCheck = true;
    disp(MSG_DEBUG,"Done enablePointCheck()");

}

bool NIBR::Surface::isPointInside_basedOnWindingNumber(float* p) {

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    return (igl::fast_winding_number(fwn_bvh,2,point) > 0.5);
}

// Returns true is a point is inside a surface mesh
// For cases where the point is on exactly on the surface, this function will return false (in most cases).
// For such possible edge cases, it might be useful to get squaredDistToPoint(float *p) to check the distance from the surface.
bool NIBR::Surface::isPointInside(float* p) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return false;
    }

    float ijk[3];
    maskAndBoundary.to_ijk(p,ijk);
    int64_t A[3] = {int64_t(std::round(ijk[0])), int64_t(std::round(ijk[1])), int64_t(std::round(ijk[2]))};

    int8_t vox = (maskAndBoundary)(A[0],A[1],A[2]);

    if (vox == OUTSIDE) return false;
    if (vox == INSIDE)  return true;

    // Only check the closed component
    // Points can't be inside open surfaces
    return (compClosedAndOpen[0].nv > 0) ? this->compClosedAndOpen[0].isPointInside_basedOnWindingNumber(p) : false;

}
