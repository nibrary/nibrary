#include "core.h"

using namespace NIBR;

// A: voxel indices, i.e. center of voxel in image space
// p0: origin of ray
// dir: direction of ray
// 
// Returns true if ray intersects the voxel, and t is the distance from the ray origin to the intersection point.
// Returns false if ray does not intersect the voxel, and t is set to NAN.
//
//
// Note: This function handles the edge cases. 
//    In image space if A=[0,0,0] the -0.5 is outside the voxel
//    For other voxels x-0.5 is inside the voxel, and x+0.5 is outside the voxel, i.e. x+0.5 is inside voxel index x+1, not x

// Handle edge case
inline bool isInsideVoxel(int index, double value) {
    if (index == 0) {
        return value > -0.5 && value < 0.5;
    }
    return value >= -0.5 && value < 0.5;
}


bool NIBR::rayTraceVoxel(int* A, double* p0, double* dir, double &t) {
    double T[3], u, v;

    // Calculate the relative position of p0 to the cube's center A
    for (int m = 0; m < 3; m++) {
        T[m] = p0[m] - A[m];
    }

    for (int axis = 0; axis < 3; axis++) {
        int axis1 = (axis + 1) % 3;
        int axis2 = (axis + 2) % 3;

        if (dir[axis] < 0.0) {
            t = (-0.5 - T[axis]) / dir[axis];
            bool validT = A[axis] > 0 ? t >= 0.0 : t > 0.0;
            if (validT) {
                u = T[axis1] + t * dir[axis1];
                v = T[axis2] + t * dir[axis2];
                if (isInsideVoxel(A[axis1], u) && isInsideVoxel(A[axis2], v)) {
                    return true;
                }
            }
        } else {
            t = (0.5 - T[axis]) / dir[axis];
            if (t > 0.0) {
                u = T[axis1] + t * dir[axis1];
                v = T[axis2] + t * dir[axis2];
                if (isInsideVoxel(A[axis1], u) && isInsideVoxel(A[axis2], v)) {
                    return true;
                }
            }
        }
    }

    t = NAN;
    return false;
}

