#pragma once

#include "dMRI/tractography/io/tractogramReader.h"

namespace NIBR 
{

// Returns the length of the segment that is within the sphere, which perfectly fits the voxel
double segmentSphereIntersectionLength(NIBR::Segment& seg, float* sphere);

template<typename T>
double segmentSphereIntersectionLength(NIBR::Image<T>* img, int* gridPos, NIBR::Segment& seg) 
{

    float sphere[4]     = {0,0,0,img->pixDims[0]*0.5f};
    float sphere_ijk[3] = {float(gridPos[0]),float(gridPos[1]),float(gridPos[2])};

    img->to_xyz(sphere_ijk,sphere);

    return segmentSphereIntersectionLength(seg, &sphere[0]);

}

}