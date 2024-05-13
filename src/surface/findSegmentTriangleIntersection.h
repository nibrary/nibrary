#pragma once

#include "surface.h"
#include "dMRI/tractography/io/tractogramReader.h"

namespace NIBR {

    // Find intersection of a line segment with a triangle using Möller–Trumbore algorithm
    // This functions returns the angle of intersection in degrees [0,90]
    // 0 means there is no intersection.

    // In its most general form. The following takes a surface and index of a face.
    // Given a point (p), direction (dir), length -> the function finds the intersection position on the triangle (pos) and its distance (dist) to p.
    // The function returns 0 if there is no intersection.
    float findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, float* p, float* dir, float length, float* pos, float* dist);
    void  findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, float* p, const float* dir, float* dist);
    float findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, NIBR::Segment& seg);


    // This version, puts and "extent" amount of length to the ends of the segment then checks the intersection
    float findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, NIBR::Segment& seg, float extent);

    float findSignedSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, NIBR::Segment& seg, float* pos);

}