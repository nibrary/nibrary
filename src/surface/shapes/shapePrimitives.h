#pragma once

#include "surface/surface.h"
#include "surface/surface_operators.h"
#include "dMRI/tractography/tractogram.h"

namespace NIBR
{
    // Makes a sphere Surface object with given center and radius
    Surface makeSphere(float radius, int radialSegments);

    // Makes a sphere Surface object with given center and radius
    Surface makeSphere(float* center, float radius, int radialSegments);

    // Makes a cylinder Surface object between p1 and p2 with given radius and radial segments
    Surface makeCylinder(float* p1, float* p2, float radius, int radialSegments);

    // Makes a cone Surface object with given tip, base center, base radius and radial segments
    Surface makeCone(float* tip, float* baseCenter, float baseRadius, int radialSegments);

    // Makes a tube Surface object along the given streamline with specified radius and radial segments
    Surface streamline2tube(const Streamline& streamline, float radius, int radialSegments = 8, bool sphericalCaps = true, int threadId = 0);
    Surface tractogram2tube(const Tractogram& tractogram, float radius, int radialSegments = 8, bool sphericalCaps = true);
}