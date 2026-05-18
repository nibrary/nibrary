#pragma once

#include "shapePrimitives.h"

Surface makePointer_AntNeuro();
Surface make3DAxis(float length, float radius);
Surface make3DAxis(float* origin, float* xAxis, float* yAxis, float* zAxis, float length, float radius, float coneLength);
Surface makeSphereWithRibbons(float* origin, float sphereRadius, float ribbonThickness);
