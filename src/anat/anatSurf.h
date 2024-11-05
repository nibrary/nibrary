#pragma once

#include "image/image_marchingCubes.h"
#include "surface/surface_operators.h"
#include "surface/surface2imageMapper.h"

namespace NIBR {

void wmSurf2wmLabel(
        Image<int8_t>& wmLabel, 
        Surface* wmSurf, 
        SurfaceField* wmSurfLabels,
        float sulcalDepth,                  // determines WM_SULCAL_FUNDUS, e.g., 5 mm
        float wallDepth,                    // determines WM_SULCAL_WALL and WM_GYRAL_CROWN, e.g., 1 mm
        float minThicknessOfSuperficialWM,  // superficialWM will be at least this thick, e.g., 1 mm
        float maxThicknessOfSuperficialWM); // superficialWM will be at most this thick, e.g., 5 mm

}