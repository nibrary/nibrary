#pragma once

#include "math/core.h"

namespace NIBR 
{

    // Returns number of spherical harmonics coefficients for given order and whether odd coefficients are ignored
    int getNumberOfSHCoeffs(int order, bool ignoreOddCoeffs);

    // Returns spherical harmonics order and whether odd coefficients are ignored
    std::tuple<int, bool> getSHOrderFromNumberOfCoeffs(int numberOfCoefficients);

    // Input:
    //  inpCoords: Nx3 list of coordinates on a unit sphere.
    //  order: spherical harmonics order
    // Output:
    //  B: NxM basis coefficients
    void  SH_basis(std::vector<std::vector<float>>& B, std::vector<Point3D>& inpCoords, int order, bool ignoreOddCoeffs);

    void computeLegendrePolynomials(double *plm, double x, int order, bool ignoreOddCoeffs); // x is cos(theta)

}




