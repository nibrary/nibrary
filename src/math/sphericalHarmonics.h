#pragma once

#include "math/reorient.h"
#include <cstddef>
#include <functional>

namespace NIBR 
{
    
class SH
{
public:
    SH(int _order, bool _noOddCoeffs = false, OrderOfDirections _orderOfDirs = XYZ, size_t _numberOfSamples = 1024);
    ~SH() { clean(); }

    SH(const SH&) = delete;
    SH& operator=(const SH&) = delete;

    int   getCoeffCount() {return coeffCount;}
    float toSF(float *sh, float *dir);

private:

    void  clean();
    void  precompute();

    int  order;
    int  coeffCount;
    bool noOddCoeffs;

    size_t numberOfSamples;
    size_t numberOfSamples_phi;
    size_t numberOfSamples_theta;

    double scalingFactor_phi;
    double scalingFactor_theta;

    float* precomputedPhiComponent   = NULL;
    float* precomputedThetaComponent = NULL;

    OrderOfDirections            orderOfDirs;
    std::function<void(float *)> orderDirections;

};

}




