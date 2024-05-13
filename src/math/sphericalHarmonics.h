#pragma once

#include "core.h"

namespace NIBR 
{

    typedef enum {
        ORDEROFDIRECTIONS_NOTSET = -1,
        XYZ,XYz,XyZ,Xyz,xYZ,xYz,xyZ,xyz,
        XZY,XZy,XzY,Xzy,xZY,xZy,xzY,xzy,
        YXZ,YXz,YxZ,Yxz,yXZ,yXz,yxZ,yxz,
        YZX,YZx,YzX,Yzx,yZX,yZx,yzX,yzx,
        ZYX,ZYx,ZyX,Zyx,zYX,zYx,zyX,zyx,
        ZXY,ZXy,ZxY,Zxy,zXY,zXy,zxY,zxy
    } OrderOfDirections;

    OrderOfDirections convertOrderOfDirections(std::string ood);

    int getNumberOfSHCoeffs(int order);
    int getSHOrderFromNumberOfCoeffs(int numberOfCoefficients);

    // Input:
    //  inpCoords: Nx3 list of coordinates on a unit sphere.
    //  order: spherical harmonics order
    // Output:
    //  B: NxM basis coefficients
    void  SH_basis(std::vector<std::vector<float>>& B, std::vector<std::vector<float>>& inpCoords, int order);

    // TODO: Below should be a separate class.
    namespace SH 
    {
        int   getNumberOfSHCoeffs();
        void  precompute(int _sphericalHarmonicOrder, OrderOfDirections _orderOfDirs, size_t _num);
        void  clean();
        float toSF(float *sh, float *dir);        
    }

}
