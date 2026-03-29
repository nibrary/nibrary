#include "fod_image.h"
#include "math/sphericalHarmonics_aux.h"

using namespace NIBR;

void NIBR::FOD_Image::setSH() {

    if (isspheresliced==false) {

        auto [order, ignoreOdd] = getSHOrderFromNumberOfCoeffs(imgDims[3]);

        shOrder         = order;
        hasNoOddCoeffs  = ignoreOdd;
        isAsym          = !ignoreOdd;
        
    } else {
        if (!isAsym) {
            shOrder         = 16;
            hasNoOddCoeffs  = true;
        } else {
            shOrder         = 13;
            hasNoOddCoeffs  = false;
        }
    }

    if (!isAsym) {
        disp(MSG_DETAIL,"(symmetric FOD with order %d)", shOrder);
    } else {
        disp(MSG_DETAIL,"(asymmetric FOD with order %d)",shOrder);
    }

    if (sh!=NULL) {
        delete sh;
        sh = NULL;
    }

    sh = new SH(shOrder, hasNoOddCoeffs, orderOfDirections, SHprecomputationResolution);

}

void NIBR::FOD_Image::fillDiscVolSph() {

    // Prep sphere parameters
    // discVolSphDim      = (!isAsym) ? 21 : 15; // For even, 1038 points on half-sphere (2076 points on full-sphere); for odd, 1004 points on full-sphere, i.e. AFOD is less densely sampled
    discVolSphDim      = (!isAsym) ? 13 : 11;
    discVolSphRadius   = (float(discVolSphDim)-1)/2.0 - 0.5;
    discVolSphShift    = discVolSphRadius + 0.5;

    float R = (float(discVolSphDim)-1)/2.0;
    float zs;
    

    if (!isAsym) {
        discVolSphInds      = new int[discVolSphDim*discVolSphDim*((discVolSphDim/2)+1)];
        zs                  = 0;
    } else {
        discVolSphInds      = new int[discVolSphDim*discVolSphDim*discVolSphDim];
        zs                  = -R;
    }
        
    int ind = 0;
    for (float x=-R; x<=R; x++)
        for (float y=-R; y<=R; y++)
            for (float z=zs; z<=R; z++) {
                float dist = std::sqrt(x*x+y*y+z*z);
                if (std::abs(dist-discVolSphRadius)<(std::sqrt(3)/2)) {
                    discVolSphInds[size_t((x+R)+((y+R)+(z-zs)*discVolSphDim)*discVolSphDim)] = ind++;                  
                    std::array<float,3> p = {x,y,z}; 
                    normalize(p);
                    discVolSphCoords.emplace_back(p);
                }
                else
                    discVolSphInds[size_t((x+R)+((y+R)+(z-zs)*discVolSphDim)*discVolSphDim)] = -1;
                
            }

}


int64_t NIBR::FOD_Image::vertexCoord2volInd(float* vertexCoord) {
    
    int x,y,z;
    
    if (!isAsym) { 
        if (vertexCoord[2]<0) {
            x = std::nearbyint(-vertexCoord[0]*discVolSphRadius) + discVolSphShift;
            y = std::nearbyint(-vertexCoord[1]*discVolSphRadius) + discVolSphShift;
            z = std::nearbyint(-vertexCoord[2]*discVolSphRadius);
        } else {
            x = std::nearbyint( vertexCoord[0]*discVolSphRadius) + discVolSphShift;
            y = std::nearbyint( vertexCoord[1]*discVolSphRadius) + discVolSphShift;
            z = std::nearbyint( vertexCoord[2]*discVolSphRadius);
        }
    } else {
        x = std::nearbyint(vertexCoord[0]*discVolSphRadius) + discVolSphShift;
        y = std::nearbyint(vertexCoord[1]*discVolSphRadius) + discVolSphShift;
        z = std::nearbyint(vertexCoord[2]*discVolSphRadius) + discVolSphShift;
    }

    int64_t volInd = discVolSphInds[x+(y+z*discVolSphDim)*discVolSphDim];
    
    return volInd;
}
