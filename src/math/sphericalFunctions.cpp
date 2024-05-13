#include "base/config.h"
#include "sphericalFunctions.h"
#include <algorithm>

namespace NIBR {
    namespace SF {
        std::vector<std::vector<float>>                      sfCoords;
        std::vector<std::vector<std::tuple<int,float>>>      sfNeighbors;
        bool                                                 sfIsEven  = false;
        bool                                                 sfHasOwnCoordinates = false;
    }
}

std::vector<std::vector<float>>&                     NIBR::SF::getSFCoords()    {return NIBR::SF::sfCoords;}
std::vector<std::vector<std::tuple<int,float> >>&    NIBR::SF::getSFNeighbors() {return NIBR::SF::sfNeighbors;}

using namespace NIBR::SF;

namespace { // private members, only SF can access
    int     sfDim     = 0;
    int*    sfInds    = NULL;
    float   sfRadius  = 0;
    float   sfShift   = 0;
}

void NIBR::SF::init(bool _sfIsEven, int _sfDim)
{
    NIBR::SF::clean();
    
    sfHasOwnCoordinates = false;
    sfIsEven            = _sfIsEven;
    sfDim               = 0;

    if (_sfDim <= 0)
        sfDim  = sfIsEven ? 13 : 11; // should be an odd number
    else if (_sfDim%2 == 0) {
        sfDim  = sfIsEven ? 13 : 11; // should be an odd number
        disp(MSG_WARN,"sfDim cannot be an even number. Switching to default value: %d",sfDim);
    } else
        sfDim  = _sfDim;

    sfRadius   = (float(sfDim)-1)/2.0 - 0.5;
    sfShift    = sfRadius + 0.5;

    float R = (float(sfDim)-1)/2.0;
    float zs;
    

    if (sfIsEven) {
        sfInds      = new int[sfDim*sfDim*((sfDim/2)+1)];
        zs          = 0;
    } else {
        sfInds      = new int[sfDim*sfDim*sfDim];
        zs          = -R;
    }
        
    int ind = 0;
    for (float x=-R; x<=R; x++)
        for (float y=-R; y<=R; y++)
            for (float z=zs; z<=R; z++) {
                float dist = std::sqrt(x*x+y*y+z*z);
                if (std::abs(dist-sfRadius)<(std::sqrt(3)/2)) {
                    sfInds[size_t((x+R)+((y+R)+(z-zs)*sfDim)*sfDim)] = ind++;                  
                    float p[3] = {x,y,z}; 
                    normalize(p);
                    std::vector<float> vertex{p[0],p[1],p[2]};
                    sfCoords.push_back(vertex);
                }
                else
                    sfInds[size_t((x+R)+((y+R)+(z-zs)*sfDim)*sfDim)] = -1;
                
            }

    sfNeighbors.resize(sfCoords.size());
    for (size_t v=0; v<sfCoords.size(); v++) {
        for (size_t u=0; u<sfCoords.size(); u++)
            sfNeighbors[v].push_back(std::make_tuple(u,distS2(sfCoords[v],sfCoords[u])));
        std::sort(sfNeighbors[v].begin(), sfNeighbors[v].end(), [](auto a, auto b) {return std::get<1>(a) < std::get<1>(b);} );
    }

}

void NIBR::SF::init(std::vector<std::vector<float>>& coordinates, bool _sfIsEven)
{
    NIBR::SF::clean();
    sfHasOwnCoordinates = true;
    sfIsEven            = _sfIsEven;
    sfCoords            = coordinates;

    sfNeighbors.resize(sfCoords.size());

    float dist = 0;
    for (size_t v=0; v<sfCoords.size(); v++) {

        for (size_t u=0; u<sfCoords.size(); u++) {
            if (sfIsEven) {
                std::vector<float> reverseDir = {-sfCoords[u][0],-sfCoords[u][1],-sfCoords[u][2]};
                dist = std::min(distS2(sfCoords[v],sfCoords[u]), distS2(sfCoords[v],reverseDir));
            } else {
                dist = distS2(sfCoords[v],sfCoords[u]);
            }
            sfNeighbors[v].push_back(std::make_tuple(u,dist));
        }

        std::sort(sfNeighbors[v].begin(), sfNeighbors[v].end(), [](auto a, auto b) {return std::get<1>(a) < std::get<1>(b);} );
    }

}


void NIBR::SF::coordinateNeighbors(std::vector<std::tuple<int,float>>& neighbors, float* coord, float distThresh) {

    std::vector<std::tuple<int,float>>().swap(neighbors);

    for (auto nei : sfNeighbors[NIBR::SF::coordinate2index(coord)])
        if (std::get<1>(nei) < distThresh)
            neighbors.push_back(nei);

}

void NIBR::SF::clean()
{
    if (sfInds!=NULL) {
        delete[] sfInds;
        sfInds = NULL;
    }
    std::vector<std::vector<float>>().swap(sfCoords);
    std::vector<std::vector<std::tuple<int,float>>>().swap(sfNeighbors);
}

int64_t NIBR::SF::coordinate2index(float* coord) {
    
    if (sfHasOwnCoordinates) {

        float minDist = FLT_MAX;
        float dist    = FLT_MAX;

        int64_t n = 0;

        for (int i = 0; i < int(sfCoords.size()); i++) {

            if (sfIsEven) {
                std::vector<float> reverseDir = {-sfCoords[i][0],-sfCoords[i][1],-sfCoords[i][2]};
                dist = std::min(distS2(coord,sfCoords[i]),distS2(coord,reverseDir));
            } else {
                dist = distS2(coord,sfCoords[i]);
            }

            if (dist < minDist) {
                minDist = dist;
                n = i;
            }

        }

        return n;

    } else {

        int x,y,z;
        
        if (sfIsEven) { 
            if (coord[2]<0) {
                x = std::nearbyint(-coord[0]*sfRadius) + sfShift;
                y = std::nearbyint(-coord[1]*sfRadius) + sfShift;
                z = std::nearbyint(-coord[2]*sfRadius);
            } else {
                x = std::nearbyint( coord[0]*sfRadius) + sfShift;
                y = std::nearbyint( coord[1]*sfRadius) + sfShift;
                z = std::nearbyint( coord[2]*sfRadius);
            }
        } else {
            x = std::nearbyint(coord[0]*sfRadius) + sfShift;
            y = std::nearbyint(coord[1]*sfRadius) + sfShift;
            z = std::nearbyint(coord[2]*sfRadius) + sfShift;
        }

        return sfInds[x+(y+z*sfDim)*sfDim];

    }

}
