#pragma once

#include "base/nibr.h"
#include "core.h"

namespace NIBR 
{

    namespace SF {

        std::vector<Point3D>&                               getSFCoords();
        std::vector<std::vector<std::tuple<int,float> >>&   getSFNeighbors();

        void    init(bool _sfIsEven, int _sfDim);
        void    init(std::vector<Point3D>& coordinates, bool _sfIsEven);

        void    clean();
        int64_t coordinate2index(float* coord);
        void    coordinateNeighbors(std::vector<std::tuple<int,float> >& neighbors, float* coord, float distThresh);

    }

}