#pragma once

#include "math/core.h"

namespace NIBR
{
    struct Segment {
        int    streamlineNo;
        float  p[3];
        float  dir[3];
        float  length;
        void*  data;
    };

    using Streamline        = std::vector<Point3D>;
    using StreamlineBatch   = std::vector<Streamline>;
    using Tractogram        = std::vector<Streamline>;
}