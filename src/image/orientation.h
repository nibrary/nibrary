#pragma once

#include <limits>
#include <vector>
#include <string>
#include <stdexcept>

namespace NIBR
{

    std::vector<std::string> aff2axcodes(const float aff[3][4]);
    // void aff2RAS(float outAff[3][4],const float inpAff[3][4], const int64_t dims[3]);

}