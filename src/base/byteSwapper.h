#pragma once

#include "base/nibr.h"
#include <cstddef>
#include <algorithm>

namespace NIBR 
{

    template<typename T>
    inline void swapByteOrder(T& a) 
    {
        unsigned char* byteArray = reinterpret_cast<unsigned char*>(&a);
        for(std::size_t i = 0; i < (sizeof(a)/2); i++)
            std::swap(byteArray[sizeof(a) - 1 - i],byteArray[i]);
    }

}
