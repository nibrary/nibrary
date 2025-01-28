#pragma once

#include "base/nibr.h"
#include <cstddef>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#undef max
#undef min
#else
#include <unistd.h>
#endif

namespace NIBR 
{

    template<typename T>
    inline void swapByteOrder(T& a) 
    {
        char* byteArray = reinterpret_cast<char*>(&a);
        for(std::size_t i = 0; i < (sizeof(a)/2); i++)
            std::swap(byteArray[sizeof(a) - 1 - i],byteArray[i]);
    }

}
