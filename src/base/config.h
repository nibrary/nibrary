#pragma once

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <set>

#ifdef BUILD_FOR_WINDOWS
#include <windows.h>
#include <io.h>
#undef max
#undef min
#else
#include <unistd.h>
#endif
#include "byteSwapper.h"
#include "multithreader.h"
#include "dataTypeHandler.h"
#include "verbose.h"

namespace NIBR 
{

    std::string&   SGNTR();
    VERBOSE_LEVEL& VERBOSE();

    void INITIALIZE();
    void TERMINATE();
    
    int  RUNTIME();
    int  MSECRUNTIME();
    void TIC();
    void TOC();
    void TOC(MESSAGE msg);

    typedef enum {
        UNKNOWNSPACEUNIT,
        METER,
        MM,
        MICRON
    } SPACEUNIT;

    typedef enum {
        UNKNOWNTIMEUNIT,
        SEC,
        MSEC,
        USEC,
        HZ,
        PPM,
        RADS
    } TIMEUNIT;

}
