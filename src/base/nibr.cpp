#include "base/nibr.h"

namespace NIBR 
{
    std::string version = std::to_string(NIBRARY_VERSION_MAJOR) + "." + 
                          std::to_string(NIBRARY_VERSION_MINOR) + "." + 
                          std::to_string(NIBRARY_VERSION_PATCH);

    std::string sgntr   = "nibrary v" + version;
    std::string&  SGNTR() {return sgntr;}
}
