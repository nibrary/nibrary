#include "stringOperations.h"
#include <regex>

using namespace NIBR;

std::vector<std::string> NIBR::splitString (const std::string &s, char delim) 
{
    
    std::vector<std::string> out;
    
    std::stringstream ss(s);
    
    std::string item;

    while (getline (ss, item, delim)) {
        out.push_back (item);
    }

    return out;
}




bool NIBR::isNumber(const std::string &s)
{
    return std::regex_match( s, std::regex( ( "((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?" ) ) );
}

std::tuple<bool,Point,float> NIBR::getCenterAndRadius(std::string inp) {
    
    Point p = {0,0,0};
    float r = 0;

    std::vector<std::string> vals = NIBR::splitString(inp, ',');

    if (vals.size()!=4) {
        return std::make_tuple(false,p,r);
    } else {
        try {
            p.x = std::stof(vals[0]);
            p.y = std::stof(vals[1]);
            p.z = std::stof(vals[2]);
            r   = std::stof(vals[3]);
            return std::make_tuple(true,p,r);
        }
        catch (const std::invalid_argument&) {
            return std::make_tuple(false,p,r);
        }
    }

}
