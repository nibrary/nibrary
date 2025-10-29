#pragma once

#include "base/nibr.h"
#include "math/core.h"
#include <sstream>

namespace NIBR 
{

    double isNumber(const std::string &s);

    std::vector<std::string> splitString (const std::string &s, char delim);
    
    std::tuple<bool,Point3D,float> getCenterAndRadius(std::string inp);

    template <typename T>
    std::string to_string_with_precision(const T val, const int n = 4)
    {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << val;
        return std::move(out).str();
    }

    // Locale guard to set numeric locale to "C" within a scope
    class CNumericLocaleGuard {
    public:
        explicit CNumericLocaleGuard() : m_old_locale_name() {
            m_old_locale_name = std::setlocale(LC_NUMERIC, NULL);
            std::setlocale(LC_NUMERIC, "C");
        }
        ~CNumericLocaleGuard() {std::setlocale(LC_NUMERIC, m_old_locale_name.c_str());}

        // A locale guard should be tied to the scope it was created in.
        // Disabling copy/move prevents accidental misuse.
        CNumericLocaleGuard(const CNumericLocaleGuard&) = delete;
        CNumericLocaleGuard& operator=(const CNumericLocaleGuard&) = delete;
        CNumericLocaleGuard(CNumericLocaleGuard&&) = delete;
        CNumericLocaleGuard& operator=(CNumericLocaleGuard&&) = delete;

    private:
        std::string m_old_locale_name;
    };
    
}