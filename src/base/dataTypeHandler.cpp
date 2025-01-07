#include "dataTypeHandler.h"

using namespace NIBR;

std::string NIBR::getTypeName(const std::type_info& T) {

    if (T == typeid(bool))                          return "BOOL";
    if (T == typeid(int8_t))                        return "INT8";
    if (T == typeid(uint8_t))                       return "UINT8";
    if (T == typeid(int16_t))                       return "INT16";
    if (T == typeid(uint16_t))                      return "UINT16";
    if (T == typeid(int32_t))                       return "INT32";
    if (T == typeid(uint32_t))                      return "UINT32";
    if (T == typeid(int64_t))                       return "INT64";
    if (T == typeid(uint64_t))                      return "UINT64";
    if (T == typeid(float))                         return "FLOAT32";
    if (T == typeid(double))                        return "FLOAT64";
    if (T == typeid(long double))                   return "FLOAT128";
    if (T == typeid(std::complex<float>))           return "COMPLEX64";
    if (T == typeid(std::complex<double>))          return "COMPLEX128";
    if (T == typeid(std::complex<long double>))     return "COMPLEX256";

    return "UNKNOWN";

}

DATATYPE NIBR::getTypeId(const std::type_info& T) {

    if (T == typeid(bool))                          return BOOL_DT;
    if (T == typeid(int8_t))                        return INT8_DT;
    if (T == typeid(uint8_t))                       return UINT8_DT;
    if (T == typeid(int16_t))                       return INT16_DT;
    if (T == typeid(uint16_t))                      return UINT16_DT;
    if (T == typeid(int32_t))                       return INT32_DT;
    if (T == typeid(uint32_t))                      return UINT32_DT;
    if (T == typeid(int64_t))                       return INT64_DT;
    if (T == typeid(uint64_t))                      return UINT64_DT;
    if (T == typeid(float))                         return FLOAT32_DT;
    if (T == typeid(double))                        return FLOAT64_DT;
    if (T == typeid(long double))                   return FLOAT128_DT;
    if (T == typeid(std::complex<float>))           return COMPLEX64_DT;
    if (T == typeid(std::complex<double>))          return COMPLEX128_DT;
    if (T == typeid(std::complex<long double>))     return COMPLEX256_DT;

    return UNKNOWN_DT;

}

DATATYPE NIBR::getTypeId(std::string type)
{
    if (type=="BOOL")         return BOOL_DT;
    if (type=="INT8")         return INT8_DT;
    if (type=="UINT8")        return UINT8_DT;
    if (type=="INT16")        return INT16_DT;
    if (type=="UINT16")       return UINT16_DT;
    if (type=="INT")          return INT32_DT;
    if (type=="INT32")        return INT32_DT;
    if (type=="UINT32")       return UINT32_DT;
    if (type=="INT64")        return INT64_DT;
    if (type=="UINT64")       return UINT64_DT;
    if (type=="FLOAT")        return FLOAT32_DT;
    if (type=="FLOAT32")      return FLOAT32_DT;
    if (type=="FLOAT64")      return FLOAT64_DT;
    if (type=="FLOAT128")     return FLOAT128_DT;
    if (type=="COMPLEX64")    return COMPLEX64_DT;
    if (type=="COMPLEX128")   return COMPLEX128_DT;
    if (type=="COMPLEX256")   return COMPLEX256_DT;

    return UNKNOWN_DT;
}