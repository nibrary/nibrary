#include "dataTypeHandler.h"

using namespace NIBR;

std::unordered_map<DATATYPE, TypeInfoRef> NIBR::TYPE = {
    {BOOL_DT,      typeid(bool)                        },
    {INT8_DT,      typeid(int8_t)                      },
    {UINT8_DT,     typeid(uint8_t)                     },
    {INT16_DT,     typeid(int16_t)                     },
    {UINT16_DT,    typeid(uint16_t)                    },
    {INT32_DT,     typeid(int32_t)                     },
    {UINT32_DT,    typeid(uint32_t)                    },
    {INT64_DT,     typeid(int64_t)                     },
    {UINT64_DT,    typeid(uint64_t)                    },
    {FLOAT32_DT,   typeid(float)                       },
    {FLOAT64_DT,   typeid(double)                      },
    {FLOAT128_DT,  typeid(long double)                 },
    {COMPLEX64_DT, typeid(std::complex<float>)         },
    {COMPLEX128_DT,typeid(std::complex<double>)        },
    {COMPLEX256_DT,typeid(std::complex<long double>)   }
};

std::unordered_map<TypeInfoRef, DATATYPE, Hasher, EqualTo> NIBR::TYPEIDS = {
    {typeid(bool),                          BOOL_DT},
    {typeid(int8_t),                        INT8_DT},
    {typeid(uint8_t),                       UINT8_DT},
    {typeid(int16_t),                       INT16_DT},
    {typeid(uint16_t),                      UINT16_DT},
    {typeid(int32_t),                       INT32_DT},
    {typeid(uint32_t),                      UINT32_DT},
    {typeid(int64_t),                       INT64_DT},
    {typeid(uint64_t),                      UINT64_DT},
    {typeid(float),                         FLOAT32_DT},
    {typeid(double),                        FLOAT64_DT},
    {typeid(long double),                   FLOAT128_DT},
    {typeid(std::complex<float>),           COMPLEX64_DT},
    {typeid(std::complex<double>),          COMPLEX128_DT},
    {typeid(std::complex<long double>),     COMPLEX256_DT}
};

std::unordered_map<TypeInfoRef, std::string, Hasher, EqualTo> NIBR::TYPENAMES = {
    {typeid(bool),                          "BOOL"},
    {typeid(int8_t),                        "INT8"},
    {typeid(uint8_t),                       "UINT8"},
    {typeid(int16_t),                       "INT16"},
    {typeid(uint16_t),                      "UINT16"},
    {typeid(int32_t),                       "INT32"},
    {typeid(uint32_t),                      "UINT32"},
    {typeid(int64_t),                       "INT64"},
    {typeid(uint64_t),                      "UINT64"},
    {typeid(float),                         "FLOAT32"},
    {typeid(double),                        "FLOAT64"},
    {typeid(long double),                   "FLOAT128"},
    {typeid(std::complex<float>),           "COMPLEX64"},
    {typeid(std::complex<double>),          "COMPLEX128"},
    {typeid(std::complex<long double>),     "COMPLEX256"}
};

DATATYPE NIBR::datatype(std::string type)
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