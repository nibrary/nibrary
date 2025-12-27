#include "base/verbose.h"
#include "proxsuite/proxqp/dense/model.hpp"
#include "streamlineAutoencoder_aux.h"

using namespace NIBR;

torch::Dtype NIBR::to_dType(std::string dataType) 
{
        
    if (dataType == "float64" || dataType == "double") {
        return torch::kFloat64;
    } else if (dataType == "float32" || dataType == "float") {
        return torch::kFloat32;
    } else if (dataType == "float16" || dataType == "half") {
        return torch::kFloat16;
    } else {
        disp(MSG_WARN, "Unknown data type '%s'. Defaulting to float32.", dataType.c_str());
        return torch::kFloat32;
    }
    
}

ModelType NIBR::to_modelType(std::string aeType) 
{
    if (aeType == "1D") {
        return SAE_1D;
    } else if (aeType == "FINSR") {
        return SAE_FINSR;
    } else {
        disp(MSG_ERROR, "Unknown model type '%s'.", aeType.c_str());
        return SAE_UNKNOWN;
    }
}