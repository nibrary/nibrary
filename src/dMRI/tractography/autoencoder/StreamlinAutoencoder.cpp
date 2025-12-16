#include "StreamlineAutoencoder.h"
#include <torch/script.h>
#include <torch/torch.h>
#include <c10/util/Half.h>

using namespace NIBR;

// Common initialization method
void StreamlineAutoencoder::init(const std::string& _moduleFile, int _inpDim, int _latDim, const std::string& dataType, double _distScaler, bool _useCPU, size_t _newGenSize) {

    moduleFile  = _moduleFile;
    inpDim      = _inpDim;
    latDim      = _latDim;
    distScaler  = _distScaler;
    useCPU      = _useCPU;
    ready       = false;
    newGenSize  = _newGenSize;

    // Set device
    device      = torch::kCPU;
    if (!useCPU && torch::cuda::is_available()) {
        device = torch::Device(torch::kCUDA);
        disp(MSG_DETAIL,"Using CUDA");
    } else {
        disp(MSG_DETAIL,"Using CPU");
    }

    // Set data type
    if (dataType == "float64" || dataType == "double") {
        this->dtype = torch::kFloat64;
        disp(MSG_DETAIL, "Model data type set to: double (float64)");
    } else if (dataType == "float32" || dataType == "float") {
        this->dtype = torch::kFloat32;
        disp(MSG_DETAIL, "Model data type set to: float (float32)");
    } else if (dataType == "float16" || dataType == "half") {
        this->dtype = torch::kFloat16;
        disp(MSG_DETAIL, "Model data type set to: half (float16)");
    } else {
        this->dtype = torch::kFloat32; // Default to float32
        if (!dataType.empty()) {
            disp(MSG_WARN, "Unknown data type '%s'. Defaulting to float32.", dataType.c_str());
        }
    }
    

    // Set other parameters
    if (moduleFile.empty()) {
        std::filesystem::path curPath = std::filesystem::absolute(__FILE__);
        std::filesystem::path parentPath = curPath.parent_path().parent_path().parent_path();
        std::filesystem::path default_model = parentPath / "models" / "conv_autoencoder_scripted.pt";
        moduleFile  = default_model.string();
        inpDim      = 256;
        latDim      = 64;
        distScaler  = 0.08;
    } else {
        if (getFileExtension(moduleFile) != "pt") {
            disp(MSG_ERROR,"Torch script module file must have extension .pt");
            ready = false;
            return;
        }
    }

    if (!std::filesystem::exists(moduleFile)) {
        disp(MSG_ERROR,"Model not found");
        ready = false;
        return;
    }

    try {
        disp(MSG_DETAIL,"Loading model %s", moduleFile.c_str());
        module = torch::jit::load(moduleFile);
        module.to(device, dtype);
        module.eval();
        ready = true;
    } catch (const c10::Error& e) {
        disp(MSG_ERROR,"Error loading the model: %s", e.what());
        ready = false;
    }
}


StreamlineAutoencoder::StreamlineAutoencoder(const std::tuple<std::string, int, int, std::string>& moduleSpec, bool _useCPU)
    : device(torch::kCPU), useCPU(_useCPU), ready(false) {
    init(std::get<0>(moduleSpec), std::get<1>(moduleSpec), std::get<2>(moduleSpec), std::get<3>(moduleSpec), 1.0f, _useCPU, 0);
}

StreamlineAutoencoder::StreamlineAutoencoder(const std::tuple<std::string, int, int, std::string>& moduleSpec, bool _useCPU, size_t _newGenSize)
    : device(torch::kCPU), useCPU(_useCPU), ready(false) {
    init(std::get<0>(moduleSpec), std::get<1>(moduleSpec), std::get<2>(moduleSpec), std::get<3>(moduleSpec), 1.0f, _useCPU, _newGenSize);
}

StreamlineAutoencoder::StreamlineAutoencoder(const std::tuple<std::string, int, int, std::string, double>& moduleSpec, bool _useCPU)
    : device(torch::kCPU), useCPU(_useCPU), ready(false) {
    init(std::get<0>(moduleSpec), std::get<1>(moduleSpec), std::get<2>(moduleSpec), std::get<3>(moduleSpec), std::get<4>(moduleSpec), _useCPU, 0);
}

StreamlineAutoencoder::StreamlineAutoencoder(const std::tuple<std::string, int, int, std::string, double>& moduleSpec, bool _useCPU, size_t _newGenSize)
    : device(torch::kCPU), useCPU(_useCPU), ready(false) {
    init(std::get<0>(moduleSpec), std::get<1>(moduleSpec), std::get<2>(moduleSpec), std::get<3>(moduleSpec), std::get<4>(moduleSpec), _useCPU, _newGenSize);
}