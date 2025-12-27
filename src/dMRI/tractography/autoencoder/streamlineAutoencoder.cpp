#include "base/fileOperations.h"
#include "base/verbose.h"
#include "streamlineAutoencoder.h"

using namespace NIBR;

NIBR::StreamlineAutoencoder::StreamlineAutoencoder(const StreamlineAutoencoderDefinition& _modelDef, bool _useCPU) : device(torch::kCPU) {

    model  = _modelDef;
    useCPU = _useCPU;
    ready  = false;

    // Set device
    if (!useCPU && torch::cuda::is_available()) {
        device = torch::Device(torch::kCUDA);
        disp(MSG_DETAIL,"Using CUDA");
    } else {
        disp(MSG_DETAIL,"Using CPU");
    }

    // Set other parameters
    if (model.moduleFile.empty()) {
        std::filesystem::path curPath = std::filesystem::absolute(__FILE__);
        std::filesystem::path parentPath = curPath.parent_path().parent_path().parent_path();
        std::filesystem::path default_model = parentPath / "models" / "conv_autoencoder_scripted.pt";
        model.moduleFile  = default_model.string();
        model.aeType      = SAE_1D;
        model.inpDim      = 256;
        model.latDim      = 64;
        model.finDim      = 0;
        model.dType       = torch::kFloat32;
        model.distScaler  = 0.08;
    } else {
        if (getFileExtension(model.moduleFile) != "pt") {
            disp(MSG_ERROR,"Torch script module file must have extension .pt");
            ready = false;
            return;
        }
    }

    if (!std::filesystem::exists(model.moduleFile)) {
        disp(MSG_ERROR,"Model not found");
        ready = false;
        return;
    }

    try {
        disp(MSG_DETAIL,"Loading model %s", model.moduleFile.c_str());
        module = torch::jit::load(model.moduleFile);
        module.to(device, model.dType);
        module.eval();
        ready = true;
    } catch (const c10::Error& e) {
        disp(MSG_ERROR,"Error loading the model: %s", e.what());
        ready = false;
    }
    
}
