#include "streamlineAutoencoder_enc_dec.h"
#include "dMRI/tractography/io/tractogramReader.h"

using namespace NIBR;

bool NIBR::encodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize)
{
    if (existsFile(out) && !force) return true;

    // Prepare input
    NIBR::TractogramReader _tractogram(inp, false);
    if(!_tractogram.isReady()) {
        disp(MSG_ERROR, "Failed opening tractogram.");
        return false;
    }

    NIBR::Tractogram tracObj = _tractogram.getTractogram();
    
    // Template that deduces the type T (float, double, at::Half) from its argument.
    auto process_with_type = [&](auto type_placeholder) -> bool {
        
        using T = decltype(type_placeholder);

        // Write to a temporary file to prevent corruption on failure
        std::string tmp_out = out + ".tmp";
        std::ofstream ofs(tmp_out, std::ios::binary);
        if (!ofs) {
            disp(MSG_ERROR, "Failed to open temporary output file for writing.");
            return false;
        }

        // Encode the streamlines
        // The encode function handles batching and resampling internally
        disp(MSG_INFO,"Encoding %d streamlines...", tracObj.size());
        auto latent = model.encode<T>(tracObj, batchSize);

        for (const auto& vec : latent) {
            ofs.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T));
        }
        ofs.close();

        // On success, rename the temporary file to the final output file
        if (std::rename(tmp_out.c_str(), out.c_str()) != 0) {
            disp(MSG_ERROR, "Failed to rename temporary file.");
            return false;
        }
        return true;
    };

    // Dispatch to the generic lambda with the correct type
    if (model.getDType() == torch::kFloat) {
        return process_with_type(float{});
    } else if (model.getDType() == torch::kDouble) {
        return process_with_type(double{});
    } else if (model.getDType() == torch::kHalf) {
        return process_with_type(at::Half{});
    } else {
        disp(MSG_ERROR, "Unsupported data type for encoding: %s", c10::toString(model.getDType()));
        return false;
    }
}

