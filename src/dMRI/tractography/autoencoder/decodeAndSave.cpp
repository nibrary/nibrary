#include "streamlineAutoencoder_enc_dec.h"
#include "dMRI/tractography/io/tractogramWriter.h"

using namespace NIBR;

bool NIBR::decodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize, float resampleStepSize)
{
    if (existsFile(out) && !force) return true;

    // Template that deduces the type T (float, double, at::Half) from its argument.
    auto read_decode_and_save = [&](auto type_placeholder) -> bool {

        using T = decltype(type_placeholder);

        int latDimMultiplier = model.getLatDimMultiplier();

        // 1. Open and read the input latent file
        std::ifstream ifs(inp, std::ios::binary);
        if (!ifs.is_open()) {
            disp(MSG_ERROR, "Failed to open file: %s", inp.c_str());
            return false;
        }

        ifs.seekg(0, std::ios::end);
        if (ifs.tellg() == 0) {
            disp(MSG_WARN, "Input file is empty: %s", inp.c_str());
            return true; // Success, nothing to do.
        }

        int N = ifs.tellg() / (sizeof(T) * latDimMultiplier * model.getLatDim());
        ifs.seekg(0, std::ios::beg);

        std::vector<std::vector<T>> latent(N, std::vector<T>(latDimMultiplier * model.getLatDim()));
        for (int i = 0; i < N; ++i) {
            ifs.read(reinterpret_cast<char*>(latent[i].data()), latDimMultiplier * model.getLatDim() * sizeof(T));
        }
        ifs.close();

        // 2. Call decode
        NIBR::Tractogram streamlines = model.decode<T>(latent, batchSize, resampleStepSize);

        if (streamlines.empty()) {
            disp(MSG_ERROR, "Decoding failed or produced no streamlines.");
            return false;
        }

        // 3. Write the final result
        return writeTractogram(out, streamlines);
    };

    // Dispatch to the generic lambda with the correct type
    if (model.getDType() == torch::kFloat) {
        return read_decode_and_save(float{});
    } else if (model.getDType() == torch::kDouble) {
        return read_decode_and_save(double{});
    } else if (model.getDType() == torch::kHalf) {
        return read_decode_and_save(at::Half{});
    } else {
        disp(MSG_ERROR, "Unsupported data type for decoding: %s", c10::toString(model.getDType()));
        return false;
    }
}