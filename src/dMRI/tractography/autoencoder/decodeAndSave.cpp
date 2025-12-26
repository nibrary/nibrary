#include "streamlineAutoencoder.h"
#include <type_traits>
#include "base/fileOperations.h"
#include <c10/util/Half.h>

using namespace NIBR;

bool NIBR::decodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize)
{
    if (existsFile(out) && !force) return true;

    // Template that deduces the type T (float, double, at::Half) from its argument.
    auto read_decode_and_save = [&](auto type_placeholder) -> bool {

        using T = decltype(type_placeholder);

        int latDimMultiplier = 2;
        if(model.newGenSize > 0) {
            latDimMultiplier = 1;
        }

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

        int N = ifs.tellg() / (sizeof(T) * latDimMultiplier * model.latDim);
        ifs.seekg(0, std::ios::beg);

        std::vector<std::vector<T>> latent(N, std::vector<T>(latDimMultiplier * model.latDim));
        for (int i = 0; i < N; ++i) {
            ifs.read(reinterpret_cast<char*>(latent[i].data()), latDimMultiplier * model.latDim * sizeof(T));
        }
        ifs.close();

        // 2. Call decodeStreamlines
        NIBR::Tractogram streamlines = decodeStreamlines<T>(latent, model, batchSize);

        if (streamlines.empty()) {
            disp(MSG_ERROR, "Decoding failed or produced no streamlines.");
            return false;
        }

        // 3. Write the final result
        return writeTractogram(out, streamlines);
    };

    // Dispatch to the generic lambda with the correct type
    if (model.dtype == torch::kFloat) {
        return read_decode_and_save(float{});
    } else if (model.dtype == torch::kDouble) {
        return read_decode_and_save(double{});
    } else if (model.dtype == torch::kHalf) {
        return read_decode_and_save(at::Half{});
    } else {
        disp(MSG_ERROR, "Unsupported data type for decoding: %s", c10::toString(model.dtype));
        return false;
    }
}