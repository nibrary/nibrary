#include "streamlineAutoencoder.h"
#include "base/fileOperations.h"
#include <c10/util/Half.h>

using namespace NIBR;

bool NIBR::encodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize)
{
    if (existsFile(out) && !force) return true;

    // Prepare input
    NIBR::TractogramReader _tractogram(inp, false);
    if(!_tractogram.isReady()){
        disp(MSG_ERROR, "Failed opening tractogram.");
        return false;
    }

    NIBR::Tractogram tracObj = _tractogram.getTractogram();

    int latDimMultiplier = 2;
    if(model.newGenSize > 0) {
        latDimMultiplier = 1;
    }

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

        int N = _tractogram.numberOfStreamlines;
        int batchCnt = (N < batchSize) ? 1 : (N + batchSize - 1) / batchSize;
        std::vector<std::vector<std::vector<T>>> results(batchCnt);

        auto run = [&](NIBR::MT::TASK task) -> void {
            int bas = std::min(batchSize, N - (int)task.no * batchSize);
            NIBR::StreamlineBatch streamlines(bas);
            for (int i = 0; i < bas; i++) {
                int idx = i + (int)task.no * batchSize;
                auto tmp = tracObj[idx];
                streamlines[i] = resampleStreamline_withStepCount(tmp, model.inpDim);
            }
            results[task.no] = encode_batch<T>(streamlines, model);
        };

        NIBR::MT::MTRUN(batchCnt, "Encoding streamlines", run);

        for (const auto& encoded_batch : results) {
            for (const auto& encoded_streamline : encoded_batch) {
                ofs.write(reinterpret_cast<const char*>(encoded_streamline.data()), latDimMultiplier * model.latDim * sizeof(T));
            }
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
    if (model.dtype == torch::kFloat) {
        return process_with_type(float{});
    } else if (model.dtype == torch::kDouble) {
        return process_with_type(double{});
    } else if (model.dtype == torch::kHalf) {
        return process_with_type(at::Half{});
    } else {
        disp(MSG_ERROR, "Unsupported data type for encoding: %s", c10::toString(model.dtype));
        return false;
    }
}

