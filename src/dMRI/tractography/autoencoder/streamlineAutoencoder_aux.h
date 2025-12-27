#pragma once

#include "base/verbose.h"
#include "dMRI/tractography/tractogram.h"
#include <torch/torch.h>
#include <torch/script.h>
#include <c10/util/Half.h>

namespace NIBR
{

    typedef enum {
        SAE_UNKNOWN,
        SAE_1D,
        SAE_FINSR,
    } ModelType;

    struct StreamlineAutoencoderDefinition {
        std::string  moduleFile; // Path to the torch script module
        ModelType    aeType;     // Model type
        int          inpDim;     // Input dimension
        int          latDim;     // Latent dimension
        int          finDim;     // FINSR dimension
        torch::Dtype dType;      // Data type
        double       distScaler; // Scaling factor to match real-space distance (this value can be obtained with the modelTest command)
    };

    // Convenience functions to get torch data type
    template<typename T>
    constexpr torch::Dtype get_dtype();

    template<>
    constexpr torch::Dtype get_dtype<float>()    { return torch::kFloat;  }

    template<>
    constexpr torch::Dtype get_dtype<double>()   { return torch::kDouble; }

    template<>
    constexpr torch::Dtype get_dtype<at::Half>() { return torch::kHalf;   }

    torch::Dtype to_dType(std::string dataType);

    ModelType    to_modelType(std::string aeType);

    // Flattens a streamline and returns it.
    // If flip is true, it also returns the flipped version.
    template <typename T>
    std::vector<T> flatten_streamlines(const NIBR::StreamlineBatch& streamlines, bool flip) {
        
        if (streamlines.empty()) return {};

        size_t number_of_streamlines = streamlines.size();
        size_t points_per_streamline = streamlines[0].size(); 

        // Safety check: ensure all streamlines have the same length
        for (const auto& sl : streamlines) {
            if (sl.size() != points_per_streamline) {
                disp(MSG_ERROR, "Streamline length mismatch. Expected %zu, got %zu.", points_per_streamline, sl.size());
                return {};
            }
        }

        size_t multiplier = flip ? 2 : 1;
        std::vector<T> flattened(multiplier * number_of_streamlines * 3 * points_per_streamline);
        size_t index = 0;

        for (size_t i = 0; i < number_of_streamlines; ++i) {

            // Append original streamline
            for (int j = 0; j < 3; ++j) {
                for (size_t k = 0; k < points_per_streamline; ++k) {
                    flattened[index++] = static_cast<T>(streamlines[i][k][j]);
                }
            }

            if (flip) {
                // Append flipped streamline
                for (int j = 0; j < 3; ++j) {
                    for (int k = points_per_streamline - 1; k != -1; --k) {
                        flattened[index++] = static_cast<T>(streamlines[i][k][j]);
                    }
                }
            }
        }
        return flattened;
    }
    


}