#pragma once
#include <Eigen/Dense>
#include "base/verbose.h"
#include "dMRI/tractography/tractogram.h"
#include "dMRI/tractography/io/tractogramWriter.h"
#include <torch/torch.h>
#include <torch/script.h>
#include <c10/util/Half.h>

namespace NIBR
{
    class StreamlineAutoencoder {

        public:
            std::string moduleFile;                 // Path to Torch script module file
            torch::jit::script::Module module;      // Torch script module
            int inpDim;                             // input model dimension < 3 x inpDim>
            int latDim;                             // latent space dimension
            torch::Dtype  dtype  = torch::kFloat32; // data type
            double distScaler    = 1.0;             // scaling factor to match real-space distance (this value can be obtained with the modelTest command)
            torch::Device device = torch::kCPU;     // CPU / GPU
            bool useCPU{false};
            bool ready{false};
            size_t newGenSize = 0;

            StreamlineAutoencoder(const std::tuple<std::string,int,int,std::string>&        _moduleSpec, bool _useCPU);
            StreamlineAutoencoder(const std::tuple<std::string,int,int,std::string>&        _moduleSpec, bool _useCPU, size_t _newGenSize);
            StreamlineAutoencoder(const std::tuple<std::string,int,int,std::string,double>& _moduleSpec, bool _useCPU);
            StreamlineAutoencoder(const std::tuple<std::string,int,int,std::string,double>& _moduleSpec, bool _useCPU, size_t _newGenSize);

            bool isReady() {return ready;}

        private:
            void init(const std::string& _moduleFile, int _inpDim, int _latDim, const std::string& dataType, double _distScaler, bool _useCPU, size_t _newGenSize);

    };

    // Convenience functions to get torch data type
    template<typename T>
    constexpr torch::Dtype get_dtype();

    template<>
    constexpr torch::Dtype get_dtype<float>()    { return torch::kFloat;  }

    template<>
    constexpr torch::Dtype get_dtype<double>()   { return torch::kDouble; }

    template<>
    constexpr torch::Dtype get_dtype<at::Half>() { return torch::kHalf; }


    // Flattens a streamline and returns it together with its flipped version
    // Output has both the original and flipped streamlines.
    template <typename T>
    std::vector<T> flatten_and_flip_streamlines(const NIBR::StreamlineBatch& streamlines) {
        
        if (streamlines.empty()) return {};

        size_t number_of_streamlines = streamlines.size();
        size_t points_per_streamline = streamlines[0].size(); // Assumes all are resampled to the same size

        std::vector<T> flattened(2 * number_of_streamlines * 3 * points_per_streamline);
        size_t index = 0;

        for (size_t i = 0; i < number_of_streamlines; ++i) {

            // Append original streamline
            for (int j = 0; j < 3; ++j) {
                for (size_t k = 0; k < points_per_streamline; ++k) {
                    flattened[index++] = static_cast<T>(streamlines[i][k][j]);
                }
            }

            // Append flipped streamline
            for (int j = 0; j < 3; ++j) {
                for (int k = points_per_streamline - 1; k != -1; --k) {
                    flattened[index++] = static_cast<T>(streamlines[i][k][j]);
                }
            }
        }
        return flattened;
    }

    template <typename T>
    std::vector<T> flatten_streamlines(const NIBR::StreamlineBatch& streamlines) {
        
        if (streamlines.empty()) return {};

        size_t number_of_streamlines = streamlines.size();
        size_t points_per_streamline = streamlines[0].size(); // Assumes all are resampled to the same size

        std::vector<T> flattened(number_of_streamlines * 3 * points_per_streamline);
        size_t index = 0;

        for (size_t i = 0; i < number_of_streamlines; ++i) {

            // Append original streamline
            for (int j = 0; j < 3; ++j) {
                for (size_t k = 0; k < points_per_streamline; ++k) {
                    flattened[index++] = static_cast<T>(streamlines[i][k][j]);
                }
            }
        }
        return flattened;
    }

    // This helper function encodes a batch of streamlines.
    template <typename T>
    std::vector<std::vector<T>> encode_batch(                               // output latent representations <number of streamlines x latDim>
        const NIBR::StreamlineBatch& streamlines,    // input tractogram <number of streamlines x (variable) number of points x 3>
        StreamlineAutoencoder& model                                        // model
    )
    {

        disp(MSG_DETAIL,"Flattening...");
        // Flatten the streamlines for this batch
        std::vector<T> flattened_data;
        if (model.newGenSize > 0) {
            disp(MSG_DETAIL,"Newgen model detected.");
            flattened_data = flatten_streamlines<T>(streamlines);
        } else {
            flattened_data = flatten_and_flip_streamlines<T>(streamlines);
        }
            
        disp(MSG_DETAIL,"Flattening completed");

        int N = streamlines.size();


        // TODO: don't use if/else. make sure this works before
        if(model.newGenSize > 0) {
            disp(MSG_DETAIL,"Creating tensor...");
            auto input   = torch::from_blob(flattened_data.data(), {(long)N, 3, model.inpDim}, model.dtype).to(model.device).clone();    
            disp(MSG_DETAIL,"Tensor created");

            // Run the model's encode method
            std::vector<torch::jit::IValue> inputs;
            inputs.push_back(input);
            disp(MSG_DETAIL,"Running encoder...");
            auto encoded = model.module.get_method("encode")(inputs).toTensor().to(torch::kCPU).contiguous();

            std::vector<std::vector<T>> latent(N, std::vector<T>(model.latDim));

            const T* src_ptr = encoded.template data_ptr<T>();
            for (int i = 0; i < N; i++) {
                T* dst_ptr = latent[i].data();
                std::memcpy(dst_ptr, src_ptr + (i * model.latDim), model.latDim * sizeof(T));
            }

            disp(MSG_DETAIL,"Encoder completed.");
            return latent;

        } else {
            // Create tensor with the correct type and move to the model's device
            disp(MSG_DETAIL,"Creating tensor...");
            auto input   = torch::from_blob(flattened_data.data(), {2 * (long)N, 3, model.inpDim}, model.dtype).to(model.device).clone();    
            disp(MSG_DETAIL,"Tensor created");

            // Run the model's encode method
            std::vector<torch::jit::IValue> inputs;
            inputs.push_back(input);
            disp(MSG_DETAIL,"Running encoder...");
            auto encoded = model.module.get_method("encode")(inputs).toTensor().to(torch::kCPU).contiguous();

            std::vector<std::vector<T>> latent(N, std::vector<T>(2 * model.latDim));

            // Copy encoded tensor data to latent vector
            const T* src_ptr = encoded.template data_ptr<T>();
            for (int i = 0; i < N; i++) {
                T* dst_ptr = latent[i].data();
                std::memcpy(dst_ptr, src_ptr + (i * 2 * model.latDim), 2 * model.latDim * sizeof(T));
            }

            disp(MSG_DETAIL,"Encoder completed.");
            return latent;
        }
        

        
    }

    // Encoder
    template <typename T>
    std::vector<std::vector<T>> encodeStreamlines(                              // output latent representations <number of streamlines x latDim>
        const NIBR::StreamlineBatch& streamlines,        // input tractogram <number of streamlines x (variable) number of points x 3>
        StreamlineAutoencoder& model,                                           // model
        int batchSize                                                           // batch-size
        )
    {
        // Safety Check: Ensure the C++ type matches the model's tensor type.
        if (model.dtype != get_dtype<T>()) {
            disp(MSG_ERROR, "Type mismatch: encodeStreamlines called with %s but model requires %s.", typeid(T).name(), c10::toString(model.dtype));
            return {};
        }

        int N = streamlines.size();
        if (N == 0) return {};

        int batchCnt = (N < batchSize) ? 1 : (N + batchSize - 1) / batchSize;

        int latDimMultiplier = 2;
        if(model.newGenSize > 0) {
            latDimMultiplier = 1;
        }
        std::vector<std::vector<T>> latent(N, std::vector<T>(latDimMultiplier * model.latDim));

        // Iterate through the whole tractogram in parallel
        auto run = [&](NIBR::MT::TASK task) -> void {

            int bas = ((int(task.no)+1)*batchSize < N) ? batchSize : (N-int(task.no)*batchSize);
            
            int idx = task.no * batchSize;

            // 1. Prepare the sub-vector of streamlines for the current batch
            NIBR::StreamlineBatch batch_streamlines(streamlines.begin() + idx, streamlines.begin() + idx + bas);

            disp(MSG_DETAIL,"Encoding batch between indices %d - %d", idx, idx + bas);

            // 2. Call the helper function to perform the encoding
            auto latent_batch = encode_batch<T>(batch_streamlines, model);

            disp(MSG_DETAIL,"Encoding completed");

            // 3. Copy the encoded tensor data to the main latent vector
            for (int i = 0; i < bas; i++) {
                std::swap(latent[i+idx],latent_batch[i]);
            }

        };
        NIBR::MT::MTRUN(batchCnt, "Encoding streamlines", run);

        return latent;

    }

    // This helper function decodes a batch of latent vectors.
    template <typename T>
    NIBR::StreamlineBatch decode_batch(              // output tractogram <number of streamlines x (variable) number of points x 3>
        const std::vector<std::vector<T>>& latent,                          // input latent representations <number of streamlines x latDim>
        StreamlineAutoencoder& model                                        // model
    )
    {
        int N = latent.size();
        if (N == 0) return {};

        // 1. Prepare input tensor from C++ vectors (ignoring the flipped representation).
        // Should also work if don't contain flipped representations
        std::vector<T> blob(N * model.latDim);
        for (int i = 0; i < N; ++i) {
            std::copy(latent[i].begin(), latent[i].begin() + model.latDim, blob.begin() + i * model.latDim);
        }

        auto input   = torch::from_blob(blob.data(), {(long)N, model.latDim}, model.dtype).to(model.device).clone();

        // 2. Decode the batch with the model
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input);
        at::Tensor decoded = model.module.get_method("decode")(inputs).toTensor().to(torch::kCPU).contiguous();

        // 3. Unpack the output tensor into the final C++ streamline format
        NIBR::StreamlineBatch streamlines(N);
        AT_DISPATCH_FLOATING_TYPES(decoded.scalar_type(), "decode_batch_dispatcher", [&] {
            for (int i = 0; i < N; ++i) {
                auto decoded_data = decoded[i].data_ptr<scalar_t>();
                NIBR::Streamline trk(model.inpDim);
                for (int j = 0; j < model.inpDim; ++j) {
                    trk[j][0] = static_cast<float>(decoded_data[j]);
                    trk[j][1] = static_cast<float>(decoded_data[j + model.inpDim]);
                    trk[j][2] = static_cast<float>(decoded_data[j + 2 * model.inpDim]);
                }
                streamlines[i] = std::move(trk);
            }
        });

        return streamlines;
    }

    // Decoder
    template <typename T>
    NIBR::Tractogram decodeStreamlines(             // output tractogram <number of streamlines x (fixed) inpDim x 3>
        const std::vector<std::vector<T>>& latent,                              // input latent representations <number of streamlines x latDim>
        StreamlineAutoencoder& model,                                           // model
        int batchSize                                                           // batch-size
        )
    {
        // Safety Check: Ensure the C++ type matches the model's tensor type.
        if (model.dtype != get_dtype<T>()) {
            disp(MSG_ERROR, "Type mismatch: decodeStreamlines called with %s but model requires %s.", typeid(T).name(), c10::toString(model.dtype));
            return {};
        }

        int N = latent.size();
        if (N == 0) return {};

        int batchCnt = (N < batchSize) ? 1 : (N + batchSize - 1) / batchSize;
        
        // 1. Pre-allocate the entire output vector to allow for parallel writes.
        NIBR::StreamlineBatch streamlines(N);

        // 2. Define the lambda for parallel execution.
        auto run = [&](NIBR::MT::TASK task) -> void {
            int idx = task.no * batchSize;
            int bas = std::min(batchSize, N - idx);

            // Prepare the sub-vector of latent data for the current batch.
            std::vector<std::vector<T>> latent_batch(latent.begin() + idx, latent.begin() + idx + bas);

            // Call the helper to decode this single batch.
            NIBR::StreamlineBatch decoded_batch = decode_batch<T>(latent_batch, model);

            // Copy the batch results into the correct slice of the pre-allocated output vector.
            for (int i = 0; i < bas; ++i) {
                int global_idx = i + idx;
                streamlines[global_idx] = std::move(decoded_batch[i]);
            }
        };
        
        // 3. Run the parallel processing.
        NIBR::MT::MTRUN(batchCnt, "Decoding streamlines", run);

        return streamlines;
    }
}
