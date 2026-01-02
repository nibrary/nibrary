#pragma once

#include "streamlineAutoencoder_aux.h"
#include "dMRI/tractography/utility/resampleStreamline.h"

namespace NIBR
{

    class StreamlineAutoencoder {

        public:

            StreamlineAutoencoder(const StreamlineAutoencoderDefinition& _modelDef, bool _useCPU = false);

            StreamlineAutoencoderDefinition model;
            torch::jit::script::Module      module;
            torch::Device                   device;
            bool                            useCPU;
            bool                            ready;

            int          getInpDim()              {return model.inpDim;}
            ModelType    getAeType()              {return model.aeType;}
            int          getLatDim()              {return model.latDim;}
            int          getFinDim()              {return model.finDim;}
            torch::Dtype getDType()               {return model.dType;}
            double       getDistScaler()          {return model.distScaler;}
            int          getLatDimMultiplier()    {return (model.aeType == SAE_FINSR) ? 1 : 2;}
            bool         isReady()                {return ready;}

            template <typename T>
            std::vector<std::vector<T>> encode(const NIBR::StreamlineBatch& streamlines,  int batchSize = 0, bool performResampling = false);

            template <typename T>
            NIBR::StreamlineBatch       decode(const std::vector<std::vector<T>>& latent, int batchSize = 0, bool useSecondHalf = false, float resampleStepSize = 0);

        private:

            template <typename T> 
            std::vector<std::vector<T>> encode_batch(const NIBR::StreamlineBatch& streamlines, bool performResampling = false);

            template <typename T>
            NIBR::StreamlineBatch       decode_batch(const std::vector<std::vector<T>>& latent, bool useSecondHalf = false, float resampleStepSize = 0);

    };


    // Encoders
    template <typename T> 
    std::vector<std::vector<T>> StreamlineAutoencoder::encode_batch(const NIBR::StreamlineBatch& streamlines, bool performResampling)
    {

        int N = streamlines.size();

        if (N == 0) return std::vector<std::vector<T>>();

        std::vector<T> flattened_data;
        int            latDimMultiplier;
        
        auto runFlatten = [&](const NIBR::StreamlineBatch& toFlatten) -> void {
            disp(MSG_DETAIL,"Flattening...");
            if (model.aeType == SAE_FINSR) {
                disp(MSG_DETAIL,"FINSR model detected.");
                flattened_data   = flatten_streamlines<T>(toFlatten, false);
                latDimMultiplier = 1;
            } else {
                flattened_data   = flatten_streamlines<T>(toFlatten, true);
                latDimMultiplier = 2;
            } 
            disp(MSG_DETAIL,"Flattening completed");
        };

        // Resample if needed
        if (performResampling) {
            disp(MSG_DETAIL, "Resampling streamlines to %d points...", model.inpDim);
            NIBR::StreamlineBatch resampled_streamlines(N);
            for (int i = 0; i < N; ++i) {
                resampled_streamlines[i] = resampleStreamline_withStepCount(streamlines[i], model.inpDim);
            }
            disp(MSG_DETAIL,"Resampling completed");
            runFlatten(resampled_streamlines);
        } else {
            runFlatten(streamlines);
        }

        disp(MSG_DETAIL,"Creating tensor...");
        auto input = torch::from_blob(flattened_data.data(), {latDimMultiplier * N, 3, model.inpDim}, model.dType).to(device).clone();    
        disp(MSG_DETAIL,"Tensor created");

        // Run the model's encode method
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input);
        disp(MSG_DETAIL,"Running encoder...");
        auto encoded = module.get_method("encode")(inputs).toTensor().to(torch::kCPU).contiguous();

        // Copy encoded tensor data to latent vector
        std::vector<std::vector<T>> latent(N, std::vector<T>(latDimMultiplier * model.latDim));

        const T* src_ptr = encoded.template data_ptr<T>();
        for (int i = 0; i < N; i++) {
            T* dst_ptr = latent[i].data();
            std::memcpy(dst_ptr, src_ptr + (i * latDimMultiplier * model.latDim), latDimMultiplier * model.latDim * sizeof(T));
        }

        disp(MSG_DETAIL,"Encoder completed.");
        return latent;

    }

    template <typename T>
    std::vector<std::vector<T>> StreamlineAutoencoder::encode(const NIBR::StreamlineBatch& streamlines, int batchSize, bool performResampling)
    {
        // Safety Check: Ensure the C++ type matches the model's tensor type.
        if (model.dType != get_dtype<T>()) {
            disp(MSG_ERROR, "Type mismatch: encodeStreamlines called with %s but model requires %s.", typeid(T).name(), c10::toString(model.dType));
            return std::vector<std::vector<T>>();
        }

        int N = streamlines.size();
        if (N == 0) return std::vector<std::vector<T>>();

        if (batchSize == 0) return encode_batch<T>(streamlines);

        int batchCnt = (N < batchSize) ? 1 : (N + batchSize - 1) / batchSize;

        int latDimMultiplier = getLatDimMultiplier();

        std::vector<std::vector<T>> latent(N, std::vector<T>(latDimMultiplier * model.latDim));

        // Iterate through the whole tractogram in parallel
        auto run = [&](NIBR::MT::TASK task) -> void {

            int bas = ((int(task.no)+1)*batchSize < N) ? batchSize : (N-int(task.no)*batchSize);
            
            int idx = task.no * batchSize;

            // 1. Prepare the sub-vector of streamlines for the current batch
            NIBR::StreamlineBatch batch_streamlines(streamlines.begin() + idx, streamlines.begin() + idx + bas);

            disp(MSG_DETAIL,"Encoding batch between indices %d - %d", idx, idx + bas);

            // 2. Call the helper function to perform the encoding
            auto latent_batch = encode_batch<T>(batch_streamlines, performResampling);

            disp(MSG_DETAIL,"Encoding completed");

            // 3. Copy the encoded tensor data to the main latent vector
            for (int i = 0; i < bas; i++) {
                std::swap(latent[i+idx],latent_batch[i]);
            }

        };
        NIBR::MT::MTRUN(batchCnt, "Encoding streamlines", run);

        return latent;

    }

    // Decoders
    template <typename T>
    NIBR::StreamlineBatch StreamlineAutoencoder::decode_batch(const std::vector<std::vector<T>>& latent, bool useSecondHalf, float resampleStepSize)
    {
        int N = latent.size();
        if (N == 0) return NIBR::StreamlineBatch();

        int latDimShift = useSecondHalf ? model.latDim : 0;

        // 1. Prepare input tensor from C++ vectors
        std::vector<T> blob(N * model.latDim);
        for (int i = 0; i < N; ++i) {
            std::copy(latent[i].begin() + latDimShift, latent[i].begin() + model.latDim + latDimShift, blob.begin() + i * model.latDim);
        }

        auto input   = torch::from_blob(blob.data(), {(long)N, model.latDim}, model.dType).to(device).clone();

        // 2. Decode the batch with the model
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input);
        at::Tensor decoded = module.get_method("decode")(inputs).toTensor().to(torch::kCPU).contiguous();

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
                if (resampleStepSize > 0) {
                    trk = resampleStreamline_withStepSize(trk, resampleStepSize);
                }
                streamlines[i] = std::move(trk);
            }
        });

        return streamlines;
    }

    template <typename T>
    NIBR::StreamlineBatch StreamlineAutoencoder::decode(const std::vector<std::vector<T>>& latent, int batchSize, bool useSecondHalf, float resampleStepSize)
    {
        // Safety Check: Ensure the C++ type matches the model's tensor type.
        if (model.dType != get_dtype<T>()) {
            disp(MSG_ERROR, "Type mismatch: decodeStreamlines called with %s but model requires %s.", typeid(T).name(), c10::toString(model.dType));
            return NIBR::StreamlineBatch();
        }

        int N = latent.size();
        if (N == 0) return NIBR::StreamlineBatch();

        if (batchSize == 0) return decode_batch<T>(latent,useSecondHalf);

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
            NIBR::StreamlineBatch decoded_batch = decode_batch<T>(latent_batch,useSecondHalf,resampleStepSize);

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
