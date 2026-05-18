#include "parallelStreamlineGenerator.h"

using namespace NIBR;

void NIBR::generateParallelStreamlineBatch(
    const Streamline& streamline, 
    StreamlineBatch& outBatch, 
    const std::vector<float>& scale_N, 
    const std::vector<float>& scale_B,
    int threadId) 
{
    int len = streamline.size();
    int par = scale_N.size();

    if (len < 2) {
        outBatch.clear();
        return;
    }

    outBatch.resize(par);
    for (int p = 0; p < par; p++) {
        outBatch[p].resize(len);
    }

    float T[3], N[3], B[3], curr_T[3];
    const int window = 6; 

    // 1. Initial tangent
    int lookahead = std::min(len - 1, window);
    vec3sub(T, streamline[lookahead], streamline[0]);
    float lenT = norm(T);

    while (lenT < EPS6 && lookahead < len - 1) {
        lookahead++;
        vec3sub(T, streamline[lookahead], streamline[0]);
        lenT = norm(T);
    }

    if (lenT < EPS6) { 
        lenT = 1.0f; 
        T[0] = 1.0f; 
        T[1] = 0.0f; 
        T[2] = 0.0f; 
    } 
    
    T[0] /= lenT; 
    T[1] /= lenT; 
    T[2] /= lenT;

    // 2. Initial normal and binormal
    NIBR::MT::RNDM()[threadId]->getAUnitRandomPerpVector(N, T);
    normalize(N);
    cross(B, T, N);

    for (int p = 0; p < par; p++)
        for (int i = 0; i < 3; i++)
            outBatch[p][0][i] = streamline[0][i] + N[i] * scale_N[p] + B[i] * scale_B[p];

    for (auto l = 1; l < (len - 1); l++) {

        // Take the symmetric chord around this point
        int prev_idx = std::max(0, l - window);
        int next_idx = std::min(len - 1, l + window);
        vec3sub(curr_T, streamline[next_idx], streamline[prev_idx]);
        float lenCT = norm(curr_T);
        
        // If the curve has repeat points or perfect 180 turn inside the window, the chord cancels out. 
        // In that case, search for the nearest true geometric neighbors.
        if (lenCT < EPS6) {

            // Scan backward for the first distinct point
            prev_idx = l - 1;
            while (prev_idx >= 0) {
                float temp_v[3];
                vec3sub(temp_v, streamline[l], streamline[prev_idx]);
                if (norm(temp_v) > EPS6) break;
                prev_idx--;
            }
            
            // Scan forward for the first distinct point
            next_idx = l + 1;
            while (next_idx < len) {
                float temp_v[3];
                vec3sub(temp_v, streamline[next_idx], streamline[l]);
                if (norm(temp_v) > EPS6) break;
                next_idx++;
            }

            prev_idx = std::max(0, prev_idx);
            next_idx = std::min(len - 1, next_idx);
            vec3sub(curr_T, streamline[next_idx], streamline[prev_idx]);
            lenCT = norm(curr_T);
        }

        if (lenCT > EPS6) {
            curr_T[0] /= lenCT; 
            curr_T[1] /= lenCT; 
            curr_T[2] /= lenCT;

            float c = dot(T, curr_T);

            // Compute parallel transport
            if (c < -0.999f) {
                N[0] = -N[0]; N[1] = -N[1]; N[2] = -N[2];
            } else {
                float v[3];
                cross(v, T, curr_T);
                float v_cross_N[3];
                cross(v_cross_N, v, N);
                float v_dot_N = dot(v, N);
                float factor = v_dot_N / (1.0f + c);

                N[0] = N[0] * c + v_cross_N[0] + v[0] * factor;
                N[1] = N[1] * c + v_cross_N[1] + v[1] * factor;
                N[2] = N[2] * c + v_cross_N[2] + v[2] * factor;
            }

            // Ensure N is orthogonal to the chord
            float n_dot_t = dot(N, curr_T);
            N[0] -= n_dot_t * curr_T[0];
            N[1] -= n_dot_t * curr_T[1];
            N[2] -= n_dot_t * curr_T[2];
            
            normalize(N);

            T[0] = curr_T[0]; T[1] = curr_T[1]; T[2] = curr_T[2];
        }

        cross(B, T, N);

        for (int p = 0; p < par; p++) {
            for (int i = 0; i < 3; i++) {
                outBatch[p][l][i] = streamline[l][i] + N[i] * scale_N[p] + B[i] * scale_B[p];
            }
        }
    }

    for (int p = 0; p < par; p++)
        for (int i = 0; i < 3; i++)
            outBatch[p][len - 1][i] = streamline[len - 1][i] + N[i] * scale_N[p] + B[i] * scale_B[p];        
}


void NIBR::getParallelStreamlines(Tractogram& out, NIBR::TractogramReader* tractogram, float radius, int ringCount, int pointCountPerRing)
{
    int numStreamlines = tractogram->numberOfStreamlines;
    int par = ringCount * pointCountPerRing + 1;
    
    if (numStreamlines < 1) return;
    out.resize(numStreamlines * par);

    // Prepare offsets
    float radStep = radius / float(ringCount);
    float angStep = TWOPI / float(pointCountPerRing);
    
    std::vector<float> scale_N(1, 0.0f);
    std::vector<float> scale_B(1, 0.0f);
    
    for (float i = 1; i < (ringCount + 1); i++) {
        for (float j = 0; j < pointCountPerRing; j++) {
            scale_N.push_back(i * radStep * std::sin(angStep * j));
            scale_B.push_back(i * radStep * std::cos(angStep * j));
        }
    }

    tractogram->reset();

    auto genParallel = [&](const NIBR::MT::TASK& task)->void {
        auto [success, streamline, streamlineId] = tractogram->getNextStreamline();
        if (!success || streamline.size() < 2) return;

        StreamlineBatch batch;
        generateParallelStreamlineBatch(streamline, batch, scale_N, scale_B, task.threadId);

        for (int p = 0; p < par; p++) {
            out[streamlineId * par + p] = std::move(batch[p]);
        }
    };

    NIBR::MT::MTRUN(numStreamlines, "Generating parallel streamlines", genParallel);
}

StreamlineBatch NIBR::getParallelStreamlines(const Streamline& streamline, float radius, int ringCount, int pointCountPerRing, int threadId)
{
    StreamlineBatch batch;
    float radStep = radius / float(ringCount);
    float angStep = TWOPI / float(pointCountPerRing);
    
    std::vector<float> scale_N(1, 0.0f);
    std::vector<float> scale_B(1, 0.0f);
    
    for (float i = 1; i < (ringCount + 1); i++) {
        for (float j = 0; j < pointCountPerRing; j++) {
            scale_N.push_back(i * radStep * std::sin(angStep * j));
            scale_B.push_back(i * radStep * std::cos(angStep * j));
        }
    }

    generateParallelStreamlineBatch(streamline, batch, scale_N, scale_B, threadId);
    return batch;
}


void NIBR::getParallelStreamlines(Tractogram& out, NIBR::TractogramReader* tractogram, float sigma, int par)
{
    int numStreamlines = tractogram->numberOfStreamlines;

    if (numStreamlines < 1) return;
    out.resize(numStreamlines * par);

    std::vector<float> scale_N(1, 0.0f);
    std::vector<float> scale_B(1, 0.0f);
    
    for (int n = 0; n < par; n++) {
        scale_N.push_back(NIBR::MT::RNDM()[0]->normal_m0_s1() * sigma);
        scale_B.push_back(NIBR::MT::RNDM()[0]->normal_m0_s1() * sigma);
    }

    tractogram->reset();

    auto genParallel = [&](const NIBR::MT::TASK& task)->void {
        auto [success, streamline, streamlineId] = tractogram->getNextStreamline();
        if (!success || streamline.size() < 2) return;

        StreamlineBatch batch;
        generateParallelStreamlineBatch(streamline, batch, scale_N, scale_B, task.threadId);

        for (int p = 0; p < par; p++) {
            out[streamlineId * par + p] = std::move(batch[p]);
        }
    };

    NIBR::MT::MTRUN(numStreamlines, "Generating parallel streamlines", genParallel);
}

StreamlineBatch NIBR::getParallelStreamlines(const Streamline& streamline, float sigma, int N, int threadId)
{
    StreamlineBatch batch;
    std::vector<float> scale_N(1, 0.0f);
    std::vector<float> scale_B(1, 0.0f);
    
    for (int n = 0; n < N; n++) {
        scale_N.push_back(NIBR::MT::RNDM()[threadId]->normal_m0_s1() * sigma);
        scale_B.push_back(NIBR::MT::RNDM()[threadId]->normal_m0_s1() * sigma);
    }

    generateParallelStreamlineBatch(streamline, batch, scale_N, scale_B, threadId);
    return batch;
}