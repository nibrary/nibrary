#include "resampleStreamline.h"
#include "streamline_operators.h"
#include "math/interp1D.h"

using namespace NIBR;

#define SPEEDINTERVAL 6
#define RESIDUALTHRES 0.05  // When resampling is with step size, there will generally be a left over piece. If this piece is very short, we append it add the last segment. Otherwise, we split into two even parts and append on both ends. This value determines what is short. 0.05 means, "short" is 5% of step size.

// This is declared only in this scope to be used for NIBR resampling functions
std::vector<std::vector<float>> runStreamlineResampler(std::vector<std::vector<float>>& streamline, std::vector<float>& lenVec, int N, float step, float residual)
{

    if (N<2)
        return std::vector<std::vector<float>>();

    std::vector<std::vector<float>> points;

    if (N==2) {
        points.push_back(streamline.front());
        points.push_back(streamline.back());
        return points;
    }

    float target = (residual>0) ? residual*0.5f : step;
    bool  isLast = false;
    float divS   = 1.0f/float(SPEEDINTERVAL);

    float T0[3], T1[3], p[3], Tn1[3], Tn2[3];

    points.push_back(streamline.front());

    vec3sub(T0, streamline[1], streamline[0]);

    int k=0;

    for (size_t l=1; l<streamline.size(); l++) {

        if (l!=(streamline.size()-1)) {
            vec3sub(Tn1, streamline[l],   streamline[l-1]);
            vec3sub(Tn2, streamline[l+1], streamline[l]  );
            vec3add(T1,Tn1,Tn2);
            T1[0] *= 0.5f;
            T1[1] *= 0.5f;
            T1[2] *= 0.5f;
        } else {
            vec3sub(T1, streamline[l],   streamline[l-1]);
        }

        do {

            while (lenVec[k]>=target) {

                float segLen = (k==0) ? lenVec[k] : lenVec[k] - lenVec[k-1];
                float s = ((1.0f - (lenVec[k]-target)/segLen)+float(k%SPEEDINTERVAL))*divS;

                hermiteInterp(p, streamline[l-1], T0, streamline[l], T1, s);
                points.push_back(std::vector<float>{p[0],p[1],p[2]});

                target += step;                

                if (points.size()==size_t(N-1)) {
                    isLast = true;
                    break;
                }

            }

            k++;

            if (isLast)
                break;

        } while (k%SPEEDINTERVAL>0);


        if (isLast)
            break;

        T0[0] = T1[0];
        T0[1] = T1[1];
        T0[2] = T1[2];

    }

    points.push_back(streamline.back());
    return points;

}


std::vector<std::vector<float>> NIBR::resampleStreamline_withStepSize(std::vector<std::vector<float>>& streamline, float step)
{

    if (step<=0)
        return std::vector<std::vector<float>>();

    auto  lenVec = getStreamlineLength_hermiteWithSpeed(streamline, SPEEDINTERVAL);
    int   N      = lenVec.back()/step;
    float res    = lenVec.back() - float(N) * step;

    N = (res>(step*RESIDUALTHRES) && (res > EPS4)) ? N+3 : N+2;

    if (N<2)
        return std::vector<std::vector<float>>();

    return runStreamlineResampler(streamline, lenVec, N, step, res);
}


std::vector<std::vector<float>> NIBR::resampleStreamline_withStepCount(std::vector<std::vector<float>>& streamline, int N)
{

    if (N<2)
        return std::vector<std::vector<float>>();
    
    auto lenVec = getStreamlineLength_hermiteWithSpeed(streamline, SPEEDINTERVAL);
    float step  = lenVec.back()/float(N-1);

    return runStreamlineResampler(streamline, lenVec, N, step, 0);

}


std::vector<std::vector<float>> NIBR::resampleStreamline_withStepSize(NIBR::TractogramReader* tractogram, int index, float step)
{
    auto streamline = tractogram->readStreamlineVector(index);
    return NIBR::resampleStreamline_withStepSize(streamline, step);
}


std::vector<std::vector<float>> NIBR::resampleStreamline_withStepCount(NIBR::TractogramReader* tractogram, int index, int N)
{
    auto streamline = tractogram->readStreamlineVector(index);
    return NIBR::resampleStreamline_withStepCount(streamline, N);
}



std::vector<std::vector<std::vector<float>>> NIBR::resampleTractogram_withStepSize(NIBR::TractogramReader* _tractogram, float step)
{

    std::vector<std::vector<std::vector<float>>> out;

    if (_tractogram->numberOfStreamlines==0)
        return out;

    out.resize(_tractogram->numberOfStreamlines);

    int numberOfThreads = (_tractogram->numberOfStreamlines<size_t(MT::MAXNUMBEROFTHREADS())) ? _tractogram->numberOfStreamlines : MT::MAXNUMBEROFTHREADS();

    NIBR::TractogramReader *tractogram = new NIBR::TractogramReader[numberOfThreads]();

    for (int t = 0; t < numberOfThreads; t++)
        tractogram[t].copyFrom(*_tractogram);

    auto resample = [&](NIBR::MT::TASK task)->void {
        out[task.no] = NIBR::resampleStreamline_withStepSize(&tractogram[task.threadId], task.no, step);
    };

    if (VERBOSE()>=VERBOSE_INFO) {
		NIBR::MT::MTRUN(_tractogram->numberOfStreamlines, "Resampling streamlines", resample);
    } else {
		NIBR::MT::MTRUN(_tractogram->numberOfStreamlines, resample);
    }


    for (int t = 0; t < numberOfThreads; t++) 
        tractogram[t].destroyCopy();
    delete[] tractogram;

    return out;

}




std::vector<std::vector<std::vector<float>>> NIBR::resampleTractogram_withStepCount(NIBR::TractogramReader* _tractogram, int N)
{

    std::vector<std::vector<std::vector<float>>> out;

    if (_tractogram->numberOfStreamlines==0)
        return out;

    out.resize(_tractogram->numberOfStreamlines);

    int numberOfThreads = (_tractogram->numberOfStreamlines<size_t(MT::MAXNUMBEROFTHREADS())) ? _tractogram->numberOfStreamlines : MT::MAXNUMBEROFTHREADS();

    NIBR::TractogramReader *tractogram = new NIBR::TractogramReader[numberOfThreads]();

    for (int t = 0; t < numberOfThreads; t++)
        tractogram[t].copyFrom(*_tractogram);

    auto resample = [&](NIBR::MT::TASK task)->void {
        out[task.no] = NIBR::resampleStreamline_withStepCount(&tractogram[task.threadId], task.no, N);
    };

    if (VERBOSE()>=VERBOSE_INFO) {
		NIBR::MT::MTRUN(_tractogram->numberOfStreamlines, "Resampling streamlines", resample);
    } else {
		NIBR::MT::MTRUN(_tractogram->numberOfStreamlines, resample);
    }


    for (int t = 0; t < numberOfThreads; t++) 
        tractogram[t].destroyCopy();
    delete[] tractogram;

    return out;

}