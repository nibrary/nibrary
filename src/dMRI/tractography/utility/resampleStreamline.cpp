#include "resampleStreamline.h"
#include "streamline_operators.h"
#include "math/interp1D.h"

using namespace NIBR;

#define SPEEDINTERVAL 6
#define RESIDUALTHRES 0.05  // When resampling is with step size, there will generally be a left over piece. If this piece is very short, we append it add the last segment. Otherwise, we split into two even parts and append on both ends. This value determines what is short. 0.05 means, "short" is 5% of step size.

// This is declared only in this scope to be used for NIBR resampling functions
Streamline runStreamlineResampler(const Streamline& streamline, std::vector<float>& lenVec, int N, float step, float residual)
{

    if (N<2)
        return Streamline();

    Streamline points;

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

    for (std::size_t l=1; l<streamline.size(); l++) {

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
                points.push_back({p[0],p[1],p[2]});

                target += step;                

                if (points.size()==std::size_t(N-1)) {
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


Streamline NIBR::resampleStreamline_withStepSize(const Streamline& streamline, float step)
{

    if (step<=0)
        return Streamline();

    auto  lenVec = getStreamlineLength_hermiteWithSpeed(streamline, SPEEDINTERVAL);
    int   N      = lenVec.back()/step;
    float res    = lenVec.back() - float(N) * step;

    N = (res>(step*RESIDUALTHRES) && (res > EPS4)) ? N+3 : N+2;

    if (N<2)
        return Streamline();

    return runStreamlineResampler(streamline, lenVec, N, step, res);
}


Streamline NIBR::resampleStreamline_withStepCount(const Streamline& streamline, int N)
{

    if (N<2)
        return Streamline();
    
    auto lenVec = getStreamlineLength_hermiteWithSpeed(streamline, SPEEDINTERVAL);
    float step  = lenVec.back()/float(N-1);

    return runStreamlineResampler(streamline, lenVec, N, step, 0);

}


StreamlineBatch NIBR::resampleTractogram_withStepSize(const StreamlineBatch& batch_in, float step)
{
    StreamlineBatch batch_out(batch_in.size());

    auto resample = [&](const NIBR::MT::TASK& task)->void {
        batch_out[task.no] = NIBR::resampleStreamline_withStepSize(batch_in[task.no], step);
    };
    
    if (VERBOSE()>=VERBOSE_INFO) {
		NIBR::MT::MTRUN(batch_in.size(), "Resampling streamlines", resample);
    } else {
		NIBR::MT::MTRUN(batch_in.size(), resample);
    }

    return batch_out;
}



StreamlineBatch NIBR::resampleTractogram_withStepCount(const StreamlineBatch& batch_in, int N)
{
    StreamlineBatch batch_out(batch_in.size());

    auto resample = [&](const NIBR::MT::TASK& task)->void {
        batch_out[task.no] = NIBR::resampleStreamline_withStepCount(batch_in[task.no], N);
    };
    
    if (VERBOSE()>=VERBOSE_INFO) {
		NIBR::MT::MTRUN(batch_in.size(), "Resampling streamlines", resample);
    } else {
		NIBR::MT::MTRUN(batch_in.size(), resample);
    }

    return batch_out;
}

StreamlineBatch NIBR::resampleBatch(const StreamlineBatch& batch_in,float stepSize, int stepCount, bool useSizeOpt) 
{
    StreamlineBatch batch_out(batch_in.size());

    auto resample = [&](const NIBR::MT::TASK& task)->void {
        if (useSizeOpt) {
            batch_out[task.no] = NIBR::resampleStreamline_withStepSize(batch_in[task.no], stepSize);
        } else {
            batch_out[task.no] = NIBR::resampleStreamline_withStepCount(batch_in[task.no], stepCount);
        }
    };
    
    if (VERBOSE()>=VERBOSE_INFO) {
		NIBR::MT::MTRUN(batch_in.size(), "Resampling streamlines", resample);
    } else {
		NIBR::MT::MTRUN(batch_in.size(), resample);
    }

    return batch_out;
}