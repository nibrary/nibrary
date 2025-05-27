#pragma once

#include "base/nibr.h"
#include "dMRI/tractography/io/tractogramReader.h"
#include "streamline_operators.h"

namespace NIBR 
{

    std::vector<std::vector<float>> resampleStreamline_withStepSize (const std::vector<std::vector<float>>& streamline, float step);
    std::vector<std::vector<float>> resampleStreamline_withStepCount(const std::vector<std::vector<float>>& streamline, int N);

    std::vector<std::vector<std::vector<float>>> resampleTractogram_withStepSize(const std::vector<std::vector<std::vector<float>>>& batch_in, float step);
    std::vector<std::vector<std::vector<float>>> resampleTractogram_withStepCount(const std::vector<std::vector<std::vector<float>>>& batch_in, int N);

    std::vector<std::vector<std::vector<float>>> resampleBatch(const std::vector<std::vector<std::vector<float>>>& batch_in,float stepSize,int stepCount,bool useSizeOpt); 

}