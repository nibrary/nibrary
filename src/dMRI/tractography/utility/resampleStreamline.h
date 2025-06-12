#pragma once

#include "base/nibr.h"
#include "dMRI/tractography/io/tractogramReader.h"
#include "streamline_operators.h"

namespace NIBR 
{

    Streamline resampleStreamline_withStepSize (const Streamline& streamline, float step);
    Streamline resampleStreamline_withStepCount(const Streamline& streamline, int N);

    StreamlineBatch resampleTractogram_withStepSize(const StreamlineBatch& batch_in, float step);
    StreamlineBatch resampleTractogram_withStepCount(const StreamlineBatch& batch_in, int N);

    StreamlineBatch resampleBatch(const StreamlineBatch& batch_in,float stepSize,int stepCount,bool useSizeOpt); 

}