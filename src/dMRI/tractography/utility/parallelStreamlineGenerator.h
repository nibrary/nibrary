#pragma once

#include "base/nibr.h"
#include "dMRI/tractography/io/tractogramReader.h"

namespace NIBR 
{
    // Compute parallel streamlines around a streamline
    // radius:    the outermost radius where streamlines will be generated
    // ringCount: inner ring count
    // pointCountPerRing: number of streamlines to be generated per ring
    // Explanation:
    //     if radius=1, ringCount=2, pointCountPerRing=5, then:
    //     at the center will be the input streamline.
    //     5 streamlines will be located on a ring at radius=0.5 
    //     5 more streamlines will be located on a ring at radius=1
    //     there will be a total of 11 streamlines at the output 
    void getParallelStreamlines(Tractogram& out, NIBR::TractogramReader* tractogram, float radius, int ringCount, int pointCountPerRing);
    StreamlineBatch getParallelStreamlines(const Streamline& streamline, float radius, int ringCount, int pointCountPerRing, int threadId);

    void getParallelStreamlines(Tractogram& out, NIBR::TractogramReader* tractogram, float sigma, int N);
    StreamlineBatch getParallelStreamlines(const Streamline& streamline, float sigma, int N, int threadId);

    // Low level function used to generate a batch of parallel streamlines around a given streamline
    void generateParallelStreamlineBatch(const Streamline& streamline, StreamlineBatch& outBatch, const std::vector<float>& scale_N, const std::vector<float>& scale_B,int threadId);

    
}
