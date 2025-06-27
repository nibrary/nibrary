#pragma once

#include "base/nibr.h"
#include "math/core.h"
#include "math/gaussian.h"
#include "dMRI/tractography/mappers/tractogram2imageMapper.h"
#include "dMRI/tractography/utility/segmentOperators.h"

using SegmentContribution       = std::pair<int,float>;                             // int: streamlineId, float: amount of contribution
using ContributionsInAVoxel     = std::pair<int,std::vector<SegmentContribution>>;  // int: voxel index,  vector of all segment contributions in the voxel
using Contributions             = std::vector<ContributionsInAVoxel>;               // Vector of all contributions

namespace NIBR 
{

    template<typename T>
    inline void allocateGrid_4streamlineIdAndGaussionWeightedLength_memVector(Tractogram2ImageMapper<T>* tim) {
        tim->template allocateGrid<std::vector<SegmentContribution>>();
    }

    template<typename T>
    inline void deallocateGrid_4streamlineIdAndGaussionWeightedLength_memVector(Tractogram2ImageMapper<T>* tim) {
        tim->template deallocateGrid<std::vector<SegmentContribution>>();
    }

    template<typename T>
    inline void processor_4streamlineIdAndGaussionWeightedLength_memVector(Tractogram2ImageMapper<T>* tim, int* gridPos, NIBR::Segment& seg) {

        int64_t ind = tim->img->sub2ind(gridPos[0],gridPos[1],gridPos[2]);
        
        float voxCenter[3], p[3], step[3];

        static float stepSize = [&]() { return tim->img->pixDims[0] * 0.05; }();

        auto G =  (Gaussian3D_evaluator_4squaredDist*)(std::get<0>(*(std::tuple<void*,Contributions*>*)(tim->data)));

        tim->img->to_xyz(gridPos, &voxCenter[0]);

        int fullSteps = static_cast<int>(seg.length / stepSize);
        float rem     = seg.length - fullSteps * stepSize;
        
        for (int i = 0; i < 3; i++) {
            step[i] = stepSize * seg.dir[i];
            p[i]    = seg.p[i];
        }

        float lineIntegral = 0.0;

        // Handle full steps
        for (int i = 0; i < fullSteps; ++i) {
            lineIntegral += G->eval(squared_dist(voxCenter, p)) * stepSize;
            p[0] += step[0];
            p[1] += step[1];
            p[2] += step[2];
        }

        // Handle the remaining length
        lineIntegral += G->eval(squared_dist(voxCenter, p)) * rem;

        {
            std::lock_guard<std::mutex> lock( (tim->useMutexGrid) ? tim->mutexGrid[ind] : tim->mutexMap[ind]);

            if (tim->grid[ind]==NULL) { // Allocate memory here
                std::vector<SegmentContribution>* streamlineIdsAndGaussionWeightedLength = new std::vector<SegmentContribution>();
                tim->grid[ind] = ((void*)streamlineIdsAndGaussionWeightedLength); // Each thread handles a single std::vector<SegmentContribution> for each voxel in the output
            }

            ((std::vector<SegmentContribution>*)(tim->grid[ind]))->push_back(std::make_pair(seg.streamlineNo,lineIntegral));
        }

    }

    template<typename T>
    inline void indexStreamlineIdAndGaussionWeightedLength_memVector(Tractogram2ImageMapper<T>* tim) {

        auto contributions = (Contributions*)(std::get<1>(*(std::tuple<void*,Contributions*>*)(tim->data)));

        std::mutex busyContrib;

        auto genOut = [&](const NIBR::MT::TASK& task)->void {

            int64_t ind = task.no;

            if (tim->grid[ind]!=NULL) {

                std::vector<SegmentContribution> tmp;
                std::swap(tmp,*(std::vector<SegmentContribution>*)(tim->grid[ind]));

                {
                    std::lock_guard lock(busyContrib);
                    contributions->emplace_back(std::make_pair(ind,std::move(tmp)));
                }

                delete ((std::vector<SegmentContribution>*)(tim->grid[ind]));
                tim->grid[ind] = NULL;

            }

        };

        NIBR::MT::MTRUN(tim->img->voxCnt,genOut);

        // Sort the contributions in ascending order based on "ind"
        std::sort(contributions->begin(), contributions->end(), 
            [](const ContributionsInAVoxel& a, const ContributionsInAVoxel& b) {
                return a.first < b.first;
            }
        );

        return;
    }

}
