#pragma once

#include "base/nibr.h"
#include "dMRI/tractography/io/tractogramReader.h"
#include "dMRI/tractography/io/tractogramField.h"

namespace NIBR 
{
    
    // Compute bounding box of a tractogram
    std::vector<float>  getTractogramBBox(NIBR::TractogramReader* tractogram);

    // Returns bool vector marking true for streamlines which are in the inpBatch but not in the refBatch
    std::vector<bool>   tractogramDiff(const StreamlineBatch& inpBatch, const StreamlineBatch& refBatch);

    StreamlineBatch     tractogramTransform(const StreamlineBatch& batch_in, float M[][4]);

    TractogramField     colorTractogram(NIBR::TractogramReader* tractogram);

    // Checks for corrupted streamlines and returns streamline lengths, number of points, mean step size and standard deviation of step size
    std::tuple<bool,std::vector<float>, std::size_t, float, float> getTractogramStats(NIBR::TractogramReader* tractogram);

}

