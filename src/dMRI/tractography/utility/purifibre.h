#pragma once

#include "base/nibr.h"
#include "dMRI/tractography/io/tractogramReader.h"
#include "dMRI/tractography/io/tractogramField.h"

namespace NIBR 
{

    // Returns FICO values for each streamline in the tractogram
    std::vector<float> getFico(TractogramReader* tractogram, 
                               float trimFactor = 10.0f, 
                               float voxDim = 4.0f, 
                               std::tuple<float,int> anisotropicSmoothing = std::make_tuple(0.0f,0), 
                               float sphericalSmoothing = 15.0f);

    // Returns the indices of streamlines to keep after purification based on FICO values
    std::vector<size_t> purify(const std::vector<float>& fico, float puriFactor);


}