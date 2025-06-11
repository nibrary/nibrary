#pragma once

#include "base/nibr.h"
#include "dMRI/tractography/io/tractogramReader.h"
#include "surface/surface.h"

namespace NIBR
{

    // Output size is of length streamlineIndices.size(). Each output element contains streamline ids, and their corresponding vertex_ids and streamline distances on the mesh
    // The output of this function is not memory efficient; std::vector<std::vector<std::pair<int, float>>> is a better choice here.
    // The below function is only used for the derivative demo
    std::vector<std::unordered_map<int,float>> index2surface(TractogramReader& tractogram, Surface& surf, float sigma); 
    
    void index2surface(TractogramReader& tractogram, Surface& surf, float sigma, std::string indexFilePrefix);
    
    std::vector<std::vector<std::pair<int, float>>>         readSurfaceIndexing(TractogramReader& tractogram, const std::string& indexFilePrefix);
    std::unordered_map<int,std::unordered_map<int, float>>  readSurfaceIndexingForSelection(const std::vector<int>& streamlineIndices, const std::string& indexFilePrefix, const std::vector<std::streampos>& positions);

}
