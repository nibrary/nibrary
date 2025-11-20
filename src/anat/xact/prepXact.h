#pragma once

#include "anat/anatSurf.h"


namespace NIBR {

// cereDistThresh: The distance threshold from the brain stem to separate cerebellar white matter [mm]
// enlargeBrainStem: Enlarge or shrink brain stem [mm]
std::vector<Surface> prepXact(std::string fsPath,std::string fslFirstFolder, float cereDistThresh, float enlargeBrainStem, float meanFaceArea, const XactOutputOption& opt);

}
