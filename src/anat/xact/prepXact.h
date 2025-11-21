#pragma once

#include "anat/anatSurf.h"


namespace NIBR {

enum XactPrepOption {
    XACT_PREP_OPT_UNSET    = 0,
    XACT_PREP_OPT_COMBINED = 1,
    XACT_PREP_OPT_L_WM     = 2,
    XACT_PREP_OPT_R_WM     = 4,
    XACT_PREP_OPT_L_GM     = 8,
    XACT_PREP_OPT_R_GM     = 16,
    XACT_PREP_OPT_L_SUB    = 32,
    XACT_PREP_OPT_R_SUB    = 64,
    XACT_PREP_OPT_CSF      = 128,
    XACT_PREP_OPT_CER_WM   = 256,
    XACT_PREP_OPT_CER_GM   = 512,
    XACT_PREP_OPT_BS       = 1024,
    XACT_PREP_OPT_I_BS     = 2048,
    XACT_PREP_OPT_BG       = 4096
};

// cereDistThresh: The distance threshold from the brain stem to separate cerebellar white matter [mm]
// enlargeBrainStem: Enlarge or shrink brain stem [mm]
// inferiorBrainStemCutLevel: The cut level to separate inferior brain stem [mm]
std::vector<Surface> prepXact(
    std::string fsPath,
    std::string fslFirstFolder, 
    float cereDistThresh, 
    float enlargeBrainStem, 
    float inferiorBrainStemCutLevel, 
    float meanFaceArea, 
    const XactPrepOption& opt
);

}
