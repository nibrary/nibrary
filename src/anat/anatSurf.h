#pragma once

#include "base/nibr.h"
#include "image/image_marchingCubes.h"

enum XactOutputOption {
    XACT_OPT_UNSET    = 0,
    XACT_OPT_COMBINED = 1,
    XACT_OPT_L_WM     = 2,
    XACT_OPT_R_WM     = 4,
    XACT_OPT_L_GM     = 8,
    XACT_OPT_R_GM     = 16,
    XACT_OPT_L_SUB    = 32,
    XACT_OPT_R_SUB    = 64,
    XACT_OPT_CSF      = 128,
    XACT_OPT_CER_WM   = 256,
    XACT_OPT_CER_GM   = 512,
    XACT_OPT_BS       = 1024,
    XACT_OPT_SPINE    = 2048,
    XACT_OPT_CRAN     = 4096
};

namespace NIBR {

Surface fsAseg2Surf(std::string fsAsegFile, int asegLabel, float meanFaceArea);
Surface fsAseg2Surf(Image<int>& fsAsegImg,  int asegLabel, float meanFaceArea);

Surface fsAseg2CranSurf(std::string fsAsegFile, std::string fslFirstFolder, float meanFaceArea);
Surface fsAseg2CranSurf(Image<int>& fsAsegImg,  Surface& brainStem, float meanFaceArea);

// {cer_mask_closed,cer_mask_open,cer_wm_closed,cer_wm_open}
std::vector<Surface> fsAseg2CerebellumSurf(std::string fsAsegFile, float distThresh, float meanFaceArea);
std::vector<Surface> fsAseg2CerebellumSurf(Image<int>& fsAsegImg,  float distThresh, float meanFaceArea);

Surface fsAseg2CSFSurf(std::string fsAsegFile, float meanFaceArea);
Surface fsAseg2CSFSurf(Image<int>& fsAsegImg,  float meanFaceArea);

Surface fsAseg2SubcortexSurf(std::string fsAsegFile, float meanFaceArea, bool excludeBrainStem);
Surface fsAseg2SubcortexSurf(Image<int>& fsAsegImg,  float meanFaceArea, bool excludeBrainStem);

// {closed_wm_surf,closed_pial_surf}, field[0] is side, field[1] is FS label
std::vector<Surface> fsReconall2CorticalSurf(std::string fsPath);

// closed_brain_mask_surf, optionally makes a thin shell of brain surface, where inside is hallow
Surface fsReconall2BrainSurf(std::string brainMaskFile, float meanFaceArea, float shellThickness);

Surface fslFirst2SubcortexSurf(std::string fslFirstFolder, float meanFaceArea, bool excludeBrainStem);
Surface fslFirst2BrainstemSurf(std::string fslFirstFolder, float meanFaceArea);

}

