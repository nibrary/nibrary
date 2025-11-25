#pragma once

#include "anat/anatSurf.h"

namespace NIBR {

enum XactLabel {
    COMBINED = 0,
    L_WM     = 1,  // Left white matter
    R_WM     = 2,  // Right white matter
    L_GM     = 3,  // Left gray matter
    R_GM     = 4,  // Right gray matter
    L_SUB    = 5,  // Left subcortical
    R_SUB    = 6,  // Right subcortical
    CSF      = 7,  // Cerebrospinal fluid
    CER_WM   = 8,  // Cerebellar white matter
    CER_GM   = 9,  // Cerebellar gray matter
    BS       = 10, // Brain stem
    I_BS     = 11, // Inferior brain stem
    ABN      = 12, // Abnormality
    BG       = 13  // Background
};

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
    XACT_PREP_OPT_ABN      = 4096,
    XACT_PREP_OPT_BG       = 8192
};

enum XactTrackOption {
    XACT_TRACK_OPT_UNSET                 = 0,
    XACT_TRACK_OPT_SEED_GM               = 1,
    XACT_TRACK_OPT_STOP_AFTER_ENTRY_GM   = 2,
    XACT_TRACK_OPT_STOP_BEFORE_EXIT_GM   = 4,
    XACT_TRACK_OPT_SEED_SUB              = 8,
    XACT_TRACK_OPT_STOP_AFTER_ENTRY_SUB  = 16,
    XACT_TRACK_OPT_STOP_BEFORE_EXIT_SUB  = 32,
    XACT_TRACK_OPT_SEED_ABN              = 64,
    XACT_TRACK_OPT_STOP_AFTER_ENTRY_ABN  = 128,
    XACT_TRACK_OPT_STOP_BEFORE_EXIT_ABN  = 256,
    XACT_TRACK_OPT_SEED_BG               = 512,
    XACT_TRACK_OPT_STOP_AFTER_ENTRY_BG   = 1024,
    XACT_TRACK_OPT_STOP_BEFORE_EXIT_BG   = 2048,
    XACT_INTRACORTICAL                   = 4096,
    XACT_CRANIAL                         = 8192,
    XACT_SUBCORTICAL_DEADEND             = 16384,
    XACT_ABNORMALITY_DEADEND             = 32768
};

// cereDistThresh: The distance threshold from the brain stem to separate cerebellar white matter [mm]
// enlargeBrainStem: Enlarge or shrink brain stem [mm]
// inferiorBrainStemCutLevel: The cut level to separate inferior brain stem [mm]
std::vector<Surface> prepXact(
    std::string fsPath,
    std::string fslFirstFolder,
    Surface* abnormality, 
    float cereDistThresh, 
    float enlargeBrainStem, 
    float inferiorBrainStemCutLevel, 
    float meanFaceArea, 
    const XactPrepOption& opt
);

}
