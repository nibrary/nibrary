#include "prepXact.h"
#include "surface/surface_make.h"

using namespace NIBR;

// XACT field labels
//
// Note: WM covers SUB
//
// 1. seed:              L_WM + R_WM + CER_WM + BS
// 2. discard_seed:      L_GM + R_GM + CER_GM + SPINE + CRAN
// 3. req_end_inside:    L_GM + R_GM + L_SUB + R_SUB  + CER_GM + BS + SPINE + CRAN
// 4. discard_if_enters: CSF
// 5. stop_before_exit:  CRAN
// 6. stop_after_entry:  L_GM + R_GM + CER_GM + SPINE
// 7. stop_before_exit:  L_GM + R_GM + CER_GM + SPINE
// 8. stop_before_exit:  L_SUB + R_SUB
//
//
// Standard settings
//
// -s ${xact} 1
// -d ${xact} 2
// -p req_end_inside_A ${xact} 3
// -p req_end_inside_B ${xact} 3
// -p discard_if_enters_A ${xact} 4
// -p discard_if_enters_B ${xact} 4
// -p stop_before_exit_A ${xact} 5
// -p stop_before_exit_B ${xact} 5
//
//
// For cortical ending, choose 1 or 2
//
// Cortical ending option 1
// -p stop_after_exit_A ${xact} 6
// -p stop_after_exit_A ${xact} 6
// 
// Cortical ending option 2
// -p stop_before_exit_A ${xact} 7
// -p stop_before_exit_B ${xact} 7
//
//
// For optional subcortical ending, stop only one side. The other side can leave the subcortex or not, so side B is not set
//
// -p stop_before_exit_A ${xact} 8
//


// Labels
#define COMBINED   0
#define L_WM       1
#define R_WM       2
#define L_GM       3
#define R_GM       4
#define L_SUB      5
#define R_SUB      6
#define CSF        7
#define CER_WM     8
#define CER_GM     9
#define BS        10
#define SPINE     11
#define CRAN      12

#define FS_BRAINSTEM    16

std::vector<Surface> NIBR::prepXact(std::string fsPath,std::string fslFirstFolder, float cereDistThresh, float enlargeBrainStem, float meanFaceArea, const XactOutputOption& opt)
{
    std::vector<Surface> out(13,Surface());

    std::vector<std::string> labelNames = {"combined","l_wm","r_wm","l_gm","r_gm","l_sub","r_sub","csf","cer_wm","cer_gm","bs","spine","cran"};
    std::vector<int> labels;
    for (int l = 0 ; l < int(labelNames.size()); l++) labels.push_back(l);

    XactOutputOption optMap[] = {XACT_OPT_COMBINED, XACT_OPT_L_WM, XACT_OPT_R_WM, XACT_OPT_L_GM, XACT_OPT_R_GM, XACT_OPT_L_SUB, XACT_OPT_R_SUB, XACT_OPT_CSF, XACT_OPT_CER_WM, XACT_OPT_CER_GM, XACT_OPT_BS, XACT_OPT_SPINE, XACT_OPT_CRAN};

    // Read aseg file
    if (!existsFile(fsPath + "/mri/aparc+aseg.mgz")) {
        disp(MSG_ERROR, "%s was not found", (fsPath + "/mri/aparc+aseg.mgz").c_str());
        return out;
    }

    Image<int> asegImg(fsPath + "/mri/aparc+aseg.mgz");
    asegImg.read();

    // Convenience function to finalize surfaces and verify their closedness
    auto finalize = [&](int i)->bool{
        
        disp(MSG_INFO,"Finalizing and verifying %s",labelNames[i].c_str());

        for (size_t f = 0; f < out[i].fields.size(); f++) {
            out[i].clearField(out[i].fields[f]);
        }
        out[i].fields.clear();

        if ( (( i == L_WM) || ( i == R_WM) || ( i == L_GM) || ( i == R_GM) || ((i==L_SUB)&&(fslFirstFolder!="")) || ((i==R_SUB)&&(fslFirstFolder!="")) || (i == SPINE) || (i == BS)) && (meanFaceArea != 0)) {
            out[i].calcArea();
            out[i] = surfRemesh(out[i],out[i].area/meanFaceArea*0.5f,1,0);
        }

        if ((i != L_SUB) && (i != R_SUB) && (i != CSF) && (i != BS) && (i != CRAN)) {
            out[i] = surfMakeItSingleClosed(out[i]);
        }

        out[i] = surfMakeItWatertight(out[i]);

        disp(MSG_DEBUG,"Finalized %s",labelNames[i].c_str());

        if (out[i].nv == 0) {
            disp(MSG_ERROR,"Failed at generating closed %s surface",labelNames[i].c_str());
            return false;
        }

        return true;
        
    };
        
    // L_WM, R_WM, L_GM, R_GM
    std::vector<Surface> cortex;
    if (opt & (XACT_OPT_COMBINED | XACT_OPT_L_WM | XACT_OPT_R_WM | XACT_OPT_L_GM | XACT_OPT_R_GM | XACT_OPT_BS | XACT_OPT_SPINE | XACT_OPT_CRAN)) {
        cortex = fsReconall2CorticalSurf(fsPath); // {closed_wm_surf,closed_pial_surf}, field[0] is side, field[1] is FS label
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_L_WM)) {
        selectVertices(&out[L_WM], &cortex[0], &cortex[0].fields[0], 1); // fields[0] = side label 1 = left
        finalize(L_WM); 
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_R_WM)) {
        selectVertices(&out[R_WM], &cortex[0], &cortex[0].fields[0], 2); // fields[0] = side label 2 = right
        finalize(R_WM); 
    }

    
    auto makeGM = [&](int side)->Surface {
        Surface open_WM, open_GM;

        selectVertices(&open_WM, &cortex[0], &cortex[0].fields[1]);         // fields[1] = fs label
        selectVertices(&open_GM, &cortex[1], &cortex[1].fields[1]);         // fields[1] = fs label

        Surface open_WM_side, open_GM_side;
        selectVertices(&open_WM_side, &open_WM, &open_WM.fields[0], side);  // fields[0] = side label 1 = left
        selectVertices(&open_GM_side, &open_GM, &open_GM.fields[0], side);  // fields[0] = side label 2 = right

        open_WM_side.flipNormalsOfFaces();
        return surfGlueBoundaries(open_GM_side,open_WM_side);
    };

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_L_GM)) {
        out[L_GM] = makeGM(1); 
        finalize(L_GM); 
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_R_GM)) {
        out[R_GM] = makeGM(2); 
        finalize(R_GM); 
    }

    
    // L_SUB, R_SUB
    Surface sub;
    if (opt & (XACT_OPT_COMBINED | XACT_OPT_L_SUB | XACT_OPT_R_SUB)) {
        sub = (fslFirstFolder=="") ? fsAseg2SubcortexSurf(asegImg, meanFaceArea, true) : fslFirst2SubcortexSurf(fslFirstFolder, meanFaceArea, true);
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_L_SUB)) {
        selectVertices(&out[L_SUB], &sub, &sub.fields[0], 1); // fields[0] = side label 1 = left
        finalize(L_SUB);
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_R_SUB)) {
        selectVertices(&out[R_SUB], &sub, &sub.fields[0], 2); // fields[0] = side label 2 = right
        finalize(R_SUB);
    }


    // CSF
    if (opt & (XACT_OPT_COMBINED | XACT_OPT_CSF)) {
        out[CSF] = fsAseg2CSFSurf(asegImg, meanFaceArea);
        finalize(CSF);
    }
    
    // CER_WM, CER_GM
    std::vector<Surface> cerebellum;

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_CER_WM | XACT_OPT_CER_GM)) {
        cerebellum = fsAseg2CerebellumSurf(asegImg, cereDistThresh, meanFaceArea); // Open and closed cerebellar surfaces: {cer_mask_closed,cer_mask_open,cer_wm_closed,cer_wm_open}        
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_CER_WM)) {
        out[CER_WM] = cerebellum[2];
        finalize(CER_WM);
    }

    if (opt & (XACT_OPT_COMBINED | XACT_OPT_CER_GM)) {
        cerebellum[3].flipNormalsOfFaces();
        out[CER_GM] = surfGlueBoundaries(cerebellum[1],cerebellum[3]);
        finalize(CER_GM);
    }
 

    // BS and SPINE
    if (opt & (XACT_OPT_COMBINED | XACT_OPT_BS | XACT_OPT_SPINE | XACT_OPT_CRAN)) {            

        // We first get the complete brainstem segmentation, we will split the lower part and label it as SPINE
        out[BS] = (fslFirstFolder=="") ? fsAseg2Surf(asegImg, FS_BRAINSTEM, 0) : fslFirst2BrainstemSurf(fslFirstFolder, 0);

        out[BS] = surfGrow(out[BS],1,enlargeBrainStem);
        
        // Compute the bounding box and find the lowest point (toward spine)
        auto bsBox = surfBbox(out[BS]);
        float xLen = bsBox[1] - bsBox[0];
        float yLen = bsBox[3] - bsBox[2];
        float zLen = bsBox[5] - bsBox[4];
        int ax;
        float planeDir[3]   = {0,0,0};
        float planePoint1[3], planePoint2[3];
        if      ((xLen > yLen) && (xLen > zLen)) {ax = 0; planePoint1[0] = bsBox[0];           planePoint1[1] = bsBox[2]+yLen*0.5f; planePoint1[2] = bsBox[4]+zLen*0.5f;}
        else if ((yLen > xLen) && (yLen > zLen)) {ax = 1; planePoint1[0] = bsBox[0]+xLen*0.5f; planePoint1[1] = bsBox[2];           planePoint1[2] = bsBox[4]+zLen*0.5f;}
        else                                     {ax = 2; planePoint1[0] = bsBox[0]+xLen*0.5f; planePoint1[1] = bsBox[2]+yLen*0.5f; planePoint1[2] = bsBox[4]; }

        planeDir[ax] = 1.0f;

        planePoint2[0] = planePoint1[0];
        planePoint2[1] = planePoint1[1];
        planePoint2[2] = planePoint1[2];

        planePoint2[ax] += (bsBox[2*ax+1] - bsBox[2*ax]);

        auto cutLength = [&](float* p)->float {
            auto [cutDown, cutUp] = splitWithPlane(cortex[0], p, &planeDir[0]);
            if (cutDown.nv == 0) return 0;
            cutDown.computeBoundaries();
            if (cutDown.boundaries.empty()) return 0;
            cutDown = surfMakeItSingleOpen(cutDown);
            return cutDown.boundaryLengths[0];
        };

        int lowPoint = (cutLength(&planePoint1[0]) < cutLength(&planePoint2[0]) ) ? 1 : 2;
        
        // We next cut the brainstem segmentation 20 mm above this level.
        float cutLevel[3] = {planePoint1[0],planePoint1[1],planePoint1[2]};
        cutLevel[ax] = (lowPoint == 1) ? (planePoint1[ax]+20) : (planePoint2[ax]-20);

        disp(MSG_DEBUG,"Applying cut level [%.2f %.2f %.2f] at axis %d", cutLevel[0],cutLevel[1],cutLevel[2],ax);

        auto [upperCut, lowerCut] = splitWithPlane(out[BS], &cutLevel[0], &planeDir[0]);
    
        upperCut = surfMakeItSingleClosed(upperCut);
        lowerCut = surfMakeItSingleClosed(lowerCut);

        upperCut.calcVolume();
        lowerCut.calcVolume();

        // out[BS]    = (upperCut.volume > lowerCut.volume) ? upperCut : lowerCut;
        out[SPINE] = (upperCut.volume > lowerCut.volume) ? lowerCut : upperCut;

        finalize(BS);
        finalize(SPINE);
    }

    // CRAN
    if (opt & (XACT_OPT_COMBINED | XACT_OPT_CRAN)) {
        out[CRAN] = fsAseg2CranSurf(asegImg,out[BS],meanFaceArea);
        finalize(CRAN);
    }

    // Create combined surface
    std::vector<int> combinedLabelField;
    std::vector<int> xactField;

    for (size_t i = 1; i < out.size(); i++) {

        if ( (out[i].nv > 0) && ((opt & optMap[i]) || (opt & XACT_OPT_COMBINED))){
            
            std::vector<int> labelField(out[i].nv,labels[i]);

            if (opt & XACT_OPT_COMBINED) {
                out[COMBINED] = surfMerge(out[COMBINED],out[i]);
                combinedLabelField.insert(combinedLabelField.end(),labelField.begin(),labelField.end());
            }

            out[i].fields.push_back(out[i].makeVertField("xact",labelField));

        }

    }


    if (opt & XACT_OPT_COMBINED) {
        disp(MSG_INFO,"Finalizing the combined surface");
        out[COMBINED].fields.push_back(out[COMBINED].makeVertField("xact",combinedLabelField));
        disp(MSG_INFO,"Done");
    }
    

    return out;

}