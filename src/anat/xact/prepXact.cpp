#include "prepXact.h"
#include "surface/surface_make.h"

using namespace NIBR;

#define FS_BRAINSTEM    16

std::vector<Surface> NIBR::prepXact(
    std::string fsPath,
    std::string fslFirstFolder, 
    Surface* abnormality,
    float cereDistThresh, 
    float enlargeBrainStem, 
    float inferiorBrainStemCutLevel, 
    float meanFaceArea, 
    const XactPrepOption& opt
)
{
    std::vector<Surface> out(14,Surface());

    std::vector<std::string> labelNames = {"combined","l_wm","r_wm","l_gm","r_gm","l_sub","r_sub","csf","cer_wm","cer_gm","bs","i_bs","abn","bg"};
    std::vector<int> labels;
    for (int l = 0 ; l < int(labelNames.size()); l++) labels.push_back(l);

    XactPrepOption optMap[] = {XACT_PREP_OPT_COMBINED, XACT_PREP_OPT_L_WM, XACT_PREP_OPT_R_WM, XACT_PREP_OPT_L_GM, XACT_PREP_OPT_R_GM, XACT_PREP_OPT_L_SUB, XACT_PREP_OPT_R_SUB, XACT_PREP_OPT_CSF, XACT_PREP_OPT_CER_WM, XACT_PREP_OPT_CER_GM, XACT_PREP_OPT_BS, XACT_PREP_OPT_I_BS, XACT_PREP_OPT_ABN, XACT_PREP_OPT_BG};
    
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

        if ( (( i == XactLabel::L_WM) || ( i == XactLabel::R_WM) || ( i == XactLabel::L_GM) || ( i == XactLabel::R_GM) || ((i==XactLabel::L_SUB)&&(fslFirstFolder!="")) || ((i==XactLabel::R_SUB)&&(fslFirstFolder!="")) || (i == XactLabel::I_BS) || (i == XactLabel::BS)) && (meanFaceArea != 0)) {
            out[i].calcArea();
            out[i] = surfRemesh(out[i],out[i].area/meanFaceArea*0.5f,1,0);
        }

        if ((i != XactLabel::L_SUB) && (i != XactLabel::R_SUB) && (i != XactLabel::CSF) && (i != XactLabel::BS) && (i != XactLabel::BG)) {
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
        
    // XactLabel::L_WM, XactLabel::R_WM, XactLabel::L_GM, XactLabel::R_GM
    std::vector<Surface> cortex;
    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_L_WM | XACT_PREP_OPT_R_WM | XACT_PREP_OPT_L_GM | XACT_PREP_OPT_R_GM | XACT_PREP_OPT_BS | XACT_PREP_OPT_I_BS | XACT_PREP_OPT_BG)) {
        cortex = fsReconall2CorticalSurf(fsPath); // {closed_wm_surf,closed_pial_surf}, field[0] is side, field[1] is FS label
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_L_WM)) {
        selectVertices(&out[XactLabel::L_WM], &cortex[0], &cortex[0].fields[0], 1); // fields[0] = side label 1 = left
        finalize(XactLabel::L_WM); 
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_R_WM)) {
        selectVertices(&out[XactLabel::R_WM], &cortex[0], &cortex[0].fields[0], 2); // fields[0] = side label 2 = right
        finalize(XactLabel::R_WM); 
    }

    
    auto makeGM = [&](int side)->Surface {
        Surface open_WM, open_GM;

        selectVertices(&open_WM, &cortex[0], &cortex[0].fields[1]);         // fields[1] = fs label
        selectVertices(&open_GM, &cortex[1], &cortex[1].fields[1]);         // fields[1] = fs label

        Surface open_WM_side, open_GM_side;
        selectVertices(&open_WM_side, &open_WM, &open_WM.fields[0], side);  // fields[0] = side label (1=left, 2=right)
        selectVertices(&open_GM_side, &open_GM, &open_GM.fields[0], side);  // fields[0] = side label (1=left, 2=right)

        open_WM_side.flipNormalsOfFaces();
        return surfGlueBoundaries(open_GM_side,open_WM_side);
    };

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_L_GM)) {
        out[XactLabel::L_GM] = makeGM(1); 
        finalize(XactLabel::L_GM); 
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_R_GM)) {
        out[XactLabel::R_GM] = makeGM(2); 
        finalize(XactLabel::R_GM); 
    }

    
    // XactLabel::L_SUB, XactLabel::R_SUB
    Surface sub;
    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_L_SUB | XACT_PREP_OPT_R_SUB)) {
        sub = (fslFirstFolder=="") ? fsAseg2SubcortexSurf(asegImg, meanFaceArea, true) : fslFirst2SubcortexSurf(fslFirstFolder, meanFaceArea, true);
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_L_SUB)) {
        selectVertices(&out[XactLabel::L_SUB], &sub, &sub.fields[0], 1); // fields[0] = side label 1 = left
        finalize(XactLabel::L_SUB);
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_R_SUB)) {
        selectVertices(&out[XactLabel::R_SUB], &sub, &sub.fields[0], 2); // fields[0] = side label 2 = right
        finalize(XactLabel::R_SUB);
    }


    // XactLabel::CSF
    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_CSF)) {
        out[XactLabel::CSF] = fsAseg2CSFSurf(asegImg, meanFaceArea);
        finalize(XactLabel::CSF);
    }
    
    // XactLabel::CER_WM, XactLabel::CER_GM
    std::vector<Surface> cerebellum;

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_CER_WM | XACT_PREP_OPT_CER_GM)) {
        cerebellum = fsAseg2CerebellumSurf(asegImg, cereDistThresh, meanFaceArea); // Open and closed cerebellar surfaces: {cer_mask_closed,cer_mask_open,cer_wm_closed,cer_wm_open}        
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_CER_WM)) {
        out[XactLabel::CER_WM] = cerebellum[2];
        finalize(XactLabel::CER_WM);
    }

    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_CER_GM)) {
        cerebellum[3].flipNormalsOfFaces();
        out[XactLabel::CER_GM] = surfGlueBoundaries(cerebellum[1],cerebellum[3]);
        finalize(XactLabel::CER_GM);
    }
 

    // XactLabel::BS and XactLabel::I_BS
    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_BS | XACT_PREP_OPT_I_BS | XACT_PREP_OPT_BG)) {            

        // We first get the complete brainstem segmentation, we will split the lower part and label it as XactLabel::I_BS
        out[XactLabel::BS] = (fslFirstFolder=="") ? fsAseg2Surf(asegImg, FS_BRAINSTEM, 0) : fslFirst2BrainstemSurf(fslFirstFolder, 0);

        out[XactLabel::BS] = surfGrow(out[XactLabel::BS],1,enlargeBrainStem);
        
        // Compute the bounding box and find the lowest point (toward spine)
        auto bsBox = surfBbox(out[XactLabel::BS]);
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
        
        // We next cut the brainstem segmentation by inferiorBrainStemCutLevel mm above this level.
        float cutLevel[3] = {planePoint1[0],planePoint1[1],planePoint1[2]};
        cutLevel[ax] = (lowPoint == 1) ? (planePoint1[ax]+inferiorBrainStemCutLevel) : (planePoint2[ax]-inferiorBrainStemCutLevel);

        disp(MSG_DEBUG,"Applying cut level [%.2f %.2f %.2f] at axis %d", cutLevel[0],cutLevel[1],cutLevel[2],ax);

        auto [upperCut, lowerCut] = splitWithPlane(out[XactLabel::BS], &cutLevel[0], &planeDir[0]);
    
        upperCut = surfMakeItSingleClosed(upperCut);
        lowerCut = surfMakeItSingleClosed(lowerCut);

        upperCut.calcVolume();
        lowerCut.calcVolume();

        // out[XactLabel::BS]    = (upperCut.volume > lowerCut.volume) ? upperCut : lowerCut;
        out[XactLabel::I_BS] = (upperCut.volume > lowerCut.volume) ? lowerCut : upperCut;

        finalize(XactLabel::BS);
        finalize(XactLabel::I_BS);
    }

    // XactLabel::ABN
    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_ABN)) {
        if (abnormality != nullptr) {
            out[XactLabel::ABN] = *abnormality;
            finalize(XactLabel::ABN);
        } else {
            disp(MSG_WARN, "Abnormality surface pointer is null, skipping abnormality surface generation.");
        }
    }

    // XactLabel::BG
    if (opt & (XACT_PREP_OPT_COMBINED | XACT_PREP_OPT_BG)) {
        out[XactLabel::BG] = fsAseg2BgSurf(asegImg,out[XactLabel::BS],meanFaceArea);
        finalize(XactLabel::BG);
    }

    // Create combined surface
    std::vector<int> combinedLabelField;

    for (size_t i = 1; i < out.size(); i++) {

        if ( (out[i].nv > 0) && ((opt & optMap[i]) || (opt & XACT_PREP_OPT_COMBINED))){
            
            std::vector<int> labelField(out[i].nv,labels[i]);

            if (opt & XACT_PREP_OPT_COMBINED) {
                out[XactLabel::COMBINED] = surfMerge(out[XactLabel::COMBINED],out[i]);
                combinedLabelField.insert(combinedLabelField.end(),labelField.begin(),labelField.end());
            }

            out[i].fields.push_back(out[i].makeVertField("xact",labelField));

        }

    }


    if (opt & XACT_PREP_OPT_COMBINED) {
        disp(MSG_INFO,"Finalizing the combined surface");
        out[XactLabel::COMBINED].fields.push_back(out[XactLabel::COMBINED].makeVertField("xact",combinedLabelField));
        disp(MSG_INFO,"Done");
    }
    

    return out;

}