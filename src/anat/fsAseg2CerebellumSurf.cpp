#include "anatSurf.h"

#define RIGHT_CER_WM 46
#define RIGHT_CER_GM 47
#define LEFT_CER_WM   7
#define LEFT_CER_GM   8
#define BRAINSTEM    16

#define DEF_D_THRSH 0.5

using namespace NIBR;

std::vector<Surface> NIBR::fsAseg2CerebellumSurf(Image<int>& asegImg, float distThresh, float meanFaceArea)
{

    // Create cerebellum white matter, cerebellum gray matter, and mask images
    Image<int8_t> cer_wm_r; imgThresh(cer_wm_r, asegImg, RIGHT_CER_WM, RIGHT_CER_WM);
    Image<int8_t> cer_wm_l; imgThresh(cer_wm_l, asegImg,  LEFT_CER_WM,  LEFT_CER_WM);
    Image<int8_t> cer_gm_r; imgThresh(cer_gm_r, asegImg, RIGHT_CER_GM, RIGHT_CER_GM);
    Image<int8_t> cer_gm_l; imgThresh(cer_gm_l, asegImg,  LEFT_CER_GM,  LEFT_CER_GM);

    Image<float>    cer_wm; imgAdd(cer_wm,cer_wm_r,cer_wm_l);                 // Cerebellum WM
    Image<float>    cer_gm; imgAdd(cer_gm,cer_gm_r,cer_gm_l);                 // Cerebellum GM
    Image<float>  cer_mask; imgAdd(cer_mask,cer_gm,cer_wm);                   // Whole cerebellum mask (WM + GM)
    Image<float> brainstem; imgThresh(brainstem, asegImg, BRAINSTEM, BRAINSTEM);  // Brainstem

    // Create closed surface meshes
    Surface cer_wm_closed  = label2surface(cer_wm,1,((meanFaceArea==0) ? 0.25 : meanFaceArea));
    if (cer_wm_closed.nv == 0) {return std::vector<Surface>();}
    cer_wm_closed.prepIglAABBTree();

    Surface cer_mask_closed = label2surface(cer_mask,1,((meanFaceArea==0) ? 0.25 : meanFaceArea));
    if (cer_mask_closed.nv == 0) {return std::vector<Surface>();}
    cer_mask_closed.prepIglAABBTree();

    Surface brainstem_closed = label2surface(brainstem,1,((meanFaceArea==0) ? 0.25 : meanFaceArea));
    if (brainstem_closed.nv == 0) {return std::vector<Surface>();}
    brainstem_closed.prepIglAABBTree();

    
    // Finding openings
    distThresh = (distThresh==0) ? DEF_D_THRSH*DEF_D_THRSH : distThresh*distThresh;

    std::vector<bool> toKeep_cer_wm;
    for (int n = 0; n < cer_wm_closed.nv; n++) {
        toKeep_cer_wm.push_back(brainstem_closed.squaredDistToPoint(cer_wm_closed.vertices[n]) > distThresh);
    }    
    
    Surface cer_wm_open = applyMask(cer_wm_closed,toKeep_cer_wm);
    if (cer_wm_open.nv == 0) {
        disp(MSG_ERROR, "Failed to generated an open cerebellum white matter surface");
        return std::vector<Surface>();
    }

    cer_wm_open = surfRemoveAllButLargestComponent(cer_wm_open);

    Surface wm_opening = surfDiff(cer_wm_closed, cer_wm_open);

    auto holes  = wm_opening.calcAreasOfConnectedComponents();
    std::sort(holes.begin(), holes.end(), std::greater<double>());
    if (holes.size() == 2) { // If there are two holes, and one is very big, then remove the other
        if (holes[1] < (holes[0]/10.0)) {
            wm_opening = surfRemoveAllButLargestComponent(wm_opening);
        }
    } else if (holes.size() > 2) { // If there are more two holes, and remove those smaller than average
        double aveh = 0; // Average hole size
        for (auto h : holes) {
            aveh += h;
        }
        aveh /= double(holes.size());
        wm_opening = surfRemoveSmallConnectedComponents(wm_opening,aveh);
    }

    cer_wm_open = surfDiff(cer_wm_closed, wm_opening);
    
    wm_opening.prepIglAABBTree();

     
    std::vector<bool> toKeep_cer_mask;
    for (int n = 0; n < cer_mask_closed.nv; n++) {
        toKeep_cer_mask.push_back( wm_opening.squaredDistToPoint(cer_mask_closed.vertices[n]) > distThresh);
    }
    
    Surface cer_mask_open = applyMask(cer_mask_closed,toKeep_cer_mask);
    if (cer_mask_open.nv == 0) {
        disp(MSG_ERROR, "Failed to generated an open cerebellum surface");
        return std::vector<Surface>();
    }

    Surface mask_opening = surfDiff(cer_mask_closed, cer_mask_open);

    // Create opening labels
    auto createOpeningLabel = [&](Surface& inp, Surface& opening)->void {

        SurfaceField f;
        f.owner     = VERTEX;
        f.datatype  = "int";
        f.dimension = 1;
        f.fdata     = NULL;
        f.idata     = NULL;
        f.name      = "opening";

        f.idata    = new int*[inp.nv];

        for (int n = 0; n < inp.nv; n++) {
            f.idata[n]    = new int[1];
            
            int matched = 0;
            for (int m = 0; m < opening.nv; m++) {
                if ((inp.vertices[n][0]==opening.vertices[m][0])&&(inp.vertices[n][1]==opening.vertices[m][1])&&(inp.vertices[n][2]==opening.vertices[m][2])){
                    matched = 1;
                    break;
                }
            }

            f.idata[n][0] = matched;
        }

        inp.fields.push_back(f);

    };

    createOpeningLabel(cer_wm_closed,wm_opening);
    createOpeningLabel(cer_mask_closed,mask_opening);

    return {cer_mask_closed,cer_mask_open,cer_wm_closed,cer_wm_open};

}

std::vector<Surface> NIBR::fsAseg2CerebellumSurf(std::string asegFile, float distThresh, float meanFaceArea)
{

    // Read aseg file
    if (!existsFile(asegFile)) {
        disp(MSG_ERROR, "%s was not found", asegFile.c_str());
        return std::vector<Surface>();
    }

    Image<int> asegImg(asegFile);
    asegImg.read();

    return fsAseg2CerebellumSurf(asegImg,distThresh,meanFaceArea);

}