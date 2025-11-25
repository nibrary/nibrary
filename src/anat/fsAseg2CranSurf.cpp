#include "anatSurf.h"
#include "surface/surface2imageMapper.h"
#include "surface/surface_operators.h"
#include "surface/surface_make.h"

using namespace NIBR;

#define FS_OPTIC_CHIASM 85

Surface NIBR::fsAseg2BgSurf(Image<int>& asegImg, Surface& brainStem, float meanFaceArea)
{

    Image<int> tmp;
    tmp.createFromTemplate(asegImg,true);

    Image<int8_t> brainStemMask;
    brainStemMask.createFromTemplate(asegImg,true);
    mapSurface2Image(&brainStem,&brainStemMask,0,NULL,NULL,MASK_WITH_BOUNDARY);

    
    Image<bool> opticChiasmMask;
    opticChiasmMask.createFromTemplate(asegImg,true);    
    for (int n = 0; n < opticChiasmMask.numel; n++) {
        opticChiasmMask.data[n] = (asegImg.data[n] == FS_OPTIC_CHIASM);
    }
    imgDilate(opticChiasmMask);
    imgDilate(opticChiasmMask);

    for (int n = 0; n < tmp.numel; n++) {
        tmp.data[n] = ((asegImg.data[n] > 0) && (opticChiasmMask.data[n] == 0)) || (brainStemMask.data[n] > 0);
    }

    auto tightBrainSurf = label2surface(tmp,1,2);

    tightBrainSurf = surfGrow(tightBrainSurf, 1,  5);
    tightBrainSurf = surfGrow(tightBrainSurf, 1, -5);

    meanFaceArea = std::max(meanFaceArea,0.25f);

    if (tightBrainSurf.nv > 0) tightBrainSurf.calcArea();
    if (tightBrainSurf.nv > 0) tightBrainSurf = surfRemesh(tightBrainSurf,tightBrainSurf.area/meanFaceArea*0.5f,1,0);

    tightBrainSurf = surfMakeItSingleClosed(tightBrainSurf);

    tightBrainSurf.flipNormalsOfFaces();

    Surface box = surfMakeBox(asegImg.getBoundingBox());

    return surfMerge(box,tightBrainSurf);

}

Surface NIBR::fsAseg2BgSurf(std::string asegFile, std::string fslFirstFolder, float meanFaceArea)
{
    // Read aseg file
    if (!existsFile(asegFile)) {
        disp(MSG_ERROR, "%s was not found", asegFile.c_str());
        return Surface();
    }

    Image<int> asegImg(asegFile);
    asegImg.read();

    Surface brainStem;

    if (fslFirstFolder != "") brainStem = fslFirst2BrainstemSurf(fslFirstFolder, meanFaceArea);

    return fsAseg2BgSurf(asegImg,brainStem,meanFaceArea);

}