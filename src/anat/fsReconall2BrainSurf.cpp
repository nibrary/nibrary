#include "anatSurf.h"

using namespace NIBR;

Surface NIBR::fsReconall2BrainSurf(std::string brainMaskFile, float meanFaceArea, float shellThickness)
{

    if (!existsFile(brainMaskFile)) {
        disp(MSG_ERROR, "%s was not found", brainMaskFile.c_str());
        return Surface();
    }

    Image<int> seg(brainMaskFile);
    seg.read();
    imgThresh(seg,0.5,FLT_MAX);

    auto innerShell = label2surface(seg,1,((meanFaceArea==0) ? 0.25 : meanFaceArea));
    
    if ( (shellThickness == 0) || (innerShell.nv == 0) ) {
        return innerShell;
    }

    auto outerShell = surfMoveVerticesAlongNormal(innerShell,shellThickness);
    innerShell.flipNormalsOfFaces();

    return surfMerge(outerShell,innerShell);
    
}
