#include "anatSurf.h"

using namespace NIBR;

Surface NIBR::fslFirst2BrainstemSurf(std::string fslFirstFolder, float meanFaceArea)
{
    // Check if folder exists
    if (!existsFolder(fslFirstFolder)) {
        disp(MSG_ERROR, "FSL first folder does not exist %s", fslFirstFolder.c_str());
        return Surface();
    }

    // Check if template image exists
    std::string tempImg = getFolderPath(fslFirstFolder) + "/first_all_fast_firstseg.nii.gz";

    if (!existsFile(tempImg)) {
        disp(MSG_ERROR, "Required reference image, %s, was not found under %s", tempImg.c_str(),fslFirstFolder.c_str());
        return Surface();
    }


    std::string fname = getFolderPath(fslFirstFolder) + "/first-BrStem_first.vtk";

    if (!existsFile(fname)) {
        disp(MSG_ERROR, "first-BrStem_first.vtk was not found under %s",fslFirstFolder.c_str());
        return Surface();
    }

    Surface surf(fname);
    surf.readMesh();

    if (meanFaceArea != 0) {
        surf.calcArea();
        surf = surfRemesh(surf,surf.area/meanFaceArea*0.5f,1,0);
    }

    surf.convertFSLsurface2RASmm(tempImg);

    return surf;

}