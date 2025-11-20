#include "anatSurf.h"

using namespace NIBR;

Surface NIBR::fslFirst2SubcortexSurf(std::string fslFirstFolder, float meanFaceArea, bool excludeBrainStem)
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

    // Convert and merge subcortical surfaces
    std::vector<std::string> label = {"Accu","Amyg","Caud","Hipp","Pall","Puta","Thal"};

    Surface subCortex;
    std::vector<int> subCortexLabels;
    std::vector<int> sideLabels;
    int labelNum = 1;

    auto addToSubCortex = [&](std::string fprefix)->bool {

        std::string fname = getFolderPath(fslFirstFolder) + "/" + fprefix;

        if (!existsFile(fname)) {
            disp(MSG_ERROR, "%s was not found under %s", fprefix.c_str(),fslFirstFolder.c_str());
            return false;
        }

        Surface surf(fname);
        surf.readMesh();

        if (meanFaceArea != 0) {
            surf.calcArea();
            surf = surfRemesh(surf,surf.area/meanFaceArea*0.5f,1,0);
        }

        surf.convertFSLsurface2RASmm(tempImg);

        subCortex = surfMerge(subCortex,surf);

        for (int n = 0; n < surf.nv; n++) {
            subCortexLabels.push_back(labelNum);
            if (fprefix == "first-BrStem_first.vtk") {
                sideLabels.push_back(0);
            } else {
                sideLabels.push_back((labelNum-1)%2+1);
            }
        }

        labelNum++;

        return true;

    };


    for (const auto& l : label) {
        for (std::string side : {"L","R"}) {
            if (!addToSubCortex("first-" + side + "_" + l + "_first.vtk")) {
                return Surface();
            }
        }
    }


    if (!excludeBrainStem) {
        if (!addToSubCortex("first-BrStem_first.vtk")) {
            return Surface();
        }
    }

    auto makeAndAddField = [&](std::string fieldName,const std::vector<int> fData) {

        SurfaceField f;
        f.owner     = VERTEX;
        f.datatype  = "int";
        f.dimension = 1;
        f.fdata     = NULL;
        f.idata     = NULL;
        f.name      = fieldName;

        f.idata    = new int*[fData.size()];

        for (int n = 0; n < int(fData.size()); n++) {
            f.idata[n]    = new int[1];
            f.idata[n][0] = fData[n];
        }

        subCortex.fields.push_back(f);

    };

    makeAndAddField("side", sideLabels);
    makeAndAddField("label",subCortexLabels);

    return subCortex;

}