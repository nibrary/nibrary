#include "anatSurf.h"

using namespace NIBR;

std::vector<Surface> NIBR::fsReconall2CorticalSurf(std::string fsPath)
{

    std::vector<Surface> out;

    if (!existsFolder(fsPath)) {
        disp(MSG_ERROR, "%s does not exist", fsPath.c_str());
        return out;
    }


    auto prepSurf = [&](std::string type)->Surface {

        Surface surfType;
        std::vector<int> sideLabel;
        std::vector<int>   fsLabel;

        int sideNo = 1;
        for (std::string side : {"lh","rh"}) {

            std::string fName = fsPath + "/surf/" + side + "." + type;
            if (!existsFile(fName)) {
                disp(MSG_ERROR, "%s does not exist", fName.c_str());
                return Surface();
            }

            std::string tName = fsPath + "/mri/orig.mgz";
            if (!existsFile(tName)) {
                disp(MSG_ERROR, "%s does not exist", tName.c_str());
                return Surface();
            }

            std::string aName = fsPath + "/label/"  + side + ".aparc.annot";
            if (!existsFile(aName)) {
                disp(MSG_ERROR, "%s does not exist", aName.c_str());
                return Surface();
            }

            Surface surf(fName);
            surf.readMesh();
            surf.convertFreesurferSurface2RASmm(tName);

            for (int n = 0; n < surf.nv; n++) {
                sideLabel.push_back(sideNo);
            }
            sideNo++;

            auto annot = surf.readFreesurferAnnot(aName, "");
            for (int n = 0; n < surf.nv; n++) {
                fsLabel.push_back(std::get<0>(annot).idata[n][0]);
            }

            surfType = surfMerge(surfType,surf);

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

            surfType.fields.push_back(f);

        };


        makeAndAddField("side", sideLabel);
        makeAndAddField("label",  fsLabel);

        return surfType;

    };

    return {prepSurf("white"),prepSurf("pial")};

}
