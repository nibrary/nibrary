#include "anatSurf.h"

#define LEFT_ACCU    26
#define RIGHT_ACCU   58
#define LEFT_AMYG    18
#define RIGHT_AMYG   54
#define LEFT_CAUD    11
#define RIGHT_CAUD   50
#define LEFT_HIPP    17
#define RIGHT_HIPP   53
#define LEFT_PALL    13
#define RIGHT_PALL   52
#define LEFT_PUTA    12
#define RIGHT_PUTA   51
#define LEFT_THAL    10
#define RIGHT_THAL   49
#define BRAINSTEM    16
#define LEFT_VDC     28
#define RIGHT_VDC    60


using namespace NIBR;

Surface NIBR::fsAseg2SubcortexSurf(Image<int>& asegImg, float meanFaceArea, bool excludeBrainStem)
{

    // Convert and merge subcortical surfaces
    std::vector<int> label = {
        LEFT_ACCU,
        RIGHT_ACCU,
        LEFT_AMYG,
        RIGHT_AMYG,
        LEFT_CAUD,
        RIGHT_CAUD,
        LEFT_HIPP,
        RIGHT_HIPP,
        LEFT_PALL,
        RIGHT_PALL,
        LEFT_PUTA,
        RIGHT_PUTA,
        LEFT_THAL,
        RIGHT_THAL,
        LEFT_VDC,
        RIGHT_VDC};

    Surface subCortex;
    std::vector<int> subCortexLabels;
    std::vector<int> sideLabels;
    int labelNum = 1;

    auto addToSubCortex = [&](int l)->bool {

        Surface surf = label2surface(asegImg,l,((meanFaceArea==0) ? 0.25 : meanFaceArea));

        subCortex = surfMerge(subCortex,surf);

        for (int n = 0; n < surf.nv; n++) {
            subCortexLabels.push_back(l);
            if (l == BRAINSTEM) {
                sideLabels.push_back(0);
            } else {
                sideLabels.push_back((labelNum-1)%2+1);
            }
        }

        labelNum++;

        return true;

    };


    for (const auto& l : label) {
        if (!addToSubCortex(l)) {
            return Surface();
        }
    }


    if (!excludeBrainStem) {
        if (!addToSubCortex(BRAINSTEM)) {
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


Surface NIBR::fsAseg2SubcortexSurf(std::string asegFile, float meanFaceArea, bool excludeBrainStem)
{
    // Read aseg file
    if (!existsFile(asegFile)) {
        disp(MSG_ERROR, "%s was not found", asegFile.c_str());
        return Surface();
    }

    Image<int> asegImg(asegFile);
    asegImg.read();

    return fsAseg2SubcortexSurf(asegImg,meanFaceArea,excludeBrainStem);

}