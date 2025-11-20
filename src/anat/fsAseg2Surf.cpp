#include "anatSurf.h"

using namespace NIBR;

Surface NIBR::fsAseg2Surf(Image<int>& asegImg, int asegLabel, float meanFaceArea)
{

    Surface surf = label2surface(asegImg,asegLabel,((meanFaceArea==0) ? 0.25 : meanFaceArea));

    std::vector<int> labelField(surf.nv,asegLabel);

    surf.fields.push_back(surf.makeVertField("label",labelField));

    return surf;

}

Surface NIBR::fsAseg2Surf(std::string asegFile, int asegLabel, float meanFaceArea)
{
    // Read aseg file
    if (!existsFile(asegFile)) {
        disp(MSG_ERROR, "%s was not found", asegFile.c_str());
        return Surface();
    }

    Image<int> asegImg(asegFile);
    asegImg.read();

    return fsAseg2Surf(asegImg,asegLabel,meanFaceArea);

}