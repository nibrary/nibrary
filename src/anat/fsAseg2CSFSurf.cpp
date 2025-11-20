#include "anatSurf.h"

#define LEFT_LATERAL_VENTRICLE      4
#define LEFT_INF_LAT_VENT           5
#define LEFT_CHOROID_PLEXUS        31
#define RIGHT_LATERAL_VENTRICLE    43
#define RIGHT_INF_LAT_VENT         44
#define RIGHT_CHOROID_PLEXUS       63
#define THIRD_VENTRICLE            14
#define FOURTH_VENTRICLE           15
#define FIFTH_VENTRICLE            72
#define CSF                        24


using namespace NIBR;

Surface NIBR::fsAseg2CSFSurf(Image<int>& asegImg, float meanFaceArea)
{

    std::vector<int> labels = {
        LEFT_LATERAL_VENTRICLE,
        LEFT_INF_LAT_VENT,
        LEFT_CHOROID_PLEXUS,
        RIGHT_LATERAL_VENTRICLE,
        RIGHT_INF_LAT_VENT,
        RIGHT_CHOROID_PLEXUS,
        THIRD_VENTRICLE,
        FOURTH_VENTRICLE,
        FIFTH_VENTRICLE,
        CSF};


    Surface csfSurf;
    std::vector<int>  csfLabels;
    std::vector<int> sideLabels;
    int i = 0;

    for (auto l : labels) {
        Surface surf = label2surface(asegImg,l,((meanFaceArea==0) ? 0.25 : meanFaceArea));
        csfSurf      = surfMerge(csfSurf,surf);

        for (int n = 0; n < surf.nv; n++) {
            csfLabels.push_back(l);
            if (i < 3) {
                sideLabels.push_back(1);
            } else if (i < 6) {
                sideLabels.push_back(2);
            } else {
                sideLabels.push_back(0);
            }
        }

        i++;
        disp(MSG_DETAIL, "Generated CSF surface, label: %d - nv: %d", l, surf.nv);
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

        csfSurf.fields.push_back(f);

    };


    makeAndAddField("side", sideLabels);
    makeAndAddField("label", csfLabels);

    return csfSurf;

}

Surface NIBR::fsAseg2CSFSurf(std::string asegFile, float meanFaceArea)
{

    // Read aseg file
    if (!existsFile(asegFile)) {
        disp(MSG_ERROR, "%s was not found", asegFile.c_str());
        return Surface();
    }

    Image<int> asegImg(asegFile);
    asegImg.read();

    return fsAseg2CSFSurf(asegImg,meanFaceArea);

}