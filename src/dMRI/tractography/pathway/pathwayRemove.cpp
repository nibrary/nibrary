#include "pathway.h"
#include "pathwayRule.h"
#include "base/vectorOperations.h"

using namespace NIBR;

bool NIBR::Pathway::remove(int ruleInd) {

    disp(MSG_DETAIL,"Deleting rule: %d of %d.", ruleInd+1, ruleCnt);

    if (ruleInd >= ruleCnt) {
        disp(MSG_ERROR,"Rule index is out of bounds");
        return false;
    }

    ruleCnt--;

    // We first reverse B_pulling. See verify() for details.
    if (B_pulled) {

        for(size_t i=0;i<prules.size(); ++i)
            if (prules[i].type != seed)
                prules[i].side  = prules[i].orig_side;

        B_pulled = false;
    }

    // Explicitly call the destructors, since these are created with "new"

    // Remove sources only if not shared
    if (( img_mask[ruleInd]!=NULL) && isUnique(img_mask,  ruleInd)) { delete  img_mask[ruleInd];  img_mask[ruleInd] = NULL;}
    if ((img_label[ruleInd]!=NULL) && isUnique(img_label, ruleInd)) { delete img_label[ruleInd]; img_label[ruleInd] = NULL;}
    if ((  img_pvf[ruleInd]!=NULL) && isUnique(img_pvf,   ruleInd)) { delete   img_pvf[ruleInd];   img_pvf[ruleInd] = NULL;}

    if (( surfData[ruleInd]!=NULL) && isUnique(surfData,  ruleInd)) { surf[ruleInd]->clearField(*surfData[ruleInd]); delete  surfData[ruleInd];  surfData[ruleInd] = NULL;}
    if ((     surf[ruleInd]!=NULL) && isUnique(surf,      ruleInd)) { delete      surf[ruleInd];      surf[ruleInd] = NULL;}
   
    if (    seeds[ruleInd] != NULL)  { delete       seeds[ruleInd];     seeds[ruleInd] = NULL; }
    if (sphCenter[ruleInd] != NULL)  { delete[] sphCenter[ruleInd]; sphCenter[ruleInd] = NULL; }
    if ( pntLists[ruleInd] != NULL)  { delete    pntLists[ruleInd];  pntLists[ruleInd] = NULL; }
    if ( dirLists[ruleInd] != NULL)  { delete    dirLists[ruleInd];  dirLists[ruleInd] = NULL; }
    if (     data[ruleInd] != NULL)  { delete        data[ruleInd];      data[ruleInd] = NULL; }

    // Erase calls destructor of objects unless they are constructed with "new"
    prules.          erase(prules.          begin() + ruleInd);
    srcType.         erase(srcType.         begin() + ruleInd);
    seeds.           erase(seeds.           begin() + ruleInd);
    img_mask.        erase(img_mask.        begin() + ruleInd);
    img_label.       erase(img_label.       begin() + ruleInd);
    img_label_val.   erase(img_label_val.   begin() + ruleInd);
    img_pvf.         erase(img_pvf.         begin() + ruleInd);
    pvf_vol.         erase(pvf_vol.         begin() + ruleInd);
    miniSegment.     erase(miniSegment.     begin() + ruleInd);
    surf.            erase(surf.            begin() + ruleInd);
    surfData.        erase(surfData.        begin() + ruleInd);
    sphCenter.       erase(sphCenter.       begin() + ruleInd);
    sphRadius.       erase(sphRadius.       begin() + ruleInd);
    sphRadiusSquared.erase(sphRadiusSquared.begin() + ruleInd);
    pntLists.        erase(pntLists.        begin() + ruleInd);
    dirLists.        erase(dirLists.        begin() + ruleInd);
    data.            erase(data.            begin() + ruleInd);

    isVerified = false;
    if (verify()) {
        disp(MSG_DETAIL, "Rule %d removed successfully", ruleInd+1);
        return true;
    } else {
        disp(MSG_ERROR, "Rule %d removal failed", ruleInd+1);
        return false;
    }

}

