#include "pathway.h"
#include "pathwayRule.h"
#include <cstdint>
#include <type_traits>
#include <vector>
#include "base/fileOperations.h"

#define SUB_VOXEL_RATIO 0.25

using namespace NIBR;

bool NIBR::Pathway::addSurface(PathwayRule prule) {

    prule.surfSrc->enablePointCheck(prule.surfaceDiscretizationRes);
    surf.back() = prule.surfSrc;

    // Prepare if type is seed
    if ((prule.type == seed) && (isTracking == true)) {

        disp(MSG_DETAIL,"Prepping seed surface for tracking");

        Seeder* seedDef;
        if (prule.surface4SeedUseInside) {
            seedDef = new SeedInsideSurface();
            if (!seedDef->setSeed(surf.back(),prule.surfaceDiscretizationRes)) {
                delete seedDef;
                return false;
            }
        } else {
            disp(MSG_ERROR,"Can't add seed surface without seeding inside yet.");
            return false;
        }
        
        seeds.back() = seedDef;

        disp(MSG_DETAIL,"Seed surface is ready for tracking");
        return true;

    }

    return true;

}