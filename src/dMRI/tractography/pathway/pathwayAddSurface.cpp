#include "pathway.h"
#include "pathwayRule.h"
#include <cstdint>
#include <type_traits>
#include <vector>
#include "base/fileOperations.h"

#define SUB_VOXEL_RATIO 0.25

using namespace NIBR;

bool NIBR::Pathway::addSurface(PathwayRule prule) {

    prule.surfSrc->isClosed();

    if (prule.surfaceUseDim == surf_useDim_unknown) {
        if (prule.surfSrc->openOrClosed == OPEN) {
            NIBR::disp(MSG_DETAIL,"Surface is open, interpreting as 2D boundary: %s", prule.surfaceSource.c_str());
            surfIs2D.back() = true;
        } else if (prule.surfSrc->openOrClosed == CLOSED) {
            NIBR::disp(MSG_DETAIL,"Surface is closed, interpreting as 3D: %s", prule.surfaceSource.c_str());
            surfIs2D.back() = false;
        } else if (prule.surfSrc->openOrClosed == OPENANDCLOSED) {
            NIBR::disp(MSG_DETAIL,"Surface has open and closed components, interpreting as 3D: %s", prule.surfaceSource.c_str());
            surfIs2D.back() = false;
        }
    }

    if (prule.surfaceUseDim == surf_useDim_3D) {
        if (prule.surfSrc->openOrClosed == OPEN) {
            NIBR::disp(MSG_DETAIL,"Surface is open, interpreting as 2D boundary: %s", prule.surfaceSource.c_str());
            surfIs2D.back() = true;
        } else {
            surfIs2D.back() = false;
        }
    }

    if (prule.surfaceUseDim == surf_useDim_2D) {
        surfIs2D.back() = true;
    }

    prule.surfSrc->enablePointCheck(prule.surfaceDiscretizationRes,surfIs2D.back());

    surf.back() = prule.surfSrc;

    // Prepare if type is seed
    if ((prule.type == seed) && (isTracking == true)) {

        disp(MSG_DETAIL,"Prepping seed surface for tracking");

        Seeder* seedDef;
        if (surfIs2D.back() == false) {
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