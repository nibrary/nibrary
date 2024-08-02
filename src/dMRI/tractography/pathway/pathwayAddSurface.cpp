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
            prule.surfSrc->make2D();
        } else if (prule.surfSrc->openOrClosed == CLOSED) {
            NIBR::disp(MSG_DETAIL,"Surface is closed, interpreting as 3D: %s", prule.surfaceSource.c_str());
            prule.surfSrc->make3D();
        } else if (prule.surfSrc->openOrClosed == OPENANDCLOSED) {
            NIBR::disp(MSG_DETAIL,"Surface has open and closed components, interpreting as 3D: %s", prule.surfaceSource.c_str());
            prule.surfSrc->make3D();
        }
    }

    if (prule.surfaceUseDim == surf_useDim_3D) {
        if (prule.surfSrc->openOrClosed == OPEN) {
            NIBR::disp(MSG_DETAIL,"Surface is open, interpreting as 2D boundary: %s", prule.surfaceSource.c_str());
            prule.surfSrc->make2D();
        } else {
            prule.surfSrc->make3D();
        }
    }

    if (prule.surfaceUseDim == surf_useDim_2D) {
        prule.surfSrc->make2D();
    }

    prule.surfSrc->enablePointCheck(prule.surfaceDiscretizationRes);

    surf.back() = prule.surfSrc;

    // Prepare if type is seed
    if ((prule.type == seed) && (isTracking == true)) {

        disp(MSG_DETAIL,"Prepping seed surface for tracking");

        Seeder* seedDef;
        
        if (prule.surfSrc->openOrClosed == OPENANDCLOSED) {
            disp(MSG_DETAIL,"Seed surfaces containing both open and closed components are not supported");
            return false;
        } else if (prule.surfSrc->interpretAs2D == false) {
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