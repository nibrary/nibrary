#include "pathway.h"
#include "pathwayRule.h"
#include <cstdint>
#include <type_traits>
#include <vector>
#include "base/fileOperations.h"

#define SUB_VOXEL_RATIO 0.25

using namespace NIBR;

bool NIBR::Pathway::addSurface(PathwayRule prule) {

    if (prule.surfaceFlipNormals)
        prule.surfSrc->flipNormalsOfFaces();

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

    if ((prule.surfSrc->interpretAs2D == false) && (prule.surfSrc->openOrClosed == OPENANDCLOSED) && (prule.type == seed)) {
        disp(MSG_DETAIL,"Seed surfaces containing both open and closed components are not supported");
    }

    if (prule.surfSrc->interpretAs2D) {
        switch (prule.type) {
            case undef_type:             {disp(MSG_ERROR,"Unacceptable rule type"); return false;}
            case seed:                   {break;}
            case discard_seed:           {break;}
            case req_entry:              {break;}
            case req_exit:               {break;}
            case req_end_inside:         {break;}
            case stop_before_entry:      {break;}
            case stop_at_entry:          {break;}
            case stop_after_entry:       {disp(MSG_ERROR,"stop_after_entry option is not allowed for 2D-interpreted surfaces"); return false;}
            case stop_before_exit:       {disp(MSG_ERROR,"stop_before_exit option is not allowed for 2D-interpreted surfaces"); return false;}
            case stop_at_exit:           {break;}
            case stop_after_exit:        {break;}
            case discard_if_enters:      {break;}
            case discard_if_exits:       {break;}
            case discard_if_ends_inside: {disp(MSG_ERROR,"Unacceptable rule type"); return false;}
        }
    }

    prule.surfSrc->enablePointCheck(prule.surfaceDiscretizationRes);
    prule.surfSrc->calcNormalsOfFaces();
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