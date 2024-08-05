#include "pathway.h"

using namespace NIBR;

bool NIBR::Pathway::verify() {

    if (isVerified)
        return true;

    isVerified        = false;
    stopFlag          = false;
    isSided           = false;
    theOneSeed        = -1;
    pathwaySeed.clear();
    order_of_prules.clear();
    order_of_side_A_prules.clear();
    order_of_side_B_prules.clear();    

    // B_pulled means that an internal modification of rules was made during the previous verify() call. 
    // The internal modification pulls the "either" sided rules to "B" side in three cases:
    //   i.   If seeded, two_sided, IN_ORDER
    //   ii.  If seeded, one_sided
    //   iii. If there is no seed, two_sided, IN_ORDER
    // We will now revert the B pulling back.
    if (B_pulled) {
        for(size_t i=0;i<prules.size(); ++i)
            if (prules[i].type != seed)
                prules[i].side  = either;

        order_of_side_B_prules.clear();
        isSided  = false;
        B_pulled = false;
    }

    for(size_t i=0; i<prules.size(); ++i) {

        disp(MSG_DETAIL,"Checking rule %d", i);

        // Check whether sources can be used
        if (prules[i].src==undef_src) {
            disp(MSG_ERROR,"Pathway rule source %d is undefined", i);
            return false;
        }

        if ((prules[i].src==res_pnt_src) && (prules[i].type!=seed) ) {
            disp(MSG_ERROR,"Pathway rule source %d can't be set to reserved point unless pathway type seed", i);
            return false;
        }

        if (prules[i].src==sph_src) {
            if ( isnan(prules[i].radius   ) ||
                 isnan(prules[i].center[0]) ||
                 isnan(prules[i].center[1]) ||
                 isnan(prules[i].center[2]) ) {
                    disp(MSG_ERROR,"Pathway rule %d has NAN in sphere source", i);
                    return false;
                 }
            
            if (prules[i].radius<=0) {
                disp(MSG_ERROR,"Pathway rule %d has negative sphere radius", i);
                return false;
            }

        } 

        if ((prules[i].src==img_mask_src) && (prules[i].imageMaskSource=="")) {
            disp(MSG_ERROR,"Pathway rule %d is missing path to image", i);
            return false;
        }

        if ((prules[i].src==img_label_src) && (prules[i].imageLabelSource=="")) {
            disp(MSG_ERROR,"Pathway rule %d  is missing path to image", i);
            return false;
        }

        if ((prules[i].src==img_label_src) && (prules[i].useLabel==false)) {
            disp(MSG_DETAIL,"Using label: 1 for image %s", prules[i].src);
        }

        if ((prules[i].src==img_pvf_src) && (prules[i].imagePvfSource=="")) {
            disp(MSG_ERROR,"Pathway rule %d is missing path to image", i);
            return false;
        }

        if ((prules[i].src==surf_src) && (prules[i].surfaceSource=="")) {
            disp(MSG_ERROR,"Pathway rule %d is missing path to surface", i);
            return false;
        }

        // Check types
        if (prules[i].type==undef_type) {
            disp(MSG_ERROR,"Pathway rule %d type can't be undefined", i);
            return false;
        }
            
        if ((prules[i].type == stop_before_entry) || (prules[i].type == stop_at_entry) || (prules[i].type == stop_after_entry) || 
            (prules[i].type == stop_before_exit ) || (prules[i].type == stop_at_exit ) || (prules[i].type == stop_after_exit ))
            stopFlag = true;

        if (prules[i].side != either)
            isSided = true;

        switch (prules[i].side) {

            case side_A: 
            {
                if ( (prules[i].type == req_entry) || (prules[i].type == req_exit) || (prules[i].type == req_end_inside))
                    order_of_side_A_prules.push_back(i);
                break;
            }

            case side_B:
            {
                if ( (prules[i].type == req_entry) || (prules[i].type == req_exit) || (prules[i].type == req_end_inside))
                    order_of_side_B_prules.push_back(i);
                break;
            }

            default:
            {
                if ( (prules[i].type == req_entry) || (prules[i].type == req_exit) || (prules[i].type == req_end_inside) )
                    order_of_prules.push_back(i); 
                break;
            }

        }

        if (prules[i].type==seed) {
            pathwaySeed.push_back(i);
        }

    }

    theOneSeed = (pathwaySeed.size()==1) ? pathwaySeed[0] : -1;

    for(size_t i=0; i<prules.size(); ++i) {
        if((prules.size()>1) && (prules[i].type != seed) && (prules[i].type != discard_seed) ) {
            if ( ( (isSided && (prules[i].side == either) ) ) ||  ( (!isSided && (prules[i].side != either) ) ) ) {
                disp(MSG_ERROR,"Pathway rule side defined wrongly. All rules must have a side, or non of them can have a side. ");
                return false;
            }
        }
        if (!hasOneSeed() && (prules[i].type==discard_seed) ) {               // discard_seed requires a seed to be defined first
            disp(MSG_ERROR,"Seed must be provided first to use \"discard_seed\".");
            return false;
        }
    }

    // if (atMaxLength == ATMAXLENGTH_STOP)    stopFlag = true;
    if (skipSeedROI)                        stopFlag = true;
    if (directionality == ONE_SIDED)        stopFlag = true;

    if (seedTrials<0) {
        disp(MSG_ERROR,"seedTrials can't be negative");
        return false;
    }

    // Check noEdgeSeeds
    // if ( noEdgeSeeds && !hasOneSeed() ) {
    //     disp(MSG_ERROR,"Seed must be provided to prevent seeding at edges.");
    //     return false;
    // }

    // Check skipSeedROI
    if ( skipSeedROI && !hasOneSeed() ) {
        disp(MSG_ERROR,"Seed must be provided to skip it at the output.");
        return false;
    }

    if ( skipSeedROI && (directionality==TWO_SIDED) ) {
        disp(MSG_ERROR,"Tracking must be \"one_sided\" to skip writing of the seed region.");
        return false;
    }

    // Check stop rules
    if (stopFlag) disp(MSG_DETAIL," stopFlag"); else disp(MSG_DETAIL," no stopFlag");
    if (isSided)  disp(MSG_DETAIL," isSided");  else disp(MSG_DETAIL," no isSided");

    if (stopFlag && !hasOneSeed()) {
        disp(MSG_ERROR,"stop options can only be used with a seed.");
        return false;
    }

    if (stopFlag && (directionality==TWO_SIDED) && !isSided) {
        disp(MSG_ERROR,"stop options can be used with seeded \"two_sided\" tracking only when rules are defined for each side.");
        return false;
    }

    // Check length rules
    if (minLength<0) {
        disp(MSG_ERROR,"minlength can't be negative");
        return false;
    }
    if (maxLength<minLength) {
        disp(MSG_ERROR,"maxlength can't be smaller than minlength");
        return false;
    }


    // Check seed, directionality, sides, and apply B_pulling if needed.
    // B_pulling basically assigns sides to all prules if they are not given by the user.
    
    // With seed
    // if ( hasOneSeed() && (directionality==TWO_SIDED) && !isSided ) {} // Allowed     - B_pulling is needed when IN_ORDER is used
    // if ( hasOneSeed() && (directionality==TWO_SIDED) &&  isSided ) {} // Allowed     - B_pulling is NOT needed
    // if ( hasOneSeed() && (directionality==ONE_SIDED) && !isSided ) {} // Allowed     - B_pulling is needed
    // if ( hasOneSeed() && (directionality==ONE_SIDED) &&  isSided ) {} // Not allowed
    
    // Without seed
    // if (!hasOneSeed() && (directionality==TWO_SIDED) && !isSided ) {} // Allowed     - B_pulling is needed when IN_ORDER is used
    // if (!hasOneSeed() && (directionality==TWO_SIDED) &&  isSided ) {} // Not allowed
    // if (!hasOneSeed() && (directionality==ONE_SIDED) && !isSided ) {} // Not allowed
    // if (!hasOneSeed() && (directionality==ONE_SIDED) &&  isSided ) {} // Not allowed

    if ( hasOneSeed() && (directionality==ONE_SIDED) &&  isSided ) {     // Not allowed
        disp(MSG_ERROR,"Rules can't have sides when \"one_sided\" option is used.");
        return false;
    }

    if (!hasOneSeed() && (directionality==TWO_SIDED) &&  isSided ) {     // Not allowed
        disp(MSG_ERROR,"Rules can't have sides without seed.");
        return false;
    }
    if (!hasOneSeed() && (directionality==ONE_SIDED) ) {                 // Not allowed (for both isSided and !isSided)
        disp(MSG_ERROR,"Seed must be provided to use \"one_sided\" option.");
        return false;
    }

    // --- Prepare B_pulled rules ---
    if ( hasOneSeed() && (directionality==TWO_SIDED) && !isSided && (satisfy_requirements_in_order==IN_ORDER)) {
    
        int i=0;
        while (prules[i].type != seed) {
            prules[i].side = side_A;
            if ( (prules[i].type == req_entry) || (prules[i].type == req_exit) || (prules[i].type == req_end_inside) )
                    order_of_side_A_prules.push_back(i);
            i++;
        }
        std::reverse(order_of_side_A_prules.begin(),order_of_side_A_prules.end());
        i++;
        while (i<int(prules.size())) {
            prules[i].side = side_B;
            if ( (prules[i].type == req_entry) || (prules[i].type == req_exit) || (prules[i].type == req_end_inside) )
                    order_of_side_B_prules.push_back(i);
            i++;
        }

        order_of_prules.clear();
        isSided                = true;
        B_pulled               = true;
    
    }

    if (( hasOneSeed() && (directionality == ONE_SIDED) ) ||
        (!hasOneSeed() && (directionality == TWO_SIDED) && (satisfy_requirements_in_order == IN_ORDER))) {

        for(size_t i=0;i<prules.size(); ++i)
            if (prules[i].type != seed)
                prules[i].side = side_B;

        order_of_side_B_prules = order_of_prules;

        order_of_side_A_prules.clear();
        order_of_prules.clear();

        isSided                = true;
        B_pulled               = true;

    }

    // All good
    isVerified = true;
    return true;
}