#include "pathway.h"

bool NIBR::Pathway::isPointInsideRule(float* p, int ruleNo) {

    // disp(MSG_DEBUG,"isPointInsideRule: [%.2f,%.2f,%.2f], ruleNo: %d", p[0],p[1],p[2],ruleNo);

    switch (srcType[ruleNo]) {

        // Undefined src
        case undef_src: {
            // disp(MSG_FATAL,"Unexpected rule type \"undef_src\"");    // This should never happen, if pathway is updated
            return false;
        }

        // Used for seeding purposes only. Reserved for a random single point. 
        case res_pnt_src: {
            // disp(MSG_WARN,"Unexpected rule type \"res_pnt_src\""); // res_pnt_src is reserved for seed rule type only
            return false;
        }

        // Sphere
        case sph_src: {
            float p2c[3];
            vec3sub(p2c, sphCenter[ruleNo], p);
            float p2c_sq_norm = p2c[0]*p2c[0] + p2c[1]*p2c[1] + p2c[2]*p2c[2];
            return (p2c_sq_norm<=sphRadiusSquared[ruleNo]) ? true : false;
        }

        // Image - mask
        case img_mask_src: {
            return ((*img_mask[ruleNo])(p) == 1);
        }

        // Image - label
        case img_label_src: {
            return ((*img_label[ruleNo])(p)==img_label_val[ruleNo]) ? true : false;
        }

        // Image - partial volume
        case img_pvf_src: {

            // PVF is 4D
            if (img_pvf[ruleNo]->getDimension() == 4) {
                
                // float fraction;
                // float maxFraction = FLT_MIN;
                // int   maxVolInd   = -1;
                
                // for (int t=0; t<img_pvf[ruleNo]->imgDims[3]; t++) {
                //     fraction = (*img_pvf[ruleNo])(p,t);
                //     if ( fraction > maxFraction) {
                //         maxFraction = fraction;
                //         maxVolInd   = t;
                //     }
                // }

                // return (maxVolInd == pvf_vol[ruleNo]) ? true : false;

                return ((*img_pvf[ruleNo])(p,pvf_vol[ruleNo]) >= pvfThresh) ? true : false;

            }

            // PVF is 3D
            // disp(MSG_DEBUG,"val: %.2f",(*img_pvf[ruleNo])(p));
            return ((*img_pvf[ruleNo])(p) > 0.0f) ? true : false;            

        }

        // Surface
        case surf_src:
        case surf_ins_src: {
            return surf[ruleNo]->isPointInside(p);
        }

    }

    disp(MSG_FATAL,"Program reached unexpected state. Unknown pathway rule source: %d", srcType[ruleNo]);
    return false;

}