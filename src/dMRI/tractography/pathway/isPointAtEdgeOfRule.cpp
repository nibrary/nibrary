#include "pathway.h"
#include "math/sphere.h"

bool NIBR::Pathway::isPointAtEdgeOfRule(float* p, int ruleNo, float distThresh) {

    switch (srcType[ruleNo]) {

        // Undefined src
        case undef_src:   {break;}

        // Point
        case res_pnt_src: {break;}

        // Sphere
        case sph_src: {

            float p2c     = squared_dist(sphCenter[ruleNo],p);
            float big_sph = (sphRadius[ruleNo]+distThresh)*(sphRadius[ruleNo]+distThresh);
            float sma_sph = (sphRadius[ruleNo]-distThresh)*(sphRadius[ruleNo]-distThresh);

            if ( (p2c < big_sph) && (p2c > sma_sph) )
                return true;

            break;
        }

        // Images
        case img_mask_src:
        case img_label_src:
        case img_pvf_src: {

            float mp[3];
            int insideCnt  = 0;
            int outsideCnt = 0;

            for (int n = 0; n < 256; n++) {

                mp[0] = p[0] + FULLSPHERE256[n][0] * distThresh;
                mp[1] = p[1] + FULLSPHERE256[n][1] * distThresh;
                mp[2] = p[2] + FULLSPHERE256[n][2] * distThresh;

                if (isPointInsideRule(mp,ruleNo)) 
                    insideCnt++;
                else
                    outsideCnt++;
            }

            if ((insideCnt > 0) && (outsideCnt > 0)) {
                return true;
            }

            break;
        }

        // Surfaces
        case surf_src: {

            if (surf[ruleNo]->squaredDistToPoint(p) < (distThresh*distThresh))
                return true;

            break;
        }

    }

    return false;

}