#include "pathway.h"
#include "pathwayRule.h"

using namespace NIBR;

// If a segment is entering a pathway rule, this function returns true
// crossDist shows the fraction of segment length to enter the pathway rule
std::tuple<bool,float> NIBR::Pathway::isSegmentEntering(const LineSegment& segment, int ruleNo) {

    // crossDist is always between [0,segment.len].
    float crossDist = segment.len;

    switch (srcType[ruleNo]) {

        // Undefined src
        case undef_src: {
            disp(MSG_FATAL,"Unexpected source type \"undef_src\"");    // This should never happen, if pathway is updated
            return std::make_tuple(false,NAN);
        }

        // Used for seeding purposes only. Reserved for a random single point.
        case res_pnt_src: {
            // disp(MSG_WARN,"Unexpected seed type \"res_pnt_src\""); // res_pnt_src is reserved for seeding purposes only
            return std::make_tuple(false,NAN);
        }

        // Sphere
        case sph_src: {

            float p2c[3];
            vec3sub(p2c, sphCenter[ruleNo], segment.beg);
            float p2c_norm = norm(p2c);

            // segment.beg is inside the sphere
            if (p2c_norm<=sphRadius[ruleNo]) {
                crossDist = 0.0f;
                return std::make_tuple(true,crossDist);
            }

            // segment can't intersect the sphere
            float proj_h = dot(p2c, segment.dir);
            if ( proj_h < 0.0 ) 
                return std::make_tuple(false,NAN);
            
            float proj_vsq = (p2c_norm*p2c_norm) - (proj_h*proj_h);
            if (proj_vsq > sphRadiusSquared[ruleNo])
                return std::make_tuple(false,NAN);

            float dist = proj_h - std::sqrt( sphRadiusSquared[ruleNo] - proj_vsq );
            if (dist>segment.len)
                return std::make_tuple(false,NAN);

            // segment intersects the sphere
            crossDist = dist;
            crossDist = std::clamp(crossDist,0.0f,segment.len);
            return std::make_tuple(true,crossDist);

        }

        // Images
        case img_mask_src:
        case img_label_src:
        case img_pvf_src: {

            // 1. Check beginning point
            if (isPointInsideRule(segment.beg, ruleNo)) {
                crossDist = 0.0f;
                return std::make_tuple(true, crossDist);
            }

            // 2. Step through the segment
            float t_prev = 0.0f;
            float current_step = miniSegment[ruleNo]; 
            
            // Use a loop that ensures we don't overshoot length
            while (current_step < segment.len) {
                
                float p_curr[3];
                p_curr[0] = segment.beg[0] + segment.dir[0] * current_step;
                p_curr[1] = segment.beg[1] + segment.dir[1] * current_step;
                p_curr[2] = segment.beg[2] + segment.dir[2] * current_step;

                if (isPointInsideRule(p_curr, ruleNo)) {
                    // We were outside at t_prev, but inside at current_step.
                    // The boundary is in between. Refine using bisection.
                
                    auto findEdge = [&](float low, float high) {
                        for(int k=0; k<12; k++) {
                            float mid = (low+high)*0.5f;
                            float pm[3] = {
                                segment.beg[0] + segment.dir[0] * mid,
                                segment.beg[1] + segment.dir[1] * mid,
                                segment.beg[2] + segment.dir[2] * mid
                            };
                            if(isPointInsideRule(pm, ruleNo)) high = mid; else low = mid;
                        }
                        return high;
                    };

                    crossDist = findEdge(t_prev, current_step);
                    return std::make_tuple(true, crossDist);
                }

                t_prev = current_step;
                current_step += miniSegment[ruleNo];
            }

            // 3. Check end point
            if (isPointInsideRule(segment.end, ruleNo)) {
                // Transition happened between t_prev and segment.len
                auto findEdge = [&](float low, float high) {
                    for(int k=0; k<12; k++) {
                        float mid = (low+high)*0.5f;
                        float pm[3] = {
                            segment.beg[0] + segment.dir[0] * mid,
                            segment.beg[1] + segment.dir[1] * mid,
                            segment.beg[2] + segment.dir[2] * mid
                        };
                        if(isPointInsideRule(pm, ruleNo)) high = mid; else low = mid;
                    }
                    return high;
                };

                crossDist = findEdge(t_prev, segment.len);
                return std::make_tuple(true, crossDist);
            }

            return std::make_tuple(false, NAN);
        }

        // Surf
        case surf_src: {

            // disp(MSG_DEBUG, "beg: [%.8f , %.8f , %.8f]", segment.beg[0],segment.beg[1],segment.beg[2]) ;
            // disp(MSG_DEBUG, "end: [%.8f , %.8f , %.8f]", segment.end[0],segment.end[1],segment.end[2]);

            auto [isBegInside,isEndInside,distance,intersectingFaceInd,towardsOutside,boundaryTransitionDist] = surf[ruleNo]->intersectSegment(&segment);

            // Segment beginning is inside
            if (isBegInside) {
                crossDist = 0.0f;
                disp(MSG_DEBUG, "Entered rule (out->in) %d at beg.", ruleNo);
                return std::make_tuple(true,crossDist);
            }
            // Segment beginning is outside
            
            // If there is intersection and segment enters inside
            if (!isnan(distance) && !towardsOutside) {
                crossDist = distance;
                crossDist = std::clamp(crossDist,0.0f,segment.len);
                disp(MSG_DEBUG, "Entered rule (out->in) %d at face %d with %.4f / %.4f at intersection.", ruleNo, intersectingFaceInd, float(distance), segment.len);
                return std::make_tuple(true,crossDist);
            }

            // Segment transitions inside the boundary without intersection
            if (!isnan(boundaryTransitionDist)) {
                crossDist = boundaryTransitionDist;
                crossDist = std::clamp(crossDist,0.0f,segment.len);
                disp(MSG_DEBUG, "Entered rule (out->in) %d boundary without intersection.", ruleNo);
                return std::make_tuple(true,crossDist);
            }

            // Segment end is inside (This case should not happen)
            if (isEndInside) {
                crossDist = segment.len;
                disp(MSG_DEBUG, "Entered rule (out->in) %d at end.", ruleNo);
                return std::make_tuple(true,crossDist);
            }

            // disp(MSG_DEBUG, "Did not enter rule %d.", ruleNo);
            return std::make_tuple(false,NAN);

        }

    }

    disp(MSG_FATAL,"Program reached unexpected state. Unknow pathway rule source: %d", srcType[ruleNo]);
    return std::make_tuple(false,NAN);
    
}

