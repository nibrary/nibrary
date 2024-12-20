#include "pathway.h"
#include "pathwayRule.h"

using namespace NIBR;

// This is used for mask and label images, not for pvf images
template <typename T>
bool checkExitUsingRayTracing(float& segCrosLength, LineSegment& segment, Image<T>* img, T label) {

    // segment.beg is outside the image
    if ((*img)(segment.beg)!=label) {
        segCrosLength = 0.0f;
        return true;
    }

    double p0[3]; // segment beg in image space
    img->to_ijk(segment.beg,p0);

    int32_t  A[3]; // voxel where segment beg is
    A[0] = std::round(p0[0]);
    A[1] = std::round(p0[1]);
    A[2] = std::round(p0[2]);

    double p1[3]; // segment end in image space
    img->to_ijk(segment.end,p1);

    // If beg and end are in the same voxel, there is nothing to do
    if ( (A[0]==std::round(p1[0])) && (A[1]==std::round(p1[1])) && (A[2]==std::round(p1[2])) ) {
        segCrosLength = 1.0;
        return false;
    }

    double dir[3], length;

    vec3sub(dir,p1,p0);   // ray direction in image space (not normalized)
    length  = norm(dir);  // length of line segment
    dir[0] /= length;     // ray direction
    dir[1] /= length;
    dir[2] /= length;
   
    double t = 0.0;

    while (length > 0.0) {

        // Check if and how long does it take for the segment to leave the current voxel
        if (rayTraceVoxel(A,p0,dir,t)) {

            // This part is written to be consistent with tractogram2imageMapper
            // Ray-tracing will cut t exactly at the current voxels edge, whose 3 faces are inside and 3 faces are outside
            // tractogram2imageMapper does the mapping based on what ray-tracing returns
            // i.e. if true, there is intersection, and the part until t is inside the current voxel
            // Therefore, if the intersection happened, then t*dir[m] will still be inside the current voxel.
            // So we add EPS4 to t to push it inside the next voxel if the segment is long enough

            if (length >= EPS4) {
                t += EPS4;
            } else {
                segCrosLength = 1.0;
                return false; // reached end of segment, there was no entry
            }

            // Cut t at length
            t = std::min(t,length);

            for (int m=0;m<3;m++) {
                p0[m] += t*dir[m];
                A[m]   = std::round(p0[m]);
            }


            if (img->isInside(A)) {
                img->to_xyz(p0,p1); // p0 is in image space and p1 is now the same point in real-space
                if ((*img)(p1)!=label) {
                    vec3sub(dir,p1,segment.beg); // dir is now in real-space
                    segCrosLength = norm(dir)/segment.len;
                    if (segCrosLength>1.0f) segCrosLength = 1.0f;
                    return true;
                }
            }

        } else {
            break;
        }

        length -= t;

    }

    // segment.end is inside the image
    if ((*img)(segment.end)!=label) {
        segCrosLength = 1.0f;
        return true;
    }

    segCrosLength = 1.0;
    return false;

}

// Explicit instantiation for mask (int8_t) and label (int) images
template bool checkExitUsingRayTracing<int8_t>(float& segCrosLength, LineSegment& segment, Image<int8_t>* img, int8_t label); // For img_mask
template bool checkExitUsingRayTracing<int>   (float& segCrosLength, LineSegment& segment, Image<int>* img,    int    label); // For img_label

// If a segment is exiting a pathway rule, this function returns true
// w->segCrosLength shows the fraction of segment length to exit the pathway rule
bool NIBR::Pathway::isSegmentExiting(NIBR::Walker* w, int ruleNo) {

    // segCrosLength is always between 0 and 1.
    // If segment does not exit the rule then segCrosLength is 1.
    w->segCrosLength = 1.0f;

    switch (srcType[ruleNo]) {

        // Undefined src
        case undef_src: {
            disp(MSG_FATAL,"Unexpected seed type \"undef_src\"");    // This should never happen, if pathway is updated
            return false;
        }

        // Used for seeding purposes only. Reserved for a random single point.
        case res_pnt_src: {
            // disp(MSG_WARN,"Unexpected seed type \"res_pnt_src\""); // res_pnt_src is reserved for seeding purposes only
            return false;
        }

        // Sphere
        case sph_src: {

            float p2c[3];
            vec3sub(p2c, sphCenter[ruleNo], w->segment.beg);
            float p2c_norm = norm(p2c);

            // segment.beg is outside the sphere
            if (p2c_norm>sphRadius[ruleNo]) {
                w->segCrosLength = 0.0f;
                return true;
            }

            // segment can't intersect the sphere
            float q2c[3];
            vec3sub(q2c, sphCenter[ruleNo], w->segment.end);
            float q2c_norm = norm(q2c);

            // segment.end is inside the sphere
            if (q2c_norm<=sphRadius[ruleNo])
                return false;

            // segment is exiting the sphere
            float proj_h = dot(p2c, w->segment.dir);
            float dd     = p2c_norm*p2c_norm - proj_h*proj_h;
            float dist   = proj_h + std::sqrt(sphRadiusSquared[ruleNo] - dd);

            // segment intersects the sphere
            w->segCrosLength = dist / w->segment.len;
            return true;

        }

        // Image - mask
        case img_mask_src: {

            return checkExitUsingRayTracing<int8_t>(w->segCrosLength, w->segment, img_mask[ruleNo], 1);

        }

        // Image - label
        case img_label_src: {

            return checkExitUsingRayTracing<int>(w->segCrosLength, w->segment, img_label[ruleNo], img_label_val[ruleNo]);

        }

        // Image - partial volume
        case img_pvf_src: {

            auto isOutsidePvf=[&](float* p)->bool {
                if (img_pvf[ruleNo]->getDimension() == 4) { // PVF is 4D
                    return ((*img_pvf[ruleNo])(p,pvf_vol[ruleNo]) < pvfThresh) ? true : false;
                } else { // PVF is 3D
                    return ((*img_pvf[ruleNo])(p) <= 0.0f) ? true : false;
                }
            };


            if (isOutsidePvf(w->segment.beg)) {
                w->segCrosLength = 0.0f;
                return true;
            }
                        
            float downsampleFactor = w->segment.len * maxSegSizeScaler[ruleNo];

            // disp(MSG_DEBUG,"segment.len: %.2f, ,maxSegSizeScaler: %.2f, f: %.2f",w->segment.len, maxSegSizeScaler[ruleNo], downsampleFactor);

            if (downsampleFactor > 1)
            {
                float s = w->segment.len / float(std::ceil(downsampleFactor));
                float end[3];

                for (float t = s; t < (w->segment.len + EPS8); t = t + s)
                {
                    end[0] = w->segment.beg[0] + w->segment.dir[0] * t;
                    end[1] = w->segment.beg[1] + w->segment.dir[1] * t;
                    end[2] = w->segment.beg[2] + w->segment.dir[2] * t;
                    
                    if (isOutsidePvf(end)) {
                        w->segCrosLength = t / w->segment.len;
                        return true;
                    }
                }
            }
            else
            {
                if (isOutsidePvf(w->segment.end)) {
                    w->segCrosLength = 1.0f;
                    return true;
                }
            }

            return false;

        }

        // Surf
        case surf_src: {

            // disp(MSG_DEBUG, "beg: [%.8f , %.8f , %.8f]", w->segment.beg[0],w->segment.beg[1],w->segment.beg[2]) ;
            // disp(MSG_DEBUG, "end: [%.8f , %.8f , %.8f]", w->segment.end[0],w->segment.end[1],w->segment.end[2]);

            auto [isBegInside,isEndInside,distance,intersectingFaceInd,towardsOutside,boundaryTransitionDist] = surf[ruleNo]->intersectSegment(&w->segment);

            // Segment beginning is outside
            if (!isBegInside) {
                w->segCrosLength = 0.0f;
                disp(MSG_DEBUG, "Exited rule (in->out) %d at beg.", ruleNo);
                return true;
            }
            // Segment beginning is inside
            
            // If there is intersection and segment enters outside
            if (!isnan(distance) && towardsOutside) {
                w->segCrosLength = distance / w->segment.len;
                if (w->segCrosLength > 1.0f) w->segCrosLength = 1.0f;
                disp(MSG_DEBUG, "Exited rule (in->out) %d at face %d with %.4f / %.4f at intersection.", ruleNo, intersectingFaceInd, float(distance), w->segment.len);
                return true;
            }

            // Segment transitions outside the boundary without intersection
            if (!isnan(boundaryTransitionDist)) {
                w->segCrosLength = boundaryTransitionDist / w->segment.len;
                disp(MSG_DEBUG, "Exited rule (in->out) %d boundary without intersection.", ruleNo);
                return true;
            }

            // Segment end is outside (This case should not happen)
            if (!isEndInside) {
                w->segCrosLength = 1.0f;
                disp(MSG_DEBUG, "Exited rule (in->out) %d at end.", ruleNo);
                return true;
            }

            // disp(MSG_DEBUG, "Did not exit rule %d.", ruleNo);
            return false;
        }

    }

    disp(MSG_FATAL,"Program reached unexpected state. Unknow pathway rule source: %d", srcType[ruleNo]);
    return false;
    
}

