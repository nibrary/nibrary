#include "pathway.h"
#include "pathwayRule.h"

using namespace NIBR;

// This is used for mask and label images, not for pvf images
template <typename T>
std::tuple<bool,float> checkEntryUsingRayTracing(float& crossDist, const LineSegment& segment, Image<T>* img, T label) {

    // segment.beg is inside the image
    if ((*img)(segment.beg)==label) {
        crossDist = 0.0f;
        disp(MSG_DEBUG,"   entered at crossDist %.12f mm", crossDist);
        return std::make_tuple(true,crossDist);
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
        return std::make_tuple(false,NAN);
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
                return std::make_tuple(false,NAN); // reached end of segment, there was no entry
            }

            // Cut t at length
            t = std::min(t,length);

            for (int m=0;m<3;m++) {
                p0[m] += t*dir[m];
                A[m]   = std::round(p0[m]);
            }


            if (img->isInside(A)) {
                img->to_xyz(p0,p1); // p0 is in image space and p1 is now the same point in real-space
                if ((*img)(p1)==label) {
                    vec3sub(dir,p1,segment.beg); // dir is now in real-space
                    crossDist = norm(dir);
                    crossDist = std::clamp(crossDist,0.0f,segment.len);
                    disp(MSG_DEBUG,"   entered at crossDist %.12f mm", crossDist);
                    return std::make_tuple(true,crossDist);
                }
            }

        } else {
            break;
        }

        length -= t;

    }

    // segment.end is inside the image
    if ((*img)(segment.end)==label) {
        crossDist = segment.len;
        disp(MSG_DEBUG,"   entered at crossDist %.12f mm", crossDist);
        return std::make_tuple(true,crossDist);
    }

    return std::make_tuple(false,NAN);

}

// Explicit instantiation for mask (int8_t) and label (int) images
template std::tuple<bool,float> checkEntryUsingRayTracing<int8_t>(float& crossDist, const LineSegment& segment, Image<int8_t>* img, int8_t label); // For img_mask
template std::tuple<bool,float> checkEntryUsingRayTracing<int>   (float& crossDist, const LineSegment& segment, Image<int>* img,    int    label); // For img_label

// If a segment is entering a pathway rule, this function returns true
// crossDist shows the fraction of segment length to enter the pathway rule
std::tuple<bool,float> NIBR::Pathway::isSegmentEntering(const LineSegment& segment, int ruleNo) {

    // crossDist is always between [0,segment.len].
    float crossDist = segment.len;

    switch (srcType[ruleNo]) {

        // Undefined src
        case undef_src: {
            disp(MSG_FATAL,"Unexpected seed type \"undef_src\"");    // This should never happen, if pathway is updated
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

        // Image - mask
        case img_mask_src: {

            return checkEntryUsingRayTracing<int8_t>(crossDist, segment, img_mask[ruleNo], 1);

        }
        
        // Image - label
        case img_label_src: {

            return checkEntryUsingRayTracing<int>(crossDist, segment, img_label[ruleNo], img_label_val[ruleNo]);

        }

        // Image - partial volume
        case img_pvf_src: {

            auto isInsidePvf=[&](float* p)->bool {
                if (img_pvf[ruleNo]->getDimension() == 4) { // PVF is 4D
                    return ((*img_pvf[ruleNo])(p,pvf_vol[ruleNo]) >= pvfThresh) ? true : false;
                } else { // PVF is 3D
                    return ((*img_pvf[ruleNo])(p) > 0.0f) ? true : false;
                }
            };

            

            if (isInsidePvf(segment.beg)) {
                crossDist = 0.0f;
                return std::make_tuple(true,crossDist);
            }
                        
            float downsampleFactor = segment.len * maxSegSizeScaler[ruleNo];

            // disp(MSG_DEBUG,"segment.len: %.2f, ,maxSegSizeScaler: %.2f, f: %.2f",segment.len, maxSegSizeScaler[ruleNo], downsampleFactor);

            if (downsampleFactor > 1)
            {
                float s = segment.len / float(std::ceil(downsampleFactor));
                float end[3];

                for (float t = s; t < (segment.len + EPS8); t = t + s)
                {
                    end[0] = segment.beg[0] + segment.dir[0] * t;
                    end[1] = segment.beg[1] + segment.dir[1] * t;
                    end[2] = segment.beg[2] + segment.dir[2] * t;
                    
                    if (isInsidePvf(end)) {
                        crossDist = t;
                        crossDist = std::clamp(crossDist,0.0f,segment.len);
                        return std::make_tuple(true,crossDist);
                    }
                }
            }
            else
            {
                if (isInsidePvf(segment.end)) {
                    crossDist = segment.len;
                    return std::make_tuple(true,crossDist);
                }
            }

            return std::make_tuple(false,NAN);

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

