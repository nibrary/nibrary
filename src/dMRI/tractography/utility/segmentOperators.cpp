#include "segmentOperators.h"

float NIBR::segmentSphereIntersectionLength(NIBR::Segment& seg, float* sphere)
{

    float r  = sphere[3];
    float rr = r*r;

    auto beg = seg.p;
    float end[3];
    vec3add(end,beg,seg.dir,seg.length);
    
    float b2c[3];
    vec3sub(b2c, sphere, beg);
    float b2c_norm = norm(b2c);

    float e2c[3];
    vec3sub(e2c, sphere, end);
    float e2c_norm = norm(e2c);

    // both ends are inside the sphere
    if ((b2c_norm<=r) && (e2c_norm<=r)) {
        return seg.length;
    }

    // beg is inside, end is outside the sphere
    if ((b2c_norm<=r) && (e2c_norm>r)) {
        
        // segment is exiting the sphere
        // return the portion that remains inside the sphere
        float proj_h = dot(&b2c[0],&seg.dir[0]);
        float dd     = b2c_norm*b2c_norm - proj_h*proj_h;
        float dist   = proj_h + std::sqrt(rr - dd);
        return dist;

    }

    // beg is is outside the sphere
    if (b2c_norm>r) {
        
        float proj_h = dot(&b2c[0],&seg.dir[0]);

        // segment can't intersect sphere
        if ( proj_h < 0.0 )
            return 0;
        
        // segment can't intersect sphere
        float proj_vsq = (b2c_norm*b2c_norm) - (proj_h*proj_h);
        if (proj_vsq > rr)
            return 0;

        // segment can't intersect sphere
        float chord = std::sqrt( rr - proj_vsq );
        float dist = proj_h - chord;
        if (dist>seg.length)
            return 0;

        // segment intersects sphere
        return ((seg.length-dist) > 2.0f*chord) ? 2.0f*chord : (seg.length-dist);

    }

    return 0;

}