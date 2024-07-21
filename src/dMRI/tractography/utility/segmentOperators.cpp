#include "segmentOperators.h"

double NIBR::segmentSphereIntersectionLength(NIBR::Segment& seg, float* sphere)
{

    double r  = sphere[3];
    double rr = r*r;

    auto beg = seg.p;
    double end[3];
    vec3add(end,beg,seg.dir,seg.length);
    
    double b2c[3];
    vec3sub(b2c, sphere, beg);
    double b2c_norm = norm(b2c);

    double e2c[3];
    vec3sub(e2c, sphere, end);
    double e2c_norm = norm(e2c);

    // both ends are inside the sphere
    if ((b2c_norm<=r) && (e2c_norm<=r)) {
        return seg.length;
    }

    // beg is inside, end is outside the sphere
    if ((b2c_norm<=r) && (e2c_norm>r)) {
        
        // segment is exiting the sphere
        // return the portion that remains inside the sphere
        double proj_h = dot(&b2c[0],&seg.dir[0]);
        double dd     = b2c_norm*b2c_norm - proj_h*proj_h;
        double dist   = proj_h + std::sqrt(rr - dd);
        return dist;

    }

    // beg is is outside the sphere
    if (b2c_norm>r) {
        
        double proj_h = dot(&b2c[0],&seg.dir[0]);

        // segment can't intersect sphere
        if ( proj_h < 0.0 )
            return 0;
        
        // segment can't intersect sphere
        double proj_vsq = (b2c_norm*b2c_norm) - (proj_h*proj_h);
        if (proj_vsq > rr)
            return 0;

        // segment can't intersect sphere
        double chord = std::sqrt( rr - proj_vsq );
        double dist = proj_h - chord;
        if (dist>seg.length)
            return 0;

        // segment intersects sphere
        return ((seg.length-dist) > 2.0f*chord) ? 2.0f*chord : (seg.length-dist);

    }

    return 0;

}