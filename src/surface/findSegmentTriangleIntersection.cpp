#include "findSegmentTriangleIntersection.h"

using namespace NIBR;

// Find intersection of a line segment with a triangle using Möller–Trumbore algorithm
// This functions returns the angle of intersection in degrees [0,90]
// 0 means there is no intersection.

double NIBR::findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, double* p, double* dir, double length, double* pos, double* dist) {

    // Set triangle
    auto ref = inpSurf->vertices[inpSurf->faces[faceIndex][0]];
    auto p1  = inpSurf->vertices[inpSurf->faces[faceIndex][1]];
    auto p2  = inpSurf->vertices[inpSurf->faces[faceIndex][2]];

    double v1[3] = {double(p1[0])-double(ref[0]),double(p1[1])-double(ref[1]),double(p1[2])-double(ref[2])};
    double v2[3] = {double(p2[0])-double(ref[0]),double(p2[1])-double(ref[1]),double(p2[2])-double(ref[2])};
    
    // Determinant
    double tmp[3];
    cross(&tmp[0],&dir[0],&v2[0]);
    double det = dot(&tmp[0],&v1[0]);
    
    // Segment is parallel to the triangle
    if ((det == 0.0) || (isnan(det))) return 0.0;
    
    double T[3] = {p[0]-double(ref[0]),p[1]-double(ref[1]),p[2]-double(ref[2])};
    
    // Check if the first barycentric coordinate is within limits
    double u = dot(&tmp[0],&T[0]);
    if (( (det < 0.0) && ((u > 0.0) || (u < det)) ) ||
        ( (det > 0.0) && ((u < 0.0) || (u > det)) ) ) return 0.0;

    // Check if the second barycentric coordinate is within limits
    cross(&tmp[0],&T[0],&v1[0]);
    double v = dot(&tmp[0],&dir[0]);    
    if (( (det < 0.0) && ((v > 0.0) || ((u+v) < det)) ) ||
        ( (det > 0.0) && ((v < 0.0) || ((u+v) > det)) ) ) return 0.0;
    
    // Check if t is within range, i.e., segment crosses the triangle
    double t  = dot(&tmp[0], &v2[0]);

    disp(MSG_DEBUG,"faceInd: %d, l: %.14f, t: %.14f, t-f:%.14f", faceIndex, length, t/det, t-length*det);

    if (( (det < 0.0) && ((t > 0.0) || (t < (length*det)))) ||
        ( (det > 0.0) && ((t < 0.0) || (t > (length*det))))) {
        *dist = NAN;
        if (pos!=NULL) {
            pos[0] = NAN;
            pos[1] = NAN;
            pos[2] = NAN;
        }
        return 0.0;
    }

    t /= det;
    t  = std::clamp(t,0.0,length);
        
    double angle  = M_PI_2 - std::acos(std::clamp(std::fabs(dot(&dir[0],&inpSurf->normalsOfFaces[faceIndex][0])),-1.0,1.0));
    angle *= (180.0/M_PI);
    
    if (pos!=NULL) {
        pos[0] = p[0] + t*dir[0];
        pos[1] = p[1] + t*dir[1];
        pos[2] = p[2] + t*dir[2];
    }

    *dist  = t;

    return angle;
        
}

double NIBR::findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, float* p, double* dir, double length, double* pos, double* dist) {
    double pd[3] = {p[0],p[1],p[2]};
    return NIBR::findSegmentTriangleIntersection(inpSurf, faceIndex, pd, dir, length, pos, dist);
}
