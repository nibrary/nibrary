#include "findSegmentTriangleIntersection.h"

using namespace NIBR;

// Find intersection of a line segment with a triangle using Möller–Trumbore algorithm
// This functions returns the angle of intersection in degrees [0,90]
// 0 means there is no intersection.
float NIBR::findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, NIBR::Segment& seg) {
    
    float *ref,*v1,*v2;
    
    // Set triangle
    ref    = inpSurf->vertices[inpSurf->faces[faceIndex][0]];
    v1     = inpSurf->triangleEdge1[faceIndex];
    v2     = inpSurf->triangleEdge2[faceIndex];
    
    float tmp[3],T[3],invDet,det,u,v,t,angle;
    
    // Determinant
    cross(tmp,seg.dir,v2);
    det    = dot(&tmp[0],&v1[0]);
    
    // Segment is parallel to the triangle
    if ( (det<EPS6) && (det>-EPS6) ) 
        return 0;
    
    invDet = 1/det;
    T[0]   = seg.p[0]-ref[0];
    T[1]   = seg.p[1]-ref[1];
    T[2]   = seg.p[2]-ref[2];
    
    // Check if the first barycentric coordinate is within limits
    u = dot(&tmp[0],&T[0])*invDet;
    if ((u<0.0f) || (u>1.0f)) 
        return 0;
    
    // Check if the second barycentric coordinate is within limits
    cross(tmp,T,v1);
    v = dot(&tmp[0],&seg.dir[0])*invDet;
    if ((v<0.0f) || ((u+v)>1.0f)) 
        return 0;
    
    // Check if t is within range, i.e., segment crosses the triangle
    t = dot(&tmp[0],&v2[0])*invDet;
    
    if ((t<0) || (t>seg.length))
        return 0;
    
    angle  = M_PI_2 - std::acos( std::clamp(std::fabs(dot(&seg.dir[0],&inpSurf->normalsOfFaces[faceIndex][0])),-1.0f,1.0f));
    angle *= (180.0/M_PI);
    
    return angle;
    
}

float NIBR::findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, NIBR::Segment& seg, float extent) {
    
    float *ref,*v1,*v2;
    
    // Set triangle
    ref    = inpSurf->vertices[inpSurf->faces[faceIndex][0]];
    v1     = inpSurf->triangleEdge1[faceIndex];
    v2     = inpSurf->triangleEdge2[faceIndex];

    float tmp[3],T[3],invDet,det,u,v,t,angle;

    // Determinant
    cross(tmp,seg.dir,v2);
    det    = dot(&tmp[0],&v1[0]);

    // Segment is parallel to the triangle
    if ( (det<EPS6) && (det>-EPS6) ) 
        return 0;
    
    invDet = 1/det;
    T[0]   = seg.p[0]-ref[0]-seg.dir[0]*extent;
    T[1]   = seg.p[1]-ref[1]-seg.dir[1]*extent;
    T[2]   = seg.p[2]-ref[2]-seg.dir[2]*extent;
    
    // Check if the first barycentric coordinate is within limits
    u = dot(&tmp[0],&T[0])*invDet;
    if ((u<0.0f) || (u>1.0f)) 
        return 0;

    // Check if the second barycentric coordinate is within limits
    cross(tmp,T,v1);
    v = dot(&tmp[0],&seg.dir[0])*invDet;
    if ((v<0.0f) || ((u+v)>1.0f))  
        return 0;
    
    // Check if t is within range, i.e., segment crosses the triangle
    t = dot(&tmp[0],&v2[0])*invDet;

    if ((t<0) || (t>(seg.length+2.0f*extent)))
        return 0;
    
    angle  = M_PI_2 - std::acos(std::clamp(std::fabs(dot(&seg.dir[0],&inpSurf->normalsOfFaces[faceIndex][0])),-1.0f,1.0f));
    angle *= (180.0/M_PI);
    
    return angle;
    
}



float NIBR::findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, float* p, float* dir, float length, float* pos, float* dist) {
    
    float *ref,*v1,*v2;

    // Set triangle
    ref    = inpSurf->vertices[inpSurf->faces[faceIndex][0]];
    v1     = inpSurf->triangleEdge1[faceIndex];
    v2     = inpSurf->triangleEdge2[faceIndex];
    
    float tmp[3],T[3],invDet,det,u,v,t,angle;
    
    // Determinant
    cross(tmp,dir,v2);
    det    = dot(&tmp[0],&v1[0]);
    
    // Segment is parallel to the triangle
    if ( (det<EPS6) && (det>-EPS6) ) 
        return 0;
    
    invDet = 1/det;
    T[0]   = p[0]-ref[0];
    T[1]   = p[1]-ref[1];
    T[2]   = p[2]-ref[2];
    
    // Check if the first barycentric coordinate is within limits
    u = dot(&tmp[0],&T[0])*invDet;
    if ((u<0.0f) || (u>1.0f)) 
        return 0;
    
    // Check if the second barycentric coordinate is within limits
    cross(tmp,T,v1);
    v = dot(&tmp[0],&dir[0])*invDet;
    if ((v<0.0f) || ((u+v)>1.0f))  
        return 0;
    
    // Check if t is within range, i.e., segment crosses the triangle
    t = dot(&tmp[0],&v2[0])*invDet;
    
    if ((t<0) || (t>length)) {
        *dist = NAN;
        if (pos!=NULL) {
            pos[0] = NAN;
            pos[1] = NAN;
            pos[2] = NAN;
        }
        return 0;
    }
        
    angle  = M_PI_2 - std::acos(std::clamp(std::fabs(dot(&dir[0],&inpSurf->normalsOfFaces[faceIndex][0])),-1.0f,1.0f));
    angle *= (180.0/M_PI);
    
    if (pos!=NULL) {
        pos[0] = p[0] + t*dir[0];
        pos[1] = p[1] + t*dir[1];
        pos[2] = p[2] + t*dir[2];
    }

    *dist  = t;

    return angle;
        
}

void NIBR::findSegmentTriangleIntersection(NIBR::Surface* inpSurf, int faceIndex, float* p, const float* dir, float* dist) {
    
    float *ref,*v1,*v2;
    
    // Set triangle
    ref    = inpSurf->vertices[inpSurf->faces[faceIndex][0]];
    v1     = inpSurf->triangleEdge1[faceIndex];
    v2     = inpSurf->triangleEdge2[faceIndex];
    
    float tmp[3],T[3],invDet,det,u,v;
    
    // Determinant
    cross(tmp,dir,v2);
    det    = dot(&tmp[0],&v1[0]);
    
    // Segment is parallel to the triangle
    if ( (det<EPS6) && (det>-EPS6) ) {
        *dist = NAN;
        return;
    }
    
    invDet = 1/det;
    T[0]   = p[0]-ref[0];
    T[1]   = p[1]-ref[1];
    T[2]   = p[2]-ref[2];
    
    // Check if the first barycentric coordinate is within limits
    u = dot(&tmp[0],&T[0])*invDet;
    if ((u<0.0f) || (u>1.0f)) {
        *dist = NAN;
        return;
    }
    
    // Check if the second barycentric coordinate is within limits
    cross(tmp,T,v1);
    v = dot(&tmp[0],&dir[0])*invDet;
    if ((v<0.0f) || ((u+v)>1.0f)) {
        *dist = NAN;
        return;
    }
    
    *dist = dot(&tmp[0],&v2[0])*invDet;

    return;
        
}