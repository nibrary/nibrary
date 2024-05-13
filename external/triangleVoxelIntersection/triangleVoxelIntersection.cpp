// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox2.txt

//******************************************************//
// AABB-triangle overlap test code                      //
// by Tomas Akenine-Möller                              //
// Function: int triBoxOverlap(float boxcenter[3],      //
//          float boxhalfsize[3],float triverts[3][3]); //
// History:                                             //
//   2001-03-05: released the code in its first version //
//   2001-06-18: changed the order of the tests, faster //
//                                                      //
// Acknowledgement: Many thanks to Pierre Terdiman for  //
// suggestions and discussions on how to optimize code. //
// Thanks to David Hunt for finding a ">="-bug!         //
//******************************************************//

// The code below is a modified version of the 
// AABB-triangle overlap test code
// by Tomas Akenine-Möller

#include "math/core.h"
#include <math.h>
#include <stdio.h>

using namespace NIBR;

#define X 0
#define Y 1
#define Z 2

#define VOXHALFSIZE      0.5
#define LARGEVOXHALFSIZE 0.5000005

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

bool planeBoxOverlap(float normal[3],float d,float halfVoxSize)
{
    
  float vmin[3],vmax[3];
  
  for(int q=0;q<3;q++)
  {
    if(normal[q]>0)
    {
      vmin[q]=-halfVoxSize;
      vmax[q]= halfVoxSize;
    }
    else
    {
      vmin[q]= halfVoxSize;
      vmax[q]=-halfVoxSize;
    }
  }
  if(  (dot(&normal[0],&vmin[0])+d) > 0 ) return false;
  if(  (dot(&normal[0],&vmax[0])+d) >=0 ) return true;

  return false;
}


//======================== X-tests ========================//
#define AXISTEST_X01(a, b, fab)             \
    p0 = a*v0[Y] - b*v0[Z];                    \
    p2 = a*v2[Y] - b*v2[Z];                    \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fab;   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fab)              \
    p0 = a*v0[Y] - b*v0[Z];                    \
    p1 = a*v1[Y] - b*v1[Z];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fab;   \
    if(min>rad || max<-rad) return 0;

//======================== Y-tests ========================//
#define AXISTEST_Y02(a, b, fab)             \
    p0 = -a*v0[X] + b*v0[Z];                   \
    p2 = -a*v2[X] + b*v2[Z];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fab;   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fab)              \
    p0 = -a*v0[X] + b*v0[Z];                   \
    p1 = -a*v1[X] + b*v1[Z];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fab;   \
    if(min>rad || max<-rad) return 0;

//======================== Z-tests ========================//
#define AXISTEST_Z12(a, b, fab)             \
    p1 = a*v1[X] - b*v1[Y];                    \
    p2 = a*v2[X] - b*v2[Y];                    \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
    rad = fab;   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fab)              \
    p0 = a*v0[X] - b*v0[Y];                \
    p1 = a*v1[X] - b*v1[Y];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fab;   \
    if(min>rad || max<-rad) return 0;

bool NIBR::triangleVoxelIntersection(int* ijk,float* c0, float* c1, float* c2)
{

   float v0[3],v1[3],v2[3];
   float min,max,d,p0,p1,p2,rad,fex,fey,fez;
   float normal[3],e0[3],e1[3],e2[3];

   vec3sub(v0,c0,ijk);
   vec3sub(v1,c1,ijk);
   vec3sub(v2,c2,ijk);

   vec3sub(e0,v1,v0);
   vec3sub(e1,v2,v1);
   vec3sub(e2,v0,v2);

   // Bullet 3:  //
   //  test the 9 tests first (this was faster) //
   fex = fabs(e0[X])*VOXHALFSIZE;
   fey = fabs(e0[Y])*VOXHALFSIZE;
   fez = fabs(e0[Z])*VOXHALFSIZE;
   AXISTEST_X01(e0[Z], e0[Y], fez+fey);
   AXISTEST_Y02(e0[Z], e0[X], fez+fex);
   AXISTEST_Z12(e0[Y], e0[X], fey+fex);

   fex = fabs(e1[X])*VOXHALFSIZE;
   fey = fabs(e1[Y])*VOXHALFSIZE;
   fez = fabs(e1[Z])*VOXHALFSIZE;
   AXISTEST_X01(e1[Z], e1[Y], fez+fey);
   AXISTEST_Y02(e1[Z], e1[X], fez+fex);
   AXISTEST_Z0(e1[Y], e1[X], fey+fex);

   fex = fabs(e2[X])*VOXHALFSIZE;
   fey = fabs(e2[Y])*VOXHALFSIZE;
   fez = fabs(e2[Z])*VOXHALFSIZE;
   AXISTEST_X2(e2[Z], e2[Y], fez+fey);
   AXISTEST_Y1(e2[Z], e2[X], fez+fex);
   AXISTEST_Z12(e2[Y], e2[X], fey+fex);

   // Bullet 1: //
   //  first test overlap in the {x,y,z}-directions //
   //  find min, max of the triangle each direction, and test for overlap in //
   //  that direction -- this is equivalent to testing a minimal AABB around //
   //  the triangle against the AABB //

   // test in X-direction //
   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if(min>VOXHALFSIZE || max<-VOXHALFSIZE) return 0;

   // test in Y-direction //
   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if(min>VOXHALFSIZE || max<-VOXHALFSIZE) return 0;

   // test in Z-direction //
   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if(min>VOXHALFSIZE || max<-VOXHALFSIZE) return 0;

   // Bullet 2: 
   //  test if the box intersects the plane of the triangle 
   //  compute plane equation of triangle: normal*x+d=0 
   cross(normal,e0,e1);
   d=-dot(&normal[0],&v0[0]);  // plane eq: normal.x+d=0 
   if(!planeBoxOverlap(normal,d,VOXHALFSIZE)) return 0;

   return 1;   // box and triangle overlaps 
}


bool NIBR::triangleLargerVoxelIntersection(int* ijk,float* c0, float* c1, float* c2)
{

   float v0[3],v1[3],v2[3];
   float min,max,d,p0,p1,p2,rad,fex,fey,fez;
   float normal[3],e0[3],e1[3],e2[3];

   vec3sub(v0,c0,ijk);
   vec3sub(v1,c1,ijk);
   vec3sub(v2,c2,ijk);

   vec3sub(e0,v1,v0);
   vec3sub(e1,v2,v1);
   vec3sub(e2,v0,v2);

   // Bullet 3:  //
   //  test the 9 tests first (this was faster) //
   fex = fabs(e0[X])*LARGEVOXHALFSIZE;
   fey = fabs(e0[Y])*LARGEVOXHALFSIZE;
   fez = fabs(e0[Z])*LARGEVOXHALFSIZE;
   AXISTEST_X01(e0[Z], e0[Y], fez+fey);
   AXISTEST_Y02(e0[Z], e0[X], fez+fex);
   AXISTEST_Z12(e0[Y], e0[X], fey+fex);

   fex = fabs(e1[X])*LARGEVOXHALFSIZE;
   fey = fabs(e1[Y])*LARGEVOXHALFSIZE;
   fez = fabs(e1[Z])*LARGEVOXHALFSIZE;
   AXISTEST_X01(e1[Z], e1[Y], fez+fey);
   AXISTEST_Y02(e1[Z], e1[X], fez+fex);
   AXISTEST_Z0(e1[Y], e1[X], fey+fex);

   fex = fabs(e2[X])*LARGEVOXHALFSIZE;
   fey = fabs(e2[Y])*LARGEVOXHALFSIZE;
   fez = fabs(e2[Z])*LARGEVOXHALFSIZE;
   AXISTEST_X2(e2[Z], e2[Y], fez+fey);
   AXISTEST_Y1(e2[Z], e2[X], fez+fex);
   AXISTEST_Z12(e2[Y], e2[X], fey+fex);

   // Bullet 1: //
   //  first test overlap in the {x,y,z}-directions //
   //  find min, max of the triangle each direction, and test for overlap in //
   //  that direction -- this is equivalent to testing a minimal AABB around //
   //  the triangle against the AABB //

   // test in X-direction //
   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if(min>LARGEVOXHALFSIZE || max<-LARGEVOXHALFSIZE) return 0;

   // test in Y-direction //
   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if(min>LARGEVOXHALFSIZE || max<-LARGEVOXHALFSIZE) return 0;

   // test in Z-direction //
   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if(min>LARGEVOXHALFSIZE || max<-LARGEVOXHALFSIZE) return 0;

   // Bullet 2: 
   //  test if the box intersects the plane of the triangle 
   //  compute plane equation of triangle: normal*x+d=0 
   cross(normal,e0,e1);
   d=-dot(&normal[0],&v0[0]);  // plane eq: normal.x+d=0 
   if(!planeBoxOverlap(normal,d,LARGEVOXHALFSIZE)) return 0;

   return 1;   // box and triangle overlaps 
}
