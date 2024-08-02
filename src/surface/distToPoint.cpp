#include "surface.h"

double NIBR::Surface::squaredDistToPoint(float *p, int& faceInd, float* closestPoint) {

    Eigen::MatrixXd point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    closestPoint[0] = C(0,0);
    closestPoint[1] = C(0,1);
    closestPoint[2] = C(0,2);

    faceInd = I(0);

    return sqrD(0);
    
}

double NIBR::Surface::squaredDistToPoint(float *p, int& faceInd) {
    float closestPoint[3];
    return squaredDistToPoint(p,faceInd,&closestPoint[0]);
}

double NIBR::Surface::squaredDistToPoint(float *p) {
    int faceInd;
    float closestPoint[3];
    return squaredDistToPoint(p,faceInd,&closestPoint[0]);
}





double NIBR::Surface::distToPoint(float *p, int& faceInd, float* closestPoint) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return NAN;
    }

    double dist;
    bool isInside = isPointInside(p,dist,faceInd,closestPoint);

    if ((dist != DBL_MAX) && (dist != DBL_MIN)) return dist;

    Eigen::MatrixXd point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    dist = isInside ? std::sqrt(sqrD(0)) : -std::sqrt(sqrD(0));
    
    closestPoint[0] = C(0,0);
    closestPoint[1] = C(0,1);
    closestPoint[2] = C(0,2);

    faceInd = I(0);

    return dist;
    
}

double NIBR::Surface::distToPoint(float *p, int& faceInd) {
    float closestPoint[3];
    return distToPoint(p,faceInd,&closestPoint[0]);
}


double NIBR::Surface::distToPoint(float *p) {
    int faceInd;
    float closestPoint[3];
    return distToPoint(p,faceInd,&closestPoint[0]);
}