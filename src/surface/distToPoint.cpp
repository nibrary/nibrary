#include "surface.h"

float NIBR::Surface::squaredDistToPoint(float *p) {

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXf sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXf C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    return sqrD(0);
    
}


float NIBR::Surface::distToPoint(float *p) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return NAN;
    }

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXf sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXf C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    float dist = isPointInside(p) ? std::sqrt(sqrD(0)) : -std::sqrt(sqrD(0));
    
    return dist;
    
}


float NIBR::Surface::squaredDistToPoint(float *p, int& faceInd) {

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXf sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXf C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    faceInd = I(0);

    return sqrD(0);
    
}

float NIBR::Surface::distToPoint(float *p, int& faceInd) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return NAN;
    }

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXf sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXf C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    float dist = isPointInside(p) ? std::sqrt(sqrD(0)) : -std::sqrt(sqrD(0));
    
    faceInd = I(0);
    return dist;
    
}

float NIBR::Surface::squaredDistToPoint(float *p, int& faceInd, float* closestPoint) {

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXf sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXf C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    closestPoint[0] = C(0,0);
    closestPoint[1] = C(0,1);
    closestPoint[2] = C(0,2);

    faceInd = I(0);

    return sqrD(0);
    
}

float NIBR::Surface::distToPoint(float *p, int& faceInd, float* closestPoint) {

    if (!enabledPointCheck) {
        disp(MSG_FATAL, "enablePointCheck is not initialized");
        return NAN;
    }

    Eigen::MatrixXf point(1,3);

    point(0,0) = p[0];
    point(0,1) = p[1];
    point(0,2) = p[2];

    Eigen::VectorXf sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXf C;
    AABB_tree.squared_distance(V,F,point,sqrD,I,C);

    float dist = isPointInside(p) ? std::sqrt(sqrD(0)) : -std::sqrt(sqrD(0));
    
    closestPoint[0] = C(0,0);
    closestPoint[1] = C(0,1);
    closestPoint[2] = C(0,2);

    faceInd = I(0);

    return dist;
    
}