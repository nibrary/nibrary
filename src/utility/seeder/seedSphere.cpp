#include "seedSphere.h"

using namespace NIBR;

SeederOutputState SeedSphere::getSeed(float* p,int t) {

    SeederOutputState state = checkSeedingLimits();

    if (state!=SEED_OK)
        return state;

    // Points are always generated within the sphere.
    // They are never placed exactly on the sphere surface.
    doRandomThings[t].getARandomPointWithinSphere(p,seed_radius);
    p[0] += seed_center[0];
    p[1] += seed_center[1];
    p[2] += seed_center[2];
    curSeed++;

    return SEED_OK;

}

SeederOutputState SeedSphere::getSeed(float* p, float*, int t) {
    return getSeed(p,t);
}

void SeedSphere::computeSeedCountAndDensity() {
    double sphVol = 1.33333333333*PI*seed_radius*seed_radius*seed_radius;
    if (hasDensity) {
        count   = density*sphVol;
    } else {
        density = (sphVol==0) ? 0 : count/sphVol;
    }
}

void SeedSphere::computeMaxPossibleSeedCount() {
    maxPossibleSeedCount = INT_MAX;
}

bool SeedSphere::setSeed(Point  xyz, float r) {
    return setSeed(xyz.x,xyz.y,xyz.z,r);
}

bool SeedSphere::setSeed(float* p, float r) {
    return setSeed(p[0],p[1],p[2],r);
}

bool SeedSphere::setSeed(float x, float y, float z, float r) {

    seed_center[0] = x;
    seed_center[1] = y;
    seed_center[2] = z;
    seed_radius    = r;

    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();

    mode = SEED_SHPERE;
    return true;
}

void SeedSphere::setNumberOfThreads(int n) {

    threadCount = n;

    if (doRandomThings!=NULL) {
        delete[] doRandomThings;
        doRandomThings = NULL;
    }
        
    doRandomThings = new RandomDoer[n];

}
