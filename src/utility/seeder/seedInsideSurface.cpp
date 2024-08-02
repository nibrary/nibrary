#include "seedInsideSurface.h"
#include "surface/surface2imageMapper.h"
#include "image/image_operators.h"

using namespace NIBR;

SeedInsideSurface::~SeedInsideSurface() {

    if (doRandomThings!=NULL) {
        delete[] doRandomThings; 
        doRandomThings=NULL;
    }
}

SeederOutputState SeedInsideSurface::getSeed(float* p, int t) {

    SeederOutputState state = checkSeedingLimits();

    if (state!=SEED_OK)
        return state;

    while (true) {
        
        seed_surf->maskAndBoundary.to_xyz(seed_indices[doRandomThings[t].uniform_int()],p);
        doRandomThings[t].randomizeWithinVoxel(p,seed_surf->maskAndBoundary.pixDims);

        // Make sure that the point is not exactly on the border
        if (seed_surf->distToPoint(p) > 0.0f)
            break;

    }

    curSeed++;

    return SEED_OK;

}

SeederOutputState SeedInsideSurface::getSeed(float* p, float*, int t) {
    return getSeed(p, t);
}

void SeedInsideSurface::computeSeedCountAndDensity() {

    double imgVol = 0;

    imgVol = seed_surf->maskAndBoundary.pixDims[0]*seed_surf->maskAndBoundary.pixDims[1]*seed_surf->maskAndBoundary.pixDims[2]*double(seed_indices.size());

    if (hasDensity) {
        count   = density*imgVol;
    } else {
        density = double(count)/imgVol;
    }
}

void SeedInsideSurface::computeMaxPossibleSeedCount() {

    maxPossibleSeedCount = (seed_indices.size()==0) ? 0 : INT_MAX;
}

bool SeedInsideSurface::setSeed(Surface *surf, float discRes) {

    disp(MSG_DEBUG,"setSeed");

    seed_surf   = surf;

    seed_surf->enablePointCheck(discRes);

    seed_indices = getNonZeroIndices(&(seed_surf->maskAndBoundary));
     
    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();

    mode = SEED_INSIDE_SURFACE;

    disp(MSG_DEBUG,"Done setSeed");
    return true;
}

void SeedInsideSurface::setNumberOfThreads(int n) {

    threadCount = n;

    if (doRandomThings!=NULL) {
        delete[] doRandomThings;
        doRandomThings = NULL;
    }
        
    doRandomThings = new RandomDoer[n];

    for (int i=0; i<n; i++) {
        doRandomThings[i].init_uniform_int(seed_indices.size()-1);
    }

}