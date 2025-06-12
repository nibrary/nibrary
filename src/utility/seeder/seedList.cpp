#include "seedList.h"

using namespace NIBR;

SeederOutputState SeedList::getSeed(float* p, int) {

    SeederOutputState state = checkSeedingLimits();

    if (state!=SEED_OK)
        return state;

    p[0] = seed_coordinates->at(curSeed)[0];
    p[1] = seed_coordinates->at(curSeed)[1];
    p[2] = seed_coordinates->at(curSeed)[2];
    
    curSeed++;

    return SEED_OK;

}

SeederOutputState SeedList::getSeed(float* p, float* dir, int) {

    SeederOutputState state = checkSeedingLimits();

    if (state!=SEED_OK)
        return state;
    
    p[0]    = seed_coordinates->at(curSeed)[0];
    p[1]    = seed_coordinates->at(curSeed)[1];
    p[2]    = seed_coordinates->at(curSeed)[2];

    if (!seed_directions->empty()) {
        dir[0]  = seed_directions->at(curSeed)[0];
        dir[1]  = seed_directions->at(curSeed)[1];
        dir[2]  = seed_directions->at(curSeed)[2];
    }

    curSeed++;

    return SEED_OK;
}

void SeedList::computeSeedCountAndDensity() {
    if (hasDensity) {
        disp(MSG_WARN,"Density is ignored when seed coordinates are explicity defined.");
        count   = maxPossibleSeedCount;
    } else {
        density = 0;
    }
}

void SeedList::computeMaxPossibleSeedCount() {
    maxPossibleSeedCount = seed_coordinates->size();
}

bool SeedList::setSeed(std::vector<Point3D>& p) {

    seed_coordinates = &p;

    threadCount = 1;
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();

    mode = SEED_LIST;
    return true;
}

bool SeedList::setSeed(std::vector<Point3D>& p, std::vector<Point3D>& dir) {

    seed_coordinates = &p;
    seed_directions  = &dir;

    if (p.size() != dir.size()) {
        disp(MSG_ERROR,"Size of seed coordinates and directions must be same");
        return false;
    }
    
    threadCount = 1;
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();

    mode = SEED_LIST_WITH_DIRECTIONS;
    return true;
}

void SeedList::setNumberOfThreads(int) {
    disp(MSG_WARN,"Can't multithread when seeding from a list.");
    threadCount = 1;
    return;
}