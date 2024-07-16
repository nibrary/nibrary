#include "seed.h"
#include "../tracker/tracker.h"
#include <climits>

using namespace NIBR;

bool Seed::setSeed(PathwayRule rule) {

    if (rule.type!=seed) {
        disp(MSG_ERROR, "Unexpected seed rule type");
        return false;
    }

    rule.surfaceFieldFile4FaceDens  = surf_faceDensity_filename;
    rule.surfaceFieldFile4VertDens  = surf_vertDensity_filename;
    rule.surfaceFieldFile4DensDtype = surf_densityFile_dataType;
    rule.surfaceFieldName4Dens      = surf_density_fieldname;

    rule.surface4SeedUseNormForDir  = useSurfNorm;
    rule.surfaceUseDim                 = useInsideSurf ? surf_useDim_3D : surf_useDim_2D;

    disp(MSG_DETAIL, "Adding seed");
    if(!TRACKER::pw.add(rule)) {
        disp(MSG_ERROR, "Failed at adding seed in pathway rules");
        return false;
    }
    disp(MSG_DETAIL, "Done");

    seeder = TRACKER::pw.getSeeder();
    
    if (seeder==NULL) {
        disp(MSG_ERROR, "Failed during setting seed");
        return false;
    }

    return true;
}

SeederOutputState Seed::getSeed(float* p, int tID) {
    if (isReady) {
        SeederOutputState seedState = seeder->getSeed(p,tID);
        return seedState;
    } else {
        disp(MSG_ERROR,"Seed is not updated");
        return SEED_ERROR;
    }
}

SeederOutputState Seed::getSeed(float* p,float *dir, int tID) {

    if (isReady) {
        SeederOutputState seedState = seeder->getSeed(p,dir,tID);
        return seedState;
    } else {
        disp(MSG_ERROR,"Seed is not updated");
        return SEED_ERROR;
    }

}

void Seed::clear() {

    surf_faceDensity_filename   = "";
    surf_vertDensity_filename   = "";
    surf_densityFile_dataType   = "";
    surf_density_fieldname      = "";

    useSurfNorm      = false;
    useInsideSurf    = true;

    sCount           = INT64_MAX;
    sDensity         = 0;
    sHasDensity      = false;
    trials           = 1;

    isReady          = false;

}

bool Seed::update() {

    isReady = false;

    if (seeder==NULL)
        return false;

    if (sHasDensity) {
        seeder->setDensity(sDensity);
    } else {
        seeder->setCount(sCount);
    }
    seeder->computeSeedCountAndDensity();
    seeder->setNumberOfThreads(TRACKER::nThreads);
    
    sCount    = seeder->count;
    sDensity  = seeder->density;

    isReady = true;
    return true;
    
}



void Seed::print() {

    if (NIBR::VERBOSE()<VERBOSE_INFO) {
        return;
    }
    
	disp(MSG_INFO,"SEEDING OPTIONS");

    std::cout << "\033[32m";

    if (sCount==INT64_MAX) {
        std::cout << "seed_count         : INF" << std::endl;
        std::cout << "seed_density       : INF" << std::endl;
    } else {
        std::cout << "seed_count         : " << sCount   << std::endl;
        std::cout << "seed_density       : " << sDensity << std::endl;
    }
    
    std::cout << "seed_trials        : " << trials   << std::endl;

    std::cout << "\033[0m";

}
