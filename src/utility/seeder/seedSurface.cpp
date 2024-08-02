#include "seedSurface.h"
#include "math/triangle.h"
#include <cmath>

using namespace NIBR;

SeederOutputState SeedSurface::getSeed(float* p, int t) {
    float tmp[3];
    return getSeed(p, tmp, t);
}

SeederOutputState SeedSurface::getSeed(float* p, float* dir, int t) {

    SeederOutputState state = checkSeedingLimits();

    if (state!=SEED_OK)
        return state;

    auto V = seed_surf->vertices;
    auto F = seed_surf->faces;
    auto Nv = seed_surf->normalsOfVertices;
    auto Nf = seed_surf->normalsOfFaces;

    double T[3];

    while (true) {

        int f = doRandomThings[t].uniform_int();

        // Always do rejection sampling because faces_vec_dens is scaled with face area
        if (doRandomThings[t].uniform_01()*max4rs > faces_vec_dens[f])
            continue;        

        float* a = V[F[f][0]];
        float* b = V[F[f][1]];
        float* c = V[F[f][2]];

        doRandomThings[t].getARandomPointWithinTriangle(p, a, b, c);

        // We will assume that the surface normals are pointing outwards.
        // So the generated points will be on side of the mesh opposite to where the normal is pointing
        // We will make sure that the points are not exactly on the mesh, i.e. dist != 0.0f,
        // but within surface border, i.e. (0.0 SURFTHICKNESS]
        vec3sub(&T[0],a,p);
        double dist = dot(Nf[f],&T[0]);

        // Try and project the point inside if needed
        if (dist < 0.0f) {
            p[0] -= 2.0 * dist * Nf[f][0];
            p[1] -= 2.0 * dist * Nf[f][1];
            p[2] -= 2.0 * dist * Nf[f][2];
            vec3sub(&T[0],p,a);
            dist = dot(Nf[f],&T[0]);
        }

        // Make sure that the point is within border.
        if ((dist <= 0.0f) && (dist > SURFTHICKNESS))
            continue;
            
        if (surfNorm) {
            barycentricInterp(dir, p, a, b, c, Nv[F[f][0]], Nv[F[f][1]], Nv[F[f][2]]);
            normalize(dir);
        }
        break;

    }

    curSeed++;

    return SEED_OK;

}


void SeedSurface::computeSeedCountAndDensity() {

    seed_surf->calcArea();
    double surfArea = seed_surf->area;
    
    if (hasDensity) {
        count   = density*surfArea;
    } else {
        density = (surfArea>0) ? double(count)/surfArea : 0;
    }

    if (!useDensInp) {
        faces_vec_dens.clear();
        max4rs = 0;
        for (int n=0; n<seed_surf->nf; n++) {
            faces_vec_dens.push_back(seed_surf->areasOfFaces[n]);
            if (seed_surf->areasOfFaces[n]>max4rs)
                max4rs = seed_surf->areasOfFaces[n];
        }
    }

}

void SeedSurface::computeMaxPossibleSeedCount() {
    if (seed_surf->nf==0) {
        maxPossibleSeedCount = 0;
    } else if (max4rs==0) {
        maxPossibleSeedCount = 0;
    } else {
        maxPossibleSeedCount = INT_MAX;
    }
}

void SeedSurface::useSurfNorm(bool use) {

    if (use==true) {

        seed_surf->calcNormalsOfVertices();

        if (mode==SEED_SURFACE_MASK)    mode=SEED_SURFACE_MASK_WITH_DIRECTIONS;
        if (mode==SEED_SURFACE_RS)      mode=SEED_SURFACE_RS_WITH_DIRECTIONS;

        surfNorm = true;

    } else {

        if (mode==SEED_SURFACE_MASK_WITH_DIRECTIONS)    mode=SEED_SURFACE_MASK;
        if (mode==SEED_SURFACE_RS_WITH_DIRECTIONS)      mode=SEED_SURFACE_RS;

        surfNorm = false;
    }

}

bool SeedSurface::setSeed(Surface *surf) {
    
    seed_surf     = surf;
    surfNorm      = false;
    useDensInp    = false;

    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeSeedCountAndDensity();
    computeMaxPossibleSeedCount();

    mode          = SEED_SURFACE_MASK;
    return true;
}

bool SeedSurface::useDensity(std::vector<float>& density_vec) {

    if ((int)density_vec.size() != seed_surf->nf) {
        disp (MSG_ERROR, "Surface sampling density size does not match the surface face count");
        return false;
    }

    // Compute relative densities per face if not already provided
    faces_vec_dens.clear();
    faces_vec_dens = density_vec;
    useDensInp = true;

    max4rs = 0;
    for (float val : faces_vec_dens) {
        if (val>max4rs)
            max4rs = val;
    }

    if (mode==SEED_SURFACE_MASK)                    mode = SEED_SURFACE_RS;
    if (mode==SEED_SURFACE_MASK_WITH_DIRECTIONS)    mode = SEED_SURFACE_RS_WITH_DIRECTIONS;

    return true;

}

void SeedSurface::setNumberOfThreads(int n) {

    threadCount = n;

    if (doRandomThings!=NULL) {
        delete[] doRandomThings;
        doRandomThings = NULL;
    }
        
    doRandomThings = new RandomDoer[n];

    for (int i=0; i<n; i++) {
        doRandomThings[i].init_uniform_int(seed_surf->nf-1);
    }

}