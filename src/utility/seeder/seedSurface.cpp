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

        // Sampling using CDF
        double rand_val = doRandomThings[t].uniform_01();
        auto it = std::lower_bound(cdf.begin(), cdf.end(), rand_val);
        int cdfInd = (it == cdf.end()) ? cdf.size() - 1 : std::distance(cdf.begin(), it);
        int f = nonZeroFaces[cdfInd];
   
        float* a = V[F[f][0]];
        float* b = V[F[f][1]];
        float* c = V[F[f][2]];

        bool seedFound = false;

        for (int trial = 0; trial < 100; trial++) {

            doRandomThings[t].getARandomPointWithinTriangle(p, a, b, c);

            // We will assume that the surface normals are pointing outwards.
            // So the generated points will be on side of the mesh opposite to where the normal is pointing
            // We will make sure that the points are not exactly on the mesh, i.e. dist != 0.0f,
            // but within surface border, i.e. (0.0 SURFTHICKNESS]
            vec3sub(T,a,p);
            double dist = dot(Nf[f],T);

            // Try and project the point inside if needed
            if (dist < 0.0f) {
                disp(MSG_DEBUG,"Fixed negative dist while seeding");
                p[0] -= (dist - 0.2 * SURFTHICKNESS) * Nf[f][0];
                p[1] -= (dist - 0.2 * SURFTHICKNESS) * Nf[f][1];
                p[2] -= (dist - 0.2 * SURFTHICKNESS) * Nf[f][2];
                vec3sub(T,a,p);
                dist = dot(Nf[f],T);
            }

            if (dist == 0.0f) {
                disp(MSG_DEBUG,"Fixed zero dist while seeding");
                p[0] -= 0.2 * SURFTHICKNESS * Nf[f][0];
                p[1] -= 0.2 * SURFTHICKNESS * Nf[f][1];
                p[2] -= 0.2 * SURFTHICKNESS * Nf[f][2];
                vec3sub(T,a,p);
                dist = dot(Nf[f],T);
            }

            // Make sure that the point is within border.
            if ((dist < (0.1*SURFTHICKNESS) ) || (dist > (0.9*SURFTHICKNESS) ))
                continue;

            if (seed_surf->isPointInside(p) == false)
                continue;
                
            disp(MSG_DEBUG,"Seed dist is: %.12f", dist);

            if (surfNorm) {
                barycentricInterp(dir, p, a, b, c, Nv[F[f][0]], Nv[F[f][1]], Nv[F[f][2]]);
                normalize(dir);
            }

            seedFound = true;
            break;

        }

        if (seedFound) break;


    }

    curSeed++;

    return SEED_OK;

}

void SeedSurface::computeCDF() {
    cdf.clear();
    nonZeroFaces.clear();
    totalDensity = 0.0;
    
    // Populate nonZeroFaces and compute cumulative densities
    for (int n = 0; n < seed_surf->nf; ++n) {
        double d = faces_vec_dens[n]*seed_surf->areasOfFaces[n];
        if (d > 0.0f) {
            totalDensity += d;
            cdf.push_back(totalDensity);
            nonZeroFaces.push_back(n);
        }
    }

    // Normalize CDF to make the last element equal to 1
    if (totalDensity > 0.0) {
        for (auto& val : cdf) {
            val /= totalDensity;
        }
    }

    if (!cdf.empty()) {
        cdf.back() = 1.0;
    }
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
        for (int n = 0; n < seed_surf->nf; n++) {
            faces_vec_dens.push_back(1.0f);
        }
    }

    computeCDF();

}

void SeedSurface::computeMaxPossibleSeedCount() {
    if (seed_surf->nf == 0) {
        maxPossibleSeedCount = 0;
    } else if (totalDensity == 0.0) {
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
    seed_surf->enablePointCheck( (surf->pointCheckGridRes > 0) ? surf->pointCheckGridRes : 1.0f );

    mode          = SEED_SURFACE_MASK;
    return true;
}

bool SeedSurface::useDensity(std::vector<float>& density_vec) {

    if ((int)density_vec.size() != seed_surf->nf) {
        disp (MSG_ERROR, "Surface sampling density size does not match the surface face count");
        return false;
    }

    faces_vec_dens.clear();
    for (int n = 0; n < seed_surf->nf; n++) {
        if (density_vec[n] < 0.0f) {
            disp(MSG_ERROR, "Negative density values are not allowed.");
            return false;
        }
        faces_vec_dens.push_back(density_vec[n]);
    }
    useDensInp = true;

    computeCDF();

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