#include "surface2imageMapper.h"
#include "surface_operators.h"
#include "image/image_math.h"
#include <cstdint>
#include <utility>
#include <atomic>

using namespace NIBR;

void NIBR::surfaceMask(NIBR::Image<bool>* img, NIBR::Surface* surf)
{
    disp(MSG_DEBUG,"surfaceMask()");

    if (surf->interpretAs2D) return; // In this case, there is only BOUNDARY and OUTSIDE

    surf->getClosedAndOpenComponents();
    Surface& closed = surf->compClosedAndOpen[0];
    closed.prepIglAABBTree();

    auto getMask = [&](NIBR::MT::TASK task)->void {

        int64_t ind = task.no;

        float p[3];
        img->to_xyz(ind,p);

        if (closed.isPointInside_basedOnWindingNumber(p)) {
            img->data[ind]=INSIDE;
        }

    };
    
    NIBR::MT::MTRUN(img->numel, "Computing mask", getMask);

}



void NIBR::surfaceMaskWithBoundary(NIBR::Image<int8_t>* img, NIBR::Surface* surf)
{

    surf->getClosedAndOpenComponents();
    Surface& closed = surf->compClosedAndOpen[0];
    closed.prepIglAABBTree();

    auto getMaskAndBoundary = [&](NIBR::MT::TASK task)->void {

        int64_t ind = task.no;

        if (img->data[ind]!=0) {
            img->data[ind]=BOUNDARY;
            return;
        }

        if (surf->interpretAs2D == false) {
            float p[3];
            img->to_xyz(ind,p);
            if (closed.isPointInside_basedOnWindingNumber(p)) {
                img->data[ind]=INSIDE;
            }
        }

    };
    
    NIBR::MT::MTRUN(img->numel, "Computing mask and boundary", getMaskAndBoundary);

}



void NIBR::surfaceEDT(NIBR::Image<float>* img, NIBR::Surface* surf)
{

    surf->enablePointCheck(img->smallestPixDim);

    auto getEDT = [&](NIBR::MT::TASK task)->void {

        int64_t ind = task.no;
        float p[3];

        img->to_xyz(ind,p);

        img->data[ind] = surf->distToPoint(p);

    };
    
    NIBR::MT::MTRUN(img->numel, "Computing EDT", getEDT);
    
}

void NIBR::surfacePVF(NIBR::Image<float>* img, NIBR::Surface* surf)
{
    float div    = 4;
    float divVol = 1.0f/(div*div*div);
    float step   = 1.0f/(div+1.0f);

    surf->enablePointCheck(img->smallestPixDim);

    auto estimPartVol = [&](NIBR::MT::TASK task)->void {
        
        int64_t ind = task.no;
        int64_t ijk[3];
        float   xyz[3];

        img->ind2sub(ind,ijk);

        // Voxel is not on the BOUNDARY
        if (img->data[ind] == 0) {

            img->to_xyz(ijk,xyz);

            if (surf->isPointInside(xyz)) { // voxel is inside
                img->data[ind] = 1.0f;
            }

            return;
        
        }

        // Voxel is on the BOUNDARY
        int partVol = 0;

        float cur_ijk[3];

        for (float i=(float(ijk[0])-0.5+step); i<(float(ijk[0])+0.5-step+EPS4); i=i+step)
            for (float j=(float(ijk[1])-0.5+step); j<(float(ijk[1])+0.5-step+EPS4); j=j+step)
                for (float k=(float(ijk[2])-0.5+step); k<(float(ijk[2])+0.5-step+EPS4); k=k+step) {

                    cur_ijk[0]=i;
                    cur_ijk[1]=j;
                    cur_ijk[2]=k;
                    
                    img->to_xyz(cur_ijk,xyz);

                    if (surf->isPointInside(xyz)) { // voxel is inside
                        partVol++;
                    }

                }


        if (partVol>0) {
            disp(MSG_DEBUG,"ijk=[%d,%d,%d], pv=%.2f",ijk[0],ijk[1],ijk[2],partVol);
            img->data[ind] = float(partVol)*divVol;
        } else {
            img->data[ind] = 0;
        }


    };

    NIBR::MT::MTRUN(img->numel, "Estimating partial volumes", estimPartVol);

}

void NIBR::surfaceMAT(NIBR::Image<float>* img, NIBR::Surface* surf)
{
    surf->enablePointCheck(img->smallestPixDim);

    struct ComparePairs {
    bool operator()(const std::pair<int64_t, float>& a, const std::pair<int64_t, float>& b) const {
        return a.second > b.second;  // Compare based on the distance in descending order
        }
    };

    std::vector<std::set<std::pair<int64_t, float>, ComparePairs>> edt;
    edt.resize(MT::MAXNUMBEROFTHREADS());

    float halfVox   = img->pixDims[0]*std::sqrt(3)+EPS4;

    auto getEDTpair = [&](NIBR::MT::TASK task)->void {

        int64_t ind = task.no;
        float p[3];

        img->to_xyz(ind,p);

        if (surf->isPointInside(p)) {

            float radius = surf->distToPoint(p);

            if (radius > halfVox) {
                edt[task.threadId].insert(std::make_pair(ind, radius));
            } else {
                // These are already the smallest values. Because in image space they will correspond to a single voxel, they can't overwrite any other voxel.
                img->data[ind] = radius;
            }
            
        }

    };    
    
    NIBR::MT::MTRUN(img->numel, "Indexing EDT", getEDTpair);

    // Merge threads
    for (int t = 1; t < MT::MAXNUMBEROFTHREADS(); t++) {
        std::set_union(
            edt[0].begin(), edt[0].end(),
            edt[t].begin(), edt[t].end(),
            std::inserter(edt[0], edt[0].end()),
            ComparePairs()
        );
        edt[t].clear();
    }

    // We will save 1/8 of the biggest sphere
    std::vector<std::vector<std::vector<float>>> s;

    // First process the biggest sphere
    auto it = edt[0].begin();
    if (!edt[0].empty()) {

        float maxRadius = it->second;
        int R           = std::ceil(maxRadius/img->pixDims[0]);    // We assume isotropic voxel size

        s.resize(R+1);
        for (int i = 0; i <= R; i++) {
            s[i].resize(R+1);
            for (int j = 0; j <= R; j++) {
                s[i][j].resize(R+1);
                for (int k = 0; k <= R; k++) {
                    float r = std::sqrt(float(i*i + j*j + k*k))*img->pixDims[0];
                    s[i][j][k] = r;
                }
            }
        }
    }

    std::vector<std::pair<int64_t, float>> vedt;
    while(it!=edt[0].end()) {
        vedt.push_back(*it);
        it++;
    }
    edt[0].clear();


    std::vector<std::pair<int64_t, float>> v;

    auto getMaximalSpheres = [&](NIBR::MT::TASK task)->void {

        float radius = vedt[task.no].second;
        int R        = std::ceil(radius/img->pixDims[0]);    // We assume isotropic voxel size

        int64_t ci,cj,ck; // i,j,k of sphere center
        img->ind2sub(vedt[task.no].first,ci,cj,ck);

        bool onMedialAxis = false;

        for (int i = ci-R; i <= ci+R; i++)
            for (int j = cj-R; j <= cj+R; j++)
                for (int k = ck-R; k <= ck+R; k++) 
            {
                if (img->isInside(i, j, k)) {

                    int64_t index = img->sub2ind(i,j,k);
                    float r = s[std::abs(ci - i)][std::abs(cj - j)][std::abs(ck - k)];

                    if (r<=radius) {
                        if (r > img->data[index]) {
                            img->data[index] = radius;
                            if (!onMedialAxis)
                                onMedialAxis = true;
                        }
                    }

                }
            }

        
        if (onMedialAxis) {
            v.push_back(vedt[task.no]);
        }

    };

    // Maybe this is possible to parallelize but I ran out of ideas.
    // The problem is the subsequent computation of medial axis.
    // This requires that img->data[index] always has the largest value.
    // When parallelized, a small value can be written on medial axis
    // which could actually be overwritten later on, but there is way
    // to keep record of that. This means keeping record of which sphere
    // center wrote which voxel!
    NIBR::MT::MTRUN(vedt.size(), 1, "Computing medial axis", getMaximalSpheres);
    
    img->deallocData();
    img->allocData();

    auto writeMAT = [&](NIBR::MT::TASK task)->void {
        img->data[v[task.no].first] = v[task.no].second;
    };
    NIBR::MT::MTRUN(v.size(), writeMAT);
}

void NIBR::surfaceMIS(NIBR::Image<float>* img, NIBR::Surface* surf)
{
    surf->enablePointCheck(img->smallestPixDim);

    struct ComparePairs {
    bool operator()(const std::pair<int64_t, float>& a, const std::pair<int64_t, float>& b) const {
        return a.second > b.second;  // Compare based on the distance in descending order
        }
    };

    std::vector<std::set<std::pair<int64_t, float>, ComparePairs>> edt;
    edt.resize(MT::MAXNUMBEROFTHREADS());

    float halfVox   = img->pixDims[0]*std::sqrt(3)+EPS4;

    auto getEDTpair = [&](NIBR::MT::TASK task)->void {

        int64_t ind = task.no;
        float p[3];

        img->to_xyz(ind,p);

        if (surf->isPointInside(p)) {

            float radius = surf->distToPoint(p);

            if (radius > halfVox) {
                edt[task.threadId].insert(std::make_pair(ind, radius));
            } else {
                // These are already the smallest values. Because in image space they will correspond to a single voxel, they can't overwrite any other voxel.
                img->data[ind] = radius;
            }
            
        }

    };    
    
    NIBR::MT::MTRUN(img->numel, "Indexing EDT", getEDTpair);

    // Merge threads
    for (int t = 1; t < MT::MAXNUMBEROFTHREADS(); t++) {
        std::set_union(
            edt[0].begin(), edt[0].end(),
            edt[t].begin(), edt[t].end(),
            std::inserter(edt[0], edt[0].end()),
            ComparePairs()
        );
        edt[t].clear();
    }

    // We will save 1/8 of the biggest sphere
    std::vector<std::vector<std::vector<float>>> s;

    // First process the biggest sphere
    auto it = edt[0].begin();
    if (!edt[0].empty()) {

        float maxRadius = it->second;
        int R           = std::ceil(maxRadius/img->pixDims[0]);    // We assume isotropic voxel size

        s.resize(R+1);
        for (int i = 0; i <= R; i++) {
            s[i].resize(R+1);
            for (int j = 0; j <= R; j++) {
                s[i][j].resize(R+1);
                for (int k = 0; k <= R; k++) {
                    float r = std::sqrt(float(i*i + j*j + k*k))*img->pixDims[0];
                    s[i][j][k] = r;
                }
            }
        }
    }

    std::vector<std::pair<int64_t, float>> vedt;
    while(it!=edt[0].end()) {
        vedt.push_back(*it);
        it++;
    }
    edt[0].clear();

    std::mutex* ijkAvailability = new std::mutex[img->voxCnt];

    auto getMaximalSpheres = [&](NIBR::MT::TASK task)->void {

        float radius = vedt[task.no].second;
        int R        = std::ceil(radius/img->pixDims[0]);    // We assume isotropic voxel size

        int64_t ci,cj,ck; // i,j,k of sphere center
        img->ind2sub(vedt[task.no].first,ci,cj,ck);

        for (int i = ci-R; i <= ci+R; i++)
            for (int j = cj-R; j <= cj+R; j++)
                for (int k = ck-R; k <= ck+R; k++) 
            {
                if (img->isInside(i, j, k)) {

                    int64_t index = img->sub2ind(i,j,k);
                    float r = s[std::abs(ci - i)][std::abs(cj - j)][std::abs(ck - k)];

                    if (r<=radius) {
                        {
                            std::lock_guard<std::mutex> lock(ijkAvailability[index]);

                            if (r > img->data[index]) {
                                img->data[index] = radius;
                            }
                            
                        }
                    }
                }
            }        

    };

    // For MIS this can be parallized.
    NIBR::MT::MTRUN(vedt.size(), "Computing maximum inscribed spheres", getMaximalSpheres);

    delete[] ijkAvailability;

}