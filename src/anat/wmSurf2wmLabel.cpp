#include "anatSurf.h"

void NIBR::wmSurf2wmLabel(
        NIBR::Image<int8_t>& wmLabel, 
        NIBR::Surface* wmSurf, 
        NIBR::SurfaceField* wmSurfLabels,
        float sulcalDepth,                  // determines WM_SULCAL_FUNDUS, e.g., 5 mm
        float wallDepth,                    // determines WM_SULCAL_WALL and WM_GYRAL_CROWN, e.g., 1 mm
        float minThicknessOfSuperficialWM,  // superficialWM will be at least this thick, e.g., 1 mm
        float maxThicknessOfSuperficialWM)  // superficialWM will be at most this thick, e.g., 5 mm
{

    disp(MSG_DEBUG,"wmSurf2wmLabel()");

    NIBR::Image<float> wmDepth;
    wmDepth.createFromTemplate(wmLabel,true);

    NIBR::Image<bool> wmMask;
    wmMask.createFromTemplate(wmLabel,true);

    disp(MSG_DEBUG,"continuing wmSurf2wmLabel()...mapSurface2Image()");
    mapSurface2Image(wmSurf,&wmMask,0,NULL,NULL,MASK);

    disp(MSG_DEBUG,"continuing wmSurf2wmLabel()...selectVertices()");
    NIBR::Surface surf;
    selectVertices(&surf, wmSurf, wmSurfLabels);


    // Below is a modified version of surfaceMIS
    disp(MSG_DEBUG,"continuing wmSurf2wmLabel()...enablePointCheck()");
    surf.enablePointCheck(wmLabel.smallestPixDim);

    struct ComparePairs {
    bool operator()(const std::pair<int64_t, float>& a, const std::pair<int64_t, float>& b) const {
        return a.second > b.second;  // Compare based on the distance in descending order
        }
    };

    std::vector<std::set<std::pair<int64_t, float>, ComparePairs>> edt;
    std::vector<std::vector<int64_t>> minSupWMInds;

    edt.resize(MT::MAXNUMBEROFTHREADS());
    minSupWMInds.resize(MT::MAXNUMBEROFTHREADS());


    float smallestSphere = wmDepth.pixDims[0]*std::sqrt(3)+EPS4;

    auto getEDTpair = [&](NIBR::MT::TASK task)->void {

        int64_t ind = task.no;

        if (wmMask.data[ind]==INSIDE) { // If inside

            float p[3];

            wmDepth.to_xyz(ind,p);

            float radius = std::sqrt(surf.squaredDistToPoint(p));

            edt[task.threadId].insert(std::make_pair(ind, radius));
            
            if (radius > smallestSphere) {
                edt[task.threadId].insert(std::make_pair(ind, radius));
            } else {
                // These are already the smallest values. Because in image space they will correspond to a single voxel, they can't overwrite any other voxel.
                wmDepth.data[ind] = radius;
            }

            if (radius <= minThicknessOfSuperficialWM) {
                minSupWMInds[task.threadId].push_back(ind);
            }
            

        }
            

    };    
    
    NIBR::MT::MTRUN(wmDepth.numel, "Indexing EDT", getEDTpair);
    

    
    // Merge threads
    for (int t = 1; t < MT::MAXNUMBEROFTHREADS(); t++) {
        std::set_union(
            edt[0].begin(), edt[0].end(),
            edt[t].begin(), edt[t].end(),
            std::inserter(edt[0], edt[0].end()),
            ComparePairs()
        );
        edt[t].clear();

        minSupWMInds[0].insert(
            minSupWMInds[0].end(),
            std::make_move_iterator(minSupWMInds[t].begin()),
            std::make_move_iterator(minSupWMInds[t].end()));
    }

    // We will save 1/8 of the biggest sphere
    std::vector<std::vector<std::vector<float>>> s;

    // First process the biggest sphere
    auto it = edt[0].begin();
    if (!edt[0].empty()) {

        float maxRadius = it->second;
        int R           = std::ceil(maxRadius/wmDepth.pixDims[0]);    // We assume isotropic voxel size

        s.resize(R+1);
        for (int i = 0; i <= R; i++) {
            s[i].resize(R+1);
            for (int j = 0; j <= R; j++) {
                s[i][j].resize(R+1);
                for (int k = 0; k <= R; k++) {
                    float r = std::sqrt(float(i*i + j*j + k*k))*wmDepth.pixDims[0];
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

    std::mutex* ijkAvailability = new std::mutex[wmDepth.voxCnt];

    auto getMaximalSpheres = [&](NIBR::MT::TASK task)->void {

        float radius = vedt[task.no].second;
        int R        = std::ceil(radius/wmDepth.pixDims[0]);    // We assume isotropic voxel size

        int64_t ci,cj,ck; // i,j,k of sphere center
        wmDepth.ind2sub(vedt[task.no].first,ci,cj,ck);

        for (int i = ci-R; i <= ci+R; i++)
            for (int j = cj-R; j <= cj+R; j++)
                for (int k = ck-R; k <= ck+R; k++) 
            {
                if (wmDepth.isInside(i, j, k)) {

                    int64_t index = wmDepth.sub2ind(i,j,k);
                    float r = s[std::abs(ci - i)][std::abs(cj - j)][std::abs(ck - k)];

                    if (r<=radius) {
                        {
                            std::lock_guard<std::mutex> lock(ijkAvailability[index]);

                            if (r > wmDepth.data[index]) {
                                wmDepth.data[index] = radius;
                            }
                            
                        }
                    }
                }
            }        

    };

    NIBR::MT::MTRUN(vedt.size(), "Computing white matter depth", getMaximalSpheres);

    delete[] ijkAvailability;

    wmLabel.createFromTemplate(wmDepth,true);
    mapSurface2Image(&surf,&wmLabel,0,NULL,NULL,ONLY_BOUNDARY); // Marks all voxels that includes triangle(s) regardless of the voxel center being inside the surface or not

    std::set<int64_t> boundary; // if inside the mask, these will be labeled as WM_GRAY_CROWN, WM_SULCAL_WALL or WM_SULCAL_FUNDUS

    for (int64_t n = 0; n < wmLabel.voxCnt; n++) {

        if (wmLabel.data[n]==1) {

            int64_t i,j,k;
            wmLabel.ind2sub(n,i,j,k);
            
            
            boundary.insert(n);
        
            if (wmMask.data[n] == 0) {
                if (wmMask(i-1,j,  k  )) boundary.insert(wmLabel.sub2ind(i-1,j,  k  ));
                if (wmMask(i+1,j,  k  )) boundary.insert(wmLabel.sub2ind(i+1,j,  k  ));
                if (wmMask(i,  j-1,k  )) boundary.insert(wmLabel.sub2ind(i,  j-1,k  ));
                if (wmMask(i,  j+1,k  )) boundary.insert(wmLabel.sub2ind(i,  j+1,k  ));
                if (wmMask(i,  j,  k-1)) boundary.insert(wmLabel.sub2ind(i,  j,  k-1));
                if (wmMask(i,  j,  k+1)) boundary.insert(wmLabel.sub2ind(i,  j,  k+1));
            }
            

        }

    }

    // First mark everything inside as WM_DEEP
    for (int64_t n = 0; n < wmLabel.voxCnt; n++) {
        wmLabel.data[n] = wmMask.data[n]; // Note that WM_DEEP = 1 so this line works
    }

    // Then mark WM_SULCAL_FUNDUS and WM_GYRAL_CROWN
    std::set<int64_t> sulciInds;
    for (auto n : boundary) {
        if (wmMask.data[n]) {
            if (wmDepth.data[n] > sulcalDepth) {
                wmLabel.data[n] = WM_SULCAL_FUNDUS;
                sulciInds.insert(n);
            } else {
                wmLabel.data[n] = WM_GYRAL_CROWN;
            }
        }
    }

    // Then region grow WM_SULCAL_FUNDUS to mark WM_SULCAL_WALL
    std::set<int64_t> walls;

    for (auto n : sulciInds) {

        int64_t i,j,k,ind;
        wmLabel.ind2sub(n,i,j,k);

        ind = wmLabel.sub2ind(i-1,j,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){walls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
        ind = wmLabel.sub2ind(i+1,j,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){walls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
        ind = wmLabel.sub2ind(i,j-1,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){walls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
        ind = wmLabel.sub2ind(i,j+1,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){walls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
        ind = wmLabel.sub2ind(i,j,k-1); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){walls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
        ind = wmLabel.sub2ind(i,j,k+1); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){walls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}

    }

    auto thereAreWalls = !walls.empty();   
    
    std::set<int64_t> newWalls;

    while (thereAreWalls) {

        for (auto n : walls) {

            int64_t i,j,k,ind;
            wmLabel.ind2sub(n,i,j,k);

            ind = wmLabel.sub2ind(i-1,j,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){newWalls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
            ind = wmLabel.sub2ind(i+1,j,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){newWalls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
            ind = wmLabel.sub2ind(i,j-1,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){newWalls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
            ind = wmLabel.sub2ind(i,j+1,k); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){newWalls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
            ind = wmLabel.sub2ind(i,j,k-1); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){newWalls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}
            ind = wmLabel.sub2ind(i,j,k+1); if ((wmLabel.data[ind] == WM_GYRAL_CROWN) && (wmDepth.data[ind] > wallDepth) ){newWalls.insert(ind); wmLabel.data[ind] = WM_SULCAL_WALL;}

        }

        walls         = newWalls;
        newWalls.clear();
        thereAreWalls = !walls.empty();   
        
    }
    // WM_SULCAL_FUNDUS, WM_SULCAL_WALL and WM_GYRAL_CROWN are now marked, and won't change

    // Next label WM_SUPERFICIAL 

    // First mark all those in WM_DEEP that are less deep than maxThicknessOfSuperficialWM as WM_SUPERFICIAL
    for (int64_t n = 0; n < wmLabel.voxCnt; n++) {
        if ( (wmLabel.data[n] == WM_DEEP) && (wmDepth.data[n] <= maxThicknessOfSuperficialWM) ) {
            wmLabel.data[n] = WM_SUPERFICIAL;
        }
    }

    // Then mark all those in WM_DEEP that are within minThicknessOfSuperficialWM
    for (auto n : minSupWMInds[0]) {
        if (wmLabel.data[n] == WM_DEEP) {
            wmLabel.data[n] = WM_SUPERFICIAL;
        }
    }

    // Lastly, make a thin line of WM_DEEP so all WM_SUPERFICIAL touches WM_DEEP
    for (int64_t n = 0; n < wmLabel.voxCnt; n++) {
        if (wmLabel.data[n] == WM_SUPERFICIAL) {
            
            int64_t i,j,k,ind;
            wmLabel.ind2sub(n,i,j,k);

            ind = wmLabel.sub2ind(i-1,j,k); if (wmLabel.data[ind] == 0) {wmLabel.data[n] = WM_DEEP; continue;}
            ind = wmLabel.sub2ind(i+1,j,k); if (wmLabel.data[ind] == 0) {wmLabel.data[n] = WM_DEEP; continue;}
            ind = wmLabel.sub2ind(i,j-1,k); if (wmLabel.data[ind] == 0) {wmLabel.data[n] = WM_DEEP; continue;} 
            ind = wmLabel.sub2ind(i,j+1,k); if (wmLabel.data[ind] == 0) {wmLabel.data[n] = WM_DEEP; continue;} 
            ind = wmLabel.sub2ind(i,j,k-1); if (wmLabel.data[ind] == 0) {wmLabel.data[n] = WM_DEEP; continue;} 
            ind = wmLabel.sub2ind(i,j,k+1); if (wmLabel.data[ind] == 0) {wmLabel.data[n] = WM_DEEP; continue;} 

        }
    }

    // Finally mark WM_BOUNDARY, which are the voxels that have triangles but the voxel centers are not within WM
    for (auto n : boundary) {
        if (wmLabel.data[n] == 0) {
            wmLabel.data[n] = WM_BOUNDARY;
        }
    }

    disp(MSG_DEBUG,"Done wmSurf2wmLabel()");

}