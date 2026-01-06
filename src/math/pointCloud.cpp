#include "base/multithreader.h"
#include "pointCloud.h"
#include "base/verbose.h"
#include "nanoflann/nanoflann.hpp"
#include <cstddef>
#include <limits>
#include <algorithm>
#include <tuple>
#include <algorithm>
#include <random>

using namespace NIBR;

#define MAXMINIBATCHSIZE 65536

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>,PointCloud,-1> KDTree;

struct NIBR::Clusterer::Impl {
    
    std::vector<Eigen::VectorXf> centers;
    std::vector<Eigen::VectorXd> sums;
    std::vector<double>          counts;
    int                          dim;
    int                          maxLeafSize;
    float                        distanceThreshold;
    float                        squaredDistanceThreshold;
    Clusterer::Mode              mode;
    bool                         shuffle;
    std::mt19937                 rng;
    size_t                       totalProcessed;
    
    // KDTrees
    PointCloud cloud,   newCloud;
    KDTree    *kdtree, *newkdtree;
    

    // Constructors
    Impl(int d, float distThresh, bool shuf) : dim(d), distanceThreshold(distThresh), squaredDistanceThreshold(distThresh*distThresh), mode(Clusterer::Mode::DISCOVERY), shuffle(shuf), kdtree(nullptr) {
        std::random_device rd;
        rng.seed(rd());
        maxLeafSize = ((256/dim) > 10) ? 10 : (256/dim);

        cloud.points = &centers;
        kdtree = new KDTree(dim, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafSize));
        kdtree->buildIndex();
        
        newCloud.points = &newCenters;
        newkdtree = new KDTree(dim, newCloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafSize));
        newkdtree->buildIndex();

        totalProcessed = 0;
    }

    // Destructors
    ~Impl() {
        if (kdtree)    delete kdtree;
        if (newkdtree) delete newkdtree;
    }

    // Reusable buffers to avoid reallocations
    std::vector<int>             assignmentsBuffer;
    std::vector<size_t>          unassignedIndicesBuffer;
    std::vector<Eigen::VectorXf> newCenters;
    std::vector<Eigen::VectorXd> newSums;
    std::vector<double>          newCounts;
    std::mutex                   newCenterMutex[MAXMINIBATCHSIZE];

};

NIBR::Clusterer::Clusterer(int dim,float distanceThreshold, bool shuffle) {
    pImpl = new Impl(dim, distanceThreshold, shuffle);
}

NIBR::Clusterer::~Clusterer() {
    delete pImpl;
}

void NIBR::Clusterer::setMode(Mode mode) {
    pImpl->mode = mode;
}

void NIBR::Clusterer::setCenters(const std::vector<Eigen::VectorXf>& centers) {
    pImpl->centers = centers;
    if (!centers.empty()) {
        pImpl->dim = centers[0].size();
        pImpl->sums.assign(centers.size(), Eigen::VectorXd::Zero(pImpl->dim));
        pImpl->counts.assign(centers.size(), 0.0);
    } else {
        pImpl->sums.clear();
        pImpl->counts.clear();
    }
    pImpl->cloud.points = &pImpl->centers;
    pImpl->kdtree->buildIndex();
}

std::vector<Eigen::VectorXf> NIBR::Clusterer::getClusterCenters() const {
    return pImpl->centers;
}

std::vector<double> NIBR::Clusterer::getClusterCounts() const {
    return pImpl->counts;
}

int NIBR::Clusterer::getClusterCenterCount() const {
    return pImpl->centers.size();
}

std::tuple<int, float> findClosestCenter(const KDTree* kdtree, const Eigen::VectorXf& point) {
    if (kdtree->size_ == 0) return {-1, std::numeric_limits<float>::quiet_NaN()};
    
    size_t closestCenterIndex;
    float  squaredDist;
    nanoflann::KNNResultSet<float> resultSet(1);
    resultSet.init(&closestCenterIndex, &squaredDist);
    kdtree->findNeighbors(resultSet, point.data(), nanoflann::SearchParameters());
    
    return std::make_tuple((int)closestCenterIndex, squaredDist);
}

void NIBR::Clusterer::updateCenters() {
    
    // Update centers
    for (size_t i = 0; i < pImpl->centers.size(); ++i) {
        if (pImpl->counts[i] > 0) {
            pImpl->centers[i] = (pImpl->sums[i] / pImpl->counts[i]).cast<float>();
        }
    }

    // Reset sums and counts
    for (auto& s : pImpl->sums) s.setZero();
    std::fill(pImpl->counts.begin(), pImpl->counts.end(), 0.0f);

    // Rebuild tree with new centers
    pImpl->cloud.points = &pImpl->centers;
    pImpl->kdtree->buildIndex();

    pImpl->totalProcessed = 0;
}

void NIBR::Clusterer::process(const std::vector<Eigen::VectorXf>& points) {
    
    if (points.empty()) return;
    pImpl->totalProcessed += points.size();
    if (pImpl->mode == Mode::DISCOVERY) {
        disp(MSG_INFO,"Current cluster count: %d. Adding %d points. Total point count: %d",pImpl->centers.size(),points.size(),pImpl->totalProcessed);
    }

    // -----------------------------------------------------
    // PHASE 1: Assign existing cluster centers in parallel
    // -----------------------------------------------------

    size_t N = points.size();
    if (pImpl->assignmentsBuffer.size() < N) pImpl->assignmentsBuffer.resize(N);
    int* assignments = pImpl->assignmentsBuffer.data();
    
    auto searchTask = [&](NIBR::MT::TASK task) -> void {

        auto [closestIdx, sqDist] = findClosestCenter(pImpl->kdtree, points[task.no]);

        if (sqDist < pImpl->squaredDistanceThreshold) {
            assignments[task.no] = closestIdx;
        } else {
            assignments[task.no] = -1;
        }

    };
    
    if (!pImpl->centers.empty()) {
        NIBR::MT::MTRUN(N, searchTask);
    } else {
        for (size_t i = 0; i < N; ++i) {
            assignments[i] = -1;
        }
    }
    

    // -------------------------------
    // PHASE 2: Accumulate assignments
    // -------------------------------
    
    pImpl->unassignedIndicesBuffer.clear();
    
    if (pImpl->mode == Mode::DISCOVERY || pImpl->mode == Mode::REFINEMENT) {
        if (pImpl->mode == Mode::REFINEMENT) {
            disp(MSG_INFO,"Assigning %d points to existing clusters.",N);
        }
        for (size_t i = 0; i < N; ++i) {
            int idx = assignments[i];
            
            if (idx != -1) {
                // Assign to existing cluster: Update stats
                pImpl->sums[idx]   += points[i].cast<double>();
                pImpl->counts[idx] += 1.0;
            } else {
                // Save for batch processing
                if (pImpl->mode != Mode::REFINEMENT) {
                    pImpl->unassignedIndicesBuffer.push_back(i); // Store absolute index
                }
            }
        }
        if (pImpl->mode == Mode::REFINEMENT) {
            disp(MSG_INFO,"Done");
        }
    } else if (pImpl->mode == Mode::COVERAGE) {
        // In coverage, we ignore updates to existing, just look for unassigned
        for (size_t i = 0; i < N; ++i) {
            if (assignments[i] == -1) {
                pImpl->unassignedIndicesBuffer.push_back(i); // Store absolute index
            }
        }
    }

    // If REFINEMENT mode, we are done
    if (pImpl->mode == Mode::REFINEMENT) return;

    disp(MSG_INFO,"%d points are far from any cluster. Computing new clusters...",pImpl->unassignedIndicesBuffer.size());

    // If there are no unassigned, we are done
    if (pImpl->unassignedIndicesBuffer.empty()) return;

    
    // ----------------------------
    // PHASE 3: Create new clusters
    // ----------------------------

    // We will process unassigned points in batches
    int unaCnt = pImpl->unassignedIndicesBuffer.size();

    if (pImpl->shuffle) {
        std::shuffle(pImpl->unassignedIndicesBuffer.begin(), pImpl->unassignedIndicesBuffer.end(), pImpl->rng);
    }

    // Reuse buffers for new centers
    pImpl->newCenters.clear();
    pImpl->newSums.clear();
    pImpl->newCounts.clear();

    int miniBatchSize;
    if      (unaCnt <= 100000)  miniBatchSize = 1000;
    else if (unaCnt <= 1000000) miniBatchSize = 10000;
    else                        miniBatchSize = MAXMINIBATCHSIZE;

    for (int start = 0; start < unaCnt; start += miniBatchSize) {

        int end              = std::min(start + miniBatchSize, unaCnt);
        int currentBatchSize = end - start;

        // ------------------------------------------------
        // Step A: Parallel assignment to known NEW centers
        // ------------------------------------------------
        // We use a simple byte array: 0 = Assigned, 1 = Miss
        std::vector<uint8_t> localMisses(currentBatchSize, 0);

        auto attemptAssign = [&](NIBR::MT::TASK task) -> void {

            if (pImpl->newCenters.empty()) {
                localMisses[task.no] = 1;
                return;
            }

            size_t realIdx = pImpl->unassignedIndicesBuffer[start + task.no];
            const auto& p  = points[realIdx];

            auto [closestIdx, sqDist] = findClosestCenter(pImpl->newkdtree, p);

            if (sqDist < pImpl->squaredDistanceThreshold) {
                // Found a match. Update stats.
                size_t mutexIdx = closestIdx % MAXMINIBATCHSIZE;
                std::lock_guard<std::mutex> lock(pImpl->newCenterMutex[mutexIdx]);

                pImpl->newSums[closestIdx]    += p.cast<double>();
                pImpl->newCounts[closestIdx]  += 1.0;
            } else {
                localMisses[task.no] = 1;
            }
        };

        NIBR::MT::MTRUN(currentBatchSize, attemptAssign);

        // Update centers based on accumulated stats (Batch Update)
        for (size_t k = 0; k < pImpl->newCenters.size(); ++k) {
             if (pImpl->newCounts[k] > 0) {
                 pImpl->newCenters[k] = (pImpl->newSums[k] / pImpl->newCounts[k]).cast<float>();
             }
        }

        // -----------------------
        // Step B: Serial creation
        // -----------------------
        // Now we take the points that failed the tree search and 
        // cluster them against EACH OTHER and add to the global list.
        // This MUST be serial to prevent duplicates.
        
        size_t centersAddedBeforeThisStep = pImpl->newCenters.size();

        for (int i = 0; i < currentBatchSize; ++i) {
            if (localMisses[i] == 1) {
                size_t realIdx = pImpl->unassignedIndicesBuffer[start + i];
                const auto& p  = points[realIdx];

                int bestCenter = -1;
                
                // Greedy linear scan of RECENTLY added centers
                for (size_t k = centersAddedBeforeThisStep; k < pImpl->newCenters.size(); ++k) {
                    float dist = (p - pImpl->newCenters[k]).squaredNorm();
                    if (dist < pImpl->squaredDistanceThreshold) {
                        bestCenter = k;
                        break;
                    }
                }

                if (bestCenter != -1) {
                    // Update recent center
                    double n = pImpl->newCounts[bestCenter] + 1.0;
                    double alpha = 1.0 / n;
                    pImpl->newCenters[bestCenter] += (float)alpha * (p - pImpl->newCenters[bestCenter]);
                    pImpl->newSums[bestCenter]    += p.cast<double>();
                    pImpl->newCounts[bestCenter]   = n;
                } else {
                    // Create BRAND NEW center
                    pImpl->newCenters.push_back(p);
                    pImpl->newSums.push_back(p.cast<double>());
                    pImpl->newCounts.push_back(1.0);
                }
            }
        }

        // --------------------
        // Step C: Rebuild Tree
        // --------------------
        if (pImpl->newCenters.size() > centersAddedBeforeThisStep) {
            pImpl->newCloud.points = &pImpl->newCenters;
            pImpl->newkdtree->buildIndex();
        }
    }

    // -------------------------------------
    // PHASE 4: Merge New Clusters to Global
    // -------------------------------------
    if (!pImpl->newCenters.empty()) {

        disp(MSG_INFO,"%d new cluster centers added.",pImpl->newCenters.size());

        pImpl->centers.insert(pImpl->centers.end(), pImpl->newCenters.begin(), pImpl->newCenters.end());
        pImpl->sums.insert(pImpl->sums.end(), pImpl->newSums.begin(), pImpl->newSums.end());
        pImpl->counts.insert(pImpl->counts.end(), pImpl->newCounts.begin(), pImpl->newCounts.end());

        pImpl->cloud.points = &pImpl->centers;
        pImpl->kdtree->buildIndex();
    }
    
}
