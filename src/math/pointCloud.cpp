#include "base/multithreader.h"
#include "pointCloud.h"
#include "nanoflann/nanoflann.hpp"
#include <limits>
#include <algorithm>
#include <utility>

using namespace NIBR;

struct NIBR::Clusterer::Impl {
    std::vector<Eigen::VectorXf> centers;
    std::vector<Eigen::VectorXd> sums;
    std::vector<double>          counts;
    float                        distanceThreshold;
    float                        squaredDistanceThreshold;
    Clusterer::Mode              mode;
    int                          dim;
    int                          batchSize;

    // Reusable buffers for processBatch to avoid reallocations
    std::vector<int>             assignmentsBuffer;
    std::vector<size_t>          unassignedIndicesBuffer;
    std::vector<Eigen::VectorXf> newCentersBuffer;
    std::vector<Eigen::VectorXd> newSumsBuffer;
    std::vector<double>          newCountsBuffer;

    // KDTree
    PointCloud cloud;
    typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>,PointCloud,-1> KDTree;
    KDTree* kdtree;

    Impl(float distThresh, int bSize) : distanceThreshold(distThresh), squaredDistanceThreshold(distThresh*distThresh), mode(Clusterer::Mode::DISCOVERY), dim(0), batchSize(bSize), kdtree(nullptr) {
    }

    ~Impl() {
        if (kdtree) delete kdtree;
    }

    void rebuildTree() {
        if (kdtree) delete kdtree;
        cloud.points = &centers;
        if (centers.empty()) {
            kdtree = nullptr;
            return;
        }
        dim = centers[0].size();
        int maxLeafSize = ((256/dim) > 10) ? 10 : (256/dim);
        kdtree = new KDTree(dim, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafSize));
        kdtree->buildIndex();
    }


};

NIBR::Clusterer::Clusterer(float distanceThreshold, int batchSize) {
    pImpl = new Impl(distanceThreshold, batchSize);
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
    pImpl->rebuildTree();
}

std::vector<Eigen::VectorXf> NIBR::Clusterer::getClusterCenters() const {
    return pImpl->centers;
}

std::vector<double> NIBR::Clusterer::getClusterCounts() const {
    return pImpl->counts;
}

std::pair<int, float> NIBR::Clusterer::assign(const Eigen::VectorXf& point) {
    if (!pImpl->kdtree) return {-1, std::numeric_limits<float>::quiet_NaN()};
    
    size_t closestCenterIndex;
    float  squaredDist;
    nanoflann::KNNResultSet<float> resultSet(1);
    resultSet.init(&closestCenterIndex, &squaredDist);
    pImpl->kdtree->findNeighbors(resultSet, point.data(), nanoflann::SearchParameters());
    
    return {(int)closestCenterIndex, squaredDist};
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
    pImpl->rebuildTree();
}

void NIBR::Clusterer::process(const std::vector<Eigen::VectorXf>& points) {
    if (points.empty()) return;

    if (pImpl->dim == 0) pImpl->dim = points[0].size();
    if (!pImpl->kdtree && !pImpl->centers.empty()) pImpl->rebuildTree();

    size_t N = points.size();
    size_t batchSize = pImpl->batchSize;
    size_t initialCenters = pImpl->centers.size();

    for (size_t start = 0; start < N; start += batchSize) {
        size_t end = std::min(start + batchSize, N);
        processBatch(points, start, end);
        
        if (pImpl->centers.size() > initialCenters) {
             pImpl->rebuildTree();
             initialCenters = pImpl->centers.size();
        }
    }
}

void NIBR::Clusterer::processBatch(const std::vector<Eigen::VectorXf>& points, size_t start, size_t end) {
    
    const size_t N = end - start;
    
    // Ensure buffer size
    if (pImpl->assignmentsBuffer.size() < N) pImpl->assignmentsBuffer.resize(N);
    int* assignments = pImpl->assignmentsBuffer.data();
    
    // ------------------------------------------
    // PHASE 1: Parallel Nearest Neighbor Search
    // ------------------------------------------
    
    auto searchTask = [&](NIBR::MT::TASK task) -> void {
        // Handle case with no centers/tree
        if (!pImpl->kdtree) {
            assignments[task.no] = -1;
            return;
        }

        size_t closestIdx;
        float  sqDist;
        nanoflann::KNNResultSet<float> resultSet(1);
        resultSet.init(&closestIdx, &sqDist);
        
        pImpl->kdtree->findNeighbors(resultSet, points[start + task.no].data(), nanoflann::SearchParameters());

        if (sqDist < pImpl->squaredDistanceThreshold) {
            assignments[task.no] = (int)closestIdx;
        } else {
            assignments[task.no] = -1;
        }
    };
    
    // Execute Parallel Search
    NIBR::MT::MTRUN(N, searchTask);

    // ----------------------------------
    // PHASE 2: Accumulation & Filtering
    // ----------------------------------
    
    pImpl->unassignedIndicesBuffer.clear();
    
    if (pImpl->mode == Mode::DISCOVERY || pImpl->mode == Mode::REFINEMENT) {
        for (size_t i = 0; i < N; ++i) {
            int idx = assignments[i];
            
            if (idx != -1) {
                // Assign to existing cluster: Update stats
                pImpl->sums[idx]   += points[start + i].cast<double>();
                pImpl->counts[idx] += 1.0;
            } else {
                // Save for Phase 3
                if (pImpl->mode != Mode::REFINEMENT) {
                    pImpl->unassignedIndicesBuffer.push_back(start + i); // Store absolute index
                }
            }
        }
    } else if (pImpl->mode == Mode::COVERAGE) {
        // In coverage, we ignore updates to existing, just look for unassigned
        for (size_t i = 0; i < N; ++i) {
            if (assignments[i] == -1) {
                pImpl->unassignedIndicesBuffer.push_back(start + i); // Store absolute index
            }
        }
    }

    // If Refinement mode, we are done (we don't create new clusters)
    if (pImpl->mode == Mode::REFINEMENT) return;
    if (pImpl->unassignedIndicesBuffer.empty()) return;

    // ------------------------------
    // PHASE 3: New Cluster Creation
    // ------------------------------
    
    // Reuse buffers for new centers
    pImpl->newCentersBuffer.clear();
    pImpl->newSumsBuffer.clear();
    pImpl->newCountsBuffer.clear();

    for (size_t idx : pImpl->unassignedIndicesBuffer) {
        const auto& p = points[idx];
        int bestCenter = -1;
        
        // Check against currently forming new centers
        size_t ncSize = pImpl->newCentersBuffer.size();
        for (size_t k = 0; k < ncSize; ++k) {
            float dist = (p - pImpl->newCentersBuffer[k]).squaredNorm();
            if (dist < pImpl->squaredDistanceThreshold) {
                bestCenter = k;
                break;
            }
        }

        if (bestCenter != -1) {
            // Update the local NEW center immediately
            double n = pImpl->newCountsBuffer[bestCenter] + 1.0;
            double alpha = 1.0 / n;
            pImpl->newCentersBuffer[bestCenter] += (float)alpha * (p - pImpl->newCentersBuffer[bestCenter]);
            pImpl->newSumsBuffer[bestCenter]    += p.cast<double>();
            pImpl->newCountsBuffer[bestCenter]   = n;
        } else {
            // Create brand new center
            pImpl->newCentersBuffer.push_back(p);
            pImpl->newSumsBuffer.push_back(p.cast<double>());
            pImpl->newCountsBuffer.push_back(1.0);
        }
    }

    // -------------------------------------
    // PHASE 4: Merge New Clusters to Global
    // -------------------------------------
    if (!pImpl->newCentersBuffer.empty()) {
        pImpl->centers.insert(pImpl->centers.end(), pImpl->newCentersBuffer.begin(), pImpl->newCentersBuffer.end());
        pImpl->sums.insert(pImpl->sums.end(), pImpl->newSumsBuffer.begin(), pImpl->newSumsBuffer.end());
        pImpl->counts.insert(pImpl->counts.end(), pImpl->newCountsBuffer.begin(), pImpl->newCountsBuffer.end());
    }
}
