#include "base/multithreader.h"
#include "pointCloud.h"
#include "nanoflann/nanoflann.hpp"
#include <limits>
#include <algorithm>
#include <utility>
#include <random>

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
    bool                         shuffle;

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
    std::mt19937 rng;

    Impl(float distThresh, int bSize, bool shuf) : distanceThreshold(distThresh), squaredDistanceThreshold(distThresh*distThresh), mode(Clusterer::Mode::DISCOVERY), dim(0), batchSize(bSize), shuffle(shuf), kdtree(nullptr) {
        std::random_device rd;
        rng.seed(rd());
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

NIBR::Clusterer::Clusterer(float distanceThreshold, int batchSize, bool shuffle) {
    pImpl = new Impl(distanceThreshold, batchSize, shuffle);
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

int NIBR::Clusterer::getClusterCenterCount() const {
    return pImpl->centers.size();
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

    size_t N                = points.size();
    size_t batchSize        = pImpl->batchSize;
    size_t initialCenters   = pImpl->centers.size();

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

    if (pImpl->shuffle) {
        std::shuffle(pImpl->unassignedIndicesBuffer.begin(), pImpl->unassignedIndicesBuffer.end(), pImpl->rng);
    }

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



/*
std::vector<Eigen::VectorXf> findClusterCenters(
    const std::vector<Eigen::VectorXf>& prevClusterCenters, 
    const std::vector<Eigen::VectorXf>& newPoints,
    float   distanceThreshold,
    int     miniBatchSize,
    bool    selectRandomMiniBatch, 
    bool    shuffleMiniBatch, 
    int     maxIteration,
    bool    updateCenters
)
{
    int dim         = newPoints[0].size(); // Number of dimensions
    int batchSize   = newPoints.size();    // Number of points
    int maxLeafSize = ((256/dim) > 10) ? 10 : (256/dim);
    float squaredDistanceThreshold = distanceThreshold * distanceThreshold;
    
    // Build KD-Tree for existing cluster centers
    PointCloud cloud;
    cloud.points = prevClusterCenters;
    typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>,PointCloud,-1> KDTree;
    KDTree kdtree(dim, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafSize) ) ;
    kdtree.buildIndex();

    // Setup mini-batch generator
    std::vector<int> indices(batchSize);
    std::iota(indices.begin(), indices.end(), 0);
    auto rand = std::mt19937(std::random_device{}());
    if (selectRandomMiniBatch) std::shuffle(indices.begin(), indices.end(), rand);

    int nextMiniBatchStart = 0;
    auto getMiniBatch = [&]() -> std::vector<Eigen::VectorXf> {
        std::vector<Eigen::VectorXf> miniBatch;
        miniBatch.reserve(miniBatchSize);
        
        for (int i = 0; i < miniBatchSize; i++) {
            if (nextMiniBatchStart + i >= batchSize) break;
            int idx = indices[nextMiniBatchStart + i];
            miniBatch.push_back(newPoints[idx]);
        }
        
        nextMiniBatchStart += miniBatch.size();

        if (shuffleMiniBatch) std::shuffle(miniBatch.begin(), miniBatch.end(), rand);

        return miniBatch;
    };

    // Compute maxIteration
    if (selectRandomMiniBatch) {
        if (maxIteration < 1) maxIteration = 1;
    } else {
        if (maxIteration < 1) maxIteration = std::ceil(float(batchSize) / float(miniBatchSize));    
    }

    // Do the clustering
    std::vector<Eigen::VectorXf> clusterCenters = prevClusterCenters;
    std::vector<float>           clusterSize(clusterCenters.size(), 1.0f);

    // Clustering
    disp(MSG_INFO,"Clustering...");

    for (int iter = 0; iter < maxIteration; iter++) {
        
        // Get a mini-batch
        std::vector<Eigen::VectorXf> miniBatch = getMiniBatch();
        disp(MSG_DETAIL,"Current mini-batch size: %d", miniBatch.size());

        if (miniBatch.empty()) break;

        // Perform clustering operations on the mini-batch
        // This will be done in two steps
        // Step 1. If a point is far from all existing cluster centers, keep it as "unassigned".
        // Step 2. Within in mini-batch, locally cluster all the "unassigned" points, and form localClusterCenters.
        // Step 3. Append localClusterCenters to clusterCenters

        // For efficiency split the unassigned points into smaller batches of:
        // 10000 points in the first iteration when there will be many new clusters,
        // and 1000 points in the other iteration when there will be many existing clusters, and less new clusters.

        std::vector<Eigen::VectorXf> localClusterCenters;       // New clusters found in this mini-batch.
        std::vector<float>           localClusterSize;
        
        // Build an empty KD-Tree for local clusters
        PointCloud localCloud;
        typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>,PointCloud,-1> KDTree;
        KDTree localKdtree(dim, localCloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafSize));
        localKdtree.buildIndex();

        std::vector<std::atomic<bool>> unassigned(miniBatch.size());    // True if a point in the mini-batch did not belong to any existing cluster
        for (std::size_t i = 0; i < miniBatch.size(); ++i) unassigned[i] = false;

        std::vector<size_t> unaInd;                              // Indices of the points, which did not belong to any existing cluster
        std::vector<Eigen::VectorXf> unassignedClusterCenters;   // Clusters of the unassigned points that were not clusered within a mini-batch
        std::vector<float>           unassignedClusterSize;
        std::vector<int>             assignments(miniBatch.size(), -1);
        size_t taskOffset = 0;

        auto addToGlobalCluster = [&](NIBR::MT::TASK task) -> void {

            if (!clusterCenters.empty()) {
                size_t closestCenterIndex;
                float  squaredDistToClosestClusterCenter;

                nanoflann::KNNResultSet<float> resultSet(1);
                resultSet.init(&closestCenterIndex, &squaredDistToClosestClusterCenter);
                kdtree.findNeighbors(resultSet, miniBatch[task.no].data(), nanoflann::SearchParameters());
                if (squaredDistToClosestClusterCenter < squaredDistanceThreshold) {
                    assignments[task.no] = closestCenterIndex;
                    return;
                }
            }

            unassigned[task.no].store(true);

        };

        std::mutex mx;

        auto addToLocalCluster = [&](NIBR::MT::TASK task) -> void {

            size_t rInd = unaInd[task.no+taskOffset];

            if (!localClusterCenters.empty()) {
                size_t closestCenterIndex;
                float  squaredDistToClosestClusterCenter;

                nanoflann::KNNResultSet<float> resultSet(1);
                resultSet.init(&closestCenterIndex, &squaredDistToClosestClusterCenter);
                localKdtree.findNeighbors(resultSet, miniBatch[rInd].data(), nanoflann::SearchParameters());
                if (squaredDistToClosestClusterCenter < squaredDistanceThreshold) return;
            }

            {
                mx.lock();

                for (size_t ind = 0; ind < unassignedClusterCenters.size(); ind++) {

                    float sum = 0;

                    for (int i = 0; i < dim; i++) {
                        float d1 = (miniBatch[rInd][i] - unassignedClusterCenters[ind][i]);
                        sum     += d1 * d1;
                    }

                    if (sum < squaredDistanceThreshold) {
                        
                        // Update the cluster center
                        float n = unassignedClusterSize[ind] + 1.0f;
                        float alpha = 1.0f / n;
                        for (int i = 0; i < dim; i++) {
                            unassignedClusterCenters[ind][i] += alpha * (miniBatch[rInd][i] - unassignedClusterCenters[ind][i]);
                        }
                        unassignedClusterSize[ind] = n;

                        mx.unlock();
                        return;
                    }
                }

                unassignedClusterCenters.push_back(miniBatch[rInd]);
                unassignedClusterSize.push_back(1.0f);
                mx.unlock();
                return;                
            }   

        };


        // Add the newly found localClusterCenters in the global clusterCenters
        auto doUnassigned = [&]() -> int {

            unaInd.clear();

            for (size_t r = 0; r < unassigned.size(); r++) {
                if (unassigned[r]) unaInd.push_back(r);
            }
            int unaCnt     = unaInd.size();

            int splitSize;
            if      (unaCnt <= 100000)  splitSize = 1000;
            else if (unaCnt <= 1000000) splitSize = 10000;
            else                        splitSize = 100000;

            taskOffset     = 0;
            int begInd     = 0;
            int endInd     = 0;
            int splitCnt   = (unaCnt / splitSize < 1) ? unaCnt : unaCnt/splitSize;
            int curSplitNo = 0;
            std::string preamble = "\033[1;32mNIBRARY::INFO: \033[0;32m";
            if (NIBR::VERBOSE()>=VERBOSE_INFO) {std::cout << preamble << "Clustering unassigned points " << ": 0%" << "\033[0m" << '\r' << std::flush;}
            float progressScaler = 100.0f/float(splitCnt);                
            while (endInd != unaCnt) {
                begInd = endInd;
                endInd = ((begInd + splitSize) <= unaCnt) ? (begInd + splitSize) : unaCnt;
                int curSplitSize = endInd - begInd;
                NIBR::MT::MTRUN(curSplitSize, addToLocalCluster);
                localClusterCenters.insert(localClusterCenters.end(), unassignedClusterCenters.begin(), unassignedClusterCenters.end());
                localClusterSize.insert(localClusterSize.end(), unassignedClusterSize.begin(), unassignedClusterSize.end());
                localCloud.points = localClusterCenters;
                localKdtree.buildIndex();
                unassignedClusterCenters.clear();
                unassignedClusterSize.clear();
                taskOffset += curSplitSize;
                if (NIBR::VERBOSE()>=VERBOSE_INFO) {std::cout << "\r\033[K" << std::flush;}
                if (NIBR::VERBOSE()>=VERBOSE_INFO) {std::cout << preamble << "Clustering unassigned points: " << std::fixed << std::setprecision(2) << (++curSplitNo)*progressScaler << "%" << "\033[0m" << std::flush;}
            }
            if (NIBR::VERBOSE()>=VERBOSE_INFO) {std::cout << "\r\033[K" << preamble << "Clustering unassigned points: 100%" << std::endl;}

            return unaCnt;

        }; 

        
        // Assign points into existing clusters, and find "unassigned" points, which were not assinged to any cluster
        NIBR::MT::MTRUN( batchSize, "Assigning clusters " + to_string_with_precision(iter+1,0) + " / " + to_string_with_precision(maxIteration,0), addToGlobalCluster);

        // Update global cluster centers
        if (updateCenters) {
            for (size_t i = 0; i < miniBatch.size(); i++) {
                int clusterIdx = assignments[i];
                if (clusterIdx != -1) {
                    float n = clusterSize[clusterIdx] + 1.0f;
                    float alpha = 1.0f / n;
                    for (int d = 0; d < dim; d++) {
                        clusterCenters[clusterIdx][d] += alpha * (miniBatch[i][d] - clusterCenters[clusterIdx][d]);
                    }
                    clusterSize[clusterIdx] = n;
                }
            }
        }

        // Find localClusterCenters that are the cluster centers of the "unassigned" points, 
        int unassignedCnt = doUnassigned();

        // Append the localClusterCenters to global clusterCenters
        clusterCenters.insert(clusterCenters.end(), localClusterCenters.begin(), localClusterCenters.end());
        clusterSize.insert(clusterSize.end(), localClusterSize.begin(), localClusterSize.end());
        
        // Update the global KD-tree
        cloud.points = clusterCenters;
        kdtree.buildIndex();

        disp(MSG_INFO,"Unassigned points: %d - New clusters found: %d - Total clusters: %d", unassignedCnt, localClusterCenters.size(), clusterCenters.size());

    }
    disp(MSG_INFO,"Done");
    
    return clusterCenters;
}
*/