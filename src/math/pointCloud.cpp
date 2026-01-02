#include "base/config.h"
#include "base/verbose.h"
#include "base/multithreader.h"
#include "base/stringOperations.h"
#include "pointCloud.h"
#include "nanoflann/nanoflann.hpp"
#include <numeric>
#include <random>

using namespace NIBR;

std::vector<Eigen::VectorXf> NIBR::findClusterCenters(
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

std::pair<int, float> NIBR::getNearestCluster(
    const std::vector<Eigen::VectorXf>& centers, 
    const Eigen::VectorXf& point
)
{
    if (centers.empty()) return {-1, FLT_MAX};

    int bestIdx = -1;
    float bestDist = FLT_MAX;

    // For small number of centers, linear search is faster than building a tree
    // But since this function might be called in a loop, we assume the user handles efficiency.
    // Here we do a simple linear search as a fallback/utility.
    // If efficiency is needed for many points, use ClusterRefiner or findClusterCenters.
    for (size_t i = 0; i < centers.size(); ++i) {
        float dist = (centers[i] - point).squaredNorm();
        if (dist < bestDist) {
            bestDist = dist;
            bestIdx = i;
        }
    }

    return {bestIdx, bestDist};
}

struct NIBR::ClusterRefiner::Impl {
    std::vector<Eigen::VectorXf> centers;
    std::vector<Eigen::VectorXf> accumulatedCenters;
    std::vector<float> counts;
    int dim;
    
    // KDTree
    PointCloud cloud;
    typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>,PointCloud,-1> KDTree;
    KDTree* kdtree;

    Impl(const std::vector<Eigen::VectorXf>& c) : centers(c) {
        if (centers.empty()) {
            dim = 0;
            kdtree = nullptr;
            return;
        }
        dim = c[0].size();
        accumulatedCenters.resize(centers.size());
        for(auto& ac : accumulatedCenters) ac = Eigen::VectorXf::Zero(dim);
        counts.assign(centers.size(), 0.0f);

        cloud.points = centers;
        int maxLeafSize = ((256/dim) > 10) ? 10 : (256/dim);
        kdtree = new KDTree(dim, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafSize));
        kdtree->buildIndex();
    }

    ~Impl() {
        if (kdtree) delete kdtree;
    }
};

NIBR::ClusterRefiner::ClusterRefiner(const std::vector<Eigen::VectorXf>& centers) {
    pImpl = new Impl(centers);
}

NIBR::ClusterRefiner::~ClusterRefiner() {
    delete pImpl;
}

void NIBR::ClusterRefiner::update(const std::vector<Eigen::VectorXf>& points) {
    if (points.empty() || !pImpl->kdtree) return;

    std::vector<int> assignments(points.size());

    // Parallel assignment
    auto assignPoints = [&](NIBR::MT::TASK task) -> void {
        size_t closestCenterIndex;
        float  squaredDist;
        nanoflann::KNNResultSet<float> resultSet(1);
        resultSet.init(&closestCenterIndex, &squaredDist);
        pImpl->kdtree->findNeighbors(resultSet, points[task.no].data(), nanoflann::SearchParameters());
        assignments[task.no] = closestCenterIndex;
    };

    NIBR::MT::MTRUN(points.size(), "Refining clusters", assignPoints);

    // Accumulate (Serial for safety, but fast since it's just addition)
    // Could be parallelized with thread-local storage if needed, but this is likely fine.
    for (size_t i = 0; i < points.size(); ++i) {
        int idx = assignments[i];
        pImpl->accumulatedCenters[idx] += points[i];
        pImpl->counts[idx] += 1.0f;
    }
}

std::vector<Eigen::VectorXf> NIBR::ClusterRefiner::getRefinedCenters() {
    std::vector<Eigen::VectorXf> refined = pImpl->centers;
    for (size_t i = 0; i < refined.size(); ++i) {
        if (pImpl->counts[i] > 0) {
            refined[i] = pImpl->accumulatedCenters[i] / pImpl->counts[i];
        }
    }
    return refined;
}

std::vector<Eigen::VectorXf> NIBR::refineClusterCenters(
    const std::vector<Eigen::VectorXf>& centers,
    const std::vector<Eigen::VectorXf>& points
)
{
    ClusterRefiner refiner(centers);
    refiner.update(points);
    return refiner.getRefinedCenters();
}