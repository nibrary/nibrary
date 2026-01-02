#pragma once

#include <Eigen/Dense>
#include <vector>

namespace NIBR
{
    // Used for KD-tree for clustering and scoring 
    struct PointCloud {
        std::vector<Eigen::VectorXf> points;

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return points.size(); }

        // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class
        inline float kdtree_distance(const float* p1, const size_t idx_p2, size_t size) const {
            float dist = 0;
            const float* p2 = points[idx_p2].data();
            for (size_t i = 0; i < size; ++i) {
                float d = p1[i] - p2[i];
                dist += d * d;
            }
            return dist;
        }

        // Returns the dim'th component of the idx'th point in the class
        inline float kdtree_get_pt(const size_t idx, int dim) const { return points[idx](dim); }

        // Optional bounding-box computation
        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const { return false; }
    };
    
    std::vector<Eigen::VectorXf> findClusterCenters(
        const std::vector<Eigen::VectorXf>& prevClusterCenters, 
        const std::vector<Eigen::VectorXf>& newPoints,
        float   distanceThreshold,
        int     miniBatchSize           = 0,
        bool    selectRandomMiniBatch   = false, 
        bool    shuffleMiniBatch        = false, 
        int     maxIteration            = 0,
        bool    updateCenters           = true
    );

    std::vector<Eigen::VectorXf> refineClusterCenters(
        const std::vector<Eigen::VectorXf>& centers,
        const std::vector<Eigen::VectorXf>& points
    );

    // Returns {index, squaredDistance} of the nearest cluster center
    std::pair<int, float> getNearestCluster(
        const std::vector<Eigen::VectorXf>& centers, 
        const Eigen::VectorXf& point
    );

    class ClusterRefiner {
    public:
        ClusterRefiner(const std::vector<Eigen::VectorXf>& centers);
        ~ClusterRefiner();

        void update(const std::vector<Eigen::VectorXf>& points);
        std::vector<Eigen::VectorXf> getRefinedCenters();

    private:
        struct Impl;
        Impl* pImpl;
    };

    std::vector<Eigen::VectorXf> refineClusterCenters(
        const std::vector<Eigen::VectorXf>& centers,
        const std::vector<Eigen::VectorXf>& points
    );

}
