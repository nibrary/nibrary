#pragma once

#include <Eigen/Dense>
#include <vector>

namespace NIBR
{
    // Used for KD-tree for clustering and scoring 
    struct PointCloud {
        const std::vector<Eigen::VectorXf>* points;

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return points->size(); }

        // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class
        inline float kdtree_distance(const float* p1, const size_t idx_p2, size_t size) const {
            float dist = 0;
            const float* p2 = (*points)[idx_p2].data();
            for (size_t i = 0; i < size; ++i) {
                float d = p1[i] - p2[i];
                dist += d * d;
            }
            return dist;
        }

        // Returns the dim'th component of the idx'th point in the class
        inline float kdtree_get_pt(const size_t idx, int dim) const { return (*points)[idx](dim); }

        // Optional bounding-box computation
        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const { return false; }
    };
    
    class Clusterer {
    public:
        enum class Mode {
            DISCOVERY,  // Pass 1: Online updates, create new clusters
            REFINEMENT, // Pass 2: Accumulate sums/counts (update at end)
            COVERAGE    // Pass 3: Fixed centers, create new only if needed
        };

        Clusterer(int dim, float distanceThreshold, bool shuffle = true);
        ~Clusterer();
        
        Clusterer(const Clusterer&) = delete;
        Clusterer& operator=(const Clusterer&) = delete;

        // Set the operation mode
        void setMode(Mode mode);

        // Process a batch of points according to the current mode
        void process(const std::vector<Eigen::VectorXf>& points);

        // Update centers based on accumulated sums and counts.
        // - centers[i] = sums[i] / counts[i]
        // - Resets sums and counts for the next pass.
        // - Rebuilds KD-Tree.
        void updateCenters();

        // Returns {index, squaredDistance} of the nearest cluster center.
        // Returns {-1, NaN} if no centers exist.
        std::pair<int, float> assign(const Eigen::VectorXf& point);

        // Getters
        std::vector<Eigen::VectorXf> getClusterCenters()     const;
        std::vector<double>          getClusterCounts()      const;
        std::vector<size_t>          getClusterMedoids()     const;
        int                          getClusterCenterCount() const;

        // Configuration
        void setCenters(const std::vector<Eigen::VectorXf>& centers); 

    private:
        struct Impl;
        Impl* pImpl;
    };

}
