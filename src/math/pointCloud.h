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
            Eigen::VectorXf vec(size);
            for (size_t i = 0; i < size; ++i)
                vec(i) = p1[i];
            return (vec - points[idx_p2]).squaredNorm();
        }

        // Returns the dim'th component of the idx'th point in the class
        inline float kdtree_get_pt(const size_t idx, int dim) const { return points[idx](dim); }

        // Optional bounding-box computation
        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const { return false; }
    };    
}
