#include "orientation.h"

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-result"
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wpedantic"
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma clang diagnostic ignored "-Wunused-local-typedefs"
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic ignored "-Warray-bounds"
#pragma clang diagnostic ignored "-Wtype-limits"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wreorder"
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wuninitialized"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-result"
#pragma clang diagnostic ignored "-Wcomment"
#endif

#include <Eigen/Dense>
#include <Eigen/SVD>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <optional>

// Below implements orientation handling functions, inspired from:
// https://github.com/nipy/nibabel/blob/d23a8336582d07d683d3af73cf097c3be2de9af8/nibabel/orientations.py#L356

using namespace Eigen;

const std::vector<std::pair<std::string, std::string>> labels = {{"L", "R"}, {"P", "A"}, {"I", "S"}};

using Orientation = std::vector<std::pair<int, int>>;

// Function to compute the orientation of input axes from a 4x4 affine matrix
Orientation io_orientation(const MatrixXd& affine) {

    MatrixXd RZS = affine.block(0, 0, 3, 3);

    // Compute norms for each column and handle zero norms to avoid division by zero
    VectorXd zooms = RZS.colwise().norm();
    for (int i = 0; i < zooms.size(); ++i) {
        zooms(i) = zooms(i) != 0 ? 1.0 / zooms(i) : 1.0;
    }

    // Normalize columns by zooms
    MatrixXd RS = RZS.array().rowwise() * zooms.transpose().array();
    JacobiSVD<MatrixXd> svd(RS, ComputeFullU | ComputeFullV);
    
    MatrixXd R = svd.matrixU() * svd.matrixV().transpose();

    Orientation ornt(3, std::make_pair(-1, 0)); // Default to -1 for dropped axes
    std::vector<bool> used(3, false);           // Track used output axes to ensure unique assignments

    for (int in_ax = 0; in_ax < 3; ++in_ax) {
        VectorXd col   = R.col(in_ax).cwiseAbs();
        int best_axis  = -1;
        double max_val = 0.0;

        for (int out_ax = 0; out_ax < 3; ++out_ax) {
            if (used[out_ax]) continue;
            if (col[out_ax] > max_val) {
                max_val     = col[out_ax];
                best_axis   = out_ax;
            }
        }

        if (best_axis != -1) {
            ornt[in_ax]     = std::make_pair(best_axis, R(best_axis, in_ax) < 0 ? -1 : 1);
            used[best_axis] = true;
        }
    }

    return ornt;
}

// Function to convert orientation to axis direction codes using fixed labels
std::vector<std::string> ornt2axcodes(const Orientation& ornt) {
    std::vector<std::string> axcodes;

    for (const auto& [axno, direction] : ornt) {
        if (axno == -1) {
            axcodes.push_back("None"); // Default value for dropped axes
            continue;
        }

        const auto& labelPair = labels[axno];
        axcodes.push_back(direction == 1 ? labelPair.second : labelPair.first);
    }

    return axcodes;
}

// Main function to compute axis direction codes from a 4x4 affine matrix
std::vector<std::string> aff2axcodesEigen(const Eigen::MatrixXd& aff) {
    auto ornt = io_orientation(aff); // Compute orientation
    return ornt2axcodes(ornt); // Convert to axis codes
}


std::vector<std::string> NIBR::aff2axcodes(const float aff[3][4]) {
    
    Eigen::MatrixXd affMatrix(4, 4);
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            affMatrix(i, j) = static_cast<double>(aff[i][j]);
        }
    }
    
    affMatrix.row(3) << 0, 0, 0, 1;

    return aff2axcodesEigen(affMatrix);
}

/*
UNUSED
void NIBR::aff2RAS(float outAff[3][4], const float inpAff[3][4], const int64_t dims[3]) {
    
    auto axcodes = aff2axcodes(inpAff);
    
    // First, copy inpAff to outAff
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            outAff[i][j] = inpAff[i][j];
        }
    }

    // Modify the last column of the matrix based on dims and axcodes
    for (int i = 0; i < 3; ++i) {
        if (axcodes[i] != "R" && axcodes[i] != "A" && axcodes[i] != "S") {
            outAff[i][3] = dims[i] - inpAff[i][3];
        }
    }

    // Determine the order of R, A, S in axcodes and reorder rows to achieve RAS
    std::vector<int> order(3, 0); // Store indices to sort axcodes as RAS
    for (int i = 0; i < 3; ++i) {
        if (axcodes[i] == "R") order[0] = i;
        else if (axcodes[i] == "A") order[1] = i;
        else if (axcodes[i] == "S") order[2] = i;
    }

    // Temporary matrix to hold the reordered affine
    float tempAff[3][4];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            tempAff[i][j] = outAff[order[i]][j];
        }
    }

    // Copy the reordered affine back to outAff
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            outAff[i][j] = tempAff[i][j];
        }
    }
}
*/



/*

// Below is a general approach for arbitrary dimensions, which we do not need in our case.
// So, we will use the code above to simplify things.

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <limits>
#include <vector>
#include <optional>
#include <string>
#include <stdexcept>

using namespace Eigen;

using Orientation = std::vector<std::pair<std::optional<int>, int>>;
using LabelPair   = std::pair<std::string, std::string>;
using Labels      = std::vector<LabelPair>;

Orientation io_orientation(const MatrixXd& affine, double tol = -1.0) {
    int q = affine.rows() - 1;
    int p = affine.cols() - 1;
    MatrixXd RZS = affine.block(0, 0, q, p);

    // Compute zooms to normalize columns
    VectorXd zooms = RZS.colwise().norm();
    for (int i = 0; i < zooms.size(); ++i) {
        if (zooms(i) == 0) zooms(i) = 1.0;
    }

    MatrixXd RS = RZS.array().rowwise() / zooms.transpose().array();
    JacobiSVD<MatrixXd> svd(RS, ComputeFullU | ComputeFullV);
    VectorXd S = svd.singularValues();

    if (tol < 0) {
        tol = S.maxCoeff() * std::max(RS.rows(), RS.cols()) * std::numeric_limits<double>::epsilon();
    }

    MatrixXd R = svd.matrixU() * svd.matrixV().transpose();

    Orientation ornt;
    std::vector<bool> used(q, false); // Track used output axes to handle ties

    for (int in_ax = 0; in_ax < p; ++in_ax) {
        VectorXd col = R.col(in_ax).cwiseAbs();
        std::optional<int> best_axis;
        double max_val = 0.0;

        for (int out_ax = 0; out_ax < q; ++out_ax) {
            if (used[out_ax]) continue; // Skip already used axes
            if (col[out_ax] > max_val) {
                max_val = col[out_ax];
                best_axis = out_ax;
            }
        }

        if (best_axis.has_value()) {
            ornt.emplace_back(best_axis, R(best_axis.value(), in_ax) < 0 ? -1 : 1);
            used[best_axis.value()] = true; // Mark this axis as used
        } else {
            // Handle case where an input axis doesn't correspond well to any output axis
            ornt.emplace_back(std::std::nullopt, 0);
        }
    }

    // Handle any remaining axes if p > q
    while (ornt.size() < static_cast<size_t>(p)) {
        ornt.emplace_back(std::std::nullopt, 0);
    }

    return ornt;
}


std::vector<std::optional<std::string>> ornt2axcodes(
    const Orientation& ornt,
    const Labels& labels = {{"L", "R"}, {"P", "A"}, {"I", "S"}}
) {
    std::vector<std::optional<std::string>> axcodes;

    for (const auto& [axno, direction] : ornt) {
        // Handle dropped axes
        if (!axno.has_value()) {
            axcodes.emplace_back(std::std::nullopt); // Append a dropped axis indicator
            continue;
        }

        // Validate direction
        if (direction != 1 && direction != -1) {
            throw std::invalid_argument("Direction should be -1 or 1.");
        }

        // Assign label based on direction
        try {
            const auto& labelPair = labels.at(axno.value()); // Use .at() for bounds checking
            axcodes.emplace_back(direction == 1 ? labelPair.second : labelPair.first);
        } catch (const std::out_of_range&) {
            throw std::invalid_argument("Axis index " + std::to_string(axno.value()) + " is out of range.");
        }
    }

    return axcodes;
}

std::vector<std::optional<std::string>> aff2axcodes(
    const Eigen::MatrixXd& aff,
    const Labels& labels = {{"L", "R"}, {"P", "A"}, {"I", "S"}},
    const std::optional<double>& tol = std::std::nullopt
) {
    // Compute the orientation from the affine matrix
    auto ornt = io_orientation(aff, tol);

    // Convert the orientation to axis codes using the provided labels
    auto axcodes = ornt2axcodes(ornt, labels);

    return axcodes;
}
*/

