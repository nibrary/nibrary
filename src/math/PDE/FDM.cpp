#include "FDM.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

using namespace NIBR;


bool NIBR::FDMsolveLaplacianWithDirichletBC(Image<double>& problem, Image<double>& solution) {

    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::Triplet<double>      Triplet;

    std::vector<Triplet>   tripletList;
    std::vector<double>    b;                // Right-hand side vector
    std::map<int64_t, int> indexMap;         // Map from grid index to equation index

    const int64_t nx = problem.imgDims[0];
    const int64_t ny = problem.imgDims[1];
    const int64_t nz = problem.imgDims[2];

    const int64_t numel = problem.numel;

    // Neighbor offsets for 6-connected grid
    auto neighborOffsets = get3DNeighbors(CONN6);

    // Build index mapping and identify boundary points
    int eqIndex = 0;
    for (int64_t n = 0; n < numel; ++n) {
        if (problem.data[n] == FDM_INTERIOR) {
            indexMap[n] = eqIndex++;
        }   
    }

    b.resize(eqIndex, 0.0);

    // Assemble the matrix A and vector b
    for (int64_t k = 0; k < nz; ++k) {
        for (int64_t j = 0; j < ny; ++j) {
            for (int64_t i = 0; i < nx; ++i) {

                int64_t idx = problem.sub2ind(i, j, k);

                double val  = problem.data[idx];

                if (val == FDM_INTERIOR) {

                    int    row           = indexMap[idx];
                    double neighborCount = 0;

                    for (const auto& offset : neighborOffsets) {

                        int ni = i + offset[0];
                        int nj = j + offset[1];
                        int nk = k + offset[2];

                        if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {

                            int64_t nidx        = problem.sub2ind(ni, nj, nk);
                            double  neighborVal = problem.data[nidx];

                            if (neighborVal == FDM_INTERIOR) {
                                int col = indexMap[nidx];
                                tripletList.emplace_back(Triplet(row, col, -1.0));
                                neighborCount++;
                            } else if (neighborVal != FDM_EXTERIOR) {
                                // Dirichlet boundary condition
                                // neighborVal is the known value at the boundary
                                b[row] += neighborVal;
                                neighborCount++;
                            }
                            // FDM_EXTERIOR points are ignored
                        }
                        // Neighbors outside the grid are treated as FDM_EXTERIOR
                    }

                    // Set the diagonal entry
                    tripletList.emplace_back(Triplet(row, row, neighborCount));
                }
                // Else, skip FDM_EXTERIOR and Dirichlet BC points
            }
        }
    }

    // Build sparse matrix A
    SpMat A(eqIndex, eqIndex);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A.makeCompressed();

    // Iterative solver with preconditioning
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> solver;

    // Adjust solver parameters
    solver.setTolerance(1e-8);
    solver.setMaxIterations(1000);

    // Add regularization term
    SpMat I(eqIndex, eqIndex);
    I.setIdentity();
    A += 1e-6 * I;

    disp(MSG_INFO,"Solving Laplace's equation");

    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        disp(MSG_ERROR,"Solver initialization failed");
        return false;
    }

    Eigen::VectorXd solution_vec = solver.solve(Eigen::Map<Eigen::VectorXd>(b.data(), b.size()));
    if (solver.info() != Eigen::Success) {
        disp(MSG_ERROR,"Solving linear system failed");
        return false;
    }

    disp(MSG_INFO, "Solver converged in %d iterations", solver.iterations());
    disp(MSG_INFO, "Estimated error: %f", solver.error());

    // Map the solution back to the solution image
    solution.createFromTemplate(problem, true);

    for (int64_t n = 0; n < numel; ++n) {
        double val = problem.data[n];

        if (val == FDM_INTERIOR) {
            int eq_idx       = indexMap[n];
            solution.data[n] = solution_vec[eq_idx];
        } else if (val != FDM_EXTERIOR) {
            // Dirichlet boundary condition
            solution.data[n] = val;
        } else {
            // FDM_EXTERIOR point
            solution.data[n] = NAN;
        }
    }

    return true;
}
