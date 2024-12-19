#include "recon_transhi2015.h"

#include "constraints.h"

#ifdef __GNUC__

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
#include <simde/hedley.h>
#include <Eigen/Dense>
#include <proxsuite/proxqp/dense/dense.hpp>
#pragma GCC diagnostic pop

#else

#include <simde/hedley.h>
#include <Eigen/Dense>
#include <proxsuite/proxqp/dense/dense.hpp>

#endif

using namespace NIBR;

using namespace proxsuite::proxqp;
using proxsuite::nullopt; // c++17 simply use std::nullopt

#define CONSTRAINT_THRESHOLD 25

float TranShi2015::checkConstraints(const Eigen::VectorXd &d) {

    Eigen::VectorXd y = Y * d.head(d.size() - 2);

    float positiveSum = y.array().max(0.0).sum();
    float negativeSum = y.array().min(0.0).sum();
    if (std::fabs(negativeSum) < EPS8) negativeSum = EPS8; // to prevent division by 0

    float ratio = std::fabs(positiveSum / negativeSum);

    // disp(MSG_DEBUG, "positivityRatio: %.2f", ratio);

    return (ratio > CONSTRAINT_THRESHOLD) ? CONSTRAINT_THRESHOLD : ratio;

}

Eigen::VectorXd TranShi2015::setConstraints(Eigen::MatrixXd &aH,  Eigen::VectorXd &g,
                                            Eigen::MatrixXd  &C,  Eigen::VectorXd &l, 
                                            int &minConstraintIndex, float& positivityRatio)
{
    // TIC();
    Eigen::VectorXd d, y;

    int constraintSearchStepSize = 20;
    
    // Constraint loop. This loop constrains most of the FODs to a positive value. 
    //constraintThreshold determines the proportion of positive values to negative values.
    for (size_t constInd = init_minConstraintIndex; constInd < TRANSHI2015CONSTRAINTS.size(); constInd += constraintSearchStepSize) {
        
        // Updating the constraints
        auto BC_new     = Eigen::MatrixXd(TRANSHI2015CONSTRAINTS[constInd].size(), Y.cols());
        Eigen::Index i  = 0;

        for(size_t index : TRANSHI2015CONSTRAINTS[constInd]) {
            BC_new.row(i) = Y_constraint.row(index - 1);
            i++;
        }
        
        // Updating the QP inequality constraints
        size_t BC_m = BC_new.rows();
        size_t BC_n = BC_new.cols();

        C = Eigen::MatrixXd::Zero(BC_m + 2, BC_n + 2);

        C.block(0, 0, BC_m, BC_n) = BC_new;
        C(BC_m, BC_n)             = 1;
        C(BC_m + 1, BC_n + 1)     = 1;
        
        l = Eigen::VectorXd::Zero(BC_m + 2);
        
        // Solving QP
        Results<double> results_dense_solver = dense::solve<double>(aH, g, qp_A, qp_b, C, l, nullopt);
        
        d = results_dense_solver.x;

        minConstraintIndex  = constInd;

        positivityRatio     = checkConstraints(d);

        // disp(MSG_DEBUG, "Constraint iteration constInd: %d", constInd);

        if (positivityRatio == CONSTRAINT_THRESHOLD) break;

    }
    
    // TOC(MSG_DEBUG);
    return d;
}

void TranShi2015::updateConstraints(Eigen::MatrixXd &aH,  Eigen::VectorXd &g,
                                    Eigen::MatrixXd  &C,  Eigen::VectorXd &l, 
                                    Eigen::VectorXd   d,  
                                    int &minConstraintIndex, float& positivityRatio)
{

    // TIC();

    int constraintSearchStepSize = (checkConstraints(d) == CONSTRAINT_THRESHOLD) ? -10 : 10;
    
    // Constraint loop. This loop constrains most of the FODs to a positive value. 
    //constraintThreshold determines the proportion of positive values to negative values.
    
    for (size_t constInd = (minConstraintIndex+constraintSearchStepSize); ((constInd < TRANSHI2015CONSTRAINTS.size()) && (constInd >= init_minConstraintIndex)); constInd += constraintSearchStepSize) {
        
        // Updating the constraints
        auto BC_new     = Eigen::MatrixXd(TRANSHI2015CONSTRAINTS[constInd].size(), Y.cols());
        Eigen::Index i  = 0;

        for(size_t index : TRANSHI2015CONSTRAINTS[constInd]) {
            BC_new.row(i) = Y_constraint.row(index - 1); // because indices in TRANSHI2015CONSTRAINTS are in Matlab format and start with 1 not 0
            i++;
        }
        
        // Updating the QP inequality constraints
        size_t BC_m = BC_new.rows();
        size_t BC_n = BC_new.cols();

        Eigen::MatrixXd C_new = Eigen::MatrixXd::Zero(BC_m + 2, BC_n + 2);
        C_new.block(0, 0, BC_m, BC_n) = BC_new;
        C_new(BC_m, BC_n)             = 1;
        C_new(BC_m + 1, BC_n + 1)     = 1;
        
        Eigen::VectorXd l_new = Eigen::VectorXd::Zero(BC_m + 2);
        
        // Solving QP
        Results<double> results_dense_solver = dense::solve<double>(aH, g, qp_A, qp_b, C_new, l_new, nullopt);
        
        auto d_new = results_dense_solver.x;

        positivityRatio = checkConstraints(d_new);

        if (constraintSearchStepSize > 0) {            

            if (positivityRatio == CONSTRAINT_THRESHOLD) { 
                minConstraintIndex = constInd;
                C    = C_new;
                l    = l_new;
                d    = d_new;
                // TOC(MSG_DEBUG);
                break;
            }
            
        } else {

            if (positivityRatio == CONSTRAINT_THRESHOLD) { 
                minConstraintIndex = constInd;
                C    = C_new;
                l    = l_new;
                d    = d_new;
            } else {
                // TOC(MSG_DEBUG);
                break;
            }

        }
        
    }

    // TOC(MSG_DEBUG);

}