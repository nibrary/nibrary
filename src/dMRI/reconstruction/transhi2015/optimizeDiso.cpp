#include "recon_transhi2015.h"

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

#include <simde/hedley.h>
#include <Eigen/Dense>
#include <proxsuite/proxqp/dense/dense.hpp>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

using namespace NIBR;

using namespace proxsuite::proxqp;

// Compartment optimization function
REASON_OF_TERMINATION TranShi2015::optimizeDiso(Eigen::MatrixXd &aH,  Eigen::VectorXd &g, 
                                      Eigen::MatrixXd &C,   Eigen::VectorXd &l, 
                                      Eigen::VectorXd &s_, 
                                      Eigen::VectorXd &d,
                                      Eigen::MatrixXd &A_plus, 
                                      double &D_exAx_iso, 
                                      const double xi,
                                      int &iterCnt)
{
    // TIC();
    // Proxsuite api is used for the optimization loop so we can more efficiently update the aH and q values.
    proxsuite::proxqp::dense::QP<double> qp(aH.rows(), qp_A.rows(), C.rows(), false, proxsuite::proxqp::DenseBackend::Automatic);
    
    qp.init(aH, g, qp_A, qp_b, C, l, proxsuite::nullopt, true);

    qp.settings.initial_guess   = proxsuite::proxqp::InitialGuessStatus::WARM_START_WITH_PREVIOUS_RESULT;
    qp.settings.max_iter        = 100; //ProxSuite default: 10000
    qp.settings.max_iter_in     = 500; //ProxSuite default: 1500
    
    auto solveQP=[&]()->void {
        
        // Updating aH and g
        aH = A_plus.transpose() * A_plus;
        g = -A_plus.transpose() * s_ + xi * I;

        // Solving QP
        qp.update(aH, g, proxsuite::nullopt, proxsuite::nullopt, proxsuite::nullopt, proxsuite::nullopt, proxsuite::nullopt);
        qp.solve();
        d = qp.results.x;
    };
    
    bool solved = true;

    Eigen::VectorXd grad_vec1, grad_vec2;

    for(iterCnt = 0; iterCnt < (maxIter/constraintUpdateCount); iterCnt++) {
                
        // Calculates gradient expressed in eq. 16 in the paper
        grad_vec1 =  (((-bval).array() * (-bval * D_exAx_iso).array().exp()) * d(d.size() - 2));
        grad_vec2 = -(s_-A_plus*d);

        double gradient = grad_vec2.transpose() * grad_vec1;

        // Calculate gradient step
        double step = gradient*deltaStep;
        step = std::clamp(step,-init_D_exAx_iso*0.1,init_D_exAx_iso*0.1);

        // Checking if convergence is reached.
        
        if(std::abs(step) < EPS7) {
            if (!solved) {solveQP();}
            // TOC(MSG_DEBUG);
            return TOO_SMALL_GRADIENT;
        }
        
        if( (D_exAx_iso-step) < EPS4) {
            if (!solved) {solveQP();}
            // TOC(MSG_DEBUG);
            return TOO_SMALL_DISO;
        }
        

        // Descent
        D_exAx_iso -= step;
        A_plus.col(A_plus.cols() - 2) = (-bval * D_exAx_iso).array().exp();

        if (fastOptimization) {

            if (iterCnt < 20) {
                solveQP();
                solved = true;
            } else if ((iterCnt  <  50) && (iterCnt %  2 == 0)) {
                solveQP();
                solved = true;
            } else if ((iterCnt  < 100) && (iterCnt %  4 == 0)) {
                solveQP();
                solved = true;
            } else if ((iterCnt >= 100) && (iterCnt % 10 == 0)) {
                solveQP();
                solved = true;
            } else if ((iterCnt+1) == ((maxIter/constraintUpdateCount))) {
                solveQP();
                solved = true;
            } else {
                solved = false;
            }

        } else {
            solveQP();
        }

    }

    // TOC(MSG_DEBUG);
    return END_OF_ITERATION;

}
