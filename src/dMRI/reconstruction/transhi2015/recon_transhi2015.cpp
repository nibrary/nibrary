#include "recon_transhi2015.h"
#include "constraints.h"
#include "../../gradient/dMRI_grad.h"
#include "math/sphere.h"
#include "image/image_operators.h"
#include "surface/surface.h"
#include "surface/surface.h"
#include <optional>
#include <cstdint>
#include <chrono>

using namespace NIBR;

void NIBR::recon_transhi2015(Image<float>& FOD, Image<float>& TM, 
                             Image<float>& dMRI,
                             std::vector<std::vector<float>> grad, 
                             Image<bool>& mask, 
                             int     shOrder                    = TRANSHI2015_SH_ORDER,
                             double  deltaStep                  = TRANSHI2015_DELTASTEP,
                             double  init_D_inAx                = TRANSHI2015_INIT_D_INAX,
                             double  init_D_trapped             = TRANSHI2015_INIT_D_TRAPPED,
                             double  init_D_exAx_iso            = TRANSHI2015_INIT_D_EXAX,
                             int     init_minNumConstraint      = TRANSHI2015_INIT_MINNUMCONSTRAINTS,
                             int     constraintUpdateCount      = TRANSHI2015_CONSTRAINT_UPDATE_COUNT,
                             bool    fastOptimization           = TRANSHI2015_FAST_OPTIMIZATION,
                             double  bValLow                    = TRANSHI2015_BVAL_LOW,
                             double  bValHigh                   = TRANSHI2015_BVAL_HIGH,
                             int     maxIter                    = TRANSHI2015_MAXITER, 
                             double  xi_init                    = TRANSHI2015_XI_INIT, 
                             double  xi_step                    = TRANSHI2015_XI_STEP, 
                             int     xi_stepCount               = TRANSHI2015_XI_STEPCOUNT,
                             int     maxCrossings               = TRANSHI2015_MAXCROSSINGS,
                             double  noiseFloor                 = TRANSHI2015_NOISEFLOOR)
{
    TIC();
    TranShi2015 model(dMRI,grad);
    TOC(MSG_DETAIL);

    TIC();
    model.setMask(mask);
    model.shOrder                   = shOrder;
    model.deltaStep                 = deltaStep;
    model.init_D_inAx               = init_D_inAx;
    model.init_D_trapped            = init_D_trapped;
    model.init_D_exAx_iso           = init_D_exAx_iso;
    model.init_minNumConstraint     = init_minNumConstraint;
    model.constraintUpdateCount     = constraintUpdateCount;
    model.fastOptimization          = fastOptimization;
    model.bValLow                   = bValLow;
    model.bValHigh                  = bValHigh;
    model.maxIter                   = maxIter;
    model.xi_init                   = xi_init;
    model.xi_step                   = xi_step;
    model.xi_stepCount              = xi_stepCount;
    model.maxCrossings              = maxCrossings;
    model.noiseFloor                = noiseFloor;
    TOC(MSG_DETAIL);

    model.print();
    
    model.recon();

    FOD = model.getFOD();
    TM  = model.getTM();
}

TranShi2015::TranShi2015(Image<float>& _dMRI,std::vector<std::vector<float>>& _grad)
{
    dMRI        = &_dMRI; 
    gradTable   = _grad;

    // rotateGradTable(&grad, dMRI.ijk2xyz);
    for (std::vector<float>& g : gradTable) {
        g[0] = -g[0];
        // std::swap(g[1],g[2]);
    }
    

    // Create the sphere
    std::vector<std::vector<float>> vertices;
    for (int n = 0; n < 1024; n++) {
        vertices.push_back({FULLSPHERE1024_VERTS[n][0],FULLSPHERE1024_VERTS[n][1],FULLSPHERE1024_VERTS[n][2]});
    }
    std::vector<std::vector<int>> faces;
    for (int n = 0; n < 2044; n++) {
        faces.push_back({FULLSPHERE1024_FACES[n][0],FULLSPHERE1024_FACES[n][1],FULLSPHERE1024_FACES[n][2]});
    }

    Surface tmp(vertices,faces);
    sphere = tmp;
    sphere.getNeighboringVertices();
    
}


void TranShi2015::recon()
{

    TIC();
    prep();
    TOC(MSG_DETAIL);

    // std::ofstream file("data.csv");      // Used with the gradient_csv_save function

    auto reconstructVoxel = [&](NIBR::MT::TASK task)->void {

        int64_t voxelIndex = voxelIndicesToReconstruct[task.no]; // This is the index of the voxel to reconstruct

        int64_t i,j,k; // i,j,k are the 3D coordinates of the voxel in the image, which are computed in the line below
        mask->ind2sub(voxelIndex,i,j,k);

        float* s    = dMRI->at(i,j,k,0);  // s is pointing to the input diffusion signal measured for this voxel, i.e., the values along the 4th dimension. s has D elements where D = dMRI.imgDims[3]. So the signal is allocated between s[0] to s[dMRI.imgDims[3]-1].
        float* fod  =  FOD->at(i,j,k,0);  // fod is pointing to the output FOD coefficients. These need to be computed below. The data is already allocated when FOD.create was called above. fod has N elements. So we can write on fod[0] to fod[N-1].
        float* tm   =   TM->at(i,j,k,0);  // tm is pointing to the output TM values. These need to be computed below. The data is already allocated when TM.create was called above. tm has 5 elements. So we can write on tm[0] to tm[4].

        // Calculating the averages of Bval < TranShi2015_default.bValLow indices for the normalization
        double S0 = 0.0;
        for(size_t b0 : B0_ind) {
            S0 += (s[b0]>0.0) ? s[b0] : 0.0;
        }

        if (std::fabs(S0)<EPS8) return;

        S0 = S0 / double(B0_ind.size());

        // Normalizing and saving the accepted signal values to a new vector
        Eigen::VectorXd s_NF(DTI_ind.size());
        Eigen::Index ind = 0;
        for(size_t b : DTI_ind) {
            s_NF(ind) = (s[b] / (S0 + EPS8)) - noiseFloor;
            if (s_NF(ind)>1.0) s_NF(ind) = 1.0;
            ind++;
        }

        float  targetPeakThresh = 0.025;
        int    targetPeakCount  = maxCrossings;

        Eigen::VectorXd d(init_A_plus.rows());  // d: spherical harmonic coeff of FOD
        Eigen::VectorXd y(Y.rows());            // y: spherical function of FOD

        double D_exAx_iso       = 0;
        int minConstraintIndex  = init_minConstraintIndex;
        float positivityRatio   = 0;
        int iterCnt             = 0;
        int totIterCnt          = 0;
        double residual         = 1;        

        // Regularization is increased starting from xi by steps of xi_stepsize,
        // until the desired number of peaks is reached or no peaks are found within the range of xi_NumSteps.
        // Otherwise, max regularization of (xi_init + xi_NumSteps*xi_stepsize)

        double xi;
        // double final_gradient;   // Used with the optimize_save_final_gradient function
        
        REASON_OF_TERMINATION res_opt = OPTIM_END_UNSET;

        for(xi = xi_init; xi <= (xi_init+xi_step*xi_stepCount); xi += xi_step) {

            disp(MSG_DEBUG,"task.no: %d, voxel index: %d", task.no, voxelIndex);

            D_exAx_iso = init_D_exAx_iso;

            Eigen::MatrixXd A_plus = init_A_plus;            
            
            // Initialization constraints for the quadratic problem.
            // Finds a solution to the problem:
            // min 1/2 x'Hx + g'x
            // s.t. Ax = b, l <= Cx <= u
            
            Eigen::MatrixXd aH =  A_plus.transpose() * A_plus;
            Eigen::VectorXd g  = -A_plus.transpose() * s_NF + xi * I;

            // "setConstraints" function constrains most of the FODs to a positive value. 
            // It also suggests an initial solution to the optimization problem by solving
            // the problem for the fixed (initial) D_exAx_iso value.
            // The output "d", basically contains the FOD coefficients, as well as the signal fractions
            // for the extra-axonal and the trapped compartments.
            // "d" is used as the initial condition for the main compartment optimization step, "optimize", below.
            // "setConstraints" also sets the inequality constraint l <= Cx, which goes into "optimize".
            
            Eigen::MatrixXd C;
            Eigen::VectorXd l;
            d = setConstraints(aH, g, C, l, minConstraintIndex, positivityRatio);
            
            // "optimize" solves for the D_exAx_iso value, for the previously set constraints and compartment fractions.
            totIterCnt = 0;
            for (int updtCnt = 0; updtCnt < constraintUpdateCount; updtCnt++) {
                disp(MSG_DEBUG,"D_exAx_iso optimizations starting: %d", updtCnt);
                res_opt     = optimizeDiso(aH, g, C, l, s_NF, d, A_plus, D_exAx_iso, xi, iterCnt);
                totIterCnt += iterCnt;
                disp(MSG_DEBUG,"Finished at %d iteration", iterCnt);
                if ((res_opt == TOO_SMALL_GRADIENT) || (res_opt == TOO_SMALL_DISO)) break;
                updateConstraints(aH, g, C, l, d, minConstraintIndex, positivityRatio);
                disp(MSG_DEBUG,"Constraints updated");
            }

            // d = optimize_save_final_gradient(aH, g, A, b, C, l, s_NF, d, A_plus, D_exAx_iso, xi, final_gradient);
            // d = gradient_csv_save(aH, g, A, b, C, l, s_NF, d, A_plus, D_exAx_iso, xi, file);

            // TIC();
            // FOD sparsity check. y is the spherical function of the FOD
            y = Y * d.head(d.size() - 2);

            // Returns an array consisting of the peaks in y
            std::vector<float> peakValues = getPeakValues(y, targetPeakThresh);

            if (peakValues.empty()) {
                residual = (A_plus*d-s_NF).norm() / s_NF.norm();
                // TOC(MSG_DEBUG);
                break;
            }

            // Test for how dominant the largest fibers are with respect to all the peaks.
            if (xi == xi_init) {
                targetPeakThresh = getTargetPeakCountAndThresh(peakValues, targetPeakCount);
            }

            int peakCount = 0;
            for(auto p : peakValues) {
                if (p > targetPeakThresh) peakCount++;
            }
            
            if (peakCount <= targetPeakCount) {
                residual = (A_plus*d-s_NF).norm() / s_NF.norm();
                // TOC(MSG_DEBUG);
                break;
            }

            // TOC(MSG_DEBUG);
            
        }

        // Save values to fod and tm
        for(int64_t n = 0; n < FOD->imgDims[3]; n++) {
            fod[n] = static_cast<float>(d(n));
        }

        float yPosSum = y.array().max(0.0).sum();

        tm[0]  = static_cast<float>(4.0 * M_PI * yPosSum / double(y.size()));           // signal fraction of intra axonal compartment
        tm[1]  = static_cast<float>(d(d.size() - 2));                                   // alpha: signal fraction of extra axonal compartment
        tm[2]  = static_cast<float>(d(d.size() - 1) + noiseFloor);                      // gamma: signal fraction of trapped compartment
        tm[3]  = static_cast<float>(D_exAx_iso);                                        // Isotropic diffusivity component of extra-axonal space
        tm[4]  = static_cast<float>(res_opt);
        tm[5]  = static_cast<float>(totIterCnt);
        tm[6]  = static_cast<float>(xi);                                                // regularization parameter
        tm[7]  = static_cast<float>(TRANSHI2015CONSTRAINTS[minConstraintIndex].size());
        tm[8]  = static_cast<float>(positivityRatio);
        tm[9]  = static_cast<float>(residual);                                          // Ration between the euclidean norm of the residuals and signal

    };
    NIBR::MT::MTRUN(voxelIndicesToReconstruct.size(), "Fitting compartment model", reconstructVoxel);

    // file.close();      // Used with the gradient_csv_save function

    return;
}
