#include "recon_transhi2015.h"
#include "constraints.h"

#include "../../gradient/dMRI_grad.h"
#include "math/sphere.h"
#include "image/image_operators.h"
#include "surface/surface.h"

using namespace NIBR;

// Converts nested std::vector of floats into and Eigen matrix
Eigen::MatrixXf VecToMatf(std::vector<std::vector<float>> Vec) {
    size_t rows = Vec.size();
    size_t cols = Vec[0].size();
    Eigen::MatrixXf output(rows, cols);
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            output(i, j) = Vec[i][j];
        }
    }
    return output;
}

// Converts nested std::vector of doubles into and Eigen matrix
Eigen::MatrixXd VecToMatd(std::vector<std::vector<double>> Vec) {
    size_t rows = Vec.size();
    size_t cols = Vec[0].size();
    Eigen::MatrixXd output(rows, cols);
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            output(i, j) = Vec[i][j];
        }
    }
    return output;
}

std::vector<int> MakeImageMask(int i, int j, int k, int kernel_size, Image<bool> mask) {
    
    std::vector<int> newMask;
    int i_ker = i + 2;
    int j_ker = j + 1;
    int k_ker = k;
    for(int i_ind = i_ker - kernel_size; i_ind <= i_ker + kernel_size; i_ind++) {
        for(int j_ind = j_ker - kernel_size; j_ind <= j_ker + kernel_size; j_ind++) {
            for(int k_ind = k_ker - kernel_size; k_ind <= k_ker + kernel_size; k_ind++) {
                int vox_ind = mask.sub2ind(i_ind, j_ind, k_ind);
                newMask.push_back(vox_ind);
            }
        }
    }
    return newMask;
}

void TranShi2015::prep()
{

    // sphere.printInfo();
    // disp(MSG_INFO,"Sphere vert 1: %.2f, %.2f, %.2f", sphere.vertices[0][0], sphere.vertices[0][1], sphere.vertices[0][2]);


    std::vector<float> bvals;
    std::vector<Point3D> bvecs;
    for(size_t i = 0; i < gradTable.size(); i++) {
        float b = gradTable[i][3];
        if(b == 0) gradTable[i] = {0, 0, 0, 0};
        if(b < bValLow) B0_ind.push_back(i);
        if(b >= bValLow && b <= bValHigh) {
            DTI_ind.push_back(i);
            bvecs.push_back({gradTable[i][0], gradTable[i][1], gradTable[i][2]});
            bvals.push_back(std::round(b / 50) * 50); // Assume all bvals are a multiple of 50. To speed up the getG.
        }
    }
    std::vector<double> bval_double(bvals.begin(), bvals.end());
    Eigen::Map<Eigen::VectorXd> bval_conversion(bval_double.data(), bval_double.size());
    bval = bval_conversion; // bvals in Eigen format


    // Finding the minimum constraint index
    init_minConstraintIndex = 0;
    for(size_t i = 0; i < TRANSHI2015CONSTRAINTS.size(); i++) {
        if(TRANSHI2015CONSTRAINTS[i].size() >= static_cast<size_t>(init_minNumConstraint)) {
            init_minConstraintIndex = i;
            break;
        }
    }

    // Constraint points are selected from a sphere.
    std::vector<Point3D> vertices;
    for (int n = 0; n < 2562; n++) {
        vertices.push_back({DENSESPHEREVERT[n][0],DENSESPHEREVERT[n][1],DENSESPHEREVERT[n][2]});
    }
    std::vector<std::vector<float>> BS;
    SH_basis(BS, vertices, shOrder);
    Y_constraint = VecToMatf(BS).cast<double>();

    // The spherical function can be synthesized on a different sphere.
    // This is used to check the positive/negative energy ratio and the peaks.
    vertices.clear();
    for (int n = 0; n < sphere.nv; n++) {
        vertices.push_back({sphere.vertices[n][0],sphere.vertices[n][1],sphere.vertices[n][2]});
    }
    BS.clear();
    SH_basis(BS, vertices, shOrder);
    Y = VecToMatf(BS).cast<double>();


    // Setting up the Spherical harmonics basis and converting them to Eigen format
    std::vector<std::vector<float>> B;
    SH_basis(B, bvecs, shOrder);
    Eigen::MatrixXd B_mat = VecToMatf(B).cast<double>();

    Eigen::MatrixXd G = VecToMatd(G_mat(bvals, init_D_inAx, init_D_trapped, shOrder));

    if (B_mat.rows() != G.rows() || B_mat.cols() != G.cols()) {
        NIBR::disp(MSG_ERROR,"B and G must have the same dimensions");
    }

    init_A_plus = B_mat.cwiseProduct(G);
    size_t m    = init_A_plus.rows();
    size_t n    = init_A_plus.cols();
    Eigen::MatrixXd BGnew(m, n + 2);
    BGnew.block(0, 0, m, n) = init_A_plus;
    init_A_plus = BGnew;

    for(size_t r = 0; r < size_t(BGnew.rows()); r++) {
        init_A_plus(r, BGnew.cols()-2) = std::exp(-bval(r) * init_D_exAx_iso);
        init_A_plus(r, BGnew.cols()-1) = 1.0;
    }
    
    // I is basically the I matrix in the paper.
    // Ix gives intra-axonal compartment signal fraction
    // Below definition add two more columns, 
    // first is for the extra-axonal (alpha), and the other for the trapped (gamma)
    // Ix + alpha + gamma = 1, which is put as an equality constraint in the quadratic programming step
    I = Eigen::VectorXd(Y_constraint.cols() + 2);
    I.setZero();
    I(0) = sqrt(4 * M_PI);
    
    auto AllSumMtx = I;
    AllSumMtx.segment(AllSumMtx.size() - 2, 2).setOnes();

    // Equality constraints
    qp_A = AllSumMtx.transpose().replicate(n+2,1);
    qp_b = Eigen::VectorXd(AllSumMtx.size());
    qp_b.fill(1 - noiseFloor);

    // Create output images
    int indexOrder[7] = {3,0,1,2,4,5,6};
    FOD = new Image<float>(indexOrder);
    int64_t dims[4] = {dMRI->imgDims[0],dMRI->imgDims[1],dMRI->imgDims[2],int64_t(B[0].size())};
    FOD->create(4, dims, dMRI->pixDims, dMRI->ijk2xyz, true);

    TM = new Image<float>(indexOrder);
    dims[3] = 10;
    TM->create(4, dims, dMRI->pixDims, dMRI->ijk2xyz, true);

    voxelIndicesToReconstruct = getNonZeroIndices(mask);
    // voxelIndicesToReconstruct = MakeImageMask(72, 67, 52, 10, *mask);   //72, 59, 57, 10,
    
    // disp(MSG_DEBUG, "Y size:            [%d , %d]", Y.rows(),            Y.cols());
    // disp(MSG_DEBUG, "Y_constraint size: [%d , %d]", Y_constraint.rows(), Y_constraint.cols());

    NIBR::disp(MSG_INFO, "Preprocessing completed");

}
