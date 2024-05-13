#include "image_operators.h"

using namespace NIBR;
using namespace SF;
using namespace SH;

void NIBR::sf2sh(NIBR::Image<float>* out, NIBR::Image<float>* inp, std::vector<std::vector<float>>& coords, int shOrder) {
    
    std::vector<std::vector<float>> Ylm;
    SH_basis(Ylm, coords, shOrder);

    int coeffCount = Ylm[0].size();

    // Create image if needed
    if ( (out == NULL) || (out->data == NULL) ) {
        int64_t imgDims[4] = {inp->imgDims[0], inp->imgDims[1], inp->imgDims[2], coeffCount};
        out->create(4, &imgDims[0], inp->pixDims, inp->ijk2xyz, true);
    }

    // We will find and only process those voxels which have non-zero values
    std::vector<std::vector<int64_t>> nnzVoxelSubs;
    auto findNonZeroVoxels = [&](NIBR::MT::TASK task)->void {
        
        int64_t i   = task.no;
        
        for (int64_t j=0; j<inp->imgDims[1]; j++)
            for (int64_t k=0; k<inp->imgDims[2]; k++)
                for (int64_t t=0; t<inp->imgDims[3]; t++) {

                    if (inp->data[inp->sub2ind(i,j,k,t)]!=0){
                        std::vector<int64_t> tmp{i,j,k};
                        NIBR::MT::PROC_MX().lock();
                        nnzVoxelSubs.push_back(tmp);
                        NIBR::MT::PROC_MX().unlock();
                        break;
                    }
                    
                }
        
    };
    NIBR::MT::MTRUN(inp->imgDims[0],findNonZeroVoxels);

    // Apply spherical harmonics expansion
    auto shExpansion = [&](NIBR::MT::TASK task)->void {
        
        std::vector<int64_t> sub = nnzVoxelSubs[task.no];
        
        for (int n=0; n<coeffCount; n++) {
            out->data[out->sub2ind(sub[0],sub[1],sub[2],n)] = 0;
            for (int64_t t=0; t<inp->imgDims[3]; t++)
                out->data[out->sub2ind(sub[0],sub[1],sub[2],n)] += Ylm[t][n]*inp->data[inp->sub2ind(sub[0],sub[1],sub[2],t)];
        }
        
    };
    NIBR::MT::MTRUN(nnzVoxelSubs.size(),"Applying spherical harmonics expansion",shExpansion);

    SH::clean();

}


void NIBR::sh2sf(NIBR::Image<float>* out, NIBR::Image<float>* inp, std::vector<std::vector<float>>& coords) {
    
    int shOrder = getSHOrderFromNumberOfCoeffs(inp->imgDims[3]);

    std::vector<std::vector<float>> Ylm;
    SH_basis(Ylm, coords, shOrder);

    // Create image if needed
    if ( (out == NULL) || (out->data == NULL) ) {
        int64_t imgDims[4] = {inp->imgDims[0], inp->imgDims[1], inp->imgDims[2], int(coords.size())};
        out->create(4, &imgDims[0], inp->pixDims, inp->ijk2xyz, true);
    }

    // We will find and only process those voxels which have non-zero values
    std::vector<std::vector<int64_t>> nnzVoxelSubs;
    auto findNonZeroVoxels = [&](NIBR::MT::TASK task)->void {
        
        int64_t i   = task.no;
        
        for (int64_t j=0; j<inp->imgDims[1]; j++)
            for (int64_t k=0; k<inp->imgDims[2]; k++)
                for (int64_t t=0; t<inp->imgDims[3]; t++) {

                    if (inp->data[inp->sub2ind(i,j,k,t)]!=0){
                        std::vector<int64_t> tmp{i,j,k};
                        NIBR::MT::PROC_MX().lock();
                        nnzVoxelSubs.push_back(tmp);
                        NIBR::MT::PROC_MX().unlock();
                        break;
                    }
                    
                }
        
    };
    NIBR::MT::MTRUN(inp->imgDims[0],findNonZeroVoxels);

    float scale = 4.0 * PI / float(coords.size());

    // Apply spherical harmonics synthesis
    auto shSynthesis = [&](NIBR::MT::TASK task)->void {
        
        std::vector<int64_t> sub = nnzVoxelSubs[task.no];
        
        for (size_t n = 0; n < coords.size(); n++) {
            out->data[out->sub2ind(sub[0],sub[1],sub[2],n)] = 0;
            for (int64_t t = 0; t < inp->imgDims[3]; t++)
                out->data[out->sub2ind(sub[0],sub[1],sub[2],n)] += scale * Ylm[n][t]*inp->data[inp->sub2ind(sub[0],sub[1],sub[2],t)];
        }
        
    };
    NIBR::MT::MTRUN(nnzVoxelSubs.size(),"Applying spherical harmonics synthesis",shSynthesis);

}


/*
// Here is the alternative approach for sh2sf conversion using precomputed coefficients
void NIBR::sh2sf(NIBR::Image<float>* out, NIBR::Image<float>* inp, std::vector<std::vector<float>>& coords) {
    
    int shOrder = getSHOrderFromNumberOfCoeffs(inp->imgDims[3]);

    // Precompute SH constants if not already computed
    SH::precompute(shOrder,XYZ,1024);


    // Create image if needed
    if ( (out == NULL) || (out->data == NULL) ) {
        int64_t imgDims[4] = {inp->imgDims[0], inp->imgDims[1], inp->imgDims[2], int(coords.size())};
        out->create(4, &imgDims[0], inp->pixDims, inp->ijk2xyz, true);
    }

    // We will find and only process those voxels which have non-zero values
    std::vector<std::vector<int64_t>> nnzVoxelSubs;
    auto findNonZeroVoxels = [&](NIBR::MT::TASK task)->void {
        
        int64_t i   = task.no;
        
        for (int64_t j=0; j<inp->imgDims[1]; j++)
            for (int64_t k=0; k<inp->imgDims[2]; k++)
                for (int64_t t=0; t<inp->imgDims[3]; t++) {

                    if (inp->data[inp->sub2ind(i,j,k,t)]!=0){
                        std::vector<int64_t> tmp{i,j,k};
                        NIBR::MT::PROC_MX().lock();
                        nnzVoxelSubs.push_back(tmp);
                        NIBR::MT::PROC_MX().unlock();
                        break;
                    }
                    
                }
        
    };
    NIBR::MT::MTRUN(inp->imgDims[0],findNonZeroVoxels);

    float scale = 4.0 * PI / float(coords.size());

    // Apply spherical harmonics synthesis
    auto shSynthesis = [&](NIBR::MT::TASK task)->void {
        
        std::vector<int64_t> sub = nnzVoxelSubs[task.no];
        
        for (int n = 0; n < coords.size(); n++) {
            float dir[3] = {coords[n][0],coords[n][1],coords[n][2]};
            float shCoeffs[inp->imgDims[3]];
            for (int t = 0; t < inp->imgDims[3]; t++) {
                shCoeffs[t] = inp->data[inp->sub2ind(sub[0],sub[1],sub[2],t)];
            }
            out->data[out->sub2ind(sub[0],sub[1],sub[2],n)] = scale*SH::toSF(shCoeffs,dir);
        }
        
    };
    NIBR::MT::MTRUN(nnzVoxelSubs.size(),"Applying spherical harmonics synthesis",shSynthesis);

}
*/



DATATYPE NIBR::getImageDataType(std::string imgFname) {
    Image<int8_t> tmpImg(imgFname);
    return tmpImg.inputDataType;
}
