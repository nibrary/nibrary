#include "seedImage.h"
#include "seeder.h"
#include "image/image_operators.h"
#include "image/image_morphological.h"

using namespace NIBR;

SeederOutputState SeedImage::getSeed(float* p, int t) {

    SeederOutputState state = checkSeedingLimits();

    if (state!=SEED_OK) {
        return state;
    }

    switch (mode) {
    
    case SEED_IMAGE_MASK:

        seed_img_mask.to_xyz(seed_indices[doRandomThings[t].uniform_int()],p);
        doRandomThings[t].randomizeWithinVoxel(p,seed_img_mask.pixDims);
        break;

    case SEED_IMAGE_RS:

        while (true) {

            seed_img_mask.to_xyz(seed_indices[doRandomThings[t].uniform_int()],p);
            doRandomThings[t].randomizeWithinVoxel(p,seed_img_pvf->pixDims);

            float val = (volInd<0) ? (*seed_img_pvf)(p) : (*seed_img_pvf)(p,volInd);

            if (val>=(doRandomThings[t].uniform_01()*max4rs)) {
                break;
            }

        }

        break;

    default:
        disp(MSG_ERROR,"Wrong seeding mode. Expected image mode.");
        return SEED_ERROR;
    }


    MT::PROC_MX().lock();
    state = checkSeedingLimits();
    if (state == SEED_OK)
        curSeed++;
    MT::PROC_MX().unlock();

    return state;

}

SeederOutputState SeedImage::getSeed(float* p, float*, int t) {
    return getSeed(p,t);
}

void SeedImage::computeSeedCountAndDensity() {

    double imgVol = 0;

    switch (mode) {
    
    case SEED_IMAGE_MASK:
        imgVol = seed_img_mask.pixDims[0]*seed_img_mask.pixDims[1]*seed_img_mask.pixDims[2]*double(seed_indices.size());
        break;

    case SEED_IMAGE_RS:
        imgVol = seed_img_pvf->pixDims[0]*seed_img_pvf->pixDims[1]*seed_img_pvf->pixDims[2]*double(seed_indices.size());
        break;

    default:
        break;
    }


    if (hasDensity) {
        count   = density*imgVol;
    } else {
        density = double(count)/imgVol;
    }
}

void SeedImage::computeMaxPossibleSeedCount() {
    if (seed_indices.size()==0) {
        maxPossibleSeedCount = 0;
    } else if ((mode==SEED_IMAGE_RS) && (max4rs==0))  {
        maxPossibleSeedCount = 0;
    } else {
        maxPossibleSeedCount = INT_MAX;
    }
}


// Seeder will random samples from voxel that have value 1. This image has to be binary, i.e. 0 and 1s only.
bool SeedImage::setSeed(Image<int8_t>* img) {

    auto minMax = imgMinMax(img);

    if ( (std::get<0>(minMax)<0) || (std::get<1>(minMax)>1) ) {
        disp(MSG_ERROR,"Seed image mask must be binary.");
        return false;
    }

    seed_indices = getIndices(img,1);
    disp(MSG_DEBUG,"Found %d non-zero voxels in %s.", seed_indices.size(), img->filePath.c_str());
    seed_img_mask.createFromTemplate(*img,false);

    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();

    mode          = SEED_IMAGE_MASK;
    return true;

}

// Seeder will random samples from voxel that have value label.
bool SeedImage::setSeed(Image<int>* img,int label) {

    seed_indices = getIndices(img,label);
    disp(MSG_DEBUG,"Found %d non-zero voxels in %s.", seed_indices.size(), img->filePath.c_str());
    seed_img_mask.createFromTemplate(*img,false);

    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();

    mode          = SEED_IMAGE_MASK;
    return true;

}

// Seeder will do rejection sampling
bool SeedImage::setSeed(Image<float>* img) {

    seed_img_pvf = img;

    auto minMax  = imgMinMax(seed_img_pvf);
            
    if (std::get<0>(minMax)<0) {
        disp(MSG_ERROR,"Seed image can't have negative values.");                
        return false;
    }
    
    max4rs          = std::get<1>(minMax);

    imgThresh(seed_img_mask, *img, 0, FLT_MAX);
    imgPad(seed_img_mask, 1);
    imgDilate(seed_img_mask);
    seed_indices    = getNonZeroIndices(&seed_img_mask);
    seed_img_mask.deallocData(); // We only need the indices and not the data

    disp(MSG_DEBUG,"Found %d non-zero voxels in %s. Max value is %.2f.", seed_indices.size(), img->filePath.c_str(), max4rs);

    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();
    
    mode = SEED_IMAGE_RS;
    return true;

}


// Seeder will do rejection sampling
bool SeedImage::setSeed(Image<float>* img, int _volInd) {

    seed_img_pvf = img;
    volInd       = _volInd;

    auto minMax  = imgMinMax(seed_img_pvf, volInd);
            
    if (std::get<0>(minMax)<0) {
        disp(MSG_ERROR,"Seed image can't have negative values.");                
        return false;
    }
    
    max4rs          = std::get<1>(minMax);

    Image<float> tmp;
    getImageSlice(&tmp,img,volInd);
    imgThresh(seed_img_mask, tmp, 0, FLT_MAX);
      
    imgPad(seed_img_mask, 1);
    imgDilate(seed_img_mask);
    seed_indices    = getNonZeroIndices(&seed_img_mask);
    seed_img_mask.deallocData(); // We only need the indices and not the data

    disp(MSG_DEBUG,"Found %d non-zero voxels in %s. Max value is %.2f.", seed_indices.size(), img->filePath.c_str(), max4rs);

    threadCount = 1;
    setNumberOfThreads(threadCount);
    computeMaxPossibleSeedCount();
    computeSeedCountAndDensity();
    
    mode = SEED_IMAGE_RS;
    return true;

}

void SeedImage::setNumberOfThreads(int n) {

    threadCount = n;

    if (doRandomThings!=NULL) {
        delete[] doRandomThings;
        doRandomThings = NULL;
    }
        
    doRandomThings = new RandomDoer[n];

    for (int i=0; i<n; i++) {
        doRandomThings[i].init_uniform_int(seed_indices.size()-1);
    }

}

