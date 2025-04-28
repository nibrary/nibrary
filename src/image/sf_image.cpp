#include "sf_image.h"

using namespace NIBR;

NIBR::SF_Image::SF_Image() {
    sphericalDomainResolution = 0;
    hasDirections             = false;
    isEven                    = false;
}

NIBR::SF_Image::SF_Image(std::string _filePath, bool _isEven) {

    sphericalDomainResolution = 0;
    hasDirections             = false;
    isEven                    = _isEven;

    setFilePath(_filePath);
    readHeader();

    if (numberOfDimensions!=4) {
        disp(MSG_FATAL, "Spherical function image must have 4 dimensions");
    }

    if (isEven) {
        if (imgDims[3]==1)      sphericalDomainResolution = 1;
        if (imgDims[3]==6)      sphericalDomainResolution = 3;
        if (imgDims[3]==38)     sphericalDomainResolution = 5;
        if (imgDims[3]==90)     sphericalDomainResolution = 7;
        if (imgDims[3]==162)    sphericalDomainResolution = 9;
        if (imgDims[3]==230)    sphericalDomainResolution = 11;
        if (imgDims[3]==334)    sphericalDomainResolution = 13;
        if (imgDims[3]==542)    sphericalDomainResolution = 15;
        if (imgDims[3]==658)    sphericalDomainResolution = 17;
        if (imgDims[3]==818)    sphericalDomainResolution = 19;
        if (imgDims[3]==1038)   sphericalDomainResolution = 21;
        if (imgDims[3]==1278)   sphericalDomainResolution = 23;
        if (imgDims[3]==1482)   sphericalDomainResolution = 25;
        if (imgDims[3]==1730)   sphericalDomainResolution = 27;
        if (imgDims[3]==2154)   sphericalDomainResolution = 29;
        if (imgDims[3]==2398)   sphericalDomainResolution = 31;
    } else {
        if (imgDims[3]==1)      sphericalDomainResolution = 1;
        if (imgDims[3]==7)      sphericalDomainResolution = 3;
        if (imgDims[3]==56)     sphericalDomainResolution = 5;
        if (imgDims[3]==152)    sphericalDomainResolution = 7;
        if (imgDims[3]==284)    sphericalDomainResolution = 9;
        if (imgDims[3]==416)    sphericalDomainResolution = 11;
        if (imgDims[3]==608)    sphericalDomainResolution = 13;
        if (imgDims[3]==1004)   sphericalDomainResolution = 15;
        if (imgDims[3]==1232)   sphericalDomainResolution = 17;
        if (imgDims[3]==1544)   sphericalDomainResolution = 19;
        if (imgDims[3]==1976)   sphericalDomainResolution = 21;
        if (imgDims[3]==2444)   sphericalDomainResolution = 23;
        if (imgDims[3]==2840)   sphericalDomainResolution = 25;
        if (imgDims[3]==3320)   sphericalDomainResolution = 27;
        if (imgDims[3]==4148)   sphericalDomainResolution = 29;
        if (imgDims[3]==4640)   sphericalDomainResolution = 31;
    }

    SF::init(isEven,sphericalDomainResolution);
}

NIBR::SF_Image::SF_Image(std::string _filePath, std::string _dirPath, bool _isEven) {

    sphericalDomainResolution = 0;
    hasDirections             = true;
    isEven                    = _isEven;

    setFilePath(_filePath);
    readHeader();

    if (numberOfDimensions!=4) {
        disp(MSG_FATAL, "Spherical function image must have 4 dimensions");
    }

    auto tmp = NIBR::readTripletsFromTextFileTo2DVector(_dirPath);
    if (int(std::get<1>(tmp).size()) != int(this->imgDims[3])) {
        disp(MSG_ERROR,"Input image dimensions (%d) and number of directions (%d) do not match.",this->imgDims[3],std::get<1>(tmp).size());
        return;
    }

    dirs = std::get<1>(tmp);

    SF::init(dirs,isEven);

}

NIBR::SF_Image::~SF_Image() { }

bool NIBR::SF_Image::read() {
    
    if (Image<float>::read()==false) {
        disp(MSG_ERROR, "Can't read spherical function image.");
        return false;
    }

    return true;
}

void NIBR::SF_Image::smooth(float angle) {

    if (angle<=0)
        return;

    float dist = deg2rad(angle);
    
    // We will find and only process those voxels which have non-zero values
    std::vector<std::vector<int64_t>> nnzVoxelSubs;
    auto findNonZeroVoxels = [&](NIBR::MT::TASK task)->void {
        
        int64_t i   = task.no;
        
        for (int64_t j=0; j<imgDims[1]; j++)
            for (int64_t k=0; k<imgDims[2]; k++)
                for (int64_t t=0; t<imgDims[3]; t++) {

                    if (data[sub2ind(i,j,k,t)]!=0){
                        std::vector<int64_t> tmp{i,j,k};
                        NIBR::MT::PROC_MX().lock();
                        nnzVoxelSubs.push_back(tmp);
                        NIBR::MT::PROC_MX().unlock();
                        break;
                    }
                    
                }
        
    };
    NIBR::MT::MTRUN(imgDims[0],findNonZeroVoxels);

    float* sdata = new float[numel];
    for (int n=0; n<numel; n++)
        sdata[n] = 0;

    // Apply smoothing
    auto applySmoothing = [&](NIBR::MT::TASK task)->void {
        
        std::vector<int64_t> sub = nnzVoxelSubs[task.no];    
        
        float viCnt = 0;
        for (int t=0; t<valCnt; t++) {
            viCnt = 0;

            for (auto vi : SF::getSFNeighbors()[t]) {
                if (std::get<1>(vi)<dist) {
                    sdata[sub2ind(sub[0],sub[1],sub[2],t)] += data[sub2ind(sub[0],sub[1],sub[2],std::get<0>(vi))];
                    viCnt++;
                } else break; // distances are sorted in sfNeighbors
            }
            if (viCnt>0)
                sdata[sub2ind(sub[0],sub[1],sub[2],t)] /= viCnt;
            
        }
        
    };
    NIBR::MT::MTRUN(nnzVoxelSubs.size(),"Applying spherical smoothing",applySmoothing);

    delete[] data;
    data = sdata;

}
