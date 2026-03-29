#include "fod_image.h"
#include "base/fileOperations.h"
#include "image/image_operators.h"
#include <cstdint>

using namespace NIBR;

NIBR::FOD_Image::FOD_Image() {
    sh                          = NULL;
    isAsym                      = false;
    hasNoOddCoeffs              = true;
    discretizationFlag          = false;
    isspheresliced              = false;
    sphereFileName              = "";
    discVolSphInds              = NULL;
    SHprecomputationResolution  = 1024;
    orderOfDirections           = NIBR::XYZ;
    indexOrder[0] = 3;
    indexOrder[1] = 0;
    indexOrder[2] = 1;
    indexOrder[3] = 2;
    indexOrder[4] = 4;
    indexOrder[5] = 5;
    indexOrder[6] = 6;
}

NIBR::FOD_Image::FOD_Image(std::string _filePath) {
    sh                          = NULL;
    isAsym                      = false;
    hasNoOddCoeffs              = true;
    discretizationFlag          = false;
    isspheresliced              = false;
    sphereFileName              = "";
    discVolSphInds              = NULL;
    SHprecomputationResolution  = 1024;
    orderOfDirections           = NIBR::XYZ;
    indexOrder[0] = 3;
    indexOrder[1] = 0;
    indexOrder[2] = 1;
    indexOrder[3] = 2;
    indexOrder[4] = 4;
    indexOrder[5] = 5;
    indexOrder[6] = 6;
    setFilePath(_filePath);
    readHeader();
}

NIBR::FOD_Image::FOD_Image(std::string _filePath,std::string _sphereFileName) {
    sh                          = NULL;
    isAsym                      = false;
    hasNoOddCoeffs              = true;
    discretizationFlag          = false;
    isspheresliced              = true;
    sphereFileName              = _sphereFileName;
    discVolSphInds              = NULL;
    SHprecomputationResolution  = 1024;
    orderOfDirections           = NIBR::XYZ;
    indexOrder[0] = 3;
    indexOrder[1] = 0;
    indexOrder[2] = 1;
    indexOrder[3] = 2;
    indexOrder[4] = 4;
    indexOrder[5] = 5;
    indexOrder[6] = 6;
    setFilePath(_filePath);
    readHeader();
}



NIBR::FOD_Image::~FOD_Image() {
    if (discVolSphInds!=NULL) {
         delete[] discVolSphInds;
         discVolSphInds = NULL;
    }
    if (sh!=NULL) {
         delete sh;
         sh = NULL;
    }
}

bool NIBR::FOD_Image::read() {
    
    if (numberOfDimensions!=4) {
        disp(MSG_FATAL,"FOD image must have 4 dimensions");
        return false;
    }
    
    if (Image<float>::read()==false) {
        disp(MSG_FATAL,"Can't read FOD image");
        return false;
    }
    disp(MSG_DETAIL,"FOD image is read");
    if (NIBR::VERBOSE() == VERBOSE_DEBUG) {
        this->printInfo();
    }

    if (isspheresliced) {

        auto sphere = readTripletsFromTextFile(sphereFileName);

        if (!std::get<0>(sphere)) {
            disp(MSG_FATAL,"Can't read spherical domain coordinates");
            return false;
        }

        for (auto p : std::get<1>(sphere)) {
            sphereCoords.push_back(p);
        }

        if (int64_t(sphereCoords.size())!=this->valCnt) {
            disp(MSG_FATAL,"Number of spherical domain coordinates does not match number of volumes in FOD image");
            return false;
        }

    }
    
    // Initialize spherical harmonic features
    setSH();

    disp(MSG_DETAIL,"Number of spherical harmonic coefficients: %d", sh->getCoeffCount());
    disp(MSG_DETAIL,"Number of volumes in FOD image: %d", imgDims[3]);
    
    if ((discretizationFlag==false) && (isspheresliced==false)) {
        if (imgDims[3] == sh->getCoeffCount()) {
            fodAmp = &FOD_Image::ampWithDiscretizationOFF;
            disp(MSG_DETAIL,"Discretization is not used");
        } else {
            disp(MSG_FATAL,"Number of volumes in FOD image does not match number of spherical harmonic coefficients. Can't continue without discretization.");
            return false;   
        }
    }
    
    
    // We will find and only process those voxels which have non-zero values
    disp(MSG_DETAIL,"Counting non-zero voxels");
    std::vector<std::vector<int64_t>> nnzVoxelSubs = getNonZero3DVoxelSubs(this);
    disp(MSG_DETAIL,"Number of non-zero voxels: %d", nnzVoxelSubs.size());
    
    // We will make a new data array and load the data there
    int64_t nnt;
    if (discretizationFlag) {
        fillDiscVolSph();
        nnt = discVolSphCoords.size();
    } else {
        nnt = sh->getCoeffCount();
    }
    
    std::vector<std::vector<float>> Ylm;

    if (isspheresliced) {
        SH_basis(Ylm,sphereCoords,shOrder,hasNoOddCoeffs);
    }

    float* ddata = new float[nnt*voxCnt];
    int    shNum = sh->getCoeffCount();

    auto loadingTask = [&](const NIBR::MT::TASK& task)->void {
        
        // Get voxel subscripts
        std::vector<int64_t> sub = nnzVoxelSubs[task.no];

        // Get voxel values
        float *FOD = new float[shNum];
        
        if (isspheresliced==false) {
            for (int t=0; t<shNum; t++)
                FOD[t] = data[sub2ind(sub[0],sub[1],sub[2],t)];
        } else {
            for (int n=0; n<shNum; n++) {
                FOD[n] = 0.0f;
                for (int64_t t=0; t<imgDims[3]; t++)
                    FOD[n] += Ylm[t][n]*data[sub2ind(sub[0],sub[1],sub[2],t)];
            }
        }    
        
        if (discretizationFlag==true) {
            for (int64_t t=0; t<nnt; t++) {
                float dir[3] = {discVolSphCoords[t][0],discVolSphCoords[t][1],discVolSphCoords[t][2]};
                ddata[t + sub[0]*nnt + sub[1]*nnt*imgDims[0] + sub[2]*nnt*imgDims[0]*imgDims[1]] = sh->toSF(FOD,dir);
            }
        } else {
            for (int64_t t=0; t<nnt; t++) {
                ddata[t + sub[0]*nnt + sub[1]*nnt*imgDims[0] + sub[2]*nnt*imgDims[0]*imgDims[1]] = FOD[t];
            }
        }
 
        delete[] FOD;
        
    };
    NIBR::MT::MTRUN(nnzVoxelSubs.size(),NIBR::MT::MAXNUMBEROFTHREADS(),"Loading FOD",loadingTask);
    
    
    delete[] data;
    
    data            = ddata;
    imgDims[3]      = nnt;
    valCnt          = nnt;
    
    for (int i=0; i<7; i++) 
        s2i[i] = 1;
    
    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            s2i[indexOrder[i]] *= imgDims[indexOrder[j]];
        
    setInterpolationMethod(interpMethod);
    
    if (discretizationFlag) {
        discVolSphCoords.clear();
        delete sh; sh = NULL;
        fodAmp = &FOD_Image::ampWithDiscretizationON;
    } else {
        fodAmp = &FOD_Image::ampWithDiscretizationOFF;
    }

    return true;
}

float NIBR::FOD_Image::ampWithDiscretizationON(float *p, float* tan) {
    return (*this)(p,vertexCoord2volInd(tan));
}

float NIBR::FOD_Image::ampWithDiscretizationOFF(float *p, float* tan) {
    float* tmp = (float*) malloc(sh->getCoeffCount()*sizeof(float));
    (*this)(p,tmp);
    float val  = sh->toSF(tmp,tan);
    free(tmp);
    return val;
}
