#include "fod_image.h"
#include "base/fileOperations.h"
#include "image/image_operators.h"
#include <cstdint>

using namespace NIBR;

NIBR::FOD_Image::FOD_Image() {
	iseven                      = true;
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
    iseven                      = true;
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
    iseven                      = true;
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
            std::vector<float> tmp = {p.x,p.y,p.z};
            sphereCoords.push_back(tmp);
        }

        if (int64_t(sphereCoords.size())!=this->valCnt) {
            disp(MSG_FATAL,"Number of spherical domain coordinates does not match number of volumes in FOD image");
            return false;
        }

    }
    
    // We need spherical harmonic basis functions in any case
    setSHorder();
    SH::precompute(shOrder, orderOfDirections, SHprecomputationResolution);
    
    if ((discretizationFlag==false) && (isspheresliced==false)) {
        fodAmp = &FOD_Image::ampWithDiscretizationOFF;
        return true;
    }
    
    
    // We will find and only process those voxels which have non-zero values
    std::vector<std::vector<int64_t>> nnzVoxelSubs = getNonZero3DVoxelSubs(this);
    
    // We will make a new data array and load the data there
    int64_t nnt;
    if (discretizationFlag) {
        fillDiscVolSph();
        nnt = discVolSphCoords.size();
    } else {
        nnt = getNumberOfSHCoeffs(shOrder);
    }
    
    std::vector<std::vector<float>> Ylm;

    if (isspheresliced) {
        SH_basis(Ylm,sphereCoords,shOrder);
    }

    float* ddata = new float[nnt*voxCnt];
    
    auto loadingTask = [&](NIBR::MT::TASK task)->void {

        int shNum = getNumberOfSHCoeffs(shOrder);
        
        std::vector<int64_t> sub = nnzVoxelSubs[task.no];
        float *FOD               = new float[shNum];
        
        if (isspheresliced==false) {
            for (int t=0; t<shNum; t++)
                FOD[t] = data[sub2ind(sub[0],sub[1],sub[2],t)];
        } else {
            for (int n=0; n<shNum; n++) {
                for (int64_t t=0; t<imgDims[3]; t++)
                    FOD[n] += Ylm[t][n]*data[sub2ind(sub[0],sub[1],sub[2],t)];
            }
        }
    
        
        if (discretizationFlag==true) {
            for (int64_t t=0; t<nnt; t++) {
                float dir[3]          = {discVolSphCoords[t][0],discVolSphCoords[t][1],discVolSphCoords[t][2]};
                ddata[t + sub[0]*nnt + sub[1]*nnt*imgDims[0] + sub[2]*nnt*imgDims[0]*imgDims[1]] = SH::toSF(FOD,dir);
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
        SH::clean();
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
    float* tmp = (float*) malloc(SH::getNumberOfSHCoeffs()*sizeof(float));
    (*this)(p,tmp);
    float val  = SH::toSF(tmp,tan);
    free(tmp);
    return val;
}
