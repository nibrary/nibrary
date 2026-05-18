#include "dMRI/tractography/mappers/tractogram2imageMapper.h"
#include "dMRI/tractography/utility/parallelStreamlineGenerator.h"
#include <atomic>
#include <set>
#include <cfenv>

using namespace NIBR;

template<typename T>
NIBR::Tractogram2ImageMapper<T>::Tractogram2ImageMapper(NIBR::TractogramReader* _tractogram, NIBR::Image<T>* _img) {

    tractogram      = _tractogram;
    mask            = NULL;
    weightFile      = NULL;
    weightType      = NO_WEIGHT;
    mapOnce         = false;
    img             = _img;
    mutexGrid       = NULL;
    useMutexGrid    = true;

    img->readHeader();

    smoothing = std::make_tuple(0,0);

    maskFromImage = false;

}

template<typename T>
NIBR::Tractogram2ImageMapper<T>::Tractogram2ImageMapper::~Tractogram2ImageMapper() {

    if ((maskFromImage) && (mask!=NULL)) {
        for (int i = 0; i < img->imgDims[0]; i++) {
            for (int j = 0; j < img->imgDims[1]; j++) {
                delete[] mask[i][j];
            }
            delete[] mask[i];
        }
        delete[] mask;
    }

    if (weightFile!=NULL) {
        for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
            fclose(weightFile[t]);
        }
        delete[] weightFile;
    }

    if (mutexGrid!=NULL)
        delete[] mutexGrid;
    
    mutexMap.clear();

}

template<typename T>
bool NIBR::Tractogram2ImageMapper<T>::setMask(NIBR::Image<int>* maskImg) {
        
    maskImg->readHeader();

    for (int i=0; i<3; i++) {
        if (img->pixDims[i]!=maskImg->pixDims[i]) {
            disp(MSG_ERROR, "Mask and template do not have same voxel dimensions.");
            return false;
        }

        if (img->imgDims[i]!=maskImg->imgDims[i]) {
            disp(MSG_ERROR, "Mask and template do not have same image dimensions.");
            return false;
        }

        for (int j=0; j<4; j++)
            if (img->xyz2ijk[i][j]!=maskImg->xyz2ijk[i][j]) {
                disp(MSG_ERROR, "Mask and template images do not have the same coordinate space.");
                return false;
            }
    }

    maskImg->read();
    
    mask = new bool**[static_cast<uint32_t>(maskImg->imgDims[0])];
    for (int i = 0; i < maskImg->imgDims[0]; i++) {
        mask[i] = new bool*[static_cast<uint32_t>(maskImg->imgDims[1])];
        for (int j = 0; j < maskImg->imgDims[1]; j++) {
            mask[i][j] = new bool[static_cast<uint32_t>(maskImg->imgDims[2])];
            for (int k = 0; k < maskImg->imgDims[2]; k++) {
                if ((*maskImg)(i,j,k)>0) {
                    mask[i][j][k] = true;
                    if (!useMutexGrid) {
                        mutexMap.try_emplace(static_cast<uint32_t>(maskImg->sub2ind(i, j, k)));
                    }
                } else {
                    mask[i][j][k] = false;
                }                
            }
        }
    }

    maskFromImage = true;

    return true;

}


template<typename T>
bool NIBR::Tractogram2ImageMapper<T>::setMask(NIBR::Image<int>* maskImg, int selectedLabel) {
        
    maskImg->readHeader();

    for (int i=0; i<3; i++) {
        if (img->pixDims[i]!=maskImg->pixDims[i]) {
            disp(MSG_ERROR, "Mask and template do not have same voxel dimensions.");
            return false;
        }

        if (img->imgDims[i]!=maskImg->imgDims[i]) {
            disp(MSG_ERROR, "Mask and template do not have same image dimensions.");
            return false;
        }

        for (int j=0; j<4; j++)
            if (img->xyz2ijk[i][j]!=maskImg->xyz2ijk[i][j]) {
                disp(MSG_ERROR, "Mask and template images do not have the same coordinate space.");
                return false;
            }
    }

    maskImg->read();

    mask = new bool**[static_cast<uint32_t>(maskImg->imgDims[0])];
    for (int i = 0; i < maskImg->imgDims[0]; i++) {
        mask[i] = new bool*[static_cast<uint32_t>(maskImg->imgDims[1])];
        for (int j = 0; j < maskImg->imgDims[1]; j++) {
            mask[i][j] = new bool[static_cast<uint32_t>(maskImg->imgDims[2])];
            for (int k = 0; k < maskImg->imgDims[2]; k++) {
                if ((*maskImg)(i,j,k)==selectedLabel) {
                    mask[i][j][k] = true;
                    if (!useMutexGrid) {
                        mutexMap.try_emplace(static_cast<uint32_t>(maskImg->sub2ind(i, j, k)));
                    }
                } else {
                    mask[i][j][k] = false;
                }
            }
        }
    }

    maskFromImage = true;

    return true;

}


template<typename T>
void NIBR::Tractogram2ImageMapper<T>::setMask(bool*** _mask) {

    mask = _mask;

    if (!useMutexGrid) {
        for (int i = 0; i < img->imgDims[0]; i++) {
            for (int j = 0; j < img->imgDims[1]; j++) {
                for (int k = 0; k < img->imgDims[2]; k++) {
                    if (mask[i][j][k]) {
                        mutexMap.try_emplace(static_cast<uint32_t>(img->sub2ind(i, j, k)));
                    }
                }
            }
        }
    }

    // NIBR::disp(MSG_DETAIL,"mutexMap size: %d", mutexMap.size());

}

template<typename T>
void NIBR::Tractogram2ImageMapper<T>::setWeights(std::string _weightFile, WEIGHTTYPE _weightType) {

    CNumericLocaleGuard lockScopeForNumericReading;

    weightType = _weightType;

    // Initialize weights and make copies for multithreader
    weightFile = new FILE*[NIBR::MT::MAXNUMBEROFTHREADS()]();
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        weightFile[t] = fopen(_weightFile.c_str(), "rb");
    }

    tractogram->getNumberOfPoints(); // Precomputed the cumulative lengths for future use

}

template<typename T>
void NIBR::Tractogram2ImageMapper<T>::setWeights(std::vector<float> _weights, WEIGHTTYPE _weightType) {
    weightType = _weightType;
    weights    = _weights;

    tractogram->getNumberOfPoints(); // Precomputed the cumulative lengths for future use
}

template<typename T>
void NIBR::Tractogram2ImageMapper<T>::run(
        std::function<void(Tractogram2ImageMapper<T>* tim, int* _gridPos, NIBR::Segment& _seg)> processor_f,
        std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f
        ) 
{
    // NIBR::disp(MSG_DETAIL,"Starting complete run");

    if (useMutexGrid) {
        mutexGrid = new std::mutex[img->voxCnt];
    }

    tractogram->reset();

    NIBR::disp(MSG_DEBUG,"Calling fetchAndProcess");

    auto fetchAndProcess = [&](const NIBR::MT::TASK& task)->void {
        StreamlineBatch kernel(std::get<1>(smoothing)+1);
        
        auto [success,streamline,streamlineId] = tractogram->getNextStreamline();

        NIBR::disp(MSG_DEBUG,"Read %d",streamlineId);

        kernel[0] = std::move(streamline);

        processStreamline(kernel,streamlineId,task.threadId, processor_f);
    };

    // Process the tractogram and fill
    NIBR::MT::MTRUN(tractogram->numberOfStreamlines, "Tractogram to image mapping", fetchAndProcess);

    // Compile output
    // NIBR::disp(MSG_DETAIL,"Compiling output all");
    outputCompiler_f(this);

}

template<typename T>
void NIBR::Tractogram2ImageMapper<T>::run(
        std::function<void(Tractogram2ImageMapper<T>* tim, int* _gridPos, NIBR::Segment& _seg)> processor_f,
        std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f,
        int beginInd,
        int endInd
        ) 
{
    if (!tractogram->isPreloaded()) {
        NIBR::disp(MSG_DETAIL,"This function requires preloaded tractograms");
    }

    // NIBR::disp(MSG_DETAIL,"Starting partial run");

    if (useMutexGrid) {
        mutexGrid = new std::mutex[img->voxCnt];
    }

    tractogram->reset();

    // Process the tractogram and fill
    NIBR::MT::MTRUN(endInd-beginInd+1, "Tractogram part to image mapping", [&](const NIBR::MT::TASK& task)->void {

            StreamlineBatch kernel(std::get<1>(smoothing)+1);

            // NIBR::disp(MSG_DETAIL,"Reading streamline %d", int(task.no+beginInd));
            kernel[0] = tractogram->getStreamline(task.no+beginInd);

            // NIBR::disp(MSG_DETAIL,"Processing streamline %d", int(task.no+beginInd));
            processStreamline(kernel,task.no+beginInd,task.threadId, processor_f);

        });

    // Compile output
    // NIBR::disp(MSG_DETAIL,"Compiling output part");
    outputCompiler_f(this);
    
}

template<typename T1>
bool NIBR::Tractogram2ImageMapper<T1>::processStreamline(StreamlineBatch& kernel, int streamlineId, uint16_t threadNo, std::function<void(Tractogram2ImageMapper<T1>* tim, int* gridPos, NIBR::Segment& seg)> f) {

    // If streamline is empty, exit.
    NIBR::disp(MSG_DEBUG,"Starting processing");

    int len = kernel[0].size();

    NIBR::disp(MSG_DEBUG,"Len: %d",len);

    if (len==0) return true;
    
    double p0[3], p1[3], dir[3], length, lengthR, lengthScale, t;

    int32_t A[3], B[3];

    // Apply anisotropic smoothing
    if (std::get<1>(smoothing)>1) {
        const Streamline tmp = kernel[0];
        kernel = getParallelStreamlines(tmp, std::get<0>(smoothing), std::get<1>(smoothing), threadNo);
    }

    for (const auto& streamline : kernel) {
    
        NIBR::Segment seg;
        seg.streamlineNo = streamlineId;
        // NIBR::disp(MSG_DETAIL,"seg.streamlineNo: %d",seg.streamlineNo);

        if (weightType!=NO_WEIGHT) {

            // NIBR::disp(MSG_DETAIL,"Making segment data.");
            seg.data = new float;

            if (weights.empty()) {
                if (weightType==STREAMLINE_WEIGHT) {
                    fseek(weightFile[threadNo], sizeof(float)*streamlineId, SEEK_SET);
                    // NIBR::disp(MSG_DETAIL,"Reading weight.");
                    std::fread((float*)seg.data, sizeof(float), 1, weightFile[threadNo]);
                    // NIBR::disp(MSG_DETAIL,"Weight: %.2f", *(float*)seg.data);
                } else {
                    const auto& cumLen = tractogram->getNumberOfPoints();
                    fseek(weightFile[threadNo], sizeof(float)*cumLen[streamlineId], SEEK_SET);
                }
            } else {
                if (weightType==STREAMLINE_WEIGHT) {
                    *((float*)(seg.data)) = weights[streamlineId];
                }
            }

        }
        
        // End of segment in real and its corner in image space
        img->to_ijk(streamline[0].data(),p0);
        A[0]  = std::round(p0[0]);
        A[1]  = std::round(p0[1]);
        A[2]  = std::round(p0[2]);
        
        // If streamline has 1 point and no segment
        if (len==1) {

            if ( img->isInside(A) && ((mask==NULL) || mask[A[0]][A[1]][A[2]]) ) {

                seg.p[0]   = streamline[0][0];
                seg.p[1]   = streamline[0][1];
                seg.p[2]   = streamline[0][2];
                seg.dir[0] = 1.0f;
                seg.dir[1] = 0.0f;
                seg.dir[2] = 0.0f;
                seg.length = 0.0f;
                
                f(this, A, seg);
                
            }
            
            continue;
        }
        
        std::set<size_t> voxelInds;
        std::pair<std::set<size_t>::iterator,bool> voxelChecker;

        auto sub2ind = [&] () -> size_t {return A[0] + A[1]*img->imgDims[1] + A[2]*img->imgDims[1]*img->imgDims[2];};
        
        auto runFun = [&]() {

            if ( img->isInside(A) && ((mask==NULL) || mask[A[0]][A[1]][A[2]]) ) {

                // NIBR::disp(MSG_DETAIL,"seg.p:      [%.2f, %.2f, %.2f]",  seg.p[0],  seg.p[1],  seg.p[2]);
                // NIBR::disp(MSG_DETAIL,"seg.dir:    [%.2f, %.2f, %.2f]",seg.dir[0],seg.dir[1],seg.dir[2]);
                // NIBR::disp(MSG_DETAIL,"seg.length: %.2f",seg.length);

                if (mapOnce) {
                    voxelChecker = voxelInds.insert(sub2ind());
                    if (voxelChecker.second) {
                        f(this, A, seg);
                    }
                } else {
                    f(this, A, seg);
                }
                
            }

        };        
        
        // If streamline has many points and segments
        for (int i=0; i<len-1; i++) {

            // Read if the segment weight is given
            if (weightType==SEGMENT_WEIGHT) {
                if (weights.empty()) {
                    std::fread((float*)seg.data, sizeof(float), 1, weightFile[threadNo]);
                } else {
                    const auto& cumLen = tractogram->getNumberOfPoints();
                    *((float*)(seg.data)) = weights[cumLen[streamlineId]+i];
                }
            }

            // End of segment in real and its corner in image space
            img->to_ijk(streamline[i+1].data(),p1);
            for (int m=0;m<3;m++) {
                seg.p[m]   = streamline[i][m];
                seg.dir[m] = streamline[i+1][m] - streamline[i][m];
                B[m]       = std::round(p1[m]);
            }    
            
            // Find segment length and direction in real space
            lengthR = norm(seg.dir);
            for (int m=0;m<3;m++)
                seg.dir[m] /= lengthR;
            
            
            // Segment does not leave the voxel, add this segment and continue with the next one
            if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {

                seg.length = lengthR;
                runFun();
                
                for (int m=0;m<3;m++) {
                    p0[m] = p1[m];
                }
                
                continue;
                
            }
            
            // Find length and direction of segment in image space
            for (int m=0;m<3;m++)
                dir[m] = p1[m] - p0[m];
            length = norm(dir);
            for (int m=0;m<3;m++)
                dir[m] /= length;
            
            
            // Grid lengthScale
            lengthScale = std::sqrt(dir[0]*img->pixDims[0]*dir[0]*img->pixDims[0]+dir[1]*img->pixDims[1]*dir[1]*img->pixDims[1]+dir[2]*img->pixDims[2]*dir[2]*img->pixDims[2]);

            auto pushSegment = [&](float pushAmount)->bool {

                // Update segment length and run the function
                seg.length = pushAmount*lengthScale;
                runFun();

                length    -= pushAmount;               // Update remaining length in image space
                lengthR   -= pushAmount*lengthScale;   // Update remaining length in real space
                
                for (int m=0;m<3;m++) {     
                    p0[m]    += pushAmount*dir[m];               // Move the initial point of the segment in image space
                    A[m]      = std::round(p0[m]);               // Update the current voxel
                    seg.p[m] += seg.length*seg.dir[m];           // Move the initial point of the segment in real space
                }

                if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {

                    // Update segment length and run the function
                    seg.length = lengthR;
                    runFun();
                    
                    return true;
                }

                return false;

            };

            
            t = 0.0;

            while (length > 0.0) {
            
                if (rayTraceVoxel(A,p0,dir,t)) {
                    // i.e. if true, there is intersection, and the part until t is inside the current voxel
                    // Therefore, if the intersection happened, then t*dir[m] will still be inside the current voxel

                    // Cut t at length
                    t = std::min(t,length);

                    if(pushSegment(t))
                        break;

                }
                
                // Push the segment by EPS4 or by length. This does not introduce any errors. 
                // It is only manually moving the point instead of the ray-tracer since ray-tracing stops exactly at the voxel edge
                if (pushSegment(std::min(EPS4,length)))
                    break;
                            
            }

            for (int m=0;m<3;m++) {
                p0[m]    = p1[m];
                A[m]     = B[m];
            }
            

        }
        
        if (weightType!=NO_WEIGHT) {
            delete (float*)(seg.data);
        }        

    }
    
    return true;
}

// Explicit instantiations
template class NIBR::Tractogram2ImageMapper<bool>;
template class NIBR::Tractogram2ImageMapper<uint8_t>;
template class NIBR::Tractogram2ImageMapper<int8_t>;
template class NIBR::Tractogram2ImageMapper<uint16_t>;
template class NIBR::Tractogram2ImageMapper<int16_t>;
template class NIBR::Tractogram2ImageMapper<uint32_t>;
template class NIBR::Tractogram2ImageMapper<int32_t>;
template class NIBR::Tractogram2ImageMapper<uint64_t>;
template class NIBR::Tractogram2ImageMapper<int64_t>;
template class NIBR::Tractogram2ImageMapper<float>;
template class NIBR::Tractogram2ImageMapper<double>;
template class NIBR::Tractogram2ImageMapper<long double>;