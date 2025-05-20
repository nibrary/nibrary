#pragma once

#include "base/nibr.h"
#include "math/core.h"
#include "dMRI/tractography/io/tractogramReader.h"
#include "image/image.h"
#include <atomic>
#include <map>
#include <mutex>

typedef enum {
    NO_WEIGHT,
    SEGMENT_WEIGHT,
    STREAMLINE_WEIGHT,
} WEIGHTTYPE;

namespace NIBR
{

    template<typename T>
    class Tractogram2ImageMapper {

    public:

        Tractogram2ImageMapper(NIBR::TractogramReader* _tractogram, NIBR::Image<T>* _img);
        ~Tractogram2ImageMapper();
        
        void run (
            std::function<void(Tractogram2ImageMapper<T>* tim, int* _gridPos, NIBR::Segment& _seg)> processor_f,
            std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f
        );

        void setMapOnce(bool val) {mapOnce=val;}
        bool setMask(NIBR::Image<int>* maskImg);
        bool setMask(NIBR::Image<int>* maskImg, int selectedLabel);
        void setMask(bool*** _mask);
        void setMaskVector(bool*** _mask);
        void anisotropicSmoothing(std::tuple<float,int> _smoothing) {smoothing = _smoothing;}
        void setWeights(std::string _weightFile, WEIGHTTYPE _weightType);
        void setWeights(std::vector<float> _weights, WEIGHTTYPE _weightType);
        void setData(void* _data) {data = _data;}
        void setBatchSize(std::size_t _batchSize) {batchSize = _batchSize;}
        void optimizeForSmallMask(bool val) {useMutexGrid = !val;}

        bool***                                 mask;
        std::vector<void*>                      grid;          // Grid holds a void* for each img->voxCnt

        bool                                    useMutexGrid;
        std::mutex*                             mutexGrid;     // Keeps track of with voxels are available to be processed in the grid
        std::unordered_map<uint32_t,std::mutex> mutexMap;
        std::mutex                              contribMutex;                       

        template<class GRIDTYPE>
        void allocateGrid();

        template<class GRIDTYPE>
        void deallocateGrid();

        // These are handled by the user
        NIBR::Image<T>*         img;
        void*                   data;

    private:

        void runAndDeleteStreamlines (
            std::function<void(Tractogram2ImageMapper<T>* tim, int* _gridPos, NIBR::Segment& _seg)> processor_f,
            std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f
        );

        void runAndKeepStreamlines (
            std::function<void(Tractogram2ImageMapper<T>* tim, int* _gridPos, NIBR::Segment& _seg)> processor_f,
            std::function<void(Tractogram2ImageMapper<T>* tim)> outputCompiler_f
        );

        bool processStreamline(
            std::vector<float**>& _kernel, 
            int _streamlineId, 
            uint16_t _threadNo, 
            std::function<void(Tractogram2ImageMapper<T>* tim, int* _gridPos, NIBR::Segment& _seg)> f, 
            bool deleteStreamline 
        );
        
        NIBR::TractogramReader* tractogram;

        bool                    mapOnce;
        std::tuple<float,int>   smoothing;

        std::vector<float>      weights;
        FILE**                  weightFile;
        WEIGHTTYPE              weightType;

        std::size_t             batchSize;

        int*                    cumLen;
        
        bool                    maskFromImage;
    };

    // Explicit instantiations
    extern template class Tractogram2ImageMapper<bool>;
    extern template class Tractogram2ImageMapper<uint8_t>;
    extern template class Tractogram2ImageMapper<int8_t>;
    extern template class Tractogram2ImageMapper<uint16_t>;
    extern template class Tractogram2ImageMapper<int16_t>;
    extern template class Tractogram2ImageMapper<uint32_t>;
    extern template class Tractogram2ImageMapper<int32_t>;
    extern template class Tractogram2ImageMapper<uint64_t>;
    extern template class Tractogram2ImageMapper<int64_t>;
    extern template class Tractogram2ImageMapper<float>;
    extern template class Tractogram2ImageMapper<double>;
    extern template class Tractogram2ImageMapper<long double>;

    // TODO: These are not yet supported because the Image<T> is missing their definitions
    // extern template class Tractogram2ImageMapper<std::complex<float>>;
    // extern template class Tractogram2ImageMapper<std::complex<double>>;
    // extern template class Tractogram2ImageMapper<std::complex<long double>>;

}


template<typename T> template<class GRIDTYPE>
void NIBR::Tractogram2ImageMapper<T>::allocateGrid() 
{
    deallocateGrid<GRIDTYPE>();
    grid.resize(img->voxCnt,NULL);
}

template<typename T> template<class GRIDTYPE>
void NIBR::Tractogram2ImageMapper<T>::deallocateGrid() 
{
    if (!grid.empty()) {
        auto resetGrid = [&](const NIBR::MT::TASK& task)->void{
            const int64_t& ind = task.no;
            if (grid[ind]!=NULL) {
                delete ((GRIDTYPE*)(grid[ind]));
                grid[ind] = NULL;
            }
        };
        NIBR::MT::MTRUN(grid.size(), resetGrid);
    }
    grid.clear();
    grid.shrink_to_fit();
}

