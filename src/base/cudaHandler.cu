#include "config.h"

#ifdef HAVE_CUDA

#include "cudaHandler.cuh"

using namespace NIBR;

bool                NIBR::CUDAHANDLER::cuda_initialized = false;
cudaDeviceProp      NIBR::CUDAHANDLER::cuda_dev_prop;
int                 NIBR::CUDAHANDLER::cuda_maxThreadsPerBlock;
int                 NIBR::CUDAHANDLER::cuda_maxThreadsDim[3];
int                 NIBR::CUDAHANDLER::cuda_maxGridSize[3];

void NIBR::CUDAHANDLER::CUDA_PRINT_ERROR(cudaError_t cudaFunction)
{
    if (cudaFunction != cudaSuccess){
         disp(MSG_ERROR,"%s",cudaGetErrorString(cudaFunction));
    }
}

void NIBR::CUDAHANDLER::CUDA_INIT() 
{
    if (cuda_initialized==false) {
        cudaDeviceReset();
        int devNo;
        cudaGetDevice(&devNo);
        CUDA_PRINT_ERROR(cudaGetDeviceProperties (&cuda_dev_prop, devNo)); 

        cuda_maxThreadsPerBlock = cuda_dev_prop.maxThreadsPerBlock;
        cuda_maxThreadsDim[0]   = cuda_dev_prop.maxThreadsDim[0];
        cuda_maxThreadsDim[1]   = cuda_dev_prop.maxThreadsDim[1];
        cuda_maxThreadsDim[2]   = cuda_dev_prop.maxThreadsDim[2];

        cuda_maxGridSize[0]     = cuda_dev_prop.maxGridSize[0];
        cuda_maxGridSize[1]     = cuda_dev_prop.maxGridSize[1];
        cuda_maxGridSize[2]     = cuda_dev_prop.maxGridSize[2];

        cuda_initialized = true;
    }
}

void NIBR::CUDAHANDLER::CUDA_EXIT() 
{
    cudaDeviceReset();
}

void NIBR::CUDAHANDLER::CUDA_PRINT_INFO() 
{
    if (NIBR::VERBOSE()<VERBOSE_INFO)
        return;

    int devCount, devNo;
    cudaGetDeviceCount(&devCount);
    cudaGetDevice(&devNo);

    disp(MSG_INFO,"CUDA INFO");
    std::cout << "\033[32m";
    std::cout << "Device count:                 " << devCount << " (Using device #" << devNo << ")" <<std::endl;
    std::cout << "name (major,minor):           " << cuda_dev_prop.name << " (" << cuda_dev_prop.major << "." << cuda_dev_prop.minor << ")" << std::endl;
    std::cout << "computeMode:                  " << cuda_dev_prop.computeMode << std::endl;
    std::cout << "totalGlobalMem:               " << cuda_dev_prop.totalGlobalMem/1024/1024/1024 << " GB"<< std::endl;
    std::cout << "sharedMemPerBlock:            " << cuda_dev_prop.sharedMemPerBlock/1024 << " KB"<< std::endl;
    std::cout << "regsPerBlock:                 " << cuda_dev_prop.regsPerBlock << " 32-bit registers"<< std::endl;
    std::cout << "maxThreadsPerBlock:           " << cuda_maxThreadsPerBlock << std::endl;
    std::cout << "maxThreadsDim:                " << "[" << cuda_maxThreadsDim[0] << "," << cuda_maxThreadsDim[1] << "," << cuda_maxThreadsDim[2] << "]"<< std::endl;
    std::cout << "maxGridSize:                  " << "[" << cuda_maxGridSize[0] << "," << cuda_maxGridSize[1] << "," << cuda_maxGridSize[2] << "]"<< std::endl;
    std::cout << "multiProcessorCount:          " << cuda_dev_prop.multiProcessorCount << std::endl;
    std::cout << "maxThreadsPerMultiProcessor:  " << cuda_dev_prop.maxThreadsPerMultiProcessor << std::endl;
    std::cout << "\033[0m";
    

}

#endif