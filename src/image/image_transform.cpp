#include "image_transform.h"
#include "fast-marching-method/fast_marching_method.hpp"

using namespace NIBR;

namespace fmm = thinks::fast_marching_method;

template <typename T>
std::pair<std::vector<std::array<int32_t, 3>>,std::vector<float>> img2fmmArrays(NIBR::Image<T>* img) 
{
    std::vector<std::array<int32_t, 3>> subs;

    int64_t i,j,k;

    for (int n=0; n<img->numel; n++) {
        if (img->data[n]==0) {
            img->ind2sub(n,i,j,k);
            subs.push_back({int32_t(i),int32_t(j),int32_t(k)});
        }
    }

    std::vector<float> dists(subs.size(),0.0f);
    return std::make_pair(subs,dists);
}

template std::pair<std::vector<std::array<int32_t, 3>>,std::vector<float>> img2fmmArrays<int> (NIBR::Image<int>*  img);  // explicit instantiation for int
template std::pair<std::vector<std::array<int32_t, 3>>,std::vector<float>> img2fmmArrays<bool>(NIBR::Image<bool>* img);  // explicit instantiation for bool

std::pair<std::vector<std::array<int32_t, 3>>,std::vector<float>> img2fmmArrays(NIBR::Image<float>* img) 
{
    std::vector<std::array<int32_t, 3>> subs;
    std::vector<float> dists;

    int64_t i,j,k;

    for (int n=0; n<img->numel; n++) {
        if (img->data[n]!=0) {
            img->ind2sub(n,i,j,k);
            subs.push_back({int32_t(i),int32_t(j),int32_t(k)});
            dists.push_back(img->data[n]);
        }
    }

    return std::make_pair(subs,dists);
}

std::vector<float> getArrivalTimes(
    std::pair<std::vector<std::array<int32_t, 3>>,std::vector<float>>& fmmArrays,
    int64_t* imgDims,
    float* pixDims)
{
    auto grid_size      = std::array<size_t, 3>{{size_t(imgDims[0]), size_t(imgDims[1]),size_t(imgDims[2])}};
    auto grid_spacing   = std::array<float,  3>{{       pixDims[0],         pixDims[1],        pixDims[2] }};
    auto uniform_speed  = pixDims[0];

    disp(MSG_DETAIL,"Computing EDT");

    auto out = fmm::SignedArrivalTime(
        grid_size,
        fmmArrays.first,
        fmmArrays.second,
        fmm::UniformSpeedEikonalSolver<float, 3>(grid_spacing, uniform_speed)
    );

    disp(MSG_DETAIL,"Done");

    return out;
}

template <typename T>
void wrapImgEDT(NIBR::Image<T>* inp, NIBR::Image<float>* out)
{
    auto fmmArrays     = img2fmmArrays(inp);
    auto arrival_times = getArrivalTimes(fmmArrays,inp->imgDims,inp->pixDims);

    out->createFromTemplate(*inp,true);

    for (int n=0; n<out->numel; n++) { 
        out->data[n] = -arrival_times[n]*inp->pixDims[0];
    }
}

template <typename T>
void wrapImgEDT(NIBR::Image<T>* inp)
{
    auto fmmArrays     = img2fmmArrays(inp);
    auto arrival_times = getArrivalTimes(fmmArrays,inp->imgDims,inp->pixDims);

    for (int n=0; n<inp->numel; n++) { 
        inp->data[n] = -arrival_times[n]*inp->pixDims[0];
    }
}

template void wrapImgEDT<int>  (NIBR::Image<int>*   inp,  NIBR::Image<float>* out);  // explicit instantiation for int
template void wrapImgEDT<bool> (NIBR::Image<bool>*  inp,  NIBR::Image<float>* out);  // explicit instantiation for bool
template void wrapImgEDT<float>(NIBR::Image<float>* inp,  NIBR::Image<float>* out);  // explicit instantiation for float
template void wrapImgEDT<float>(NIBR::Image<float>* inp);

void NIBR::imgEDT(NIBR::Image<bool>*  inp, NIBR::Image<float>* out) {wrapImgEDT(inp, out);}
void NIBR::imgEDT(NIBR::Image<int>*   inp, NIBR::Image<float>* out) {wrapImgEDT(inp, out);}
void NIBR::imgEDT(NIBR::Image<float>* inp, NIBR::Image<float>* out) {wrapImgEDT(inp, out);}
void NIBR::imgEDT(NIBR::Image<float>* inp) {wrapImgEDT(inp);}
