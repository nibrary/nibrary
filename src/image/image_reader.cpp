#include "image.h"
#include "zlib.h"
#include <cstdint>
#include <ostream>
#include <vector>

using namespace NIBR;

template<typename T>
bool NIBR::Image<T>::read() {

    if (headerIsRead==false) {
        if (readHeader()==false)
            return false;
    }

    if ((fileExtension=="nii.gz") || (fileExtension=="nii"))
        return read_nii();

    if ((fileExtension=="mgh") || (fileExtension=="mgz"))
        return read_mghz();

    disp(MSG_ERROR,"Can't read image data with this extension yet: %s",fileExtension.c_str());
    return false;
}

template <typename OUT_T,typename INP_T>
void convert(OUT_T* out, void* inp, int64_t* imgDims, float dataScaler, float dataOffset, bool swapByte)
{

    // disp(MSG_DEBUG,"OUT_T: %s", typeid(OUT_T).name());
    // disp(MSG_DEBUG,"INP_T: %s", typeid(INP_T).name());
    // disp(MSG_DEBUG,"dataScaler: %.6f", dataScaler);
    // disp(MSG_DEBUG,"dataOffset: %.6f", dataOffset);
    // disp(MSG_DEBUG,"swapByte:   %d",   int(swapByte));

    int64_t numel = 1;
    for (auto i=0; i<7; i++) numel *= imgDims[i];

    for (auto i=0; i<numel; i++) {
        if (swapByte) swapByteOrder(*((INP_T*)inp+i));
        out[i] = *((INP_T*)inp+i);
    }

    bool scaleData = ( (dataScaler==0) || ((dataScaler==1) && (dataOffset==0)) ) ? false : true;
    if (scaleData) {
        for (auto i=0; i<numel; i++)
            out[i] = dataScaler*out[i] + dataOffset;
    }

}

template <typename OUT_T,typename INP_T>
void convert(OUT_T* out, void* inp, int64_t* imgDims, int* outIndexOrder, float dataScaler, float dataOffset, bool swapByte)
{

    // disp(MSG_DEBUG,"OUT_T: %s", typeid(OUT_T).name());
    // disp(MSG_DEBUG,"INP_T: %s", typeid(INP_T).name());
    // disp(MSG_DEBUG,"dataScaler: %.6f", dataScaler);
    // disp(MSG_DEBUG,"dataOffset: %.6f", dataOffset);
    // disp(MSG_DEBUG,"swapByte:   %d",   int(swapByte));

    int64_t outS2i[7];
    int64_t numel = 1;

    for (int i=0; i<7; i++) {
        outS2i[i]      = 1;
        numel         *= imgDims[i];
    }

    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++) {
            outS2i[outIndexOrder[i]] *= imgDims[outIndexOrder[j]];
        }

    
    int64_t ind = 0;
    for (auto o=0; o<imgDims[6]; o++)
        for (auto n=0; n<imgDims[5]; n++)
            for (auto m=0; m<imgDims[4]; m++)
                for (auto l=0; l<imgDims[3]; l++)
                    for (auto k=0; k<imgDims[2]; k++)
                        for (auto j=0; j<imgDims[1]; j++)
                            for (auto i=0; i<imgDims[0]; i++) {
                                if (swapByte) swapByteOrder(*((INP_T*)inp+ind));
                                out[i*outS2i[0]  + 
                                    j*outS2i[1]  + 
                                    k*outS2i[2]  + 
                                    l*outS2i[3]  + 
                                    m*outS2i[4]  + 
                                    n*outS2i[5]  + 
                                    o*outS2i[6]] = 
                                *((INP_T*)inp+ind++);
                            }



    bool scaleData = ( (dataScaler==0) || ((dataScaler==1) && (dataOffset==0)) ) ? false : true;
    if (scaleData) {
        for (int64_t i=0; i<numel; i++)
            out[i] = dataScaler*out[i] + dataOffset;
    }

}

template<typename T>
bool NIBR::Image<T>::read_nii() {

    nifti_image* nim = nifti_image_read(filePath.c_str(),0);

    if (nifti_image_load(nim)==-1) {
        disp(MSG_FATAL,"Cannot read nifti image: %s",filePath.c_str());
        return false;
    }


    data = new T[numel]();
    // data = (T*) malloc(numel*sizeof(T));

    bool reIdx = false;
    for (int i=1; i<7; i++)
        if (indexOrder[i] != i)
            reIdx = true;

    switch (inputDataType) {
 
        case BOOL_DT:          reIdx ? convert<T,bool>       (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,bool>       (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case UINT8_DT:         reIdx ? convert<T,uint8_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,uint8_t>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case INT8_DT:          reIdx ? convert<T,int8_t>     (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,int8_t>     (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case UINT16_DT:        reIdx ? convert<T,uint16_t>   (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,uint16_t>   (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case INT16_DT:         reIdx ? convert<T,int16_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,int16_t>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case UINT32_DT:        reIdx ? convert<T,uint32_t>   (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,uint32_t>   (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case INT32_DT:         reIdx ? convert<T,int32_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,int32_t>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case UINT64_DT:        reIdx ? convert<T,uint64_t>   (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,uint64_t>   (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case INT64_DT:         reIdx ? convert<T,int64_t>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,int64_t>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case FLOAT32_DT:       reIdx ? convert<T,float>      (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,float>      (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case FLOAT64_DT:       reIdx ? convert<T,double>     (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,double>     (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        case FLOAT128_DT:      reIdx ? convert<T,long double>(data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,long double>(data,nim->data,imgDims,dataScaler,dataOffset,false);    break;

        // TODO: Implement converters for complex data types
        // case COMPLEX64:     reIdx ? convert<T,std::complex<double>)>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,std::complex<double>>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        // case COMPLEX128:    reIdx ? convert<T,std::complex<long double>)>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,std::complex<long double>>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;
        // case COMPLEX256:    reIdx ? convert<T,std::complex<long long double>)>    (data,nim->data,imgDims,indexOrder,dataScaler,dataOffset,false) : convert<T,std::complex<long long double>>    (data,nim->data,imgDims,dataScaler,dataOffset,false);    break;

        default:
            disp(MSG_FATAL,"Can't read nifti file. Unknown datatype");
            break;
    }

    nifti_image_free(nim);

    return true;
}

template<typename T>
bool NIBR::Image<T>::read_mghz() {

    std::function<size_t(void*, size_t, size_t)> readFunc;

    gzFile gzfp = NULL;
    FILE*  fp   = NULL;

    if (fileExtension == "mgz") {
        gzfp = gzopen(filePath.c_str(), "rb");
        if (!gzfp) {
            disp(MSG_ERROR, "Cannot open .mgz file %s", filePath.c_str());
            return false;
        }
        readFunc = [gzfp](void* ptr, size_t size, size_t nmemb) -> size_t {auto byteCnt = gzread(gzfp, ptr, size * nmemb); return byteCnt/size;};
    } else { // Assuming "mgh" for simplicity
        fp = fopen(filePath.c_str(), "rb");
        if (!fp) {
            disp(MSG_ERROR, "Cannot open .mgh file %s", filePath.c_str());
            return false;
        }
        readFunc = [fp](void* ptr, size_t size, size_t nmemb) -> size_t {return fread(ptr, size, nmemb, fp);};
    }

    // Seek to the start of the data
    if (gzfp!=NULL) gzseek(gzfp, 284, SEEK_SET);
    if (fp  !=NULL) fseek (fp,   284, SEEK_SET);

    // Set true if reindexing needs to be done
    bool reIdx = false;
    for (int i=1; i<7; i++)
        if (indexOrder[i] != i)
            reIdx = true;

    // Read image data
    void* inpBuffer = NULL;

    auto run_reader = [&](auto t)->bool {

        inpBuffer = (void*)(new decltype(t)[numel]);

        if (inpBuffer==NULL) {
            disp(MSG_ERROR, "Failed to allocate memory for raw image data");
            if (gzfp!=NULL) gzclose(gzfp);
            if (fp!=NULL)   fclose(fp);
            return false;
        }

        if (readFunc(inpBuffer, sizeof(decltype(t)), numel) != size_t(numel)) {
            disp(MSG_ERROR, "Failed to read the correct amount of image data");
            delete[] reinterpret_cast<decltype(t)*>(inpBuffer);
            if (gzfp!=NULL) gzclose(gzfp);
            if (fp!=NULL)   fclose(fp);
            return false;
        }

        data = new T[numel]();
        // data = (T*) malloc(numel*sizeof(T));

        if (data==NULL) {
            disp(MSG_ERROR, "Failed to allocate memory for image");
            delete[] reinterpret_cast<decltype(t)*>(inpBuffer);
            if (gzfp!=NULL) gzclose(gzfp);
            if (fp!=NULL)   fclose(fp);
            return false;
        }

        reIdx ? convert<T,decltype(t)>    (data,inpBuffer,imgDims,indexOrder,dataScaler,dataOffset,true) : convert<T,decltype(t)>    (data,inpBuffer,imgDims,dataScaler,dataOffset,true);

        delete[] reinterpret_cast<decltype(t)*>(inpBuffer);
        if (gzfp!=NULL) gzclose(gzfp);
        if (fp!=NULL)   fclose(fp);
        return true;

    };

    
    // mgh/mgz files do not support all datatypes
    switch (inputDataType) {
 
        // case BOOL_DT:          { bool        tmp=0; if(!run_reader(tmp)) return false; break;}
        case UINT8_DT:         { uint8_t     tmp=0; if(!run_reader(tmp)) return false; break;}
        // case INT8_DT:          { int8_t      tmp=0; if(!run_reader(tmp)) return false; break;}
        // case UINT16_DT:        { uint16_t    tmp=0; if(!run_reader(tmp)) return false; break;}
        case INT16_DT:         { int16_t     tmp=0; if(!run_reader(tmp)) return false; break;}
        // case UINT32_DT:        { uint32_t    tmp=0; if(!run_reader(tmp)) return false; break;}
        case INT32_DT:         { int32_t     tmp=0; if(!run_reader(tmp)) return false; break;}
        // case UINT64_DT:        { uint64_t    tmp=0; if(!run_reader(tmp)) return false; break;}
        // case INT64_DT:         { int64_t     tmp=0; if(!run_reader(tmp)) return false; break;}
        case FLOAT32_DT:       { float       tmp=0; if(!run_reader(tmp)) return false; break;}
        // case FLOAT64_DT:       { double      tmp=0; if(!run_reader(tmp)) return false; break;}
        // case FLOAT128_DT:      { long double tmp=0; if(!run_reader(tmp)) return false; break;}

        // TODO: Implement converters for complex data types
        // case COMPLEX64:     { std::complex<double>             tmp; if(!run_reader(tmp)) return false; break;}
        // case COMPLEX128:    { std::complex<long double>        tmp; if(!run_reader(tmp)) return false; break;}
        // case COMPLEX256:    { std::complex<long long double>   tmp; if(!run_reader(tmp)) return false; break;}

        default:
            disp(MSG_FATAL,"Can't read file. Unknown datatype");
            break;
    }

    return true;

}

// Explicit instantiations
template class NIBR::Image<bool>;
template class NIBR::Image<uint8_t>;
template class NIBR::Image<int8_t>;
template class NIBR::Image<uint16_t>;
template class NIBR::Image<int16_t>;
template class NIBR::Image<uint32_t>;
template class NIBR::Image<int32_t>;
template class NIBR::Image<uint64_t>;
template class NIBR::Image<int64_t>;
template class NIBR::Image<float>;
template class NIBR::Image<double>;
template class NIBR::Image<long double>;