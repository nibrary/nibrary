#include "image.h"
#include "zlib.h"

using namespace NIBR;

template<typename T>
bool NIBR::Image<T>::write(std::string filePath_) {

    std::string ext = getFileExtension(filePath_);

    if ((ext=="nii") || (ext=="nii.gz"))
        return write_nii(filePath_);

    if ((ext=="mgh") || (ext=="mgz"))
        return write_mghz(filePath_);

    return true;

}

template<typename T>
nifti_1_header* NIBR::Image<T>::getNiftiHeader() {

    nifti_1_header* header;

    const int64_t dims[8] = {numberOfDimensions,imgDims[0],imgDims[1],imgDims[2],imgDims[3],imgDims[4],imgDims[5],imgDims[6]};

    switch (dataType) {
    case BOOL_DT:          header = nifti_make_new_n1_header(dims,DT_UINT8    ); break;
    case UINT8_DT:         header = nifti_make_new_n1_header(dims,DT_UINT8    ); break;
    case INT8_DT:          header = nifti_make_new_n1_header(dims,DT_INT8     ); break;
    case UINT16_DT:        header = nifti_make_new_n1_header(dims,DT_UINT16   ); break;
    case INT16_DT:         header = nifti_make_new_n1_header(dims,DT_INT16    ); break;
    case UINT32_DT:        header = nifti_make_new_n1_header(dims,DT_UINT32   ); break;
    case INT32_DT:         header = nifti_make_new_n1_header(dims,DT_INT32    ); break;
    case UINT64_DT:        header = nifti_make_new_n1_header(dims,DT_UINT64   ); break;
    case INT64_DT:         header = nifti_make_new_n1_header(dims,DT_INT64    ); break;
    case FLOAT32_DT:       header = nifti_make_new_n1_header(dims,DT_FLOAT32  ); break;
    case FLOAT64_DT:       header = nifti_make_new_n1_header(dims,DT_FLOAT64  ); break;
    case FLOAT128_DT:      header = nifti_make_new_n1_header(dims,DT_FLOAT128 ); break;
    default:   disp(MSG_ERROR,"Unknown output datatype."); return NULL; break;
    }

    header->scl_slope       = 1;
    header->scl_inter       = 0;


    header->qform_code      = 1;
    header->sform_code      = 1;


    nifti_dmat44 R;
    for (int i=0; i<3; i++)
        for (int j=0; j<4; j++) {
            R.m[i][j] = ijk2xyz[i][j];
        }
    R.m[3][0] = R.m[3][1] = R.m[3][2] = 0;
    R.m[3][3] = 1;


    double qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac; // not used
    nifti_dmat44_to_quatern(R,&qb,&qc,&qd,&qx,&qy,&qz,&dx,&dy,&dz,&qfac);

    header->quatern_b = qb;
    header->quatern_c = qc;
    header->quatern_d = qd;
    header->qoffset_x = qx;
    header->qoffset_y = qy;
    header->qoffset_z = qz;
    header->pixdim[0] = qfac;

    header->dim[0] = numberOfDimensions;
    for (int i=0; i<numberOfDimensions; i++) {
        header->pixdim[i+1] = pixDims[i];
        header->dim[i+1]    = imgDims[i];
    }

    header->srow_x[0] = R.m[0][0];
    header->srow_x[1] = R.m[0][1];
    header->srow_x[2] = R.m[0][2];
    header->srow_x[3] = R.m[0][3];
    header->srow_y[0] = R.m[1][0];
    header->srow_y[1] = R.m[1][1];
    header->srow_y[2] = R.m[1][2];
    header->srow_y[3] = R.m[1][3];
    header->srow_z[0] = R.m[2][0];
    header->srow_z[1] = R.m[2][1];
    header->srow_z[2] = R.m[2][2];
    header->srow_z[3] = R.m[2][3];

    switch (spaceUnit) {
    case METER:     header->xyzt_units = NIFTI_UNITS_METER;      break;
    case MM:        header->xyzt_units = NIFTI_UNITS_MM;         break;
    case MICRON:    header->xyzt_units = NIFTI_UNITS_MICRON;     break;
    default:        header->xyzt_units = NIFTI_UNITS_UNKNOWN;    break;
    }

    for (size_t i=0; i<SGNTR().length() && i<80; i++) {
        header->descrip[i] = SGNTR()[i];
    }

    return header;

}


template<typename T>
void* resetIndexingAndCopyData(T* inp, int64_t* imgDims, int* indexOrder) {

    int64_t s2i[7];

    int64_t numel = 1;

    for (int i=0; i<7; i++) {
        s2i[i] = 1;
        numel *= imgDims[i];
    }

    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            s2i[i] *= imgDims[j];

    // void* out = (T*) malloc(numel*sizeof(T));
    void* out = (void*)(new T[numel]());

    auto run = [&](NIBR::MT::TASK task) {
        int64_t sub[7];
        int64_t offset;
        int64_t ind = task.no;

        for (int i = 0; i < 7; i++) {
            offset             = ind % imgDims[indexOrder[i]];
            ind               -= offset;
            ind               /= imgDims[indexOrder[i]];
            sub[indexOrder[i]] = offset;
        }

        *((T*)(out)+sub[0]*s2i[0] + sub[1]*s2i[1] + sub[2]*s2i[2] + sub[3]*s2i[3] + sub[4]*s2i[4] + sub[5]*s2i[5] + sub[6]*s2i[6]) = inp[task.no];

    };
    NIBR::MT::MTRUN(numel,NIBR::MT::MAXNUMBEROFTHREADS(),run);

    return out;

}

template<typename T>
bool NIBR::Image<T>::write_nii(std::string filePath_) {

    nifti_1_header* header = getNiftiHeader();
    nifti_image* nim       = nifti_convert_n1hdr2nim(*header,filePath_.c_str());

    bool resetIndexing = false;

    for (int i=0; i<7; i++)
        if (indexOrder[i]!=i)
            resetIndexing = true;

    if (resetIndexing)
        nim->data = resetIndexingAndCopyData(data,imgDims,indexOrder);
    else
        nim->data = data;

    nifti_image_write(nim);

    if (resetIndexing) {
        nifti_image_free(nim);
    } else {
        if( nim->fname != NULL ) free(nim->fname) ;
        if( nim->iname != NULL ) free(nim->iname) ;
        (void)nifti_free_extensions( nim ) ;
        free(nim);
    }

    free(header);
    return true;

}

template<typename T>
bool NIBR::Image<T>::write_mghz(std::string filePath_) {

    // First handle the data buffer
    void*  buffer    = data;
    float* fltBuffer = NULL;

    // Reset indexing if needed
    bool resetIndexing = false;
    for (int i=0; i<7; i++)
        if (indexOrder[i]!=i)
            resetIndexing = true;
        
    if (resetIndexing) buffer = resetIndexingAndCopyData(data,imgDims,indexOrder);

    // Apply datatype conversion if needed
    int  type      = 0;
    switch (dataType) {

        case BOOL_DT:
        case UINT8_DT:   {type = 0; break;}
        case INT16_DT:   {type = 4; break;}
        case INT32_DT:   {type = 1; break;}
        case FLOAT32_DT: {type = 3; break;}
        default: {
            type = 3; 
            disp(MSG_WARN, "mgh/mgz format don't support this datatype. Output will be in float32.");
            fltBuffer = (float*)malloc(numel*sizeof(float));
            for (int i = 0; i < numel; i++) fltBuffer[i] = ((T*)(buffer))[i];
            break;
        }
    }

    


    // Handle the file
    bool isMGZ = (getFileExtension(filePath_) == "mgz");

    std::function<void(const void*, size_t, size_t)> writeFunc;

    FILE* fp    = NULL;
    gzFile gzfp = NULL;
    
    if (isMGZ) {
        gzfp = gzopen(filePath_.c_str(), "wb");
        if (!gzfp) {
            disp(MSG_ERROR, "Cannot open .mgz file for writing: %s", filePath_.c_str());
            return false;
        }
        writeFunc = [gzfp](const void* ptr, size_t size, size_t nmemb) {gzwrite(gzfp, ptr, size * nmemb);};
    } else {
        fp = fopen(filePath_.c_str(), "wb");
        if (!fp) {
            disp(MSG_ERROR, "Cannot open .mgh file for writing: %s", filePath_.c_str());
            return false;
        }
        writeFunc = [fp](const void* ptr, size_t size, size_t nmemb) {fwrite(ptr, size, nmemb, fp);};
    }

    int w = 0;

    // Define a lambda for writing and byte-swapping
    auto writeAndSwap = [&](const auto& value) {
        auto temp = value;
        swapByteOrder(temp);
        writeFunc(&temp, sizeof(temp), 1);
        w += sizeof(temp);
    };

    
    // Write the header
    writeAndSwap(int(1)); // version
    writeAndSwap(int(imgDims[0]));
    writeAndSwap(int(imgDims[1]));
    writeAndSwap(int(imgDims[2]));
    writeAndSwap(int(imgDims[3]));
    writeAndSwap(type);     // type
    writeAndSwap(int(1));   // dof
    writeAndSwap(short(1)); // goodRASFlag

    float ci    = imgDims[0]/2.0f;
    float cj    = imgDims[1]/2.0f;
    float ck    = imgDims[2]/2.0f;

    float xr = ijk2xyz[0][0] / pixDims[0];
    float yr = ijk2xyz[0][1] / pixDims[1];
    float zr = ijk2xyz[0][2] / pixDims[2];
    float cr = ijk2xyz[0][3] + (ijk2xyz[0][0]*ci + ijk2xyz[0][1]*cj + ijk2xyz[0][2]*ck);
    float xa = ijk2xyz[1][0] / pixDims[0];
    float ya = ijk2xyz[1][1] / pixDims[1];
    float za = ijk2xyz[1][2] / pixDims[2];
    float ca = ijk2xyz[1][3] + (ijk2xyz[1][0]*ci + ijk2xyz[1][1]*cj + ijk2xyz[1][2]*ck);
    float xs = ijk2xyz[2][0] / pixDims[0];
    float ys = ijk2xyz[2][1] / pixDims[1];
    float zs = ijk2xyz[2][2] / pixDims[2];
    float cs = ijk2xyz[2][3] + (ijk2xyz[2][0]*ci + ijk2xyz[2][1]*cj + ijk2xyz[2][2]*ck);

    writeAndSwap(pixDims[0]);
    writeAndSwap(pixDims[1]);
    writeAndSwap(pixDims[2]);
    writeAndSwap(xr);
    writeAndSwap(xa);
    writeAndSwap(xs);
    writeAndSwap(yr);
    writeAndSwap(ya);
    writeAndSwap(ys);
    writeAndSwap(zr);
    writeAndSwap(za);
    writeAndSwap(zs);
    writeAndSwap(cr);
    writeAndSwap(ca);
    writeAndSwap(cs);

    short filler[97];
    writeFunc(filler, sizeof(short), 97);    

    // Write the data
    if (fltBuffer!=NULL) {
        for (int i = 0; i < numel; i++) {
            writeAndSwap(fltBuffer[i]);
        }

        // writeFunc(fltBuffer, sizeof(float), numel); 
        free(fltBuffer);
    } else {
        for (int i = 0; i < numel; i++) {
            writeAndSwap(((T*)(buffer))[i]);
        }
        // writeFunc(buffer, sizeof(T), numel);
    }

    if (resetIndexing) delete[] reinterpret_cast<T*>(buffer);
    
    if (isMGZ) gzclose(gzfp);
    else fclose(fp);

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