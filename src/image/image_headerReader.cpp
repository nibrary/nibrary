#include "image.h"
#include "zlib.h"

using namespace NIBR;

template<typename T>
bool NIBR::Image<T>::readHeader() {

    if (headerIsRead) {

        return true;

    } else {

        if ((fileExtension=="nii.gz") || (fileExtension=="nii"))
            return readHeader_nii();

        if ((fileExtension=="mgh") || (fileExtension=="mgz"))
            return readHeader_mghz();

        #if defined(HAS_DCM2NIIX)
        if ((fileExtension=="dcm") || (fileExtension=="")) {
            return readHeader_dcm();
        }
        #endif

        disp(MSG_ERROR,"Unknown file extension: %s", fileExtension.c_str());
        return false;

    }

}


template<typename T>
void NIBR::Image<T>::parseHeader() {

    voxCnt   = imgDims[0]*imgDims[1]*imgDims[2];
    valCnt   = imgDims[3]*imgDims[4]*imgDims[5]*imgDims[6];
    numel    = voxCnt*valCnt;

    if (valCnt<1) valCnt = 1;

    smallestPixDim = pixDims[0];
    smallestPixDim = (smallestPixDim < pixDims[1]) ? smallestPixDim : pixDims[1];
    smallestPixDim = (smallestPixDim < pixDims[2]) ? smallestPixDim : pixDims[2];

    // indexOrder and s2i are used for indexing the data in any desired order and for sub2ind and ind2sub conversion
    // These do not change the image dimensions.

    for (int i=0; i<7; i++)
        s2i[i] = 1;

    for (int i=1; i<7; i++)
        for (int j=0; j<i; j++)
            s2i[indexOrder[i]] *= imgDims[indexOrder[j]];

    headerIsRead   = true;

}

template<typename T>
bool NIBR::Image<T>::readHeader_nii_wrapper(nifti_image* nim) {

    if (nim==NULL) return false;

    description = std::string(nim->descrip);

    switch (nim->xyz_units) {
        case NIFTI_UNITS_METER:     spaceUnit = METER;              break;
        case NIFTI_UNITS_MM:        spaceUnit = MM;                 break;
        case NIFTI_UNITS_MICRON:    spaceUnit = MICRON;             break;
        default:                    spaceUnit = UNKNOWNSPACEUNIT;   break;
    }

    switch (nim->time_units) {
        case NIFTI_UNITS_SEC:       timeUnit = SEC;             break;
        case NIFTI_UNITS_MSEC:      timeUnit = MSEC;            break;
        case NIFTI_UNITS_USEC:      timeUnit = USEC;            break;
        case NIFTI_UNITS_HZ:        timeUnit = HZ;              break;
        case NIFTI_UNITS_PPM:       timeUnit = PPM;             break;
        case NIFTI_UNITS_RADS:      timeUnit = RADS;            break;
        default:                    timeUnit = UNKNOWNTIMEUNIT; break;
    }

    dataScaler = nim->scl_slope;
    dataOffset = nim->scl_inter;

    switch (nim->datatype) {

        case DT_UINT8:         inputDataType=UINT8_DT;      break;
        case DT_INT8:          inputDataType=INT8_DT;       break;
        case DT_UINT16:        inputDataType=UINT16_DT;     break;
        case DT_INT16:         inputDataType=INT16_DT;      break;
        case DT_UINT32:        inputDataType=UINT32_DT;     break;
        case DT_INT32:         inputDataType=INT32_DT;      break;
        case DT_UINT64:        inputDataType=UINT64_DT;     break;
        case DT_INT64:         inputDataType=INT64_DT;      break;
        case DT_FLOAT32:       inputDataType=FLOAT32_DT;    break;
        case DT_FLOAT64:       inputDataType=FLOAT64_DT;    break;
        case DT_FLOAT128:      inputDataType=FLOAT128_DT;   break;

        case DT_COMPLEX64:
            inputDataType=COMPLEX64_DT;
            disp(MSG_ERROR,"Nifti datatype: complex64 is not an accepted datatype");
            return false;

        case DT_COMPLEX128:
            inputDataType=COMPLEX128_DT;
            disp(MSG_ERROR,"Nifti datatype: complex128 is not an accepted datatype");
            return false;

        case DT_COMPLEX256:
            inputDataType=COMPLEX256_DT;
            disp(MSG_ERROR,"Nifti datatype: complex256 is not an accepted datatype");
            return false;

        case DT_BINARY:
            inputDataType=UNKNOWN_DT;
            disp(MSG_ERROR,"Nifti datatype: binary is not an accepted datatype");
            return false;

        case DT_RGB:
            inputDataType=UNKNOWN_DT;
            disp(MSG_ERROR,"Nifti datatype: rgb24 is not an accepted datatype");
            return false;

        case DT_RGBA32:
            inputDataType=UNKNOWN_DT;
            disp(MSG_ERROR,"Nifti datatype: rgba32 is not an accepted datatype");
            return false;

        case DT_ALL:
            inputDataType=UNKNOWN_DT;
            disp(MSG_ERROR,"Nifti datatype: all is not an accepted datatype");
            return false;

        default:
            inputDataType=UNKNOWN_DT; 
            disp(MSG_ERROR,"Nifti datatype: unknown or not applicable is not an accepted datatype");
            return false;
    }

    // Get dims and pixDims
    numberOfDimensions = nim->dim[0];

    if (numberOfDimensions>0) {
        for (int i=0; i<7; i++) {
            imgDims[i] = nim->dim[i+1];
            pixDims[i] = (imgDims[i]==0) ? 1 : nim->pixdim[i+1];
            imgDims[i] = (imgDims[i]==0) ? 1 : imgDims[i];
        }
    }

    // Use sform if possible otherwise use qform
    if (nim->sform_code>0) {
        for (int i=0; i<3; i++)
            for (int j=0; j<4; j++) {
                xyz2ijk[i][j] = nim->sto_ijk.m[i][j];
                ijk2xyz[i][j] = nim->sto_xyz.m[i][j];
            }
    }
    else {
        for (int i=0; i<3; i++)
            for (int j=0; j<4; j++) {
                xyz2ijk[i][j] = nim->qto_ijk.m[i][j];
                ijk2xyz[i][j] = nim->qto_xyz.m[i][j];
            }
    }

    return true;

}



template<typename T>
bool NIBR::Image<T>::readHeader_nii() {

    nifti_image* nim = nifti_image_read(filePath.c_str(),0);

    if(!readHeader_nii_wrapper(nim)) {
        nifti_image_free(nim);
        return false;
    }

    nifti_image_free(nim);

    parseHeader();

    return true;

}

template<typename T>
bool NIBR::Image<T>::readHeader_mghz() {  

    int   version, type, dof, dim;
    short goodRASFlag;
    float xr, xa, xs, yr, ya, ys, zr, za, zs, cr, ca, cs;

    std::function<void(void*, size_t, size_t)> readFunc;

    gzFile gzfp = NULL;
    FILE* fp    = NULL;

    if (fileExtension == "mgz") {
        gzfp = gzopen(filePath.c_str(), "rb");
        if (!gzfp) {
            disp(MSG_ERROR, "Cannot open .mgz file %s", filePath.c_str());
            return false;
        }
        readFunc = [gzfp](void* ptr, size_t size, size_t nmemb) {gzread(gzfp, ptr, size * nmemb);};
    } else { // Assuming "mgh" for simplicity
        fp = fopen(filePath.c_str(), "rb");
        if (!fp) {
            disp(MSG_ERROR, "Cannot open .mgh file %s", filePath.c_str());
            return false;
        }
        readFunc = [fp](void* ptr, size_t size, size_t nmemb) {fread(ptr, size, nmemb, fp);};
    }

    auto readAndSwap = [&readFunc, this](auto& value) {
        readFunc(&value, sizeof(value), 1);
        swapByteOrder(value);
    };

    readAndSwap(version);
    readAndSwap(dim); imgDims[0] = dim;
    readAndSwap(dim); imgDims[1] = dim;
    readAndSwap(dim); imgDims[2] = dim;
    readAndSwap(dim); imgDims[3] = dim;
    readAndSwap(type);
    readAndSwap(dof);
    readAndSwap(goodRASFlag);
    readAndSwap(pixDims[0]);
    readAndSwap(pixDims[1]);
    readAndSwap(pixDims[2]);
    readAndSwap(xr);
    readAndSwap(xa);
    readAndSwap(xs);
    readAndSwap(yr);
    readAndSwap(ya);
    readAndSwap(ys);
    readAndSwap(zr);
    readAndSwap(za);
    readAndSwap(zs);
    readAndSwap(cr);
    readAndSwap(ca);
    readAndSwap(cs);

    // std::cout << "Version    : " << version     << std::endl;
    // std::cout << "Width      : " << imgDims[0]  << std::endl;
    // std::cout << "Height     : " << imgDims[1]  << std::endl;
    // std::cout << "Depth      : " << imgDims[2]  << std::endl;
    // std::cout << "Volumes    : " << imgDims[3]  << std::endl;

    // std::cout << "type       : " << type        << std::endl;
    // std::cout << "dof        : " << dof         << std::endl;
    // std::cout << "goodRASFlag: " << goodRASFlag << std::endl;
    // std::cout << "xdim       : " << pixDims[0]  << std::endl;
    // std::cout << "ydim       : " << pixDims[1]  << std::endl;
    // std::cout << "zdim       : " << pixDims[2]  << std::endl;

    // std::cout << "xr         : " << xr          << std::endl;
    // std::cout << "xa         : " << xa          << std::endl;
    // std::cout << "xs         : " << xs          << std::endl;
    // std::cout << "yr         : " << yr          << std::endl;
    // std::cout << "ya         : " << ya          << std::endl;
    // std::cout << "ys         : " << ys          << std::endl;
    // std::cout << "zr         : " << zr          << std::endl;
    // std::cout << "za         : " << za          << std::endl;
    // std::cout << "zs         : " << zs          << std::endl;
    // std::cout << "cr         : " << cr          << std::endl;
    // std::cout << "ca         : " << ca          << std::endl;
    // std::cout << "cs         : " << cs          << std::endl;

    if (gzfp!=NULL) gzclose(gzfp);
    if (fp  !=NULL) fclose(fp);

    // Handle image dimensions
    numberOfDimensions = (imgDims[3] > 1) ? 4 : 3;
    if (imgDims[3] <= 1) imgDims[3] = 1;
    imgDims[4] = imgDims[5] = imgDims[6] = 1;


    // Handle data type: UCHAR, SHORT, INT, or FLOAT (specified as 0, 4, 1, or 3, respectively)
    switch (type) {
        case 0:  {inputDataType = UINT8_DT;     break;}
        case 4:  {inputDataType = INT16_DT;     break;}
        case 1:  {inputDataType = INT32_DT;     break;}
        case 3:  {inputDataType = FLOAT32_DT;   break;}
        default: {
            inputDataType = UNKNOWN_DT; 
            disp(MSG_FATAL, "Unknown .mgz file datatype");
            return false;
        }
    }

    // Handle goodRASFlag
    if (goodRASFlag!=1) {
        disp(MSG_WARN,".mgz file is not good for RAS mm conversion.");
        return false;
    }

    // Handle pixDims
    pixDims[3] = pixDims[4] = pixDims[5] = pixDims[6] = 1.0f;
    

    // Choose between sform or qform
    float ci    = imgDims[0]/2.0f;
    float cj    = imgDims[1]/2.0f;
    float ck    = imgDims[2]/2.0f;

    mat44 vox2ras;
    vox2ras.m[0][0] = pixDims[0]*xr;
    vox2ras.m[0][1] = pixDims[1]*yr;
    vox2ras.m[0][2] = pixDims[2]*zr;
    vox2ras.m[0][3] = cr - (vox2ras.m[0][0]*ci + vox2ras.m[0][1]*cj + vox2ras.m[0][2]*ck);
    vox2ras.m[1][0] = pixDims[0]*xa;
    vox2ras.m[1][1] = pixDims[1]*ya;
    vox2ras.m[1][2] = pixDims[2]*za;
    vox2ras.m[1][3] = ca - (vox2ras.m[1][0]*ci + vox2ras.m[1][1]*cj + vox2ras.m[1][2]*ck);
    vox2ras.m[2][0] = pixDims[0]*xs;
    vox2ras.m[2][1] = pixDims[1]*ys;
    vox2ras.m[2][2] = pixDims[2]*zs;
    vox2ras.m[2][3] = cs - (vox2ras.m[2][0]*ci + vox2ras.m[2][1]*cj + vox2ras.m[2][2]*ck);
    vox2ras.m[3][0] = 0;
    vox2ras.m[3][1] = 0;
    vox2ras.m[3][2] = 0;
    vox2ras.m[3][3] = 0;

    mat44 vox2rastkr;
    vox2rastkr.m[0][0] = pixDims[0]*xr;
    vox2rastkr.m[0][1] = pixDims[1]*yr;
    vox2rastkr.m[0][2] = pixDims[2]*zr;
    vox2rastkr.m[0][3] = -(vox2ras.m[0][0]*ci + vox2ras.m[0][1]*cj + vox2ras.m[0][2]*ck);
    vox2rastkr.m[1][0] = pixDims[0]*xa;
    vox2rastkr.m[1][1] = pixDims[1]*ya;
    vox2rastkr.m[1][2] = pixDims[2]*za;
    vox2rastkr.m[1][3] = -(vox2ras.m[1][0]*ci + vox2ras.m[1][1]*cj + vox2ras.m[1][2]*ck);
    vox2rastkr.m[2][0] = pixDims[0]*xs;
    vox2rastkr.m[2][1] = pixDims[1]*ys;
    vox2rastkr.m[2][2] = pixDims[2]*zs;
    vox2rastkr.m[2][3] = -(vox2ras.m[2][0]*ci + vox2ras.m[2][1]*cj + vox2ras.m[2][2]*ck);
    vox2rastkr.m[3][0] = 0;
    vox2rastkr.m[3][1] = 0;
    vox2rastkr.m[3][2] = 0;
    vox2rastkr.m[3][3] = 0;

    mat44 rastkr2vox;
    rastkr2vox = nifti_mat44_inverse(vox2rastkr);


    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            rastkr2ras[i][j] =
            vox2ras.m[i][0]*rastkr2vox.m[0][j] +
            vox2ras.m[i][1]*rastkr2vox.m[1][j] +
            vox2ras.m[i][2]*rastkr2vox.m[2][j] +
            vox2ras.m[i][3]*rastkr2vox.m[3][j];
        }
    }    

    // Make ijk2xyz and xyz2ijk
    vox2ras.m[3][3]   = 1;
    mat44 xyz2ijk_m44 = nifti_mat44_inverse(vox2ras);
    for (int i=0; i<3; i++) {
        for (int j=0; j<4; j++) {
            ijk2xyz[i][j] = vox2ras.m[i][j];
            xyz2ijk[i][j] = xyz2ijk_m44.m[i][j];
        }
    }
    
    parseHeader();

    return true;

}

#if defined(HAS_DCM2NIIX)
template<typename T>
bool NIBR::Image<T>::readHeader_dcm() {
    // dcm2niix wrapper does not allow to read the dcm header alone yet
    return true;
}
#endif

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
