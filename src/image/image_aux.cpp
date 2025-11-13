#include "image.h"

using namespace NIBR;

template<typename T>
std::string NIBR::Image<T>::getSpaceUnit() {
    std::string unit;
    switch (spaceUnit) {
    case METER:         unit = "meter";        break;
    case MM:            unit = "mm";           break;
    case MICRON:        unit = "micron";       break;
    default:            unit = "unknown unit"; break;
    }
    return unit;
}

template<typename T>
void NIBR::Image<T>::printInfo() {

    if (NIBR::VERBOSE()<VERBOSE_INFO)
        return;

    disp(MSG_INFO,"Image info");

    std::cout << "\033[32m";
    std::cout << "File name:              "   << filePath << std::endl;
    std::cout << "File type:              "   << fileExtension << std::endl;
    std::cout << "Data description:       "   << description << std::endl;
    std::cout << "Number of dimensions:   "   << numberOfDimensions << std::endl;
    std::cout << "Dimensions:             [";
    for (int i=0; i<numberOfDimensions; i++) {
        std::cout << imgDims[i];
        if (i!=(numberOfDimensions-1))
            std::cout << " x ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Pixdim:                 [" << pixDims[0] << " x " << pixDims[1] << " x " << pixDims[2] << "]" << std::endl;
    std::cout << "Index order:            [" << indexOrder[0] << "," << indexOrder[1] << "," << indexOrder[2] << "," << indexOrder[3] << "," << indexOrder[4] << "," << indexOrder[5] << "," << indexOrder[6] << "]" << std::endl;
    // std::cout << "voxCnt, valCnt, numel:  [" << voxCnt << "," << valCnt << "," << numel << "]" << std::endl;
    // std::cout << "dataScaler, dataOffset: [" << dataScaler << "," << dataOffset << "]" << std::endl;

    std::cout << "Unit of pixdim:         " << getSpaceUnit() << std::endl;

    auto imgOr = getOrientation();
    std::cout << "Orientation:            " << imgOr[0] << imgOr[1] << imgOr[2] << std::endl;

    std::cout << "Datatype:               ";

    switch (inputDataType) {

        case BOOL_DT:          std::cout << "BOOL";       break;
        case UINT8_DT:         std::cout << "UINT8";      break;
        case INT8_DT:          std::cout << "INT8";       break;
        case UINT16_DT:        std::cout << "UINT16";     break;
        case INT16_DT:         std::cout << "INT16";      break;
        case UINT32_DT:        std::cout << "UINT32";     break;
        case INT32_DT:         std::cout << "INT32";      break;
        case UINT64_DT:        std::cout << "UINT64";     break;
        case INT64_DT:         std::cout << "INT64";      break;
        case FLOAT32_DT:       std::cout << "FLOAT32";    break;
        case FLOAT64_DT:       std::cout << "FLOAT64";    break;
        case FLOAT128_DT:      std::cout << "FLOAT128";   break;

        default:
            std::cout<<"Unknown data type";
        break;
    }

    std::cout << std::endl;

    std::cout << "Datatype (internal use):" << getTypeName(typeid(T)) << std::endl;

    std::cout << std::setprecision(4) << "ijk2xyz:                " << std::endl;
    std::cout << std::setprecision(4) << "                        " << ijk2xyz[0][0] << " " << ijk2xyz[0][1] << " " << ijk2xyz[0][2] << " " << ijk2xyz[0][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << ijk2xyz[1][0] << " " << ijk2xyz[1][1] << " " << ijk2xyz[1][2] << " " << ijk2xyz[1][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << ijk2xyz[2][0] << " " << ijk2xyz[2][1] << " " << ijk2xyz[2][2] << " " << ijk2xyz[2][3] << std::endl;

    std::cout << std::setprecision(4) << "xyz2ijk:                " << std::endl;
    std::cout << std::setprecision(4) << "                        " << xyz2ijk[0][0] << " " << xyz2ijk[0][1] << " " << xyz2ijk[0][2] << " " << xyz2ijk[0][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << xyz2ijk[1][0] << " " << xyz2ijk[1][1] << " " << xyz2ijk[1][2] << " " << xyz2ijk[1][3] << std::endl;
    std::cout << std::setprecision(4) << "                        " << xyz2ijk[2][0] << " " << xyz2ijk[2][1] << " " << xyz2ijk[2][2] << " " << xyz2ijk[2][3] << std::endl;

    std::cout << "\033[0m";

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