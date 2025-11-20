#pragma once

#include "base/nibr.h"
#include "base/dataTypeHandler.h"
#include "dMRI/tractography/io/tractogramReader.h"

using namespace NIBR;

namespace NIBR
{
    typedef enum {
        STREAMLINE_OWNER,
        POINT_OWNER,
        OWNER_NOTSET
    } TractogramOwnerType;

    struct TractogramField {
        TractogramOwnerType     owner{OWNER_NOTSET};
        std::string             name{""};
        DATATYPE                datatype{UNKNOWN_DT};
        int                     dimension{0};
        void*                   data{NULL};
    };

    void clearField(TractogramField& field, TractogramReader& tractogram);
    void clearField(TractogramField& field, Tractogram& tractogram);

    std::vector<NIBR::TractogramField> findTractogramFields(TractogramReader& tractogram);
    std::vector<NIBR::TractogramField> readTractogramFields(TractogramReader& tractogram);

    TractogramField readTractogramField(TractogramReader& tractogram,std::string name);

    TractogramField makeTractogramFieldFromFile(TractogramReader& tractogram, std::string filePath, std::string name, std::string owner, std::string dataType, int dimension, bool isASCII);
    
    // Make tractogramfield from vector
    template<typename T>
    TractogramField makeTractogramFieldFromVector(TractogramReader& tractogram, std::string name, const std::vector<T>& dataVec);

}


template<class T>
void clearFieldWrapper(TractogramField& field,TractogramReader& tractogram) {    

    const auto& cumLen  = tractogram.getNumberOfPoints();

    if (field.data != NULL) {

        if (field.owner == POINT_OWNER) {

            T*** toDel = reinterpret_cast<T***>(field.data);

            for (std::size_t s = 0; s < tractogram.numberOfStreamlines; s++) {

                auto len = cumLen[s+1] - cumLen[s];

                for (uint32_t l = 0; l < len; l++) {
                    delete[] toDel[s][l];
                }
                delete[] toDel[s];
            }

            delete[] toDel; 

        }

        if (field.owner == STREAMLINE_OWNER) {

            T** toDel = reinterpret_cast<T**>(field.data);

            for (std::size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
                delete[] toDel[s];
            }
            
            delete[] toDel;

        }

    }

    field.owner     = POINT_OWNER;
    field.name      = "";
    field.datatype  = UNKNOWN_DT;
    field.dimension = 0;
    field.data      = NULL;

}

template<class T>
void clearFieldWrapper(TractogramField& field, Tractogram& tractogram) {    

    if (field.data != NULL) {

        if (field.owner == POINT_OWNER) {

            T*** toDel = reinterpret_cast<T***>(field.data);

            for (std::size_t s = 0; s < tractogram.size(); s++) {
                for (uint32_t l = 0; l < tractogram[s].size(); l++) {
                    delete[] toDel[s][l];
                }
                delete[] toDel[s];
            }

            delete[] toDel; 

        }

        if (field.owner == STREAMLINE_OWNER) {

            T** toDel = reinterpret_cast<T**>(field.data);

            for (std::size_t s = 0; s < tractogram.size(); s++) {
                delete[] toDel[s];
            }
            
            delete[] toDel;

        }

    }

    field.owner     = POINT_OWNER;
    field.name      = "";
    field.datatype  = UNKNOWN_DT;
    field.dimension = 0;
    field.data      = NULL;

}

template<typename T>
TractogramField NIBR::makeTractogramFieldFromVector(TractogramReader& tractogram, std::string name, const std::vector<T>& dataVec) {

    TractogramField field;
    field.name = name;

    // Check for empty data
    if (dataVec.empty()) return field;

    // Infer owner
    if (dataVec.size() == tractogram.numberOfStreamlines) {
        field.owner = STREAMLINE_OWNER;
    } else {
        const auto& cumLen = tractogram.getNumberOfPoints();
        std::size_t totalPoints = cumLen.back();
        if (dataVec.size() == totalPoints) {
            field.owner = POINT_OWNER;
        } else {
            disp(MSG_ERROR, "Data vector size does not match number of streamlines or number of points.");
            return field;
        }
    }

    // Infer dimension and type
    if constexpr (std::is_same<T, float>::value) {
        field.dimension = 1;
        field.datatype  = FLOAT32_DT;
    } else if constexpr (std::is_same<T, int>::value) {
        field.dimension = 1;
        field.datatype  = INT32_DT;
    } else if constexpr (std::is_same<T, std::vector<float>>::value) {
        field.dimension = static_cast<int>(dataVec[0].size());
        field.datatype  = FLOAT32_DT;
    } else if constexpr (std::is_same<T, std::vector<int>>::value) {
        field.dimension = static_cast<int>(dataVec[0].size());
        field.datatype  = INT32_DT;
    } else if constexpr (std::is_same<T, std::array<float,3>>::value) {
        field.dimension = 3;
        field.datatype  = FLOAT32_DT;
    } else if constexpr (std::is_same<T, std::array<int,3>>::value) {
        field.dimension = 3;
        field.datatype  = INT32_DT;
    } else {
        disp(MSG_ERROR, "Unsupported data type for tractogram field.");
        return field;
    }

    // Lambda to handle data access for both scalars and containers
    auto getDataVal = [&](size_t idx, int dim) {
        if constexpr (std::is_same<T, float>::value || std::is_same<T, int>::value) {
            return dataVec[idx];
        } else {
            return dataVec[idx][dim];
        }
    };

    // Allocate and Copy
    if (field.owner == POINT_OWNER) {
        // Structure: T*** (Streamline -> Point -> Dimension)
        const auto& cumLen = tractogram.getNumberOfPoints();
        
        if (field.datatype == FLOAT32_DT) {
            float*** data = new float**[tractogram.numberOfStreamlines];
            std::size_t globalPointIdx = 0;
            
            for (std::size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
                std::size_t nPoints = cumLen[s+1] - cumLen[s];
                data[s] = new float*[nPoints];
                for (std::size_t p = 0; p < nPoints; p++) {
                    data[s][p] = new float[field.dimension];
                    for(int d=0; d<field.dimension; d++) {
                        data[s][p][d] = getDataVal(globalPointIdx, d);
                    }
                    globalPointIdx++;
                }
            }
            field.data = reinterpret_cast<void*>(data);

        } else if (field.datatype == INT32_DT) {
            int*** data = new int**[tractogram.numberOfStreamlines];
            std::size_t globalPointIdx = 0;

            for (std::size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
                std::size_t nPoints = cumLen[s+1] - cumLen[s];
                data[s] = new int*[nPoints];
                for (std::size_t p = 0; p < nPoints; p++) {
                    data[s][p] = new int[field.dimension];
                    for(int d=0; d<field.dimension; d++) {
                        data[s][p][d] = getDataVal(globalPointIdx, d);
                    }
                    globalPointIdx++;
                }
            }
            field.data = reinterpret_cast<void*>(data);
        }

    } else { // STREAMLINE_OWNER
        // Structure: T** (Streamline -> Dimension)
        if (field.datatype == FLOAT32_DT) {
            float** data = new float*[tractogram.numberOfStreamlines];
            for (std::size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
                data[s] = new float[field.dimension];
                for(int d=0; d<field.dimension; d++) {
                    data[s][d] = getDataVal(s, d);
                }
            }
            field.data = reinterpret_cast<void*>(data);

        } else if (field.datatype == INT32_DT) {
            int** data = new int*[tractogram.numberOfStreamlines];
            for (std::size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
                data[s] = new int[field.dimension];
                for(int d=0; d<field.dimension; d++) {
                    data[s][d] = getDataVal(s, d);
                }
            }
            field.data = reinterpret_cast<void*>(data);
        }
    }

    return field;
}