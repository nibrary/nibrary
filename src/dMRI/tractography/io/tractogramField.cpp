#include "tractogramField.h"
#include <type_traits>
#include <cctype> // For toupper

using namespace NIBR;

auto initFieldReader(TractogramReader& tractogram) {
    
    auto input = tractogram.file;
    const size_t strLength = 256;
	char dummy[strLength];

    std::fseek(input, tractogram.streamlinePos[tractogram.numberOfStreamlines-1], SEEK_SET);
    std::fseek(input, sizeof(float)*tractogram.len[tractogram.numberOfStreamlines-1]*3,SEEK_CUR);
    int tmp = std::fgetc(input);
    if (tmp != '\n') std::ungetc(tmp,input);
    std::fgets(dummy,strLength,input); // Skip the line about line number
    std::fseek(input, sizeof(int)*(tractogram.numberOfStreamlines+tractogram.numberOfPoints),SEEK_CUR);
    tmp = std::fgetc(input);
    if (tmp != '\n') std::ungetc(tmp,input);

    return input;
}

std::string toUpperCase(const std::string& input) {
    std::string result = input;
    for (char &c : result) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
    return result;
}

std::vector<NIBR::TractogramField> NIBR::findTractogramFields(TractogramReader& tractogram)
{
    std::vector<NIBR::TractogramField> fieldList;
    if (tractogram.fileFormat == VTK_ASCII) {
        disp(MSG_ERROR,"Can't read fields from ASCII files");
        return fieldList;
    }

    auto input = initFieldReader(tractogram);

    const size_t strLength = 256;
	char  dummy[strLength];
    char* name = new char[128];
    char* type = new char[128];
    int   dimension;
    bool  cellDataFound  = false;
    bool  pointDataFound = false;
    int   tmpi;
    
    while(feof(input) == 0) {
        
        std::fgets(dummy,strLength,input);
        
        if (std::string(dummy).find("CELL_DATA")!=std::string::npos) {
            disp(MSG_DEBUG,"Cell data found");
            cellDataFound   = true;
            pointDataFound  = false;
            continue;
        }
    
        if (std::string(dummy).find("POINT_DATA")!=std::string::npos)  {
            disp(MSG_DEBUG,"Point data found");
            cellDataFound  = false;
            pointDataFound = true;
            continue;
        }

        if (std::string(dummy).find("SCALARS")!=std::string::npos) {
            std::sscanf(dummy,"SCALARS %s %s %d\n", name, type, &dimension);
            std::fgets(dummy,strLength,input);

            disp(MSG_DEBUG,"Found field %s", name);

            if (cellDataFound==true) {
                TractogramField f = {STREAMLINE_OWNER,std::string(name),NIBR::datatype(toUpperCase(std::string(type))),dimension,NULL};
                fieldList.push_back(f);
                if (std::string(type)=="float") std::fseek(input,sizeof(float)*tractogram.numberOfStreamlines*dimension,SEEK_CUR);
                if (std::string(type)=="int")   std::fseek(input,sizeof(int)*tractogram.numberOfStreamlines*dimension,SEEK_CUR);
            }
            if (pointDataFound==true) {
                TractogramField f = {POINT_OWNER,std::string(name),NIBR::datatype(toUpperCase(std::string(type))),dimension,NULL};
                fieldList.push_back(f);
                if (std::string(type)=="float") std::fseek(input,sizeof(float)*tractogram.numberOfPoints*dimension,SEEK_CUR);
                if (std::string(type)=="int")   std::fseek(input,sizeof(int)*tractogram.numberOfPoints*dimension,SEEK_CUR);
            }
            tmpi = std::fgetc(input); if (tmpi != '\n') std::ungetc(tmpi,input); // Make sure to go end of the line
        }

    }

    delete[] name;
    delete[] type;

    disp(MSG_DEBUG,"Found %d fields", fieldList.size());
    
    return fieldList;

}

std::vector<NIBR::TractogramField> NIBR::readTractogramFields(TractogramReader& tractogram)
{

    std::vector<NIBR::TractogramField> fieldList;
    if (tractogram.fileFormat == VTK_ASCII) {
        disp(MSG_ERROR,"Can't read fields from ASCII files");
        return fieldList;
    }

    auto input = initFieldReader(tractogram);

    const size_t strLength = 256;
	char  dummy[strLength];
    char* name = new char[128];
    char* type = new char[128];
    int   dimension;
    bool  cellDataFound  = false;
    bool  pointDataFound = false;
    float tmpf;
    int   tmpi;
    
    while(feof(input) == 0) {    
        
        std::fgets(dummy,strLength,input);
        
        if (std::string(dummy).find("CELL_DATA")!=std::string::npos) {
            disp(MSG_DEBUG,"Found CELL_DATA");
            cellDataFound   = true;
            pointDataFound  = false;
            continue;
        }
    
        if (std::string(dummy).find("POINT_DATA")!=std::string::npos)  {
            disp(MSG_DEBUG,"Found POINT_DATA");
            cellDataFound  = false;
            pointDataFound = true;
            continue;
        }

        if (std::string(dummy).find("SCALARS")!=std::string::npos) {

            std::sscanf(dummy,"SCALARS %s %s %d\n", name, type, &dimension);
            std::fgets(dummy,strLength,input);
            
            disp(MSG_DEBUG,"Reading %s", name);
            
            if (cellDataFound==true) {

                TractogramField field;
                field.name = name;

                float** fdata = NULL;
                int**   idata = NULL;
                
                if (std::string(type)=="float") fdata = new float*[tractogram.numberOfStreamlines];
                if (std::string(type)=="int")   idata = new int*[tractogram.numberOfStreamlines];
                
                for (size_t n=0; n<tractogram.numberOfStreamlines; n++) {
                    
                    if (std::string(type)=="float") fdata[n] = new float[dimension];
                    if (std::string(type)=="int")   idata[n] = new int[dimension];
                    
                    for (int d=0; d<dimension; d++) {
                        if (std::string(type)=="float") {
                            std::fread(&tmpf,sizeof(float),1,input); 
                            swapByteOrder(tmpf);
                            fdata[n][d] = tmpf;
                        }
                        if (std::string(type)=="int") {
                            std::fread(&tmpi,sizeof(int),1,input); 
                            swapByteOrder(tmpi);
                            idata[n][d] = tmpi;
                        }
                    }
                }
                
                field.owner     = STREAMLINE_OWNER;
                field.dimension = dimension;
                if (std::string(type)=="float") {
                    field.datatype = FLOAT32_DT;
                    field.data     = reinterpret_cast<void*>(fdata);
                }
                if (std::string(type)=="int") {
                    field.datatype = INT32_DT;
                    field.data     = reinterpret_cast<void*>(idata);
                }

                fieldList.push_back(field);
                

            }

            if (pointDataFound==true) {

                TractogramField field;
                field.name = name;

                float*** fdata = NULL;
                int***   idata = NULL;
                
                if (std::string(type)=="float") fdata = new float**[tractogram.numberOfStreamlines];
                if (std::string(type)=="int")   idata = new   int**[tractogram.numberOfStreamlines];
                
                for (size_t s=0; s<tractogram.numberOfStreamlines; s++) {

                    if (std::string(type)=="float") fdata[s] = new float*[tractogram.len[s]];
                    if (std::string(type)=="int")   idata[s] = new   int*[tractogram.len[s]];

                    for (uint32_t l=0; l<tractogram.len[s]; l++) {
                    
                        if (std::string(type)=="float") fdata[s][l] = new float[dimension];
                        if (std::string(type)=="int")   idata[s][l] = new   int[dimension];
                        
                        for (int d=0; d<dimension; d++) {
                            if (std::string(type)=="float") {
                                std::fread(&tmpf,sizeof(float),1,input); 
                                swapByteOrder(tmpf);
                                fdata[s][l][d] = tmpf;
                            }
                            if (std::string(type)=="int") {
                                std::fread(&tmpi,sizeof(int),1,input); 
                                swapByteOrder(tmpi);
                                fdata[s][l][d] = tmpi;
                            }
                        }

                    }

                }
                
                field.owner     = POINT_OWNER;
                field.dimension = dimension;
                if (std::string(type)=="float") {
                    field.datatype = FLOAT32_DT;
                    field.data     = reinterpret_cast<void*>(fdata);
                }
                if (std::string(type)=="int") {
                    field.datatype = INT32_DT;
                    field.data     = reinterpret_cast<void*>(idata);
                }

                fieldList.push_back(field);
                
            }
            tmpi = std::fgetc(input); if (tmpi != '\n') std::ungetc(tmpi,input); // Make sure to go end of the line
            
        }
        
    }

    delete[] name;
    delete[] type;

    return fieldList;

}

TractogramField NIBR::readTractogramField(TractogramReader& tractogram,std::string fieldName) {
    
    TractogramField field;
    if (tractogram.fileFormat == VTK_ASCII) {
        disp(MSG_ERROR,"Can't read fields from ASCII files");
        return field;
    }
    
    field.name = fieldName;

    auto input = initFieldReader(tractogram);

    const size_t strLength = 256;
	char  dummy[strLength];
    char* name = new char[128];
    char* type = new char[128];
    int   dimension;
    bool  cellDataFound  = false;
    bool  pointDataFound = false;
    float tmpf;
    int   tmpi;
    
    while(feof(input) == 0) {    
        
        std::fgets(dummy,strLength,input);
        
        if (std::string(dummy).find("CELL_DATA")!=std::string::npos) {
            cellDataFound   = true;
            pointDataFound  = false;
            continue;
        }
    
        if (std::string(dummy).find("POINT_DATA")!=std::string::npos)  {
            cellDataFound  = false;
            pointDataFound = true;
            continue;
        }

        if (std::string(dummy).find("SCALARS")!=std::string::npos) {
            std::sscanf(dummy,"SCALARS %s %s %d\n", name, type, &dimension);
            std::fgets(dummy,strLength,input);
            
            if (name==fieldName) {
            
                if (cellDataFound==true) {

                    float** fdata = NULL;
                    int**   idata = NULL;
                    
                    if (std::string(type)=="float") fdata = new float*[tractogram.numberOfStreamlines];
                    if (std::string(type)=="int")   idata = new int*[tractogram.numberOfStreamlines];
                    
                    for (size_t n=0; n<tractogram.numberOfStreamlines; n++) {
                        
                        if (std::string(type)=="float") fdata[n] = new float[dimension];
                        if (std::string(type)=="int")   idata[n] = new int[dimension];
                        
                        for (int d=0; d<dimension; d++) {
                            if (std::string(type)=="float") {
                                std::fread(&tmpf,sizeof(float),1,input); 
                                swapByteOrder(tmpf);
                                fdata[n][d] = tmpf;
                            }
                            if (std::string(type)=="int") {
                                std::fread(&tmpi,sizeof(int),1,input); 
                                swapByteOrder(tmpi);
                                idata[n][d] = tmpi;
                            }
                        }
                    }
                    
                    field.owner     = STREAMLINE_OWNER;
                    field.dimension = dimension;
                    if (std::string(type)=="float") {
                        field.datatype = FLOAT32_DT;
                        field.data     = reinterpret_cast<void*>(fdata);
                    }
                    if (std::string(type)=="int") {
                        field.datatype = INT32_DT;
                        field.data     = reinterpret_cast<void*>(idata);
                    }
                    

                }

                if (pointDataFound==true) {

                    float*** fdata = NULL;
                    int***   idata = NULL;
                    
                    if (std::string(type)=="float") fdata = new float**[tractogram.numberOfStreamlines];
                    if (std::string(type)=="int")   idata = new   int**[tractogram.numberOfStreamlines];
                    
                    for (size_t s=0; s<tractogram.numberOfStreamlines; s++) {

                        if (std::string(type)=="float") fdata[s] = new float*[tractogram.len[s]];
                        if (std::string(type)=="int")   idata[s] = new   int*[tractogram.len[s]];

                        for (uint32_t l=0; l<tractogram.len[s]; l++) {
                        
                            if (std::string(type)=="float") fdata[s][l] = new float[dimension];
                            if (std::string(type)=="int")   idata[s][l] = new   int[dimension];
                            
                            for (int d=0; d<dimension; d++) {
                                if (std::string(type)=="float") {
                                    std::fread(&tmpf,sizeof(float),1,input); 
                                    swapByteOrder(tmpf);
                                    fdata[s][l][d] = tmpf;
                                }
                                if (std::string(type)=="int") {
                                    std::fread(&tmpi,sizeof(int),1,input); 
                                    swapByteOrder(tmpi);
                                    fdata[s][l][d] = tmpi;
                                }
                            }

                        }

                    }
                    
                    field.owner     = POINT_OWNER;
                    field.dimension = dimension;
                    if (std::string(type)=="float") {
                        field.datatype = FLOAT32_DT;
                        field.data     = reinterpret_cast<void*>(fdata);
                    }
                    if (std::string(type)=="int") {
                        field.datatype = INT32_DT;
                        field.data     = reinterpret_cast<void*>(idata);
                    }
                    
                }
                tmpi = std::fgetc(input); if (tmpi != '\n') std::ungetc(tmpi,input); // Make sure to go end of the line
                break;
            
            } else {
                
                if (cellDataFound==true) {
                    if (std::string(type)=="float") std::fseek(input,sizeof(float)*tractogram.numberOfStreamlines*dimension,SEEK_CUR);
                    if (std::string(type)=="int")   std::fseek(input,sizeof(int)*tractogram.numberOfStreamlines*dimension,SEEK_CUR);
                }
                if (pointDataFound==true) {
                    if (std::string(type)=="float") std::fseek(input,sizeof(float)*tractogram.numberOfPoints*dimension,SEEK_CUR);
                    if (std::string(type)=="int")   std::fseek(input,sizeof(int)*tractogram.numberOfPoints*dimension,SEEK_CUR);
                }
                tmpi = std::fgetc(input); if (tmpi != '\n') std::ungetc(tmpi,input); // Make sure to go end of the line
                
            }
            
            
        }
        
    }

    delete[] name;
    delete[] type;
    
    return field;
}

void NIBR::clearField(TractogramField& field,TractogramReader& tractogram)
{

    switch (field.datatype) {

        case BOOL_DT:          clearFieldWrapper<bool>(field,tractogram);    break;
        case UINT8_DT:         clearFieldWrapper<uint8_t>(field,tractogram);    break;
        case INT8_DT:          clearFieldWrapper< int8_t>(field,tractogram);    break;
        case UINT16_DT:        clearFieldWrapper<uint16_t>(field,tractogram);    break;
        case INT16_DT:         clearFieldWrapper< int16_t>(field,tractogram);    break;
        case UINT32_DT:        clearFieldWrapper<uint32_t>(field,tractogram);    break;
        case INT32_DT:         clearFieldWrapper< int32_t>(field,tractogram);    break;
        case UINT64_DT:        clearFieldWrapper<uint64_t>(field,tractogram);    break;
        case INT64_DT:         clearFieldWrapper< int64_t>(field,tractogram);    break;
        case FLOAT32_DT:       clearFieldWrapper<float>(field,tractogram);    break;
        case FLOAT64_DT:       clearFieldWrapper<double>(field,tractogram);    break;
        case FLOAT128_DT:      clearFieldWrapper<long double>(field,tractogram);    break;
        // case COMPLEX64_DT:     clearFieldWrapper<std::complex<double>>(field,tractogram);    break;
        // case COMPLEX128_DT:    clearFieldWrapper<std::complex<long double>>(field,tractogram);    break;
        // case COMPLEX256_DT:    clearFieldWrapper<std::complex<long long double>>(field,tractogram);    break;

        default:
            disp(MSG_FATAL,"Unknown datatype");
            break;
    }

    

}

TractogramField NIBR::makeTractogramFieldFromFile(TractogramReader& tractogram, std::string filePath, std::string name, std::string owner, std::string dataType, int dimension, bool isASCII) {

    // Create fields for vtk output
    TractogramField field;
    
    if (owner == "POINT") {
        field.owner = POINT_OWNER;
    } else if (owner == "STREAMLINE") {
        field.owner = STREAMLINE_OWNER;
    } else {
        disp(MSG_ERROR,"Unknown owner type. Owner can be either \"POINT\" or \"STREAMLINE\"");
        return field;
    }

    if (datatype(toUpperCase(dataType)) == UNKNOWN_DT) {
        disp(MSG_ERROR,"Unknown data type. Data type can be either \"float\" or \"int\"");
        return field;
    }

    field.name      = name;
    field.dimension = dimension;
    field.datatype  = datatype(toUpperCase(dataType));

    // Read field values
    FILE *input;
	input = fopen(filePath.c_str(),"r");

    auto readStreamlineData = [&](auto data, auto t) {
        for (size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
            data[s] = new decltype(t)[field.dimension];

            if (!isASCII) {
                std::fread(data[s], sizeof(decltype(t)), field.dimension, input);
            } else {
                for (int d = 0; d < field.dimension; d++) {
                    if (field.datatype == FLOAT32_DT) {
                        float tmp;
                        std::fscanf(input, "%f", &tmp);
                        data[s][d] = tmp;
                    } else if (field.datatype == INT32_DT) {
                        int tmp;
                        std::fscanf(input, "%i", &tmp);
                        data[s][d] = tmp;
                    }
                }
                std::fscanf(input, "\n");
            }
        }
    };


    auto readPointData = [&](auto data, auto t) {
        for (size_t s = 0; s < tractogram.numberOfStreamlines; s++) {
            data[s] = new decltype(&t)[tractogram.len[s]];

            for (uint32_t l = 0; l < tractogram.len[s]; l++) {
                data[s][l] = new decltype(t)[field.dimension];
                
                if (!isASCII) {
                    std::fread(data[s][l], sizeof(decltype(t)), field.dimension, input);
                } else {
                    for (int d = 0; d < field.dimension; d++) {
                        if (field.datatype == FLOAT32_DT) {
                            float tmp;
                            std::fscanf(input, "%f", &tmp);
                            data[s][l][d] = tmp;
                        } else if (field.datatype == INT32_DT) {
                            int tmp;
                            std::fscanf(input, "%i", &tmp);
                            data[s][l][d] = tmp;
                        }
                    }
                    std::fscanf(input, "\n");
                }
            }
        }
    };


    if (field.datatype == FLOAT32_DT) {

        if (field.owner == STREAMLINE_OWNER) {
            float** data = new float*[tractogram.numberOfStreamlines];
            readStreamlineData(data,float(0));
            field.data = (void*)data;
        }

        if (field.owner == POINT_OWNER) {
            float*** data = new float**[tractogram.numberOfStreamlines];
            readPointData(data,float(0));
            field.data = (void*)data;
        }

    }

    if (field.datatype == INT32_DT) {

        if (field.owner == STREAMLINE_OWNER) {
            int** data = new int*[tractogram.numberOfStreamlines];
            readStreamlineData(data,int(0));
            field.data = (void*)data;
        }

        if (field.owner == POINT_OWNER) {
            int*** data = new int**[tractogram.numberOfStreamlines];
            readPointData(data,int(0));
            field.data = (void*)data;
        }

    }
    
    fclose(input);


    disp(MSG_DEBUG,"field.name  %s", field.name.c_str());
    if (field.owner == POINT_OWNER)         disp(MSG_DEBUG,"field.owner POINT");
    if (field.owner == STREAMLINE_OWNER)    disp(MSG_DEBUG,"field.owner STREAMLINE");
    disp(MSG_DEBUG,"field.dimension %d", field.dimension);
    disp(MSG_DEBUG,"field.data %d", field.data);

    disp(MSG_DEBUG, "Read field");

    return field;
    
}