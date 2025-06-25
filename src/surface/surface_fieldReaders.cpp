#include "surface.h"

using namespace NIBR;

std::tuple<NIBR::SurfaceField,std::vector<std::tuple<std::string,int,int,int,int>>> NIBR::Surface::readFreesurferAnnot(std::string filePath, std::string LUTfname) {

    std::string extension = filePath.substr(filePath.find_last_of(".") + 1);

    SurfaceField annotation;
    std::vector<std::tuple<std::string,int,int,int,int>> ctable;
    
    if (extension != "annot") {
        disp(MSG_ERROR,"Not a .annot file");
        return std::make_tuple(annotation,ctable);
    } 
    
    FILE *input;
	input = fopen(filePath.c_str(),"rb");
    
    int tmpi;
        
    std::fread(&tmpi,sizeof(int),1,input); 
    swapByteOrder(tmpi);
    int nv = tmpi;
    
    int** data = new int*[nv]; // annotation data
    for (int n=0; n<nv; n++) {
        data[n] = new int[1];
        std::fread(&tmpi,sizeof(int),1,input); 
        swapByteOrder(tmpi);
        int vertexId = tmpi;
        std::fread(&tmpi,sizeof(int),1,input); 
        swapByteOrder(tmpi);
        data[vertexId][0] = tmpi;
    }

    if(!feof(input)) { // There is a colortable
    
        std::fread(&tmpi,sizeof(int),1,input); // Not used
        swapByteOrder(tmpi);

        std::fread(&tmpi,sizeof(int),1,input); 
        swapByteOrder(tmpi);
        int num = tmpi;
        char* buffer;
        int R,G,B,code;

        if (num>0) { // Original version


            // TODO: The case for original version is not tested
            ctable.resize(num);

            std::fread(&tmpi,sizeof(int),1,input); 
            swapByteOrder(tmpi);
            int len = tmpi;
            buffer = new char[len];
            std::fread(buffer,sizeof(char),len,input); // Name of original colortable
            delete[] buffer;

            for (int i=0; i<num; i++) {

                std::fread(&tmpi,sizeof(int),1,input); 
                swapByteOrder(tmpi);
                len = tmpi;
                buffer = new char[len];
                std::fread(buffer,sizeof(char),len,input); // Name of anatomical structure

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi);
                R = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi);
                G = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi);
                B = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi); // Always 0

                code = R + G*256 + B*65536;

                ctable[i] = std::make_tuple(std::string(buffer),R,G,B,code);

                delete[] buffer;

            }

            std::sort( ctable.begin(), ctable.end(), [](auto t1, auto t2){return (std::get<4>(t1) < std::get<4>(t2)) ? 1 : 0;} );
            ctable.erase( std::unique( ctable.begin(), ctable.end(),[](auto t1, auto t2){return (std::get<4>(t1) == std::get<4>(t2)) ? 1 : 0;} ), ctable.end() );

        } else {     // Other version

            int version=-num;
            if (version!=2) {
                disp(MSG_ERROR,"Unknown .annot version: %d", version);
                return std::make_tuple(annotation,ctable);
            }

            std::fread(&tmpi,sizeof(int),1,input); 
            swapByteOrder(tmpi);
            num = tmpi;

            std::fread(&tmpi,sizeof(int),1,input); 
            swapByteOrder(tmpi);
            int len = tmpi;

            buffer = new char[len];
            std::fread(buffer,sizeof(char),len,input); // Name of original colortable
            delete[] buffer;

            std::fread(&tmpi,sizeof(int),1,input); 
            swapByteOrder(tmpi);
            num = tmpi; // number of entries to read
            ctable.resize(num);

            int structureNo;

            for (int i=0; i<num; i++) {

                std::fread(&tmpi,sizeof(int),1,input); 
                swapByteOrder(tmpi);
                structureNo = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); 
                swapByteOrder(tmpi);
                len = tmpi;

                buffer = new char[len];
                std::fread(buffer,sizeof(char),len,input); // Name of anatomical structure
                

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi);
                R = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi);
                G = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi);
                B = tmpi;

                std::fread(&tmpi,sizeof(int),1,input); swapByteOrder(tmpi); // Always 0

                code = R + G*256 + B*65536;

                ctable[structureNo] = std::make_tuple(std::string(buffer),R,G,B,code);

                delete[] buffer;

            }

            std::sort( ctable.begin(), ctable.end(), [](auto t1, auto t2){return (std::get<4>(t1) < std::get<4>(t2)) ? 1 : 0;} );

        }

    }

    fclose(input);

    if (!LUTfname.empty()) {

        FILE *fLUT;
	    fLUT = fopen(LUTfname.c_str(),"rb");

        std::vector<std::vector<int>> LUT;

        while (!feof(fLUT)) {
            std::vector<int> map;
            map.resize(2);
            fscanf(fLUT, "%d %d%*[^\n]\n", &map[0], &map[1]);
            LUT.push_back(map);
        }

        fclose(fLUT);

        for (int n=0; n<nv; n++) {
            for (auto map : LUT) {
                if (data[n][0] == map[0] ) {
                    data[n][0] = map[1];
                    break;
                }
            }
        }

    }


    annotation.owner     = VERTEX;
    annotation.name      = filePath;
    annotation.datatype  = "int";
    annotation.dimension = 1;
    annotation.fdata     = NULL;
    annotation.idata     = data;
    
    return std::make_tuple(annotation,ctable);
    
}

NIBR::SurfaceField NIBR::Surface::readFreesurferLabel(std::string filePath) {

    std::string extension = filePath.substr(filePath.find_last_of(".") + 1);

    SurfaceField label;
    
    if (extension != "label") {
        disp(MSG_ERROR,"Not a .label file");
        return label;
    } 
    
    FILE *input;
	input = fopen(filePath.c_str(),"rb");
    
    const int strLength = 128;
    char dummy[strLength];
    fgets(dummy,strLength,input); // Skip first line

    int N; // number of labeled vertices to read
    std::fscanf(input,"%d\n",&N);

    int** data = new int*[nv]; // label data
    for (int n=0; n<nv; n++) data[n]=new int();

    int index;  
    
    for (int n=0; n<N; n++) {
        std::fscanf(input,"%d  %*f  %*f  %*f %*f\n",&index);
        data[index][0] = 1;
    }
    
    fclose(input);
        
    label.owner     = VERTEX;
    label.name      = filePath;
    label.datatype  = "int";
    label.dimension = 1;
    label.fdata     = NULL;
    label.idata     = data;
    
    return label;
    
}

NIBR::SurfaceField NIBR::Surface::makeFieldFromFile(std::string filePath, std::string name, std::string owner, std::string dataType, int dimension, std::string LUTfname, bool isASCII) {

    // Create fields for vtk output
    SurfaceField field;
    field.name      = name;
    field.dimension = dimension;
    
    int N;
    
    if (owner == "VERTEX") {
        field.owner     = VERTEX;
        N               = nv;
    }
    else if (owner == "FACE") {
        field.owner     = FACE;
        N               = nf;
    }
    else {
        disp(MSG_ERROR,"Unknown owner type. Owner can be either \"VERTEX\" or \"FACE\"");
        return field;
    }
    
    if (dataType == "int")
        field.datatype = "int";
    else if (dataType == "float") 
        field.datatype = "float";
    else {
        disp(MSG_ERROR,"Unknown data type. Data type can be either \"float\" or \"int\"");
        return field;
    }
    
    // Read field values
    FILE *input;
	input = fopen(filePath.c_str(),"rb");    
    
    float** fdata = NULL;
    int**   idata = NULL;

    
    if (field.datatype == "float") {
        fdata = new float*[N];
        for (int n=0;n<N;n++) {
            fdata[n]     = new float[field.dimension];
            if (isASCII==false) {
                std::fread(fdata[n],sizeof(float),field.dimension,input);
            } else {
                for (int d=0; d<field.dimension; d++) {
                    std::fscanf(input,"%f",&fdata[n][d]);
                }
                std::fscanf(input,"\n");
            }
        }
    }
    if (field.datatype == "int") {
        idata = new int*[N];
        for (int n=0;n<N;n++) {
            idata[n]     = new int[field.dimension];
            if (isASCII==false) {
                std::fread(idata[n],sizeof(int),field.dimension,input);
            } else {
                for (int d=0; d<field.dimension; d++) {
                    std::fscanf(input,"%d",&idata[n][d]);
                }
                std::fscanf(input,"\n");
            }
        }
    }
    
    fclose(input);

    if (!LUTfname.empty()) {

        FILE *fLUT;
	    fLUT = fopen(LUTfname.c_str(),"rb");

        std::vector<std::vector<int>> LUT;

        while (!feof(fLUT)) {
            std::vector<int> map;
            map.resize(2);
            fscanf(fLUT, "%d %d%*[^\n]\n", &map[0], &map[1]);
            LUT.push_back(map);
        }

        fclose(fLUT);

        if (idata!=NULL) {
            for (int n=0; n<N; n++) {
                for (auto map : LUT) {
                    if (idata[n][0] == map[0] ) {
                        idata[n][0] = map[1];
                        break;
                    }
                }
            }
        }

        if (fdata!=NULL) {
            for (int n=0; n<N; n++) {
                for (auto map : LUT) {
                    if (fdata[n][0] == map[0] ) {
                        fdata[n][0] = map[1];
                        break;
                    }
                }
            }
        }

    }
    
    field.fdata = fdata;
    field.idata = idata;

    return field;
    
}