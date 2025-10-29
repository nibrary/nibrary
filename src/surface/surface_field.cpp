#include "surface.h"
#include "surface_operators.h"
#include <vector>

using namespace NIBR;

std::vector<NIBR::SurfaceField> NIBR::Surface::findFields()
{

    CNumericLocaleGuard lockScopeForNumericReading();

	FILE *input;
	input = fopen(filePath.c_str(),"rb");

    std::vector<SurfaceField> fieldList;
    
	if (input==NULL) return fieldList;

	const size_t strLength = 256;
	char dummy[strLength];
    
	for (int i=0; i<4; i++) fgets(dummy,strLength,input);
    
    int numberOfPoints;
    std::fscanf(input,"%*s %d %*s\n",&numberOfPoints);
    
	if (numberOfPoints>0) {

        // Skip points
        std::fseek(input,sizeof(float)*numberOfPoints*3,SEEK_CUR);
        int tmp = std::fgetc(input);
        if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
        
        // Skip cells
        int numberOfFaces;
        std::fscanf(input,"%*s %d %*d\n",&numberOfFaces);
        std::fseek(input,sizeof(int)*(4*numberOfFaces),SEEK_CUR);
        tmp = std::fgetc(input);
        if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
        
        
        char* name = new char[128];
        char* type = new char[128];
        int   dimension;
        bool  cellDataFound  = false;
        bool  pointDataFound = false;
        
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
                if (cellDataFound==true) {
                    
                    SurfaceField f = {FACE,std::string(name),std::string(type),dimension,NULL,NULL};
                    fieldList.push_back(f);
                    if (std::string(type)=="float") std::fseek(input,sizeof(float)*numberOfFaces*dimension,SEEK_CUR);
                    if (std::string(type)=="int")   std::fseek(input,sizeof(int)*numberOfFaces*dimension,SEEK_CUR);
                }
                if (pointDataFound==true) {
                    SurfaceField f = {VERTEX,std::string(name),std::string(type),dimension,NULL,NULL};
                    fieldList.push_back(f);
                    if (std::string(type)=="float") std::fseek(input,sizeof(float)*numberOfPoints*dimension,SEEK_CUR);
                    if (std::string(type)=="int")   std::fseek(input,sizeof(int)*numberOfPoints*dimension,SEEK_CUR);
                }
                tmp = std::fgetc(input); if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
            }
            
            
        }
        
        delete[] name;
        delete[] type;
        
	}

	fclose(input);
	return fieldList;
}


bool NIBR::Surface::readFields() {

    CNumericLocaleGuard lockScopeForNumericReading();

	FILE *input;
	input = fopen(filePath.c_str(),"rb");
    
    if (input==NULL) return false;

	const size_t strLength = 256;
	char dummy[strLength];
    
	for (int i=0; i<4; i++) fgets(dummy,strLength,input);
    
    int numberOfPoints;
    std::fscanf(input,"%*s %d %*s\n",&numberOfPoints);
    
	if (numberOfPoints>0) {

        // Skip points
        std::fseek(input,sizeof(float)*numberOfPoints*3,SEEK_CUR);
        int tmp = std::fgetc(input);
        if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
        
        // Skip cells
        int numberOfFaces;
        std::fscanf(input,"%*s %d %*d\n",&numberOfFaces);
        std::fseek(input,sizeof(int)*(4*numberOfFaces),SEEK_CUR);
        tmp = std::fgetc(input);
        if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
        
        
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
                
                if (cellDataFound==true) {
                    
                    float** fdata = NULL;
                    int**   idata = NULL;
                    
                    if (std::string(type)=="float") fdata = new float*[nf];
                    if (std::string(type)=="int")   idata = new int*[nf];
                        
                    for (int n=0; n<nf; n++) {
                        
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
                    
                    SurfaceField f = {FACE,std::string(name),std::string(type),dimension,fdata,idata};
                    fields.push_back(f);

                }
                if (pointDataFound==true) {
                    
                    float** fdata = NULL;
                    int**   idata = NULL;
                    
                    if (std::string(type)=="float") fdata = new float*[nv];
                    if (std::string(type)=="int")   idata = new int*[nv];
                    
                    for (int n=0; n<nv; n++) {

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
                    
                    SurfaceField f = {VERTEX,std::string(name),std::string(type),dimension,fdata,idata};
                    fields.push_back(f);
                    
                }
                tmp = std::fgetc(input); if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
            }
            
            
        }
        
        delete[] name;
        delete[] type;
        
	}

	fclose(input);
    
    return true;
}

NIBR::SurfaceField NIBR::Surface::readField(std::string fieldName) {

    CNumericLocaleGuard lockScopeForNumericReading();

    float** fdata = NULL;
    int**   idata = NULL;
    
    SurfaceField field;
    field.name  = fieldName;
    field.fdata = fdata;
    field.idata = idata;
    
	FILE *input;
	input = fopen(filePath.c_str(),"rb");
    
    if (input==NULL) return field;

	const size_t strLength = 256;
	char dummy[strLength];
    
	for (int i=0; i<4; i++) fgets(dummy,strLength,input);
    
    int numberOfPoints;
    std::fscanf(input,"%*s %d %*s\n",&numberOfPoints);
    
	if (numberOfPoints>0) {

        // Skip points
        std::fseek(input,sizeof(float)*numberOfPoints*3,SEEK_CUR);
        int tmp = std::fgetc(input);
        if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
        
        // Skip cells
        int numberOfFaces;
        std::fscanf(input,"%*s %d %*d\n",&numberOfFaces);
        std::fseek(input,sizeof(int)*(4*numberOfFaces),SEEK_CUR);
        tmp = std::fgetc(input);
        if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
        
        
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
                        
                        if (std::string(type)=="float") fdata = new float*[nf];
                        if (std::string(type)=="int")   idata = new int*[nf];
                        
                        for (int n=0; n<nf; n++) {
                            
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
                        
                        field.fdata     = fdata;
                        field.idata     = idata;
                        field.owner     = FACE;
                        field.dimension = dimension;
                        field.datatype  = std::string(type);

                    }
                    if (pointDataFound==true) {
                        
                        if (std::string(type)=="float") fdata = new float*[nv];
                        if (std::string(type)=="int")   idata = new int*[nv];
                        
                        for (int n=0; n<nv; n++) {
                            
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
                        
                        field.fdata     = fdata;
                        field.idata     = idata;
                        field.owner     = VERTEX;
                        field.dimension = dimension;
                        field.datatype  = std::string(type);
                        
                    }
                    tmp = std::fgetc(input); if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
                    break;
                
                } else {
                    
                    if (cellDataFound==true) {
                        if (std::string(type)=="float") std::fseek(input,sizeof(float)*numberOfFaces*dimension,SEEK_CUR);
                        if (std::string(type)=="int")   std::fseek(input,sizeof(int)*numberOfFaces*dimension,SEEK_CUR);
                    }
                    if (pointDataFound==true) {
                        if (std::string(type)=="float") std::fseek(input,sizeof(float)*numberOfPoints*dimension,SEEK_CUR);
                        if (std::string(type)=="int")   std::fseek(input,sizeof(int)*numberOfPoints*dimension,SEEK_CUR);
                    }
                    tmp = std::fgetc(input); if (tmp != '\n') std::ungetc(tmp,input); // Make sure to go end of the line
                    
                }
                
                
            }
            
            
        }
        
        delete[] name;
        delete[] type;
        
	}

	fclose(input);

    if ((fdata == NULL) && (idata == NULL)) {
        disp(MSG_WARN, "Field not found: %s", fieldName.c_str());
    }

    return field;
}

void NIBR::Surface::deleteField(std::string _fieldName) {

    bool reCheck = true;

    while (reCheck) {

        reCheck = false;

        for (size_t i=0; i<fields.size(); i++) {
            if (fields[i].name == _fieldName ) {
                clearField(fields[i]);
                fields.erase(fields.begin()+i);
                reCheck = true;
                break;
            }
        }
        
    }

}


SurfaceField NIBR::Surface::copyField(std::string _fieldName) {

    for (size_t i=0; i<fields.size(); i++) {
        if (fields[i].name == _fieldName ) {
            return this->copyField(fields[i]);
        }
    }

    SurfaceField out;
    return out;

}


SurfaceField NIBR::Surface::copyField(SurfaceField& field) {
    
    SurfaceField out;

    if (field.owner == NOTSET)
        return out;

    out.name      = field.name;
    out.owner     = field.owner;
    out.datatype  = field.datatype;
    out.dimension = field.dimension;

    int N = (field.owner == VERTEX) ? nv : nf;

    if (out.datatype=="int") {
        out.idata = new int*[N];
        for (int n=0; n<N; n++) {
            out.idata[n] = new int[out.dimension];
            memcpy(out.idata[n], field.idata[n], out.dimension*sizeof(int));
        }
    }

    if (out.datatype=="float") {
        out.fdata = new float*[N];
        for (int n=0; n<N; n++) {
            out.fdata[n] = new float[out.dimension];
            memcpy(out.fdata[n], field.fdata[n], out.dimension*sizeof(float));
        }
    }

    return out;
    
}

void NIBR::Surface::clearField(SurfaceField& field) {

    if (field.owner==NOTSET)
        return;
    
    if (field.idata!=NULL) {
        if (field.owner == VERTEX) {
            for (int i=0; i<nv; i++) {
                delete[] field.idata[i];
            }
            delete[] field.idata;
        }
        
        if (field.owner == FACE) {
            for (int i=0; i<nf; i++) {
                delete[] field.idata[i];
            }
            delete[] field.idata;
        }
    }
    
    if (field.fdata!=NULL) {
        if (field.owner == VERTEX) {
            for (int i=0; i<nv; i++) {
                delete[] field.fdata[i];
            }
            delete[] field.fdata;
        }
        
        if (field.owner == FACE) {
            for (int i=0; i<nf; i++) {
                delete[] field.fdata[i];
            }
            delete[] field.fdata;
        }
    }
    
    field.owner = NOTSET;
    field.name.clear();
    field.datatype.clear();
    field.dimension = 0;
    field.idata = NULL;
    field.fdata = NULL;
    
}

void NIBR::Surface::convert2FaceField(SurfaceField& field) {
    SurfaceField out = NIBR::convert2FaceField(this, &field);
    this->clearField(field);
    field = this->copyField(out);
    this->clearField(out);
}

void NIBR::Surface::convert2VertField(SurfaceField& field) {
    SurfaceField out = NIBR::convert2VertField(this, &field);
    this->clearField(field);
    field = this->copyField(out);
    this->clearField(out);
}

void NIBR::Surface::convert2IntField(SurfaceField& field) {

    if (field.owner == NOTSET)  return;
    if (field.fdata == NULL)    return;

    int N = (field.owner == VERTEX) ? nv : nf;
    int D = field.dimension;

    field.idata = new int*[N];
    for (int n=0; n<N; n++) {
        field.idata[n] = new int[D];
        for (int d=0; d<D; d++) {
            field.idata[n][d] = field.fdata[n][d];
        }
        delete[] field.fdata[n];
    }
    delete[] field.fdata;
    field.fdata = NULL;

}


void NIBR::Surface::convert2FloatField(SurfaceField& field) {

    if (field.owner == NOTSET)  return;
    if (field.idata == NULL)    return;

    int N = (field.owner == VERTEX) ? nv : nf;
    int D = field.dimension;

    field.fdata = new float*[N];
    for (int n=0; n<N; n++) {
        field.fdata[n] = new float[D];
        for (int d=0; d<D; d++) {
            field.fdata[n][d] = field.idata[n][d];
        }
        delete[] field.idata[n];
    }
    delete[] field.idata;
    field.idata = NULL;
    
}

std::vector<float> NIBR::Surface::readFloatFieldData(SurfaceField& field,int d) {

    std::vector<float> out;

    if (field.owner == NOTSET)  return out;
    if (field.fdata == NULL)    return out;
    if (d>field.dimension)      return out;

    int N = (field.owner == VERTEX) ? nv : nf;      
    out.reserve(N);

    for (int n=0; n<N; n++)
        out.push_back(field.fdata[n][d]);

    return out;

}


std::vector<int> NIBR::Surface::readIntFieldData(SurfaceField& field,int d) {

    std::vector<int> out;

    if (field.owner == NOTSET)  return out;
    if (field.idata == NULL)    return out;
    if (d>field.dimension)      return out;

    int N = (field.owner == VERTEX) ? nv : nf;      
    out.reserve(N);

    for (int n=0; n<N; n++)
        out.push_back(field.idata[n][d]);

    return out;

}

NIBR::SurfaceField NIBR::Surface::makeVertField(const std::string& name, const std::vector<std::vector<int>>& idata) {
    
    SurfaceField f;
    f.name      = name;
    f.owner     = VERTEX;
    f.datatype  = "int";
    f.dimension = idata.size();
    f.fdata     = NULL;
    
    f.idata = new int*[nv];

    for (int n = 0; n < nv; n++) {
        f.idata[n] = new int[f.dimension];
        for (int d = 0; d < f.dimension; d++) {
            f.idata[n][d] = idata[d][n];
        }
    }

    return f;

}

NIBR::SurfaceField NIBR::Surface::makeVertField(const std::string& name, const std::vector<std::vector<float>>& fdata) {
    
    SurfaceField f;
    f.name      = name;
    f.owner     = VERTEX;
    f.datatype  = "float";
    f.dimension = fdata.size();
    f.idata     = NULL;
    
    f.fdata = new float*[nv];

    for (int n = 0; n < nv; n++) {
        f.fdata[n] = new float[f.dimension];
        for (int d = 0; d < f.dimension; d++) {
            f.fdata[n][d] = fdata[d][n];
        }
    }

    return f;

}

NIBR::SurfaceField NIBR::Surface::makeFaceField(const std::string& name, const std::vector<std::vector<int>>& idata) {
    
    SurfaceField f;
    f.name      = name;
    f.owner     = FACE;
    f.datatype  = "int";
    f.dimension = idata.size();
    f.fdata     = NULL;
    
    f.idata = new int*[nf];

    for (int n = 0; n < nf; n++) {
        f.idata[n] = new int[f.dimension];
        for (int d = 0; d < f.dimension; d++) {
            f.idata[n][d] = idata[d][n];
        }
    }

    return f;

}

NIBR::SurfaceField NIBR::Surface::makeFaceField(const std::string& name, const std::vector<std::vector<float>>& fdata) {
    
    SurfaceField f;
    f.name      = name;
    f.owner     = FACE;
    f.datatype  = "float";
    f.dimension = fdata.size();
    f.idata     = NULL;
    
    f.fdata = new float*[nf];

    for (int n = 0; n < nf; n++) {
        f.fdata[n] = new float[f.dimension];
        for (int d = 0; d < f.dimension; d++) {
            f.fdata[n][d] = fdata[d][n];
        }
    }

    return f;

}

NIBR::SurfaceField NIBR::Surface::makeVertField(const std::string& name, const std::vector<int>& data) {
    std::vector<std::vector<int>> mdata;
    mdata.push_back(data);
    return makeVertField(name,mdata);
}

NIBR::SurfaceField NIBR::Surface::makeVertField(const std::string& name, const std::vector<float>& data) {
    std::vector<std::vector<float>> mdata;
    mdata.push_back(data);
    return makeVertField(name,mdata);
}

NIBR::SurfaceField NIBR::Surface::makeFaceField(const std::string& name, const std::vector<int>& data) {
    std::vector<std::vector<int>> mdata;
    mdata.push_back(data);
    return makeFaceField(name,mdata);
}

NIBR::SurfaceField NIBR::Surface::makeFaceField(const std::string& name, const std::vector<float>& data) {
    std::vector<std::vector<float>> mdata;
    mdata.push_back(data);
    return makeFaceField(name,mdata);
}