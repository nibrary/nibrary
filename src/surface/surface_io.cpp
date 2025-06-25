#include "surface.h"
#include <fstream>
#include <stdio.h>

using namespace NIBR;

// Header readers
bool NIBR::Surface::readHeader(std::string _filePath)
{

    filePath = _filePath;
    extension = getFileExtension(filePath);

    disp(MSG_DEBUG,"File extension: %s", extension.c_str());

    if (!existsFile(filePath)) {
        disp(MSG_ERROR,"File not found: %s", filePath.c_str());
        return false;
    }

    if (extension == "vtk") return readVTKMeshHeader();

    if (extension == "gii") return readGIIMeshHeader();
    
    if ((extension == "orig") || (extension == "pial") || (extension == "white") || (extension == "inflated") || (extension == "sphere") || (extension == "smoothwm")) {
        return readFreesurferMeshHeader();
    }

    disp(MSG_ERROR,"Unknown file format.");
    return false;
}

bool NIBR::Surface::readVTKMeshHeader()
{

    FILE *input;
    input = fopen(filePath.c_str(), "rb");

    const size_t strLength = 256;
    char dummy[strLength];
    char type[strLength];

    fgets(dummy, strLength, input); // version

    if ( (std::string(dummy).find("3.0") != std::string::npos) && 
         (std::string(dummy).find("4.2") != std::string::npos) ) {
            // Replace the newline character with a null terminator
        size_t len = strlen(dummy);
        if (len > 0 && dummy[len - 1] == '\n') dummy[len - 1] = '\0';
        disp(MSG_WARN,"Vtk version is not supported: %s", dummy);
    }

    disp(MSG_DEBUG,"Vtk version is %s", dummy);

    fgets(dummy, strLength, input); // Custom one line metadata
    fgets(type, strLength, input);  // ASCII or BINARY
    fgets(dummy, strLength, input); // Skip DATASET POLYDATA line
    std::fscanf(input, "%*s %d %*s\n", &nv);

    if (std::string(type) == "ASCII\n")
    {
        for (int n = 0; n < nv; n++)
        {
            fgets(dummy, strLength, input);
        }
    }
    else
    {
        std::fseek(input, sizeof(float) * nv * 3, SEEK_CUR);
    }

    std::fscanf(input, "%*s %d %*s\n", &nf);

    fclose(input);

    return true;
}

bool NIBR::Surface::readGIIMeshHeader()
{
    gifti_image *gifti = gifti_read_image(filePath.c_str(), 0);
    bool trianDataFound = false;
    bool pointDataFound = false;

    for (int i = 0; i < gifti->numDA; i++) {

        if (gifti->darray[i]->intent == NIFTI_INTENT_POINTSET)
        {
            nv = gifti->darray[i]->dims[0];
            pointDataFound = true;
            continue;
        }
        if (gifti->darray[i]->intent == NIFTI_INTENT_TRIANGLE)
        {
            nf = gifti->darray[i]->dims[0];
            trianDataFound = true;
            continue;
        }

    }
    delete gifti;

    if (!trianDataFound || !pointDataFound) {
        disp(MSG_ERROR,"XML file must include NIFTI_INTENT_POINTSET and NIFTI_INTENT_TRIANGLE intents.");
        return false;
    }

    return true;
}

bool NIBR::Surface::readFreesurferMeshHeader()
{

    FILE *input;
    input = fopen(filePath.c_str(), "rb");

    // Check magic number
    unsigned char m1, m2, m3;
    std::fread(&m1, sizeof(unsigned char), 1, input);
    std::fread(&m2, sizeof(unsigned char), 1, input);
    std::fread(&m3, sizeof(unsigned char), 1, input);
    int magic = 65536 * int(m1) + 256 * int(m2) + int(m3);

    if (magic != 16777214)
    {
        disp(MSG_ERROR,"%s is not a valid Freesurfer surface file.");
        return false;
    }

    const size_t strLength = 256;
    char dummy[strLength];
    for (int i = 0; i < 2; i++)
        fgets(dummy, strLength, input);

    std::fread(&nv, sizeof(int), 1, input);
    swapByteOrder(nv);
    std::fread(&nf, 4, 1, input);
    swapByteOrder(nf);

    fclose(input);

    return true;
}

// Mesh readers
bool NIBR::Surface::readMesh() {

    if (extension == "vtk") return readVTKMesh();

    if (extension == "gii") return readGIIMesh();
    
    if ((extension == "orig") ||(extension == "pial") || (extension == "white") || (extension == "inflated") || (extension == "sphere") || (extension == "smoothwm")) {
        return readFreesurferMesh();
    }
    
    disp(MSG_ERROR,"Unknown file format.");
    return false;
}

bool NIBR::Surface::readVTKMesh()
{

    FILE *input;
    input = fopen(filePath.c_str(), "rb");

    const size_t strLength = 256;
    char dummy[strLength];
    char type[strLength];

    fgets(dummy, strLength, input); // Skip "vtk DataFile Version 3.0" line
    fgets(dummy, strLength, input); // Custom one line metadata
    fgets(type, strLength, input);  // ASCII or BINARY
    fgets(dummy, strLength, input); // Skip DATASET POLYDATA line
    std::fscanf(input, "%*s %d %*s\n", &nv);

    vertices = new float *[nv];
    if (std::string(type) == "ASCII\n")
    {
        for (int n = 0; n < nv; n++)
        {
            vertices[n] = new float[3];
            std::fscanf(input, "%f %f %f\n", &vertices[n][0], &vertices[n][1], &vertices[n][2]);
        }
    }
    else
    {
        float tmpf;
        for (int n = 0; n < nv; n++)
        {
            vertices[n] = new float[3];
            for (int i = 0; i < 3; i++)
            {
                std::fread(&tmpf, sizeof(float), 1, input);
                swapByteOrder(tmpf);
                vertices[n][i] = tmpf;
            }
        }
    }

    std::fscanf(input, "%*s %d %*s\n", &nf);

    faces = new int *[nf];
    if (std::string(type) == "ASCII\n")
    {
        for (int n = 0; n < nf; n++)
        {
            faces[n] = new int[3];
            std::fscanf(input, "%*d %d %d %d\n", &faces[n][0], &faces[n][1], &faces[n][2]);
        }
    }
    else
    {
        int tmpi;
        for (int n = 0; n < nf; n++)
        {
            faces[n] = new int[3];
            std::fread(&tmpi, sizeof(int), 1, input); // Skip the 3
            for (int i = 0; i < 3; i++)
            {
                std::fread(&tmpi, sizeof(int), 1, input);
                swapByteOrder(tmpi);
                faces[n][i] = tmpi;
            }
        }
    }

    fclose(input);

    return true;
}

bool NIBR::Surface::readGIIMesh()
{
    gifti_image *gifti = gifti_read_image(filePath.c_str(), 1);
    for (int i = 0; i < gifti->numDA; i++)
    {
        if (gifti->darray[i]->intent == NIFTI_INTENT_POINTSET)
        {
            vertices = new float *[nv];
            for (int j = 0; j < nv; j++)
            {
                vertices[j] = new float[3];
                for (int k = 0; k < 3; k++)
                    vertices[j][k] = ((float *)gifti->darray[i]->data)[j * 3 + k];
            }
        }
        if (gifti->darray[i]->intent == NIFTI_INTENT_TRIANGLE)
        {
            faces = new int *[nf];
            for (int j = 0; j < nf; j++)
            {
                faces[j] = new int[3];
                for (int k = 0; k < 3; k++)
                    faces[j][k] = ((int *)gifti->darray[i]->data)[j * 3 + k];
            }
        }
    }
    
    delete gifti;
    return true;
}

bool NIBR::Surface::readFreesurferMesh()
{

    FILE *input;
    input = fopen(filePath.c_str(), "rb");

    // Skip header since all is checked by the point
    const size_t strLength = 256;
    char dummy[strLength];
    std::fseek(input, sizeof(unsigned char) * 3, SEEK_CUR); // Skip checking magic number
    for (int i = 0; i < 2; i++)
        fgets(dummy, strLength, input);           // Skip the metadata
    std::fseek(input, sizeof(int) * 2, SEEK_CUR); // Skip getting number of vertices and faces, these are already done before

    vertices = new float *[nv];
    float tmpf;
    for (int n = 0; n < nv; n++)
    {
        vertices[n] = new float[3];
        for (int i = 0; i < 3; i++)
        {
            std::fread(&tmpf, sizeof(float), 1, input);
            swapByteOrder(tmpf);
            vertices[n][i] = tmpf;
        }
    }

    faces = new int *[nf];
    int tmpi;
    for (int n = 0; n < nf; n++)
    {
        faces[n] = new int[3];
        for (int i = 0; i < 3; i++)
        {
            std::fread(&tmpi, sizeof(int), 1, input);
            swapByteOrder(tmpi);
            faces[n][i] = tmpi;
        }
    }

    fclose(input);

    return true;
}

// Writers
bool NIBR::Surface::write(std::string _filename){

    std::string outExtension = _filename.substr(_filename.find_last_of(".") + 1);

    if(outExtension == "vtk")
        return writeVTK(_filename.c_str());
    if(outExtension == "gii")
        return writeGII(_filename.c_str());


    disp(MSG_ERROR,"Unknown file format.");
    return false;
}

bool NIBR::Surface::writeVTK(std::string _filePath)
{

    FILE *out;

    out = fopen(_filePath.c_str(), "wb");

    char buffer[256];
    sprintf(buffer, "# vtk DataFile Version 3.0\n");
    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "Generated by ");
    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, SGNTR().c_str());
    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "\n");
    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "BINARY\n");
    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "DATASET POLYDATA\n");
    fwrite(buffer, sizeof(char), strlen(buffer), out);

    sprintf(buffer, "POINTS %lu float\n", (size_t)(nv));
    fwrite(buffer, sizeof(char), strlen(buffer), out);
    float tmpf;

    for (int n = 0; n < nv; n++)
    {
        tmpf = vertices[n][0];
        swapByteOrder(tmpf);
        fwrite(&tmpf, sizeof(float), 1, out);
        tmpf = vertices[n][1];
        swapByteOrder(tmpf);
        fwrite(&tmpf, sizeof(float), 1, out);
        tmpf = vertices[n][2];
        swapByteOrder(tmpf);
        fwrite(&tmpf, sizeof(float), 1, out);
    }

    sprintf(buffer, "POLYGONS %lu %lu\n", (size_t)(nf), (size_t)(nf * 4));
    fwrite(buffer, sizeof(char), strlen(buffer), out);

    int tmpi;
    for (int n = 0; n < nf; n++)
    {
        tmpi = 3;
        swapByteOrder(tmpi);
        fwrite(&tmpi, sizeof(int), 1, out);
        tmpi = faces[n][0];
        swapByteOrder(tmpi);
        fwrite(&tmpi, sizeof(int), 1, out);
        tmpi = faces[n][1];
        swapByteOrder(tmpi);
        fwrite(&tmpi, sizeof(int), 1, out);
        tmpi = faces[n][2];
        swapByteOrder(tmpi);
        fwrite(&tmpi, sizeof(int), 1, out);
    }

    std::vector<int> cellDataInd;
    std::vector<int> pointDataInd;

    for (size_t i = 0; i < fields.size(); i++)
    {
        if (fields[i].owner == FACE)
            cellDataInd.push_back(i);
        if (fields[i].owner == VERTEX)
            pointDataInd.push_back(i);
    }

    if (cellDataInd.size() > 0)
    {

        sprintf(buffer, "CELL_DATA %lu\n", (size_t)(nf));
        fwrite(buffer, sizeof(char), strlen(buffer), out);

        for (size_t i = 0; i < cellDataInd.size(); i++)
        {

            sprintf(buffer, "SCALARS %s %s %d\n", fields[cellDataInd[i]].name.c_str(), fields[cellDataInd[i]].datatype.c_str(), fields[cellDataInd[i]].dimension);
            fwrite(buffer, sizeof(char), strlen(buffer), out);
            sprintf(buffer, "LOOKUP_TABLE default\n");
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            if (fields[cellDataInd[i]].datatype == "float")
            {
                for (size_t s = 0; s < (size_t)(nf); s++)
                {
                    for (int d = 0; d < fields[cellDataInd[i]].dimension; d++)
                    {
                        tmpf = fields[cellDataInd[i]].fdata[s][d];
                        swapByteOrder(tmpf);
                        fwrite(&tmpf, sizeof(float), 1, out);
                    }
                }
            }
            if (fields[cellDataInd[i]].datatype == "int")
            {
                for (size_t s = 0; s < (size_t)(nf); s++)
                {
                    for (int d = 0; d < fields[cellDataInd[i]].dimension; d++)
                    {
                        tmpi = fields[cellDataInd[i]].idata[s][d];
                        swapByteOrder(tmpi);
                        fwrite(&tmpi, sizeof(int), 1, out);
                    }
                }
            }
        }
    }

    if (pointDataInd.size() > 0)
    {

        sprintf(buffer, "POINT_DATA %lu\n", (size_t)(nv));
        fwrite(buffer, sizeof(char), strlen(buffer), out);

        for (size_t i = 0; i < pointDataInd.size(); i++)
        {

            sprintf(buffer, "SCALARS %s %s %d\n", fields[pointDataInd[i]].name.c_str(), fields[pointDataInd[i]].datatype.c_str(), fields[pointDataInd[i]].dimension);
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            sprintf(buffer, "LOOKUP_TABLE default\n");
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            if (fields[pointDataInd[i]].datatype == "float")
            {
                for (size_t s = 0; s < (size_t)(nv); s++)
                {
                    for (int d = 0; d < fields[pointDataInd[i]].dimension; d++)
                    {
                        tmpf = fields[pointDataInd[i]].fdata[s][d];
                        swapByteOrder(tmpf);
                        fwrite(&tmpf, sizeof(float), 1, out);
                    }
                }
            }
            if (fields[pointDataInd[i]].datatype == "int")
            {
                for (size_t s = 0; s < (size_t)(nv); s++)
                {
                    for (int d = 0; d < fields[pointDataInd[i]].dimension; d++)
                    {
                        tmpi = fields[pointDataInd[i]].idata[s][d];
                        swapByteOrder(tmpi);
                        fwrite(&tmpi, sizeof(int), 1, out);
                    }
                }
            }
        }
    }

    fclose(out);
    return true;
}

bool NIBR::Surface::writeGII(std::string _filename)
{
    // Create empty gifti
    gifti_image *gifti = new gifti_image();
    gifti_add_to_meta(&(gifti->meta),"Created by",SGNTR().c_str(),0);

    // Add point data array
    gifti_add_empty_darray(gifti, 1);
    gifti_set_DA_defaults(gifti->darray[0]);
    gifti->darray[0]->intent   = NIFTI_INTENT_POINTSET;
    gifti->darray[0]->datatype = NIFTI_TYPE_FLOAT32;
    gifti->darray[0]->num_dim  = 2;
    gifti->darray[0]->dims[0]  = nv;
    gifti->darray[0]->dims[1]  = 3;
    gifti->darray[0]->nvals    = nv * 3;
    gifti->darray[0]->encoding = GIFTI_ENCODING_B64GZ;
    gifti->darray[0]->endian   = GIFTI_ENDIAN_LITTLE;
    gifti->darray[0]->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;

    // Add triangle data array
    gifti_add_empty_darray(gifti, 1);
    gifti_set_DA_defaults(gifti->darray[1]);
    gifti->darray[1]->intent   = NIFTI_INTENT_TRIANGLE;
    gifti->darray[1]->datatype = NIFTI_TYPE_INT32;
    gifti->darray[1]->num_dim  = 2;
    gifti->darray[1]->dims[0]  = nf;
    gifti->darray[1]->dims[1]  = 3;
    gifti->darray[1]->nvals    = nf * 3;
    gifti->darray[1]->encoding = GIFTI_ENCODING_B64GZ;
    gifti->darray[1]->endian   = GIFTI_ENDIAN_LITTLE;
    gifti->darray[1]->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;

    // Allocate memory
    gifti_update_nbyper(gifti);
    gifti_alloc_DA_data(gifti, NULL, 0);

    // Write gifti file
    
    for (int i = 0; i < gifti->numDA; i++)
    {
        if (gifti->darray[i]->intent == NIFTI_INTENT_POINTSET)
            for (int j = 0; j < nv; j++)
                memcpy(static_cast<char *>(gifti->darray[i]->data) + 3 * sizeof(float) * j, vertices[j], 3 * sizeof(float));
        if (gifti->darray[i]->intent == NIFTI_INTENT_TRIANGLE)
            for (int j = 0; j < nf; j++)
                memcpy(static_cast<char *>(gifti->darray[i]->data) + 3 * sizeof(int) * j, faces[j], 3 * sizeof(int));
    }
    gifti_write_image(gifti, _filename.c_str(), 1);
    
    // Clean up
    gifti_free_image(gifti);
    return true;
}

// This conversion is for surfaces generated with FSL's first.
// It has its own convention. Surface coordinates are in image space but they are scaled by voxel size.
// Additionally, if the ijk2xyz mapping has positive determinant then direction of x is swapped.
void NIBR::Surface::convertFSLsurface2RASmm(std::string inp_template)
{

    Image<float> *img = new Image<float>(inp_template);
    img->readHeader();

    float det = img->ijk2xyz[0][0] * img->ijk2xyz[1][1] * img->ijk2xyz[2][2] - img->ijk2xyz[0][0] * img->ijk2xyz[2][1] * img->ijk2xyz[1][2] - img->ijk2xyz[1][0] * img->ijk2xyz[0][1] * img->ijk2xyz[2][2] + img->ijk2xyz[1][0] * img->ijk2xyz[2][1] * img->ijk2xyz[0][2] + img->ijk2xyz[2][0] * img->ijk2xyz[0][1] * img->ijk2xyz[1][2] - img->ijk2xyz[2][0] * img->ijk2xyz[1][1] * img->ijk2xyz[0][2];

    float p[3];
    for (int n = 0; n < nv; n++)
    {

        if (det > 0)
            p[0] = (img->imgDims[0] - vertices[n][0]) / img->pixDims[0];
        else
            p[0] = vertices[n][0] / img->pixDims[0];

        p[1] = vertices[n][1] / img->pixDims[1];
        p[2] = vertices[n][2] / img->pixDims[2];

        img->to_xyz(p, vertices[n]);
    }

    delete img;
}

void NIBR::Surface::convertFreesurferSurface2RASmm(std::string inp_template)
{

    Image<float> *img = new Image<float>(inp_template);
    img->readHeader();

    float p[3];
    for (int n = 0; n < nv; n++)
    {

        p[0] = vertices[n][0];
        p[1] = vertices[n][1];
        p[2] = vertices[n][2];

        vertices[n][0] = p[0] * img->rastkr2ras[0][0] + p[1] * img->rastkr2ras[0][1] + p[2] * img->rastkr2ras[0][2] + img->rastkr2ras[0][3];
        vertices[n][1] = p[0] * img->rastkr2ras[1][0] + p[1] * img->rastkr2ras[1][1] + p[2] * img->rastkr2ras[1][2] + img->rastkr2ras[1][3];
        vertices[n][2] = p[0] * img->rastkr2ras[2][0] + p[1] * img->rastkr2ras[2][1] + p[2] * img->rastkr2ras[2][2] + img->rastkr2ras[2][3];
    }

    delete img;
}