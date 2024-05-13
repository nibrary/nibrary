#include "surface_make.h"

using namespace NIBR;

Surface NIBR::surfMakeBox(const std::vector<float>& bbox) 
{

    float** vertices = new float*[8];

    for (int i = 0; i < 8; ++i) {
        vertices[i] = new float[3];
    }

    vertices[0][0] = bbox[0]; vertices[0][1] = bbox[2]; vertices[0][2] = bbox[4];
    vertices[1][0] = bbox[1]; vertices[1][1] = bbox[2]; vertices[1][2] = bbox[4];
    vertices[2][0] = bbox[1]; vertices[2][1] = bbox[3]; vertices[2][2] = bbox[4];
    vertices[3][0] = bbox[0]; vertices[3][1] = bbox[3]; vertices[3][2] = bbox[4];
    vertices[4][0] = bbox[0]; vertices[4][1] = bbox[2]; vertices[4][2] = bbox[5];
    vertices[5][0] = bbox[1]; vertices[5][1] = bbox[2]; vertices[5][2] = bbox[5];
    vertices[6][0] = bbox[1]; vertices[6][1] = bbox[3]; vertices[6][2] = bbox[5];
    vertices[7][0] = bbox[0]; vertices[7][1] = bbox[3]; vertices[7][2] = bbox[5];

    // Define the 12 triangles
    int faceIndices[12][3] = {
        {2, 1, 0}, {0, 3, 2},   // Bottom face
        {4, 5, 6}, {6, 7, 4},   // Top face
        {0, 1, 5}, {5, 4, 0},   // Front face
        {1, 2, 6}, {6, 5, 1},   // Right face
        {2, 3, 7}, {7, 6, 2},   // Back face
        {3, 0, 4}, {4, 7, 3}    // Left face
    };

    // Allocate memory for faces
    int** faces = new int*[12];
    for(int i = 0; i < 12; ++i) {
        faces[i] = new int[3];
        for(int j = 0; j < 3; ++j) {
            faces[i][j] = faceIndices[i][j];
        }
    }

    Surface box;
    box.nv = 8;
    box.nf = 12;
    box.vertices = vertices;
    box.faces    = faces;

    return box;

}