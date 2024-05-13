#include "image_marchingCubes.h"

using namespace NIBR;

bool NIBR::isosurface(Image<float> *img, float isoValue, Surface *surf)
{
    if (surf==NULL)
        surf = new Surface();

    auto mc = MarchingCubes (img->imgDims[0], img->imgDims[1], img->imgDims[2]);
    mc.set_ext_data(img->data);
    mc.init_all();
    mc.run(isoValue);

    auto v = mc.vertices();
    auto f = mc.triangles();

    surf->nv = mc.nverts();
    surf->nf = mc.ntrigs();

    surf->vertices = new float *[surf->nv];
    for (int n = 0; n < surf->nv; n++) {
        surf->vertices[n] = new float[3];
        applyTransform(surf->vertices[n],&v[n].x,img->ijk2xyz);
    }

    // Flip face orientation so they point outwards
    surf->faces = new int *[surf->nf];
    for (int n = 0; n < surf->nf; n++) {
        surf->faces[n] = new int[3];
        surf->faces[n][0] = f[n].v3;
        surf->faces[n][1] = f[n].v2;
        surf->faces[n][2] = f[n].v1;
    }

    mc.clean_all();

    return true;

}