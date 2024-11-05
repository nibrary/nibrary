#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Warray-bounds="
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wstringop-overflow="
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-result"
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#include <igl/AABB.h>
#include <igl/fast_winding_number.h>
#pragma GCC diagnostic pop

#else

#include <igl/AABB.h>
#include <igl/fast_winding_number.h>

#endif

#include <Eigen/Core>
#include <Eigen/Dense>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/mesh_fill_holes.cpp>
#include <geogram/mesh/mesh_repair.h>

#include "surface_operators.h"
// This file contains all the geogram function.
// Geogram has its own thread handling mechanism.
// In certain cases, it might delete all the running geogram threads, e.g., for remeshing
// That is why it is not possible to run any application in parallel
// if geogram routines are used, especially involving remeshing.
// 
// TODO: Make a task queue for geogram function. Add geogram functions in a queue, so they
// run sequentially.

using namespace NIBR;

void surf2Geo(const NIBR::Surface& surf, GEO::Mesh& M) {

    // Add vertices
    M.vertices.create_vertices(surf.nv);  
    for(int i=0; i<surf.nv; ++i) {      
        M.vertices.point(i) = GEO::vec3(surf.vertices[i][0], surf.vertices[i][1], surf.vertices[i][2]);
    }

    // Add faces
    M.facets.create_facets(surf.nf,3);
    for(int i=0; i<surf.nf; ++i) {
        M.facets.set_vertex(i, 0, surf.faces[i][0]);
        M.facets.set_vertex(i, 1, surf.faces[i][1]);
        M.facets.set_vertex(i, 2, surf.faces[i][2]);
    }

    // Many geogram functions require this 
    // GEO::mesh_connect_and_reorient_facets_no_check(M);
    M.facets.connect();
    GEO::compute_normals(M);

    return;

}

Surface geo2surf(const GEO::Mesh& M) {

    Surface surf;
    
    // Add vertices
    surf.nv       = M.vertices.nb();
    surf.vertices = new float*[surf.nv];

    for (int i = 0; i < surf.nv; ++i) {
        surf.vertices[i]    = new float[3];
        surf.vertices[i][0] = M.vertices.point(i).x;
        surf.vertices[i][1] = M.vertices.point(i).y;
        surf.vertices[i][2] = M.vertices.point(i).z;
    }

    // Add faces
    surf.nf    = M.facets.nb();
    surf.faces = new int*[surf.nf];

    for(int i = 0; i < surf.nf; ++i) {
        surf.faces[i]    = new int[3];
        surf.faces[i][0] = M.facets.vertex(i,0);
        surf.faces[i][1] = M.facets.vertex(i,1);
        surf.faces[i][2] = M.facets.vertex(i,2);
    }

    return surf;

}

NIBR::Surface NIBR::surfRepair(const NIBR::Surface& surf) {

    if (surf.nv == 0) {return Surface();}

    Surface out(surf,true);

    bool allFixed = false;

    int iterationLimit = 100;
    int iterationCnt   = 0;

    while (!allFixed && (iterationCnt<iterationLimit)) {
        int initVertexCnt = out.nv;
        int initFaceCnt   = out.nf;

        disableTerminalOutput();
        GEO::Mesh M;
        surf2Geo(out,M);
 
        GEO::remove_degree2_vertices(M);
        GEO::mesh_repair(M,GEO::MeshRepairMode(GEO::MESH_REPAIR_DEFAULT | GEO::MESH_REPAIR_QUIET), EPS8);

        out = geo2surf(M);
        enableTerminalOutput();

        out = surfRemoveVerticesWithNAN(out);
        out = surfRemoveSingularVertices(out);
        out = surfRemoveOverconnectedVertices(out);
        out = surfRemoveBadBoundaryVertices(out);

        out.categorizeEdges();
        out.categorizeVertices();

        if ((out.nv == initVertexCnt) && (out.nf == initFaceCnt) && (out.boundaryEdges.size() == out.boundaryVertices.size())){
            allFixed = true;
        } else if (out.nv == 0) {
            allFixed = true;
            disp(MSG_WARN, "Surface repair removed all the vertices.");
        }

        iterationCnt++;
    }

    if (!allFixed) {
        disp(MSG_WARN, "Surface repair failed.");
    }

    return out;

}


NIBR::Surface NIBR::surfRemoveSmallFaces(const NIBR::Surface& surf, double minArea) {

    if (surf.nv == 0) {return Surface();}
    
    disableTerminalOutput();
    GEO::Mesh M;
    surf2Geo(surfRepair(surf),M);
    
    GEO::remove_small_facets(M,minArea);
    enableTerminalOutput();

    return surfRepair(geo2surf(M));
    
}

NIBR::Surface NIBR::surfRemoveSmallConnectedComponents(const NIBR::Surface& surf,double minArea) {

    if (surf.nv == 0) {return Surface();}

    disp(MSG_DEBUG,"surfRemoveSmallConnectedComponents");

    Surface out(surf,true);

    out = surfRemoveVerticesWithNAN(out);
    out = surfRemoveSingularVertices(out);

    auto area   = out.calcAreasOfConnectedComponents();

    // out.printInfo();

    std::vector<Surface> comp;
    std::swap(comp,out.comp);

    for (int n = 0; n < int(area.size()); n++) {
        if (area[n] < minArea) {
            out = surfDiff(out,comp[n]);
        }
    }

    disp(MSG_DEBUG,"Done surfRemoveSmallConnectedComponents");

    return out;

}

NIBR::Surface NIBR::surfFillHoles(const NIBR::Surface& surf, double maxArea) {

    if (surf.nv == 0) {return Surface();}

    disableTerminalOutput();
    GEO::Mesh M;
    surf2Geo(surfRepair(surf),M);   

    GEO::fill_holes(M,maxArea);
    enableTerminalOutput();

    return surfRepair(geo2surf(M));
    
}

NIBR::Surface NIBR::surfSmooth(const NIBR::Surface& surf,int smoothing_iterations) {

    if (surf.nv == 0) {return Surface();}
    
    disableTerminalOutput();
    GEO::Mesh M;
    surf2Geo(surfRepair(surf),M);

    GEO::simple_Laplacian_smooth(M, smoothing_iterations, false);
    enableTerminalOutput();

    return surfRepair(geo2surf(M));
    
}

NIBR::Surface NIBR::surfRemesh(const NIBR::Surface& surf,int newVertexCount, int smoothing_iterations, float anisotropy) {

    if (surf.nv == 0) {return Surface();}
    
    GEO::Mesh M;
    surf2Geo(surfRepair(surf),M);
    GEO::Mesh remeshedM;

    disableTerminalOutput();
    GEO::simple_Laplacian_smooth(M, smoothing_iterations, true);
    GEO::set_anisotropy(M, anisotropy);
    GEO::remesh_smooth(M, remeshedM, newVertexCount);
    enableTerminalOutput();

    return surfRepair(geo2surf(remeshedM));
    
}

std::vector<double> NIBR::surfCalcAreasOfHoles(const NIBR::Surface& surf) {

    if (surf.nv == 0) {return std::vector<double>();}
    
    GEO::Mesh M;
    surf2Geo(surf,M);

    GEO::MeshHalfedges MH(M);
    GEO::vector<Hole> holes;

    // Identification of holes 
    {

        GEO::vector<bool> corner_is_visited(M.facet_corners.nb(), false);
        for(GEO::index_t f: M.facets) {
            for(GEO::index_t c: M.facets.corners(f)) {
                if(
                    M.facet_corners.adjacent_facet(c) == GEO::NO_FACET && !corner_is_visited[c]
                ) {
                    holes.push_back(Hole());
                    GEO::MeshHalfedges::Halfedge first(f, c);
                    GEO::MeshHalfedges::Halfedge H(f, c);
                    do {
                        holes.rbegin()->push_back(H);
                        corner_is_visited[H.corner] = true;
                        MH.move_to_next_around_border(H);
                    } while(H != first);

                }
            }
        }
    }


    // Estimation of hole areas using LOOP_SPLIT, and preparing output
    std::vector<double> hole_areas(holes.size(), 0.0);
    
    for(GEO::index_t i = 0; i < holes.size(); i++) {
        GEO::vector<GEO::trindex> triangles;
        bool ok = triangulate_hole_loop_splitting(MH, holes[i], triangles, false);
        hole_areas[i] = ok ? hole_area(M, triangles) : 0;
    }    

    // Sort comp_area in descending order
    std::sort(hole_areas.begin(), hole_areas.end(), std::greater<double>());

    return hole_areas;

}