#include "surface_operators.h"
#include "surface.h"
#include "surface2imageMapper.h"
#include "image/image_math.h"
#include "image/image_operators.h"
#include "image/image_marchingCubes.h"
#include "base/vectorOperations.h"
#include "math/core.h"
#include "math/triangle.h"
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <tuple>
#include <queue>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseCholesky>

using namespace NIBR;

Surface NIBR::surfRemoveVerticesWithNAN(Surface& surf)
{

    if (surf.nv == 0) {return Surface();}

    // disp(MSG_DEBUG,"surfRemoveVerticesWithNAN");

    std::vector<bool> vertexMask(surf.nv, false);
    
    bool runApplyMask = false;
    
    for (int n = 0; n < surf.nv; n++) {
        vertexMask[n] = !(isnan(surf.vertices[n][0]) || isnan(surf.vertices[n][1]) || isnan(surf.vertices[n][2]));
        if (!runApplyMask && !vertexMask[n]) {
            runApplyMask = true;
        }
    }

    // disp(MSG_DEBUG,"Done..calling surfRemoveVerticesWithNAN");

    return (runApplyMask) ? applyMask(surf, vertexMask) : surf;

}

Surface NIBR::surfRemoveSingularVertices(Surface& surf)
{

    if (surf.nv == 0) {return Surface();}

    // disp(MSG_DEBUG,"surfRemoveSingularVertices");

    surf.categorizeVertices();

    std::vector<bool> vertexMask(surf.nv, true);
    for (int v : surf.singularVertices) {
        vertexMask[v] = false;
    }

    // disp(MSG_DEBUG,"Done..calling surfRemoveSingularVertices");

    return (!surf.singularVertices.empty()) ? applyMask(surf, vertexMask) : surf;
}

Surface NIBR::surfRemoveOverconnectedVertices(Surface& surf)
{

    if (surf.nv == 0) {return Surface();}

    // disp(MSG_DEBUG,"surfRemoveOverconnectedVertices");

    surf.categorizeVertices();

    std::vector<bool> vertexMask(surf.nv, true);
    for (int v : surf.overconnectedVertices) {
        vertexMask[v] = false;
    }

    // disp(MSG_DEBUG,"Done..calling surfRemoveOverconnectedVertices");

    return (!surf.overconnectedVertices.empty()) ? applyMask(surf, vertexMask) : surf;

}

Surface NIBR::surfRemoveBadBoundaryVertices(Surface& surf)
{

    if (surf.nv == 0) {return Surface();}

    // disp(MSG_DEBUG,"surfRemoveBadBoundaryVertices");

    Surface out(surf,true);
    out.computeBoundaries();

    std::vector<bool> vertexMask(out.nv, true);
    bool runApplyMask = false;

    for (const auto& boundary : out.boundaries) {
        for (const auto& v : boundary) {
            if (std::find(out.boundaryVertices.begin(),out.boundaryVertices.end(),v) == out.boundaryVertices.end()) {
                vertexMask[v] = false;
                if (!runApplyMask) {
                    runApplyMask = true;
                }
            }
        }
    }

    // disp(MSG_DEBUG,"Done..calling surfRemoveBadBoundaryVertices");

    return (runApplyMask) ? applyMask(surf, vertexMask) : surf;

}

Surface NIBR::surfRemoveAllBoundaryVertices(Surface& surf)
{

    if (surf.nv == 0) {return Surface();}

    // disp(MSG_DEBUG,"surfRemoveAllBoundaryVertices");

    Surface out(surf,true);
    out.categorizeVertices();

    std::vector<bool> vertexMask(surf.nv, true);
    for (int v : surf.boundaryVertices) {
        vertexMask[v] = false;
    }

    // disp(MSG_DEBUG,"Done..calling surfRemoveAllBoundaryVertices");

    return (!surf.boundaryVertices.empty()) ? applyMask(surf, vertexMask) : surf;

}

Surface NIBR::surfMerge(const Surface& s1, const Surface& s2) 
{

    // disp(MSG_DEBUG,"surfMerge");

    if (s1.nv == 0) return s2;
    if (s2.nv == 0) return s1;

    // Allocate and copy vertices
    Surface out;

	out.nv       = s1.nv + s2.nv;
    out.vertices = new float*[out.nv];

    for (int i=0; i<s1.nv; i++) {
        out.vertices[i] = new float[3];
        memcpy(out.vertices[i], s1.vertices[i], 3*sizeof(float));
    }

    for (int i=0; i<s2.nv; i++) {
        out.vertices[i+s1.nv] = new float[3];
        memcpy(out.vertices[i+s1.nv], s2.vertices[i], 3*sizeof(float));
    }
    
    // Allocate and copy faces
    out.nf     = s1.nf + s2.nf;
    out.faces  = new int*[out.nf];
    for (int i=0; i<out.nf; i++) out.faces[i] = new int[3];

    for (int i=0; i<s1.nf; i++) {
        memcpy(out.faces[i], s1.faces[i], 3*sizeof(int));
    }

    for (int i=0; i<s2.nf; i++) {
        out.faces[i+s1.nf][0] = s2.faces[i][0] + s1.nv;
        out.faces[i+s1.nf][1] = s2.faces[i][1] + s1.nv;
        out.faces[i+s1.nf][2] = s2.faces[i][2] + s1.nv;
    }

    // disp(MSG_DEBUG,"Done surfMerge");
    return out;

}

Surface NIBR::surfDiff(const Surface& s1, const Surface& s2) 
{

    if (s2.nv == 0) {
        Surface out = s1;
        return out;
    }

    // disp(MSG_DEBUG,"surfDiff");

    std::vector<bool> vertexMask(s1.nv);
    for (int i = 0; i < s1.nv; i++) {
        for (int j = 0; j < s2.nv; j++) {
            vertexMask[i] = ((std::fabs(s1.vertices[i][0] - s2.vertices[j][0]) > EPS6 ) || (std::fabs(s1.vertices[i][1] - s2.vertices[j][1]) > EPS6 ) || (std::fabs(s1.vertices[i][2] - s2.vertices[j][2]) > EPS6 ));
            if (vertexMask[i]==0)
                break;
        }
    }

    // disp(MSG_DEBUG,"Done surfDiff..calling applyMask");

    return applyMask(s1, vertexMask);
}

Surface NIBR::surfRemoveAllButLargestComponent(const Surface& surf) {

    if (surf.nv == 0) return Surface();

    Surface out(surf,true);

    out = surfRemoveVerticesWithNAN(out);
    out = surfRemoveSingularVertices(out);

    auto area   = out.calcAreasOfConnectedComponents();

    double maxArea = 0;
    int    maxInd  = 0;

    for (int n = 0; n < int(area.size()); n++) {
        if (area[n] > maxArea) {
            maxArea = area[n];
            maxInd  = n;
        }
    }

    return out.comp[maxInd];
}

Surface NIBR::surfFillAllButLargestHole(const Surface& surf) {
    Surface out(surf,true);
    auto areas = surfCalcAreasOfHoles(out);
    if (areas.size()>1) {
        out = surfFillHoles(out,areas.front()-EPS2);
    }
    return out;
}

Surface NIBR::surfMakeItSingleOpen(const Surface& surf) {

    if (surf.nv == 0) {return Surface();}

    Surface out(surf,true);

    out.getConnectedComponents();
    out.computeBoundaries();
    if ( (out.comp.size() == 1) && (out.boundaries.size() == 1) ) 
        return out;

    out = surfRepair(out);
    
    out.getConnectedComponents();
    out.computeBoundaries();
    if ( (out.comp.size() == 1) && (out.boundaries.size() == 1) ) 
        return out;
    
    if (out.boundaries.size() == 0) return out;

    int initNumberOfVertices = surf.nv;
    int initNumberOfFaces    = surf.nf;

    bool doneIterating = false;
    
    do {

        out = surfRemoveAllButLargestComponent(out);
        out = surfFillAllButLargestHole(out);

        out.getConnectedComponents();
        out.computeBoundaries();
        
        if ( (out.comp.size() == 1) && (out.boundaries.size() == 1) ) {
            doneIterating = true;
        } else if ((initNumberOfVertices == out.nv) && (initNumberOfFaces == out.nf)) {
            out = surfRemoveAllBoundaryVertices(out);
        } else if (out.nv == 0) {
            doneIterating = true;
            disp(MSG_ERROR,"Failed at making a single open surface.");
        } else {
            disp(MSG_DETAIL, "Remove %d vertices", initNumberOfVertices - out.nv);
            initNumberOfVertices = out.nv;
            initNumberOfFaces    = out.nf;
        }

    } while (!doneIterating);

    return out;

}

Surface NIBR::surfMakeItWatertight(Surface& surf) {

    if (surf.nv == 0) {return Surface();}    

    Surface out(surf,true);

    out.categorizeVertices();
    if (out.boundaryVertices.size() == 0)
        return out;

    out = surfRepair(out);

    out.categorizeVertices();
    if (out.boundaryVertices.size() == 0)
        return out;

    int initNumberOfVertices = surf.nv;
    int initNumberOfFaces    = surf.nf;

    bool doneIterating = false;
    
    do {

        out = surfFillHoles(out, DBL_MAX);
        out = surfRemoveAllBoundaryVertices(out);
        out.categorizeVertices();

        if (out.boundaryVertices.size() == 0) {
            doneIterating = true;
        } else if ((initNumberOfVertices == out.nv) && (initNumberOfFaces == out.nf)) {
            out = surfRemoveAllBoundaryVertices(out);
        } else if (out.nv == 0) {
            doneIterating = true;
            disp(MSG_ERROR,"Failed at making a watertight surface.");
        } else {
            disp(MSG_DETAIL, "Remove %d vertices", initNumberOfVertices - out.nv);
            initNumberOfVertices = out.nv;
            initNumberOfFaces    = out.nf;
        }

    } while (!doneIterating);

    return out;

}

Surface NIBR::surfMakeItSingleClosed(const Surface& surf) {

    if (surf.nv == 0) {return Surface();}

    Surface out(surf,true);

    out.getConnectedComponents();
    out.categorizeVertices();
    if ( (out.comp.size() == 1) && (out.boundaryVertices.size() == 0) ) 
        return out;

    out = surfRepair(out);

    out.getConnectedComponents();
    out.categorizeVertices();
    if ( (out.comp.size() == 1) && (out.boundaryVertices.size() == 0) ) 
        return out;


    int initNumberOfVertices = surf.nv;
    int initNumberOfFaces    = surf.nf;

    bool doneIterating = false;
    
    do {

        out = surfRemoveAllButLargestComponent(out);
        out = surfFillHoles(out, DBL_MAX);
        out = surfRemoveAllBoundaryVertices(out);

        out.getConnectedComponents();
        out.categorizeVertices();
        
        if ( (out.comp.size() == 1) && (out.boundaryVertices.size() == 0) ) {
            doneIterating = true;
        } else if ((initNumberOfVertices == out.nv) && (initNumberOfFaces == out.nf)) {
            out = surfRemoveAllBoundaryVertices(out);
        } else if (out.nv == 0) {
            doneIterating = true;
            disp(MSG_ERROR,"Failed at making a single closed surface.");
            return Surface();
        } else {
            disp(MSG_DETAIL, "Remove %d vertices", initNumberOfVertices - out.nv);
            initNumberOfVertices = out.nv;
            initNumberOfFaces    = out.nf;
        }

    } while (!doneIterating);

    return out;

}

Surface NIBR::surfMoveVerticesAlongNormal(const Surface& surf, float shift) {

    if (surf.nv < 1) return surf;
    if (shift == 0)  return surf;

    Surface out(surf,true);
    out.calcNormalsOfVertices();

    for (int n = 0; n < out.nv; n++) {
        out.vertices[n][0] += out.normalsOfVertices[n][0]*shift;
        out.vertices[n][1] += out.normalsOfVertices[n][1]*shift;
        out.vertices[n][2] += out.normalsOfVertices[n][2]*shift;
    }

    return out;

}

Surface NIBR::meanCurvatureFlow(const Surface& surf, float dt, int iterationCount) {

    if (surf.nv < 1) return surf;
    if (dt == 0)  return surf;

    Surface out(surf, true);
    out.toEigen();

    for (int i = 0; i < iterationCount; ++i) {

        // Compute the cotangent Laplacian
        Eigen::SparseMatrix<double> L;
        igl::cotmatrix(out.V, out.F, L);

        // Compute the mass matrix
        Eigen::SparseMatrix<double> M;
        igl::massmatrix(out.V, out.F, igl::MASSMATRIX_TYPE_VORONOI, M);

        // Assemble the system matrix
        Eigen::SparseMatrix<double> A = M - dt * L;

        // Factorize the system matrix
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);

        if (solver.info() != Eigen::Success) {
            disp(MSG_FATAL,"Failed to factorize the system matrix.");
            return out; // Handle error appropriately
        }

        Eigen::MatrixXd B = M * out.V;

        Eigen::MatrixXd V_new = solver.solve(B);

        if (solver.info() != Eigen::Success) {
            disp(MSG_FATAL,"Failed to solve the linear system.");
            return out; // Handle error appropriately
        }

        out.V = V_new;
    }

    out.fromEigen();

    return out;
}

Surface surfMoveVerticesAlongNormal(const Surface& surf, float* vertexShift);

std::vector<float> NIBR::surfBbox(const Surface& surf)
{
    
    if (surf.nv<1)
        return std::vector<float>();

    std::vector<float> bbox(6,0);
    
    bbox[0] = surf.vertices[0][0];
    bbox[1] = surf.vertices[0][0];
    bbox[2] = surf.vertices[0][1];
    bbox[3] = surf.vertices[0][1];
    bbox[4] = surf.vertices[0][2];
    bbox[5] = surf.vertices[0][2];

    for (int n = 0; n < surf.nv; n++) {
        if ( surf.vertices[n][0] < bbox[0] )   bbox[0] = surf.vertices[n][0];
        if ( surf.vertices[n][0] > bbox[1] )   bbox[1] = surf.vertices[n][0];
        if ( surf.vertices[n][1] < bbox[2] )   bbox[2] = surf.vertices[n][1];
        if ( surf.vertices[n][1] > bbox[3] )   bbox[3] = surf.vertices[n][1];
        if ( surf.vertices[n][2] < bbox[4] )   bbox[4] = surf.vertices[n][2];        
        if ( surf.vertices[n][2] > bbox[5] )   bbox[5] = surf.vertices[n][2];
    }
    
    return bbox;

}


Surface NIBR::surfGrow(const Surface& surf, float discretizationResolution, float shift) {

    Image<float> bboxImg;
    surfBbox2Img(bboxImg,surf,discretizationResolution);
    bboxImg.allocData();
    
    int pad = shift*discretizationResolution*1.1f;

    if (pad > 0) {
        imgPad(bboxImg,pad); // Allow 10% extra padding
    }

    Surface tmp(surf,true);
    mapSurface2Image(&tmp, &bboxImg, 0, NULL, NULL, EDT);

    Surface out;
    if (!isosurface(&bboxImg, -shift, &out)) {
        disp(MSG_ERROR, "Failed to generate isosurface");
        return Surface();
    }

    return out;

}


float getBoundaryDist(const Surface& s1, const std::vector<int>& s1b1, const Surface& s2, const std::vector<int>& s2b2)
{
    
    auto b1 = s1b1;
    auto b2 = s2b2;

    // Find closest boundary vertices
    float minDist = FLT_MAX;
    int   b1start = 0;
    int   b2start = 0;

    for (int i = 0; i < int(b1.size()); i++) {
        for (int j = 0; j < int(b2.size()); j++) {
            float d = dist(s1.vertices[b1[i]],s2.vertices[b2[j]]);
            if (d < minDist) {
                minDist = d;
                b1start = i;
                b2start = j;
            }
        }    
    }

    // Rotate the boundary vertices so that the closest vertices are at the start
    std::rotate(b1.begin(), b1.begin() + b1start, b1.end());
    std::rotate(b2.begin(), b2.begin() + b2start, b2.end());

    // Resample boundaries so they have the same vertices
    int sampleSize = std::min(8.0f,std::min(float(b1.size()),float(b2.size())));

    auto b1res = getEvenlySeparatedSamples(b1,sampleSize);
    auto b2res = getEvenlySeparatedSamples(b2,sampleSize);

    auto calcDist = [&]()->float {
        float totDist = 0.0f;
        for (int i = 0; i < sampleSize; ++i) {
            totDist += dist(s1.vertices[b1res[i]], s2.vertices[b2res[i]]);
        }
        return totDist;
    };
        
    // Flip the direction of second boundary if needed
    float origDist = calcDist();

    std::rotate (b2res.begin(), b2res.begin()+1, b2res.end());
    std::reverse(b2res.begin(), b2res.end());

    float flipDist = calcDist();

    return (origDist < flipDist) ? origDist : flipDist;

}

Surface NIBR::surfGlueBoundaries(const Surface& inp1, const Surface& inp2) 
{

    disp(MSG_DEBUG,"surfGlueBoundaries");

    if (inp1.nv == 0) {return Surface();}
    if (inp2.nv == 0) {return Surface();}

    Surface s1(inp1,true);
    Surface s2(inp2,true);

    s1.categorizeVertices();
    if (s1.boundaryVertices.empty()) {return Surface();}

    if (s1.boundaryEdges.size() != s1.boundaryVertices.size()) {
        disp(MSG_DEBUG,"Surface contains non-manifold vertices/edges, repairing mesh...");
        s1 = surfRepair(s1);
        s1 = surfMakeItSingleOpen(s1);
        s1.categorizeVertices();
        if (s1.boundaryEdges.size() != s1.boundaryVertices.size()) {
            disp(MSG_FATAL,"Degenerate surface.");
        }
        disp(MSG_DEBUG,"Mesh repair successful.");
    }


    s2.categorizeVertices();
    if (s2.boundaryVertices.empty()) {return Surface();}

    if (s2.boundaryEdges.size() != s2.boundaryVertices.size()) {
        disp(MSG_DEBUG,"Surface contains non-manifold vertices/edges, repairing mesh...");
        s2 = surfRepair(s2);
        s2 = surfMakeItSingleOpen(s2);
        s2.categorizeVertices();
        if (s2.boundaryEdges.size() != s2.boundaryVertices.size()) {
            disp(MSG_FATAL,"Degenerate surface.");
        }
        disp(MSG_DEBUG,"Mesh repair successful.");
    }


    // disp(MSG_DEBUG,"Vertices categorizes");

    s1.computeBoundaries();
    s2.computeBoundaries();

    std::vector<std::pair<int,int>> b2b;

    if (s1.boundaries.size() <= s2.boundaries.size()) {

        for (int i = 0; i < int(s1.boundaries.size()); i++) {

            std::pair<int,int> tmp;

            float minDist = FLT_MAX;

            for (int j = 0; j < int(s2.boundaries.size()); j++) {
                float d = getBoundaryDist(s1,s1.boundaries[i],s2,s2.boundaries[j]);
                if (d < minDist) {
                    minDist    = d;
                    tmp.first  = i;
                    tmp.second = j;
                }
            }

            b2b.push_back(tmp);
            // disp(MSG_DEBUG,"Will glue s1, b%d with s2, b%d", tmp.first, tmp.second);

        }

    } else {

        for (int i = 0; i < int(s2.boundaries.size()); i++) {

            std::pair<int,int> tmp;

            float minDist = FLT_MAX;

            for (int j = 0; j < int(s1.boundaries.size()); j++) {
                float d = getBoundaryDist(s2,s2.boundaries[i],s1,s1.boundaries[j]);
                if (d < minDist) {
                    minDist    = d;
                    tmp.first  = j;
                    tmp.second = i;
                }
            }

            b2b.push_back(tmp);
            // disp(MSG_DEBUG,"Will glue s1, b%d with s2, b%d", tmp.first, tmp.second);

        }

    }

    

    auto addNewBoundaryVertices = [](Surface& surf, std::vector<int>& boundary, size_t targetSize)->bool {

        // Add new vertices and faces on the longest edges first
        int missingVertCnt = targetSize - boundary.size();

        // disp(MSG_DEBUG, "Inserting %d new vertices", missingVertCnt);

        // Allocate memory for new vertex list, and copy the previous vertices
        int     newNv       = surf.nv + missingVertCnt; // increased vertices
        float** newVertices = new float*[newNv];
        for (int n = 0; n < newNv; n++) {
            newVertices[n] = (n < surf.nv) ? surf.vertices[n] : (new float[3]());
        }
        delete[] surf.vertices;
        int newVertexIndex  = surf.nv;
        surf.nv             = newNv;
        surf.vertices       = newVertices;

        // Allocate memory for new face list, and copy the previous faces
        int     newNf       = surf.nf + missingVertCnt; // increased faces, one faces per vertex is added
        int**   newFaces    = new int*[newNf];
        for (int n = 0; n < newNf; n++) {
            if (n < surf.nf) {
                newFaces[n] = surf.faces[n];
            } else {
                newFaces[n]     = new int[3];
                newFaces[n][0]  = INT_MAX;
                newFaces[n][1]  = INT_MAX;
                newFaces[n][2]  = INT_MAX;
            }
        }
        delete[] surf.faces;
        int newFaceIndex    = surf.nf;
        surf.nf             = newNf;
        surf.faces          = newFaces;

        
        for (int i = 0; i < missingVertCnt; ++i) {

            // Find the longest edge
            float maxLen = FLT_MIN;
            int startIndex,endIndex,boundaryInd;
            startIndex=endIndex=boundaryInd=0;
            for (size_t j = 0; j < boundary.size(); j++) {
                size_t nextIndex = (j + 1) % boundary.size();
                double len = dist(surf.vertices[boundary[j]], surf.vertices[boundary[nextIndex]]);
                if (len > maxLen) {
                    startIndex  = boundary[j];
                    endIndex    = boundary[nextIndex];
                    boundaryInd = j;
                    maxLen      = len;
                }
            }

            // Calculate the midpoint
            surf.vertices[newVertexIndex][0] = (surf.vertices[startIndex][0] + surf.vertices[endIndex][0])*0.5f;
            surf.vertices[newVertexIndex][1] = (surf.vertices[startIndex][1] + surf.vertices[endIndex][1])*0.5f;
            surf.vertices[newVertexIndex][2] = (surf.vertices[startIndex][2] + surf.vertices[endIndex][2])*0.5f;
            // disp(MSG_DEBUG, "Inserted vertex: %d, between [%d - %d]", newNv-newVertexIndex,startIndex,endIndex);
            
            
            // Split the face into two
            
            // Find the face containing the boundary edge
            int faceIndex = -1;
            for (int n = 0; n < surf.nf; n++) {
                bool containsStart = ((surf.faces[n][0] == startIndex) || (surf.faces[n][1] == startIndex) || (surf.faces[n][2] == startIndex));
                bool containsEnd   = ((surf.faces[n][0] == endIndex)   || (surf.faces[n][1] == endIndex)   || (surf.faces[n][2] == endIndex));
                if (containsStart && containsEnd) {
                    faceIndex = n;
                    break;
                }
            }

            if (faceIndex == -1) {

                surf.faces[newFaceIndex][0] = startIndex;
                surf.faces[newFaceIndex][1] = newVertexIndex;
                surf.faces[newFaceIndex][2] = endIndex;
                // disp(MSG_DEBUG, "Inserted face: %d", newNf-newFaceIndex);
                // disp(MSG_ERROR, "Can't glue boundary. Boundary edge is non-manifold. Can't find face with edge [%d-%d]",startIndex,endIndex);
                // return false;

            } else {

                int* originalFace = surf.faces[faceIndex];            

                // Find the third vertex in the face that is not part of the edge
                int thirdVertexIndex = -1;
                for (int j = 0; j < 3; ++j) {
                    if (originalFace[j] != startIndex && originalFace[j] != endIndex) {
                        thirdVertexIndex = originalFace[j];
                        break;
                    }
                }

                // Create two new faces replacing the original
                float faceDir[3], tmp1[3], tmp2[3], newFaceDir[3];
                vec3sub(tmp1,surf.vertices[originalFace[1]],surf.vertices[originalFace[0]]);
                vec3sub(tmp2,surf.vertices[originalFace[2]],surf.vertices[originalFace[0]]);
                cross(faceDir,tmp1,tmp2);
                            
                vec3sub(tmp1, surf.vertices[newVertexIndex],surf.vertices[thirdVertexIndex]);
                vec3sub(tmp2, surf.vertices[startIndex], surf.vertices[thirdVertexIndex]);
                cross(newFaceDir,tmp1,tmp2);

                if (dot(faceDir,newFaceDir) > 0) {
                    originalFace[0] = thirdVertexIndex;
                    originalFace[1] = newVertexIndex;
                    originalFace[2] = startIndex;
                } else {
                    originalFace[0] = startIndex;
                    originalFace[1] = newVertexIndex;
                    originalFace[2] = thirdVertexIndex;
                }


                vec3sub(tmp1, surf.vertices[newVertexIndex], surf.vertices[thirdVertexIndex]);
                vec3sub(tmp2, surf.vertices[endIndex], surf.vertices[thirdVertexIndex]);
                cross(newFaceDir,tmp1,tmp2);

                if (dot(faceDir,newFaceDir) > 0) {
                    surf.faces[newFaceIndex][0] = thirdVertexIndex;
                    surf.faces[newFaceIndex][1] = newVertexIndex;
                    surf.faces[newFaceIndex][2] = endIndex;
                } else {
                    surf.faces[newFaceIndex][0] = endIndex;
                    surf.faces[newFaceIndex][1] = newVertexIndex;
                    surf.faces[newFaceIndex][2] = thirdVertexIndex;
                }

                // disp(MSG_DEBUG, "Inserted face: %d", newNf-newFaceIndex);
            }

            boundary.insert(boundary.begin()+boundaryInd+1,newVertexIndex);
            // disp(MSG_DEBUG, "Boundary length: %d", int(boundary.size()));

            newVertexIndex++;
            newFaceIndex++;
            
        }

        return true;

    };


    int newFacesToAddToMerge = 0;

    for (auto bpair : b2b) {

        
        // Increase number of vertices in the smaller boundary
        size_t targetSize = std::max(s1.boundaries[bpair.first].size(), s2.boundaries[bpair.second].size());
        if (s1.boundaries[bpair.first].size() < targetSize) {
            addNewBoundaryVertices(s1, s1.boundaries[bpair.first], targetSize);
        } else if (s2.boundaries[bpair.second].size() < targetSize) {
            addNewBoundaryVertices(s2, s2.boundaries[bpair.second], targetSize);
        }
        

        // For each boundary vertex-pair we will add two new faces later into the merged surface
        newFacesToAddToMerge += 2*targetSize;        

        // disp(MSG_DEBUG,"New boundary vertices added");
        
        // Find closest boundary vertices
        float minDist = FLT_MAX;
        int   b1start = 0;
        int   b2start = 0;

        for (int i = 0; i < int(s1.boundaries[bpair.first].size()); i++) {
            for (int j = 0; j < int(s2.boundaries[bpair.second].size()); j++) {
                float d = dist(s1.vertices[s1.boundaries[bpair.first][i]],s2.vertices[s2.boundaries[bpair.second][j]]);
                if (d < minDist) {
                    minDist = d;
                    b1start = i;
                    b2start = j;
                }
            }    
        }

        // Rotate the boundary vertices so that the closest vertices are at the start
        std::rotate(s1.boundaries[bpair.first].begin(),  s1.boundaries[bpair.first].begin()  + b1start, s1.boundaries[bpair.first].end());
        std::rotate(s2.boundaries[bpair.second].begin(), s2.boundaries[bpair.second].begin() + b2start, s2.boundaries[bpair.second].end());

        auto calculateBoundaryDistance = [&]()->float {
            float totDist = 0.0f;
            for (size_t i = 0; i < s1.boundaries[bpair.first].size() && i < s2.boundaries[bpair.second].size(); ++i) {
                totDist += dist(s1.vertices[s1.boundaries[bpair.first][i]], s2.vertices[s2.boundaries[bpair.second][i]]);
            }
            return totDist;
        };
            
        // Flip the direction of second boundary if needed
        float origDist = calculateBoundaryDistance();

        std::rotate (s2.boundaries[bpair.second].begin(), s2.boundaries[bpair.second].begin()+1, s2.boundaries[bpair.second].end());
        std::reverse(s2.boundaries[bpair.second].begin(), s2.boundaries[bpair.second].end());

        float flipDist = calculateBoundaryDistance();

        if (origDist < flipDist) {
            std::rotate (s2.boundaries[bpair.second].begin(), s2.boundaries[bpair.second].begin()+1, s2.boundaries[bpair.second].end());
            std::reverse(s2.boundaries[bpair.second].begin(), s2.boundaries[bpair.second].end());
        }

        // disp(MSG_DEBUG,"Aligned boundaries");

    }

    // Glue the two boundaries
    Surface out = surfMerge(s1,s2);
    // disp(MSG_DEBUG,"Surfaces merged");

    
    // Allocate memory for new face list, and copy the previous faces
    int     lastFaceInd = out.nf;
    int     newNf       = out.nf + newFacesToAddToMerge; // increased faces, one faces per vertex is added
    int**   newFaces    = new int*[newNf];
    for (int n = 0; n < newNf; n++) {
        newFaces[n] = (n < out.nf) ? out.faces[n] : (new int[3]);
    }
    delete[] out.faces;
    out.nf    = newNf;
    out.faces = newFaces;

    

    for (auto bpair : b2b) {        

        for (size_t i = 0; i < s1.boundaries[bpair.first].size(); i++) {
            
            size_t currentIndex = i;
            size_t nextIndex    = (i + 1) % s1.boundaries[bpair.first].size(); // Wrap around to the start

            // Indices of vertices in s1 and s2
            int v1_s1 = s1.boundaries[bpair.first][currentIndex];
            int v2_s1 = s1.boundaries[bpair.first][nextIndex];
            int v1_s2 = s2.boundaries[bpair.second][currentIndex];
            int v2_s2 = s2.boundaries[bpair.second][nextIndex];

            // Find the ordered edge containing the boundary edge
            bool edgeOrder = false;
            for (int n = 0; n < s1.nf; ++n) {
                auto& face = s1.faces[n];
                if      ( ((face[0] == v1_s1) && (face[1] == v2_s1)) || ((face[0] == v1_s1) && (face[2] == v2_s1)) || ((face[1] == v1_s1) && (face[2] == v2_s1)) ) {edgeOrder = true;  break;}
                else if ( ((face[1] == v1_s1) && (face[0] == v2_s1)) || ((face[2] == v1_s1) && (face[0] == v2_s1)) || ((face[2] == v1_s1) && (face[1] == v2_s1)) ) {edgeOrder = false; break;}
            }

            out.faces[lastFaceInd][0] = edgeOrder ? v1_s1 : v2_s1;
            out.faces[lastFaceInd][1] = edgeOrder ? v2_s1 : v1_s1;
            out.faces[lastFaceInd][2] = v2_s2 + s1.nv;
            lastFaceInd++;

            edgeOrder = false;
            for (int n = 0; n < s2.nf; ++n) {
                auto& face = s2.faces[n];
                if      ( ((face[0] == v1_s2) && (face[1] == v2_s2)) || ((face[0] == v1_s2) && (face[2] == v2_s2)) || ((face[1] == v1_s2) && (face[2] == v2_s2)) ) {edgeOrder = true;  break;}
                else if ( ((face[1] == v1_s2) && (face[0] == v2_s2)) || ((face[2] == v1_s2) && (face[0] == v2_s2)) || ((face[2] == v1_s2) && (face[1] == v2_s2)) ) {edgeOrder = false; break;}
            }

            out.faces[lastFaceInd][0] = v1_s1;
            out.faces[lastFaceInd][1] = edgeOrder ? (v2_s2+s1.nv) : (v1_s2+s1.nv);
            out.faces[lastFaceInd][2] = edgeOrder ? (v1_s2+s1.nv) : (v2_s2+s1.nv);
            lastFaceInd++;

            
        }

    }
    
    

    disp(MSG_DEBUG,"Done surfGlueBoundaries");
    return out;

}

std::tuple<Surface,Surface> NIBR::splitClosedAndOpenParts(Surface& surf)
{
    surf.reset();
    surf.isClosed(); // internally computed connected components and checks whether they are open or closed

    Surface closedPart;
    Surface openPart;

    for (int n = 0; n < int(surf.comp.size()); n++) {

        if (surf.compOpenOrClosed[n] == CLOSED)
            closedPart  = surfMerge(closedPart,surf.comp[n]);
        else if (surf.compOpenOrClosed[n] == OPEN)
            openPart    = surfMerge(openPart,  surf.comp[n]);
        else
            disp(MSG_ERROR,"Internal error. Surface must be either open or closed.");
        
    }

    return std::make_tuple(closedPart,openPart);

}


std::tuple<Surface,Surface> NIBR::splitWithPlane(const Surface& surf, const float planePoint[3], const float planeNormal[3])
{
    disp(MSG_DEBUG,"splitWithPlane");
    // disp(MSG_DEBUG,"planePoint=[%.2f,%.2f,%.2f], planeNormal=[%.2f,%.2f,%.2f]",planePoint[0],planePoint[1],planePoint[2],planeNormal[0],planeNormal[1],planeNormal[2]);

    if (surf.nv == 0) return std::make_tuple(Surface(),Surface());

    Surface cutSurf(surf,true);
    
    cutSurf.calcNormalsOfFaces();

    float D = -dot(planePoint, planeNormal); // Plane equation constant in Ax+By+Cz+D=0

    std::vector<std::array<float, 3>> newVertices;
    std::vector<std::array<int, 3>>   newFaces;

    int newNv = surf.nv;

    auto calculateEdgePlaneIntersection = [&](float* v1, float* v2, std::pair<std::array<float, 3>, int>& intersection)->bool {

        float edge[3]     = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
        float numerator   = -dot(planeNormal, v1) - D;
        float denominator =  dot(planeNormal, edge);

        if (denominator == 0) return false; // Edge is parallel to the plane, no intersection

        float t = numerator / denominator;
        if (t >= 0.0f && t < 1.0f) { // Intersection is within the edge segment
            for (int i = 0; i < 3; i++) {
                intersection.first[i] = v1[i] + t * edge[i];
            }
            return true;
        }
        return false;

    };

    // Find the faces that intersects the plane
    for (int n = 0; n < cutSurf.nf; n++) {

        newFaces.emplace_back(std::array<int,3>{cutSurf.faces[n][0],cutSurf.faces[n][1],cutSurf.faces[n][2]});

        // Check if face intersects the plane
        int frontOfPlane = 0;
        int backOfPlane  = 0;    

        for (int i = 0; i < 3; i++) { // for each point on the face

            float p2s[3];
            vec3sub(p2s,cutSurf.vertices[cutSurf.faces[n][i]],planePoint);

            if (NIBR::dot(p2s,planeNormal) >= 0) {
                frontOfPlane++;
            } else {
                backOfPlane++;
            }

        }

        // Face intersects plane, next find intersecting points
        if ( (frontOfPlane > 0) && (backOfPlane > 0) ) {
            
            std::vector<std::pair<std::array<float, 3>, int>> intersectionPoints; // To store intersection points and the intersecting edge no

            for (int i = 0; i < 3; i++) { // Check each edge of the face
            
                int next = (i + 1) % 3;   // Get the next vertex index to form an edge

                std::pair<std::array<float, 3>, int> intersection;
                intersection.second = i;

                if (calculateEdgePlaneIntersection(cutSurf.vertices[cutSurf.faces[n][i]], cutSurf.vertices[cutSurf.faces[n][next]], intersection)) {
                    intersectionPoints.emplace_back(std::move(intersection.first), intersection.second);
                }

            }

            
            // Remove overlapping or close intersection points, for 2 point intersection case
            if ( (intersectionPoints.size() == 2) && (dist(intersectionPoints[0].first,intersectionPoints[1].first) < EPS10) ) {
                disp(MSG_DEBUG,"Intersecting face %d does not need splitting", n);
                continue; // Face does not need spliting
            }

            // Remove overlapping or close intersection points, for 3 point intersection case
            if (intersectionPoints.size() == 3) {
                if (dist(intersectionPoints[1].first,intersectionPoints[2].first) < EPS10) intersectionPoints.pop_back();
                if (dist(intersectionPoints[0].first,intersectionPoints[1].first) < EPS10) intersectionPoints.pop_back();
                if ((intersectionPoints.size() == 2) && (dist(intersectionPoints[0].first,intersectionPoints[1].first) < EPS10)) intersectionPoints.pop_back();
                if ((intersectionPoints.size() == 3) && (dist(intersectionPoints[0].first,intersectionPoints[2].first) < EPS10)) intersectionPoints.pop_back();
                if (intersectionPoints.size() < 3) {
                    disp(MSG_DEBUG,"Removed %d overlapping vertex in face %d", 3-intersectionPoints.size());
                }
            }
                        

            // Split the face if there are two distinct intersection points
            if (intersectionPoints.size() == 2) {

                // Next we need to form new faces and add them. (After all is done, we need to remove the original face.)
                // Importantly, we will form the new faces so that their normals point to the same direction as the original face.

                // One of the two intersection points can overlap with an existing vertex of the face.
                int twoSplit = -1;
                bool edgeOnPlane = false;
                
                for (int i = 0; i < 3; i++) {
                    if (dist(intersectionPoints[0].first,surf.vertices[surf.faces[n][i]]) < EPS10) {
                        intersectionPoints[0].second = i;
                        twoSplit = 0;
                        break;
                    }
                }

                for (int i = 0; i < 3; i++) {
                    if (dist(intersectionPoints[1].first,surf.vertices[surf.faces[n][i]]) < EPS10) {
                        intersectionPoints[1].second = i;
                        if (twoSplit < 0) {
                            twoSplit = 1;
                        } else {
                            edgeOnPlane = true;
                        }
                        break;
                    }
                }

                if (edgeOnPlane) {
                    disp(MSG_DEBUG,"Intersecting face %d has an edge on plane", n);
                    continue; // Face does not need spliting
                }
                

                // We remove the current face and replace it with new split faces
                newFaces.pop_back();   

                if (twoSplit >= 0) { // We will divide the face into two new faces and add a new vertex.
                    
                    float vertex[3] = {intersectionPoints[abs(twoSplit-1)].first[0],intersectionPoints[abs(twoSplit-1)].first[1],intersectionPoints[abs(twoSplit-1)].first[2]};

                    newVertices.emplace_back(std::array<float,3>{vertex[0], vertex[1], vertex[2]});

                    newFaces.emplace_back(std::array<int,3>{surf.faces[n][intersectionPoints[twoSplit].second], surf.faces[n][(intersectionPoints[twoSplit].second+1)%3], newNv});
                    newFaces.emplace_back(std::array<int,3>{surf.faces[n][intersectionPoints[twoSplit].second], newNv, surf.faces[n][(intersectionPoints[twoSplit].second+2)%3]});
                    
                    newNv++;

                } else { // We will divide the face into three new faces and add two new vertices.

                    float vertex0[3] = {intersectionPoints[0].first[0],intersectionPoints[0].first[1],intersectionPoints[0].first[2]};
                    float vertex1[3] = {intersectionPoints[1].first[0],intersectionPoints[1].first[1],intersectionPoints[1].first[2]};

                    newVertices.emplace_back(std::array<float, 3>{vertex0[0], vertex0[1], vertex0[2]});
                    newVertices.emplace_back(std::array<float, 3>{vertex1[0], vertex1[1], vertex1[2]});
                    
                    if ( (intersectionPoints[0].second+1)%3 == intersectionPoints[1].second ) {

                        newFaces.emplace_back(std::array<int,3>{newNv,     surf.faces[n][intersectionPoints[1].second],       newNv + 1});
                        newFaces.emplace_back(std::array<int,3>{newNv + 1, surf.faces[n][(intersectionPoints[1].second+1)%3], surf.faces[n][intersectionPoints[0].second]});
                        newFaces.emplace_back(std::array<int,3>{surf.faces[n][intersectionPoints[0].second], newNv, newNv + 1});

                    } else {
                        
                        newFaces.emplace_back(std::array<int,3>{newNv + 1, surf.faces[n][intersectionPoints[0].second], newNv});
                        newFaces.emplace_back(std::array<int,3>{newNv,surf.faces[n][(intersectionPoints[0].second+1)%3], newNv + 1});
                        newFaces.emplace_back(std::array<int,3>{surf.faces[n][(intersectionPoints[0].second+1)%3], surf.faces[n][intersectionPoints[1].second], newNv + 1});

                    }

                    newNv += 2;

                }     

            } else if (intersectionPoints.size() > 2){
                disp(MSG_DEBUG,"Intersecting face %d does not need splitting", n);
            }

        }

    }
    disp(MSG_DEBUG,"Processed all faces");

    // Next we add the new vertices and faces, and also remove the faces which were split
    float** V = new float*[surf.nv + newVertices.size()];
    for (int i = 0; i < (surf.nv + int(newVertices.size())); i++) {
        
        V[i] = new float[3];

        if (i < surf.nv) {
            V[i][0] = surf.vertices[i][0];
            V[i][1] = surf.vertices[i][1];
            V[i][2] = surf.vertices[i][2];
        } else {
            V[i][0] = newVertices[i-surf.nv][0];
            V[i][1] = newVertices[i-surf.nv][1];
            V[i][2] = newVertices[i-surf.nv][2];
            // disp(MSG_DEBUG,"newPoint=[%.2f,%.2f,%.2f]",V[i][0],V[i][1],V[i][2]);
        }
         
    }

    Surface tmp;
    tmp.nv       = surf.nv + newVertices.size();
    tmp.vertices = V;

    // disp(MSG_DEBUG,"Created new vertex array");

    
    int** F = new int*[newFaces.size()];
    for (size_t i = 0; i < (newFaces.size()) ; i++) {
        F[i]    = new int[3];
        F[i][0] = newFaces[i][0];
        F[i][1] = newFaces[i][1];
        F[i][2] = newFaces[i][2];
    }
    tmp.nf    = newFaces.size();
    tmp.faces = F;

    // disp(MSG_DEBUG,"Created new face array");

    tmp = surfRepair(tmp);

    // Finally we split the surface into two
    std::vector<bool> frontMask(tmp.nv,false);
    std::vector<bool> behindMask(tmp.nv,false);

    for (int n = 0; n < tmp.nv; n++) {

        float p2s[3];
        vec3sub(p2s,tmp.vertices[n],planePoint);

        float side = NIBR::dot(p2s,planeNormal);

        if (side >= 0) frontMask[n]  = true;
        if (side <= 0) behindMask[n] = true;

    }

    auto out = std::make_tuple(surfRepair(applyMask(tmp,frontMask)),surfRepair(applyMask(tmp,behindMask)));

    disp(MSG_DEBUG,"Done splitWithPlane");

    return out;

}


void NIBR::surfaceFieldThreshold(Surface *surf, SurfaceField* field, float loVal, float hiVal) {

    if (field->fdata != NULL) {
        for (int n=0; n<surf->nv; n++) {
            if ( (field->fdata[n][0] >= loVal) && (field->fdata[n][0] <= hiVal) )
                field->fdata[n][0] = 1;
            else 
                field->fdata[n][0] = 0;
        }
    }

}

void NIBR::surfaceFieldThreshold(Surface *surf, SurfaceField* field, int loVal, int hiVal) {

    if (field->idata != NULL) {
        for (int n=0; n<surf->nv; n++) {
            if ( (field->idata[n][0] >= loVal) && (field->idata[n][0] <= hiVal) )
                field->idata[n][0] = 1;
            else 
                field->idata[n][0] = 0;
        }
    }

}


SurfaceField NIBR::makeDiscMask(Surface* surf, float x, float y, float z, float r) {

    if (surf->nv == 0) {return SurfaceField();}

    SurfaceField out;
    out.owner     = VERTEX;
    out.datatype  = "int";
    out.dimension = 1;
    out.fdata     = NULL;
    out.idata     = new int*[surf->nv];
    for (int n=0; n<surf->nv; n++) {
        out.idata[n]    = new int[1];
        out.idata[n][0] = 0;
    }

    std::set<int> pivot_set;

    auto V = surf->vertices;
    auto F = surf->faces;

    int centralVert = -1;
    for (int v=0; v<surf->nv; v++) {
        if (abs(V[v][0]-x) + abs(V[v][1]-y) + abs(V[v][2]-z) < EPS3) {
            centralVert = v;
            break;
        }
    }

    if (centralVert == -1) {
        return out;
    } else {

        surf->getNeighboringFaces();
        std::set<int> vertsHandled;
        std::set<int> vertsToHandle;
        vertsToHandle.insert(centralVert);
        
        while (!vertsToHandle.empty()) {
            int v = *vertsToHandle.begin();
            for (int f : surf->neighboringFaces[v]) {
                if (pivot_set.find(f) != pivot_set.end()) continue;

                bool addFace = true;
                // Iterate through the vertices of the face f
                for (int i=0; i<3; i++) {    
                    float* vert = V[F[f][i]];                 
                    float sqDist = (vert[0]-x)*(vert[0]-x) + (vert[1]-y)*(vert[1]-y) + (vert[2]-z)*(vert[2]-z);
                    
                    if (sqDist <= r*r) { 
                        if (vertsHandled.find(F[f][i]) == vertsHandled.end()) {
                            vertsToHandle.insert(F[f][i]);
                        }
                    } else {
                        addFace = false;
                    }
                }
                if (addFace) {
                    pivot_set.insert(f);
                }
            }
            vertsToHandle.erase(v);
            vertsHandled.insert(v);
        }
    }

    surf->getNeighboringVertices();
    for (auto faceList : pivot_set) {
        for (auto f : surf->neighboringVertices[faceList]) {
            out.idata[f][0] = 1;
        }
    }

    return out;

}

void NIBR::extractDisc(Surface* outSurf, Surface* surf, float x, float y, float z, float r) {    
    SurfaceField mask = makeDiscMask(surf,x,y,z,r);
    selectVertices(outSurf, surf, &mask);
    surf->clearField(mask);
}

// Modifies the input mask field
void NIBR::selectVertices(Surface* outSurf, Surface* surf, SurfaceField* mask, int label) {
    SurfaceField sel = surf->copyField(*mask);
    surfaceFieldThreshold(surf, &sel, label, label);
    return selectVertices(outSurf, surf, &sel);
}

// Modifies the input mask field
void NIBR::removeVertices(Surface* outSurf, Surface* surf, SurfaceField* mask, int label) {
    surfaceFieldThreshold(surf, mask, label, label);
    return removeVertices(outSurf, surf, mask);
}

void NIBR::selectVertices(Surface* outSurf, Surface* surf, SurfaceField* mask) {

    SurfaceField rem = surf->copyField(*mask);
    surf->convert2VertField(rem);

    for (int n=0; n<surf->nv; n++) {
        rem.idata[n][0] = (rem.idata[n][0] == 0) ? 1 : 0;
    }

    removeVertices(outSurf, surf, &rem);
    surf->clearField(rem);

}

void NIBR::removeVertices(Surface* outSurf, Surface* surf, SurfaceField* mask) {

    bool* includeVertex = new bool[surf->nv];
    int* matchingVertex = new int[surf->nv];

    std::vector<int> ids;
    int nv = 0;

    for (int n=0; n<surf->nv; n++) {
       
        if (int(mask->idata[n][0]) == 0) {
            includeVertex[n]  = true;
            matchingVertex[n] = nv;
            ids.push_back(n);
            nv++;
        } else {
            includeVertex[n] = false;
        }
        
    }
    
    bool* includeFace = new bool[surf->nf];
    std::vector<int> fids;
    int nf = 0;

    for (int n=0; n<surf->nf; n++) {
        if (includeVertex[surf->faces[n][0]] && includeVertex[surf->faces[n][1]] && includeVertex[surf->faces[n][2]]) {
            includeFace[n] = true;
            fids.push_back(n);
            nf++;
        }
        else
            includeFace[n] = false;
    }
   
    
    if (nv>0) {
        
        outSurf->nv         = nv;
        outSurf->vertices   = new float*[outSurf->nv];
        outSurf->nf         = nf;
        outSurf->faces      = new int*[outSurf->nf];
        
        for (int n=0; n<outSurf->nv; n++) {
            outSurf->vertices[n]    = new float[3];
            outSurf->vertices[n][0] = surf->vertices[ids[n]][0];
            outSurf->vertices[n][1] = surf->vertices[ids[n]][1];
            outSurf->vertices[n][2] = surf->vertices[ids[n]][2];
        }
        
        nf = 0;
        for (int n=0; n<surf->nf; n++) {
            if (includeFace[n]) {
                outSurf->faces[nf]    = new int[3];
                outSurf->faces[nf][0] = matchingVertex[surf->faces[n][0]];
                outSurf->faces[nf][1] = matchingVertex[surf->faces[n][1]];
                outSurf->faces[nf][2] = matchingVertex[surf->faces[n][2]];
                nf++;
            }
        }
        
        
        
        for (size_t f=0; f<surf->fields.size(); f++) {
            
            float** fdata = NULL;
            int**   idata = NULL;
            
            if (surf->fields[f].owner==VERTEX) {
                
                if (surf->fields[f].datatype=="float") {
                    fdata = new float*[outSurf->nv];
                    for (int n=0; n<outSurf->nv; n++) {
                        fdata[n] = new float[surf->fields[f].dimension];
                        for (int d=0; d<surf->fields[f].dimension; d++) {
                            fdata[n][d] = surf->fields[f].fdata[ids[n]][d];
                        }
                    }
                }
                
                if (surf->fields[f].datatype=="int") {
                    idata = new int*[outSurf->nv];
                    for (int n=0; n<outSurf->nv; n++) {
                        idata[n] = new int[surf->fields[f].dimension];
                        for (int d=0; d<surf->fields[f].dimension; d++) {
                            idata[n][d] = surf->fields[f].idata[ids[n]][d];
                        }
                    }
                }
                
            }
            
            if (surf->fields[f].owner==FACE) {
                
                if (surf->fields[f].datatype=="float") {
                    fdata = new float*[outSurf->nf];
                    for (int n=0; n<outSurf->nf; n++) {
                        fdata[n] = new float[surf->fields[f].dimension];
                        for (int d=0; d<surf->fields[f].dimension; d++) {
                            fdata[n][d] = surf->fields[f].fdata[fids[n]][d];
                        }
                    }
                }
                
                
                if (surf->fields[f].datatype=="int") {
                    idata = new int*[outSurf->nf];
                    for (int n=0; n<outSurf->nf; n++) {
                        idata[n] = new int[surf->fields[f].dimension];
                        for (int d=0; d<surf->fields[f].dimension; d++) {
                            idata[n][d] = surf->fields[f].idata[fids[n]][d];
                        }
                    }
                }
                
            }
            
            SurfaceField F = {surf->fields[f].owner,surf->fields[f].name,surf->fields[f].datatype,surf->fields[f].dimension,fdata,idata};
            outSurf->fields.push_back(F);
        }
        
    }

    delete[] includeVertex;
    delete[] includeFace;
    delete[] matchingVertex;

}

Surface NIBR::applyMask(const Surface& surf, std::vector<bool>& mask) {

    // disp(MSG_DEBUG,"applyMask");

    if (surf.nv == 0) {return Surface();}

    Surface outSurf;

    bool* includeVertex = new bool[surf.nv];
    int* matchingVertex = new int[surf.nv];

    std::vector<int> ids;
    int nv = 0;

    for (int n=0; n<surf.nv; n++) {
       
        if (mask[n] == 1) {
            includeVertex[n]  = true;
            matchingVertex[n] = nv;
            ids.push_back(n);
            nv++;
        } else {
            includeVertex[n] = false;
        }
        
    }
    
    bool* includeFace = new bool[surf.nf];
    std::vector<int> fids;
    int nf = 0;

    for (int n=0; n<surf.nf; n++) {
        if (includeVertex[surf.faces[n][0]] && includeVertex[surf.faces[n][1]] && includeVertex[surf.faces[n][2]]) {
            includeFace[n] = true;
            fids.push_back(n);
            nf++;
        }
        else
            includeFace[n] = false;
    }
   
    
    if ( (nv>0) && (nf>0) ) {
        
        outSurf.nv         = nv;
        outSurf.vertices   = new float*[outSurf.nv];
        outSurf.nf         = nf;
        outSurf.faces      = new int*[outSurf.nf];
        
        for (int n=0; n<outSurf.nv; n++) {
            outSurf.vertices[n]    = new float[3];
            outSurf.vertices[n][0] = surf.vertices[ids[n]][0];
            outSurf.vertices[n][1] = surf.vertices[ids[n]][1];
            outSurf.vertices[n][2] = surf.vertices[ids[n]][2];
        }
        
        nf = 0;
        for (int n=0; n<surf.nf; n++) {
            if (includeFace[n]) {
                outSurf.faces[nf]    = new int[3];
                outSurf.faces[nf][0] = matchingVertex[surf.faces[n][0]];
                outSurf.faces[nf][1] = matchingVertex[surf.faces[n][1]];
                outSurf.faces[nf][2] = matchingVertex[surf.faces[n][2]];
                nf++;
            }
        }
        
        
        
        for (size_t f=0; f<surf.fields.size(); f++) {
            
            float** fdata = NULL;
            int**   idata = NULL;
            
            if (surf.fields[f].owner==VERTEX) {
                
                if (surf.fields[f].datatype=="float") {
                    fdata = new float*[outSurf.nv];
                    for (int n=0; n<outSurf.nv; n++) {
                        fdata[n] = new float[surf.fields[f].dimension];
                        for (int d=0; d<surf.fields[f].dimension; d++) {
                            fdata[n][d] = surf.fields[f].fdata[ids[n]][d];
                        }
                    }
                }
                
                if (surf.fields[f].datatype=="int") {
                    idata = new int*[outSurf.nv];
                    for (int n=0; n<outSurf.nv; n++) {
                        idata[n] = new int[surf.fields[f].dimension];
                        for (int d=0; d<surf.fields[f].dimension; d++) {
                            idata[n][d] = surf.fields[f].idata[ids[n]][d];
                        }
                    }
                }
                
            }
            
            if (surf.fields[f].owner==FACE) {
                
                if (surf.fields[f].datatype=="float") {
                    fdata = new float*[outSurf.nf];
                    for (int n=0; n<outSurf.nf; n++) {
                        fdata[n] = new float[surf.fields[f].dimension];
                        for (int d=0; d<surf.fields[f].dimension; d++) {
                            fdata[n][d] = surf.fields[f].fdata[fids[n]][d];
                        }
                    }
                }
                
                
                if (surf.fields[f].datatype=="int") {
                    idata = new int*[outSurf.nf];
                    for (int n=0; n<outSurf.nf; n++) {
                        idata[n] = new int[surf.fields[f].dimension];
                        for (int d=0; d<surf.fields[f].dimension; d++) {
                            idata[n][d] = surf.fields[f].idata[fids[n]][d];
                        }
                    }
                }
                
            }
            
            SurfaceField F = {surf.fields[f].owner,surf.fields[f].name,surf.fields[f].datatype,surf.fields[f].dimension,fdata,idata};
            outSurf.fields.push_back(F);
        }
        
    }

    delete[] includeVertex;
    delete[] includeFace;
    delete[] matchingVertex;


    // disp(MSG_DEBUG,"Done");

    return outSurf;

}

// shifts the vertices of s1 towards s2. shift is between 0 and 1. 0 is on s1, 1 is on s2.
bool NIBR::shiftVerticesBetweenSurfaces(Surface* out, Surface* s1, Surface* s2, float shift) 
{

    if (s1->vertices==NULL) if (!s1->readMesh()) return false;
    if (s2->vertices==NULL) if (!s2->readMesh()) return false;

    if ((s1->nv != s2->nv) || (s1->nf != s2->nf)) return false;

    int nv = s1->nv;
    int nf = s1->nf;

    out->nv = nv;
    out->nf = nf;

    out->faces  = new int*[nf];
    for (int i=0; i<nf; i++) {
        out->faces[i] = new int[3];
        memcpy(out->faces[i], s1->faces[i], 3*sizeof(int));
    }

    float v[3];
    out->vertices = new float*[nv];
    for (int i=0; i<nv; i++) {
        out->vertices[i] = new float[3];
        vec3sub(v,s2->vertices[i],s1->vertices[i]);
        vec3add(out->vertices[i],s1->vertices[i],v,shift);
    }

    return true;

}

SurfaceField NIBR::convert2FaceField(Surface* surf, SurfaceField* field) {

    SurfaceField out;

    if (field->owner == NOTSET)
        return out;

    out.name      = field->name;
    out.owner     = FACE;
    out.datatype  = field->datatype;
    out.dimension = field->dimension;

    if (out.datatype=="int") {
        out.idata = new int*[surf->nf];
        for (int n=0; n<surf->nf; n++) {
            out.idata[n] = new int[out.dimension];
        }
    }

    if (out.datatype=="float") {
        out.fdata = new float*[surf->nf];
        for (int n=0; n<surf->nf; n++) {
            out.fdata[n] = new float[out.dimension];
        }
    }

    if (field->owner == FACE) {

        if (out.datatype=="int") {
            for (int n=0; n<surf->nf; n++)
                memcpy(out.idata[n], field->idata[n], out.dimension*sizeof(int));
        }

        if (out.datatype=="float") {
            for (int n=0; n<surf->nf; n++)
                memcpy(out.fdata[n], field->fdata[n], out.dimension*sizeof(float));
        }

    } else {

        if (out.datatype=="int") {
            for (int n=0; n<surf->nf; n++) {
                for (int d=0; d<out.dimension; d++) {

                    float val = 0.0f;

                    val += field->idata[surf->faces[n][0]][d];
                    val += field->idata[surf->faces[n][1]][d];
                    val += field->idata[surf->faces[n][2]][d];
                    
                    out.idata[n][d] = int(val / 3.0f);
                }
            }
        }

        if (out.datatype=="float") {
            for (int n=0; n<surf->nf; n++) {
                for (int d=0; d<out.dimension; d++) {

                    float val = 0.0f;

                    val += field->fdata[surf->faces[n][0]][d];
                    val += field->fdata[surf->faces[n][1]][d];
                    val += field->fdata[surf->faces[n][2]][d];
                
                    out.fdata[n][d] = float(val / 3.0f);
                }
            }
        }

    }

    return out;
}


SurfaceField NIBR::convert2VertField(Surface* surf, SurfaceField* field) {

    SurfaceField out;

    if (field->owner == NOTSET)
        return out;

    out.name      = field->name;
    out.owner     = VERTEX;
    out.datatype  = field->datatype;
    out.dimension = field->dimension;

    if (out.datatype=="int") {
        out.idata = new int*[surf->nv];
        for (int n=0; n<surf->nv; n++) {
            out.idata[n] = new int[out.dimension];
        }
    }

    if (out.datatype=="float") {
        out.fdata = new float*[surf->nv];
        for (int n=0; n<surf->nv; n++) {
            out.fdata[n] = new float[out.dimension];
        }
    }

    if (field->owner == VERTEX) {

        if (out.datatype=="int") {
            for (int n=0; n<surf->nv; n++)
                memcpy(out.idata[n], field->idata[n], out.dimension*sizeof(int));
        }

        if (out.datatype=="float") {
            for (int n=0; n<surf->nv; n++)
                memcpy(out.fdata[n], field->fdata[n], out.dimension*sizeof(float));
        }

    } else {

        surf->getNeighboringFaces();

        if (out.datatype=="int") {
            for (int n=0; n<surf->nv; n++) {
                for (int d=0; d<out.dimension; d++) {
                    int val = 0;
                    for (int v : surf->neighboringFaces[n]) {
                        val += field->idata[v][d];
                    }
                    val /= surf->neighboringFaces[n].size();
                    out.idata[n][d] = val;
                }
            }
        }

        if (out.datatype=="float") {
            for (int n=0; n<surf->nv; n++) {
                for (int d=0; d<out.dimension; d++) {
                    float val = 0;
                    for (int v : surf->neighboringFaces[n]) {
                        val += field->fdata[v][d];
                    }
                    val /= float(surf->neighboringFaces[n].size());
                    out.fdata[n][d] = val;
                }
            }
        }

    }

    return out;
}

std::unordered_map<int,float> NIBR::getVertexNeigborhoodOfAPoint(Surface* surf, float* p, int faceNo, float radius)
{
    surf->getNeighboringVertices();

    std::unordered_map<int,float> vertexNeighborhood;
    std::unordered_set<int> visitedVertices; // To store visited vertices, preventing from revisit.
    std::queue<int> vertexQueue; // A queue for BFS implementation.

    float rr = radius*radius;

    // Initially, add vertices of the provided face to the queue and mark them as visited.
    for(int i = 0; i < 3; i++){
        int vtxIndex = surf->faces[faceNo][i];
        vertexQueue.push(vtxIndex);
        visitedVertices.insert(vtxIndex);
    }

    // Iterate vertices in the queue. While the queue is not empty, more vertices within the distance may still be found.
    while (!vertexQueue.empty()){
        int currentVertex = vertexQueue.front();
        vertexQueue.pop();
        float d = squared_dist(p, surf->vertices[currentVertex]);

        if(d <= rr) {
            // If the distance is less than the radius, add it to the neighborhood
            vertexNeighborhood[currentVertex] = std::sqrt(d);

            // Check all neighboring vertices
            for(auto& neighborIndex : surf->neighboringVertices[currentVertex]){
                if(visitedVertices.insert(neighborIndex).second){
                    // If the neighboring vertex hasn't been visited yet, push it to the queue
                    vertexQueue.push(neighborIndex);
                }
            }
        }
    }

    return vertexNeighborhood;


}

