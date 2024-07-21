#include "surface.h"
#include "surface_operators.h"
#include "math/core.h"
#include "math/triangle.h"
#include <cstring>
#include <unordered_set>

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
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-result"
#pragma GCC diagnostic ignored "-Wcomment"
#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>
#pragma GCC diagnostic pop

#else

#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>

#endif

using namespace NIBR;

void NIBR::Surface::init() {
    filePath 		    = "";
    extension           = "";
	nv 			        = 0;
 	nf 			        = 0;
    ne                  = 0;
    openOrClosed        = OPENORCLOSEDUNSET;
    manifoldOrNot       = MANIFOLDORNOT;
    vertices            = NULL;
    faces               = NULL;
    centersOfFaces      = NULL;
    normalsOfFaces      = NULL;
    areasOfFaces        = NULL;
    normalsOfVertices   = NULL;
    neighboringVertices = NULL;
    neighboringFaces    = NULL;
    triangleNormal      = NULL;
    triangleEdge0       = NULL;
    triangleEdge1       = NULL;
    triangleEdge2       = NULL;
    gaussianCurvature   = NULL;
    meanCurvature       = NULL;
    area                = NAN;
    volume              = NAN;

    comp                = std::vector<Surface>();
    compOpenOrClosed    = std::vector<OpenOrClosed>();
    compClosedAndOpen   = std::vector<Surface>();
    compArea            = std::vector<double>();
    compVolume          = std::vector<double>();
    grid                = std::vector<std::vector<std::vector<std::vector<int>>>>();
    maskAndBoundary.clear();
    maskAndBoundary.setInterpolationMethod(NEAREST);
    enabledPointCheck   = false;
    pointCheckGridRes   = 0;

    edgesCategorized    = false;
    verticesCategorized = false;

    singularVertices        = std::vector<int>();
    boundaryVertices        = std::vector<int>();
    overconnectedVertices   = std::vector<int>();
    boundaryEdges           = std::vector<std::pair<int, int>>();
    overconnectedEdges      = std::vector<std::pair<int, int>>();

    boundaries              = std::vector<std::vector<int>>();
    boundaryLengths         = std::vector<double>();
    boundaryAreas           = std::vector<double>();
}

NIBR::Surface::Surface() {
	init();
}

NIBR::Surface::Surface(std::string _filePath) {
	init();
    readHeader(_filePath);
}



NIBR::Surface::Surface(const Surface& obj) {
    init();
    copyFrom(obj);
}

NIBR::Surface::Surface(const Surface& obj, bool copyMeshOnly) {
    init();

    if (copyMeshOnly) {

        nv = obj.nv;
        nf = obj.nf;

        vertices = NULL;
        if (obj.vertices!=NULL) {
            vertices = new float*[nv];
            for (int n=0; n<nv; n++) {
                vertices[n] = new float[3];
                std::memcpy(vertices[n], obj.vertices[n], 3*sizeof(float));
            }
        }

        faces = NULL;
        if (obj.faces!=NULL) {
            faces = new int*[nf];
            for (int n=0; n<nf; n++) {
                faces[n] = new int[3];
                std::memcpy(faces[n], obj.faces[n], 3*sizeof(int));
            }
        }

    } else {
        
        copyFrom(obj);
    }
}

NIBR::Surface::Surface(std::vector<std::vector<float>>& vertexList, std::vector<std::vector<int>>& faceList) {

    init();

    filePath 		    = "";
    extension           = "";
	nv 			        = vertexList.size();
 	nf 			        = faceList.size();
    openOrClosed        = OPENORCLOSEDUNSET;
    manifoldOrNot       = MANIFOLDORNOT;

    vertices = new float*[nv];
    for (int n = 0; n < nv; n++) {
        vertices[n]    = new float[3];
        vertices[n][0] = vertexList[n][0];
        vertices[n][1] = vertexList[n][1];
        vertices[n][2] = vertexList[n][2];
    }

    faces    = new int*[nf];
    for (int n = 0; n < nf; n++) {
        faces[n]    = new int[3];
        faces[n][0] = faceList[n][0];
        faces[n][1] = faceList[n][1];
        faces[n][2] = faceList[n][2];
    }

}

NIBR::Surface& NIBR::Surface::operator=(const Surface &obj) {

    if(this != &obj) { // protect against self-assignment
        copyFrom(obj);
    }
    return *this;
}

void NIBR::Surface::copyMeshFrom(const Surface& obj)
{
    if(this != &obj) { // protect against self-assignment

        this->clear();

        nv = obj.nv;
        nf = obj.nf;

        vertices = NULL;
        if (obj.vertices!=NULL) {
            vertices = new float*[nv];
            for (int n=0; n<nv; n++) {
                vertices[n] = new float[3];
                std::memcpy(vertices[n], obj.vertices[n], 3*sizeof(float));
            }
        }

        faces = NULL;
        if (obj.faces!=NULL) {
            faces = new int*[nf];
            for (int n=0; n<nf; n++) {
                faces[n] = new int[3];
                std::memcpy(faces[n], obj.faces[n], 3*sizeof(int));
            }
        }

    } else {
        this->reset();
    }

}

void NIBR::Surface::copyFrom(const Surface& obj)
{

    clear();

    nv 			         = obj.nv;
	nf 			         = obj.nf;
    ne                   = obj.ne;
    openOrClosed         = obj.openOrClosed;
    manifoldOrNot        = obj.manifoldOrNot;
    area                 = obj.area;
    volume               = obj.volume;
    filePath 		     = obj.filePath;
    extension            = obj.extension;

    vertices = NULL;
    if (obj.vertices!=NULL) {
        vertices = new float*[nv];
        for (int n=0; n<nv; n++) {
            vertices[n] = new float[3];
            std::memcpy(vertices[n], obj.vertices[n], 3*sizeof(float));
        }
    }

    faces = NULL;
    if (obj.faces!=NULL) {
        faces = new int*[nf];
        for (int n=0; n<nf; n++) {
            faces[n] = new int[3];
            std::memcpy(faces[n], obj.faces[n], 3*sizeof(int));
        }
    }

    centersOfFaces = NULL;
    if (obj.centersOfFaces!=NULL) {
        centersOfFaces = new float*[nf];
        for (int n=0; n<nf; n++) {
            centersOfFaces[n] = new float[3];
            std::memcpy(centersOfFaces[n], obj.centersOfFaces[n], 3*sizeof(float));
        }
    }

    normalsOfFaces = NULL;
    if (obj.normalsOfFaces!=NULL) {
        normalsOfFaces = new float*[nf];
        for (int n=0; n<nf; n++) {
            normalsOfFaces[n] = new float[3];
            std::memcpy(normalsOfFaces[n], obj.normalsOfFaces[n], 3*sizeof(float));
        }
    }

    areasOfFaces = NULL;
    if (obj.areasOfFaces!=NULL) {
        areasOfFaces = new float[nf];
        std::memcpy(areasOfFaces, obj.areasOfFaces, nf*sizeof(float));
    }

    normalsOfVertices = NULL;
    if (obj.normalsOfVertices!=NULL) {
        normalsOfVertices = new float*[nv];
        for (int n=0; n<nv; n++) {
            normalsOfVertices[n] = new float[3];
            std::memcpy(normalsOfVertices[n], obj.normalsOfVertices[n], 3*sizeof(float));
        }
    }

    triangleNormal = NULL;
    if (obj.triangleNormal!=NULL) {
        triangleNormal = new float*[nf];
        for (int n=0; n<nf; n++) {
            triangleNormal[n] = new float[3];
            std::memcpy(triangleNormal[n], obj.triangleNormal[n], 3*sizeof(float));
        }
    }

    triangleEdge0 = NULL;
    if (obj.triangleEdge0!=NULL) {
        triangleEdge0 = new float*[nf];
        for (int n=0; n<nf; n++) {
            triangleEdge0[n] = new float[3];
            std::memcpy(triangleEdge0[n], obj.triangleEdge0[n], 3*sizeof(float));
        }
    }

    triangleEdge1 = NULL;
    if (obj.triangleEdge1!=NULL) {
        triangleEdge1 = new float*[nf];
        for (int n=0; n<nf; n++) {
            triangleEdge1[n] = new float[3];
            std::memcpy(triangleEdge1[n], obj.triangleEdge1[n], 3*sizeof(float));
        }
    }

    triangleEdge2 = NULL;
    if (obj.triangleEdge2!=NULL) {
        triangleEdge2 = new float*[nf];
        for (int n=0; n<nf; n++) {
            triangleEdge2[n] = new float[3];
            std::memcpy(triangleEdge2[n], obj.triangleEdge2[n], 3*sizeof(float));
        }
    }

    neighboringVertices = NULL;
    if (obj.neighboringVertices!=NULL) {
        neighboringVertices = new std::set<int>[nv];
        for (int n=0; n<nv; n++) {
            neighboringVertices[n] = obj.neighboringVertices[n];
        }
    }

    neighboringFaces = NULL;
    if (obj.neighboringFaces!=NULL) {
        neighboringFaces = new std::set<int>[nv];
        for (int n=0; n<nv; n++) {
            neighboringFaces[n] = obj.neighboringFaces[n];
        }
    }

    gaussianCurvature = NULL;
    if (obj.gaussianCurvature!=NULL) {
        gaussianCurvature = new float[nv];
        for (int n=0; n<nv; n++) {
            gaussianCurvature[n] = obj.gaussianCurvature[n];
        }
    }

    meanCurvature = NULL;
    if (obj.meanCurvature!=NULL) {
        meanCurvature = new float[nv];
        for (int n=0; n<nv; n++) {
            meanCurvature[n] = obj.meanCurvature[n];
        }
    }

    for (auto field : obj.fields) {
        fields.push_back(copyField(field));
    }

    comp                = obj.comp;
    compOpenOrClosed    = obj.compOpenOrClosed;
    compClosedAndOpen   = obj.compClosedAndOpen;
    compArea            = obj.compArea;
    compVolume          = obj.compVolume;
    grid                = obj.grid;
    maskAndBoundary     = obj.maskAndBoundary;
    maskAndBoundary.setInterpolationMethod(NEAREST);
    enabledPointCheck   = obj.enabledPointCheck;
    pointCheckGridRes   = obj.pointCheckGridRes;

    edgesCategorized        = obj.edgesCategorized;
    verticesCategorized     = obj.verticesCategorized;

    singularVertices        = obj.singularVertices;
    boundaryVertices        = obj.boundaryVertices;
    overconnectedVertices   = obj.overconnectedVertices;
    boundaryEdges           = obj.boundaryEdges;
    overconnectedEdges      = obj.overconnectedEdges;

    boundaries              = obj.boundaries;
    boundaryLengths         = obj.boundaryLengths;
    boundaryAreas           = obj.boundaryAreas;

}

NIBR::Surface::~Surface() { 
    clear();
}

void NIBR::Surface::clear() {

    reset();

    if (vertices!=NULL) {
        for (int n=0; n<nv; n++) {
            delete[] vertices[n];
        }
        delete[] vertices;
    }
    
    if (faces!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] faces[n];
        }
        delete[] faces;
    }

    vertices = NULL;
    faces    = NULL;

    nv = 0;
    nf = 0;

}

void Surface::reset() {

    // disp(MSG_DEBUG,"surface reset");

    filePath.clear();
    extension.clear();
   
    if (centersOfFaces!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] centersOfFaces[n];
        }
        delete[] centersOfFaces;
    }
    
    if (normalsOfFaces!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] normalsOfFaces[n];
        }
        delete[] normalsOfFaces;
        delete[] areasOfFaces;
    }
    
    if (normalsOfVertices!=NULL) {
        for (int n=0; n<nv; n++) {
            delete[] normalsOfVertices[n];
        }
        delete[] normalsOfVertices;
    }
    
    if (neighboringVertices!=NULL) {
        delete[] neighboringVertices;
    }
    
    if (neighboringFaces!=NULL) {
        delete[] neighboringFaces;
    }
    
    if (triangleNormal!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] triangleNormal[n];
        }
        delete[] triangleNormal;
    }
    
    if (triangleEdge0!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] triangleEdge0[n];
        }
        delete[] triangleEdge0;
    }
    
    if (triangleEdge1!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] triangleEdge1[n];
        }
        delete[] triangleEdge1;
    }
    
    if (triangleEdge2!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] triangleEdge2[n];
        }
        delete[] triangleEdge2;
    }

    if (gaussianCurvature!=NULL) {
        delete[] gaussianCurvature;
    }

    if (meanCurvature!=NULL) {
        delete[] meanCurvature;
    }
    
    for (size_t i=0; i<fields.size(); i++) {
        clearField(fields[i]);
    }

    ne                  = 0;
    openOrClosed        = OPENORCLOSEDUNSET;
    manifoldOrNot       = MANIFOLDORNOT;

    V = Eigen::MatrixXf();
    F = Eigen::MatrixXi();
    AABB_tree = igl::AABB<Eigen::MatrixXf,3>();

    centersOfFaces      = NULL;
    normalsOfFaces      = NULL;
    areasOfFaces        = NULL;
    normalsOfVertices   = NULL;
    neighboringVertices = NULL;
    neighboringFaces    = NULL;
    triangleNormal      = NULL;
    triangleEdge0       = NULL;
    triangleEdge1       = NULL;
    triangleEdge2       = NULL;
    gaussianCurvature   = NULL;
    meanCurvature       = NULL;
    area                = NAN;
    volume              = NAN;

    comp.clear();
    compOpenOrClosed.clear();
    compClosedAndOpen.clear();
    compArea.clear();
    compVolume.clear();
    grid.clear();
    maskAndBoundary.clear();
    maskAndBoundary.setInterpolationMethod(NEAREST);
    enabledPointCheck = false;
    pointCheckGridRes = 0;

    edgesCategorized    = false;
    verticesCategorized = false;

    singularVertices.clear();
    boundaryVertices.clear();
    overconnectedVertices.clear();

    boundaryEdges.clear();
    overconnectedEdges.clear();

    boundaries.clear();
    boundaryLengths.clear();
    boundaryAreas.clear();

    // disp(MSG_DEBUG,"Done surface reset");

}


void NIBR::Surface::printInfo() {

    calcArea();                      
    calcVolume();
    isManifold();
    isClosed();
    auto componentAreas       = calcAreasOfConnectedComponents();
    auto componentVolumes     = calcVolumesOfConnectedComponents();

    auto box                  = surfBbox(*this);

    categorizeEdges();
    categorizeVertices();

    if (boundaryEdges.size() != boundaryVertices.size()) {
        disp(MSG_WARN,"Surface contains non-manifold vertices/edges. Boundaries won't be computed.");
    } else {
        computeBoundaries();
    }
    
    std::vector<double> holes = surfCalcAreasOfHoles(*this);

    if (NIBR::VERBOSE()<VERBOSE_INFO)
        return;
        
    disp(MSG_INFO,"Surface info");
    std::cout << "\033[32m";
    
    std::cout << "File name: " << filePath << std::endl;

    std::vector<std::string> manifoldOrNot{" Non-manifold"," Manifold"};

    if (openOrClosed == OPEN) {
        std::cout << "Open" << manifoldOrNot[int(isManifold())] << std::endl;
    } else if (openOrClosed == CLOSED) {
        std::cout << "Closed" << manifoldOrNot[int(isManifold())] << std::endl;
    } else {
        std::cout << "Mix of open and closed meshes." << manifoldOrNot[int(isManifold())] << std::endl;
    }

    std::cout << std::left << std::setw(30) << "Bounding box (xmin - xmax): " << std::right << std::setw(25) << " [ " + to_string_with_precision(box[0]) + ", " + to_string_with_precision(box[1]) + " ] " << std::endl;
    std::cout << std::left << std::setw(30) << "Bounding box (ymin - ymax): " << std::right << std::setw(25) << " [ " + to_string_with_precision(box[2]) + ", " + to_string_with_precision(box[3]) + " ] " << std::endl;
    std::cout << std::left << std::setw(30) << "Bounding box (zmin - zmax): " << std::right << std::setw(25) << " [ " + to_string_with_precision(box[4]) + ", " + to_string_with_precision(box[5]) + " ] " << std::endl;
                    
    std::cout << std::left << std::setw(30) << "Vertex count: " << std::right << std::setw(10) << nv << std::endl;
    std::cout << std::left << std::setw(30) << "Face count: "   << std::right << std::setw(10) << nf << std::endl;
    std::cout << std::left << std::setw(30) << "Average face area: "   << std::right << std::setw(10) << ((nf>0) ? (area/nf) : 0) << std::endl;

    std::cout << std::left << std::setw(30) << "Isolated vertex count: " << std::right << std::setw(10) << singularVertices.size() << std::endl;
    std::cout << std::left << std::setw(30) << "Boundary vertex count: " << std::right << std::setw(10) << boundaryVertices.size() << std::endl;
    std::cout << std::left << std::setw(30) << "Overconnected vertex count: " << std::right << std::setw(10) << overconnectedVertices.size() << std::endl;

    std::cout << std::left << std::setw(30) << "Boundary edge count: " << std::right << std::setw(10) << boundaryEdges.size() << std::endl;
    std::cout << std::left << std::setw(30) << "Overconnected edge count: " << std::right << std::setw(10) << overconnectedEdges.size() << std::endl;


    std::cout << std::left << std::setw(30) << "Euler number: "         << std::right << std::setw(10) << nv + nf - ne << std::endl;
    std::cout << std::left << std::setw(30) << "Connected components: " << std::right << std::setw(10) << comp.size() << std::endl;
    std::cout << std::left << std::setw(30) << "Number of boundaries: " << std::right << std::setw(10) << boundaries.size() << std::endl;
    std::cout << std::left << std::setw(30) << "Number of holes: "      << std::right << std::setw(10) << holes.size() << std::endl;
    
    std::cout << std::left << std::setw(30) << "Total volume / area:" << std::right << std::setw(25) << to_string_with_precision(volume) + " / " + to_string_with_precision(area) << std::endl;
    for (size_t n = 0; n < comp.size(); ++n) {
        std::cout << std::left << std::setw(30) << "   Volume / Area (#" + std::to_string(n+1) + "): " << std::right << std::setw(25) << to_string_with_precision(componentVolumes[n]) + " / " + to_string_with_precision(componentAreas[n])  << std::endl;
    }

    double totBoundaryLength = 0;
    double totBoundaryArea   = 0;
    for (auto& l : boundaryLengths) {totBoundaryLength += l;}
    for (auto& a : boundaryAreas  ) {totBoundaryArea   += a;}
    std::cout << std::left << std::setw(30) << "Total boundary len. / area:  " << std::right << std::setw(25) << to_string_with_precision(totBoundaryLength) + " / " + to_string_with_precision(totBoundaryArea) <<std::endl;

    for (size_t n = 0; n < boundaryLengths.size(); ++n) {
        std::cout << std::left << std::setw(30) << "   Length / Area (#" + std::to_string(n+1) + "):    " << std::right << std::setw(25) << to_string_with_precision(boundaryLengths[n]) + " / " + to_string_with_precision(boundaryAreas[n]) << std::endl;
    }

    double totHoleArea = 0;
    for (auto& a : holes) {totHoleArea += a;}
    std::cout << std::left << std::setw(30) << "Total hole area:  " << std::right << std::setw(25) << std::fixed << std::setprecision(4) << totHoleArea <<std::endl;

    for (size_t n = 0; n < holes.size(); ++n) {
        std::cout << std::left << std::setw(30) << "   Area (#" + std::to_string(n+1) + "):    " << std::right << std::setw(25) << std::fixed << std::setprecision(4) << holes[n] << std::endl;
    }

    std::cout << "\033[0m";
}

void NIBR::Surface::toEigen() {

    // disp(MSG_DEBUG,"toEigen()");

    if (V.size()!=0) {
        // disp(MSG_DEBUG,"Done toEigen()");
        return;
    }
        
    // Convert vertices to Eigen matrix
    V.resize(nv, 3);
    for (int i = 0; i < nv; i++) {
        V(i, 0) = vertices[i][0];
        V(i, 1) = vertices[i][1];
        V(i, 2) = vertices[i][2];
    }

    // Convert faces to Eigen matrix
    F.resize(nf, 3);
    for (int i = 0; i < nf; i++) {
        F(i, 0) = faces[i][0];
        F(i, 1) = faces[i][1];
        F(i, 2) = faces[i][2];
    }

    // disp(MSG_DEBUG,"Done toEigen()");
}

void NIBR::Surface::prepIglAABBTree() {

    // disp(MSG_DEBUG,"prepIglAABBTree()");

    if (fwn_bvh.F.size()!=0) {
        // disp(MSG_DEBUG,"Done prepIglAABBTree()");
        return;
    }
    
    toEigen();

    igl::fast_winding_number(V,F,2,fwn_bvh); // 2 is the accuracy_scale

    AABB_tree.init(V,F);
    isClosed(); // needed for distance functions

    // disp(MSG_DEBUG,"Done prepIglAABBTree()");
}

void NIBR::Surface::calcCentersOfFaces() {

    if (centersOfFaces!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] centersOfFaces[n];
        }
        delete[] centersOfFaces;
    }
    
    centersOfFaces = new float*[nf];
    for (int n=0; n<nf; n++) {
        centersOfFaces[n] = new float[3];
        for (int i=0; i<3; i++) {
            centersOfFaces[n][i] = (vertices[faces[n][0]][i] + vertices[faces[n][1]][i] + vertices[faces[n][2]][i]) / 3.0;
        }
    }
    
}

void NIBR::Surface::calcNormalsOfFaces() {
    
    float  v1[3];
    float  v2[3];
    float* p1;
    float* p2;
    float* p3;

    if (normalsOfFaces!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] normalsOfFaces[n];
        }
        delete[] normalsOfFaces;
        delete[] areasOfFaces;
    }
    
    
    normalsOfFaces = new float*[nf];
    areasOfFaces   = new float[nf];
    area           = 0;
    
    for (int n=0; n<nf; n++) {
        
        p1 = vertices[faces[n][0]];
        p2 = vertices[faces[n][1]];
        p3 = vertices[faces[n][2]];
        
        for (int i=0; i<3; i++) {
            v1[i] = p2[i] - p1[i];
            v2[i] = p3[i] - p1[i]; 
        }
        
        normalsOfFaces[n] = new float[3];
        
        cross(normalsOfFaces[n],v1,v2);
        areasOfFaces[n] = norm(normalsOfFaces[n])/2.0;
        normalize(normalsOfFaces[n]);
        area += areasOfFaces[n];
    }
    
}

void NIBR::Surface::flipNormalsOfFaces() {
    for (int n=0; n<nf; n++) {
        std::swap(faces[n][0],faces[n][2]);
    }
}


void NIBR::Surface::calcAreasOfFaces() {
    if (areasOfFaces==NULL)
        calcNormalsOfFaces();
}


void NIBR::Surface::getNeighboringVertices() {

    // disp(MSG_DEBUG,"getNeighboringVertices()");
    if (neighboringVertices!=NULL) return;
    
    neighboringVertices = new std::set<int>[nv];
    
    for (int n=0; n<nf; n++) {
        neighboringVertices[faces[n][0]].insert(faces[n][1]);
        neighboringVertices[faces[n][0]].insert(faces[n][2]);
        neighboringVertices[faces[n][1]].insert(faces[n][0]);
        neighboringVertices[faces[n][1]].insert(faces[n][2]);
        neighboringVertices[faces[n][2]].insert(faces[n][0]);
        neighboringVertices[faces[n][2]].insert(faces[n][1]);
    }
    // disp(MSG_DEBUG,"Done getNeighboringVertices()");
}

void NIBR::Surface::getNeighboringFaces() {

    if (neighboringFaces!=NULL) return;
    
    neighboringFaces = new std::set<int>[nv];
    
    for (int n=0; n<nf; n++) {
        neighboringFaces[faces[n][0]].insert(n);
        neighboringFaces[faces[n][1]].insert(n);
        neighboringFaces[faces[n][2]].insert(n);
    }
    
}

void NIBR::Surface::calcNormalsOfVertices() {
    
    if (normalsOfVertices != NULL) return;
    if (normalsOfFaces    == NULL) calcNormalsOfFaces();
    if (neighboringFaces  == NULL) getNeighboringFaces();
    
    normalsOfVertices = new float*[nv];
    
    for (int n=0; n<nv; n++) {
        
        normalsOfVertices[n]    = new float[3];
        
        for (std::set<int>::iterator it=neighboringFaces[n].begin(); it!=neighboringFaces[n].end(); ++it) {
            normalsOfVertices[n][0] += normalsOfFaces[*it][0];
            normalsOfVertices[n][1] += normalsOfFaces[*it][1];
            normalsOfVertices[n][2] += normalsOfFaces[*it][2];
        }
        
        float s = neighboringFaces[n].size();
        if (s>0) {
            normalsOfVertices[n][0] /= s;
            normalsOfVertices[n][1] /= s;
            normalsOfVertices[n][2] /= s;
        }
        
        normalize(normalsOfVertices[n]);
        
    }
    
}

void NIBR::Surface::applyAffineTransform(float** affT) {
    
    float p[3];
    for (int n=0; n<nv; n++) {
        
        p[0] = vertices[n][0];
        p[1] = vertices[n][1];
        p[2] = vertices[n][2];
        
        
        vertices[n][0] = p[0]*affT[0][0] + p[1]*affT[0][1] + p[2]*affT[0][2] + affT[0][3];
        vertices[n][1] = p[0]*affT[1][0] + p[1]*affT[1][1] + p[2]*affT[1][2] + affT[1][3];
        vertices[n][2] = p[0]*affT[2][0] + p[1]*affT[2][1] + p[2]*affT[2][2] + affT[2][3];
        
        
    }
    
}

void NIBR::Surface::calcVolume() {

    // disp(MSG_DEBUG,"Computing volume");

    volume = 0.0;

    calcVolumesOfConnectedComponents();

    for (auto& v : compVolume)
        volume += v;
   
}

void NIBR::Surface::calcArea() {
    if (areasOfFaces==NULL)
        calcNormalsOfFaces();
}

void NIBR::Surface::calcGaussianCurvature() {
    toEigen();
    Eigen::VectorXd K;
    igl::gaussian_curvature(V,F,K);

    if (gaussianCurvature!=NULL) delete[] gaussianCurvature;

    gaussianCurvature = new float[nv];

    for (int n = 0; n < nv; n++) {
        gaussianCurvature[n] = K(n);
    }

}


void NIBR::Surface::calcMeanCurvature() {

    toEigen();
    Eigen::VectorXd H;

    // Compute curvature directions via quadric fitting
    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;
    igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);

    H = 0.5*(PV1+PV2); // Mean curvature

    if (meanCurvature!=NULL) delete[] meanCurvature;

    meanCurvature = new float[nv];

    for (int n = 0; n < nv; n++) {
        meanCurvature[n] = H(n);
    }

}



void NIBR::Surface::calcTriangleVectors() {
    
    float* p0;
    float* p1;
    float* p2;

    if (triangleNormal!=NULL) {
        for (int n=0; n<nf; n++) {
            delete[] triangleNormal[n];
            delete[] triangleEdge0[n];
            delete[] triangleEdge1[n];
            delete[] triangleEdge2[n];
        }
        delete[] triangleNormal;
        delete[] triangleEdge0;
        delete[] triangleEdge1;
        delete[] triangleEdge2;
    }
    
    triangleNormal = new float*[nf];
    triangleEdge0  = new float*[nf];
    triangleEdge1  = new float*[nf];
    triangleEdge2  = new float*[nf];
    
    for (int n=0; n<nf; n++) {
        
        triangleNormal[n] = new float[3];
        triangleEdge0[n]  = new float[3];
        triangleEdge1[n]  = new float[3];
        triangleEdge2[n]  = new float[3];
        
        p0 = vertices[faces[n][0]];
        p1 = vertices[faces[n][1]];
        p2 = vertices[faces[n][2]];
        
        for (int i=0; i<3; i++) {
            triangleEdge0[n][i] = p1[i] - p2[i];
            triangleEdge1[n][i] = p1[i] - p0[i];
            triangleEdge2[n][i] = p2[i] - p0[i];
        }
        
        cross(triangleNormal[n],triangleEdge1[n],triangleEdge2[n]);
    }

    // disp(MSG_DEBUG,"Triangle vectors computed");
}

struct edge_pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2>& pair) const {
        return std::hash<T1>{}(pair.first) ^ (std::hash<T2>{}(pair.second) << 1);
    }
};

void Surface::categorizeEdges() {

    // disp(MSG_DEBUG,"Running categorizeEdges().");

    // printInfo();

    if (edgesCategorized) {
        // disp(MSG_DEBUG,"edges already categorized");
        return;
    }

    if (nv==0) {
        edgesCategorized = true;
        // disp(MSG_DEBUG,"nv = 0");
        return;
    }

    // Temporary map to store edge-to-face connectivity
    std::unordered_map<std::pair<int, int>, std::vector<int>, edge_pair_hash> edgeToFaceMap;

    // Clear existing data in sets
    boundaryEdges.clear();
    overconnectedEdges.clear();
    // manifoldEdges.clear();

    // disp(MSG_DEBUG,"Computing edgeToFaceMap");

    // Populate edgeToFaceMap with edges and their associated faces
    for (int n = 0; n < nf; ++n) {
        int* face = faces[n];
        std::pair<int, int> edges[3] = {
            {std::min(face[0], face[1]), std::max(face[0], face[1])},
            {std::min(face[1], face[2]), std::max(face[1], face[2])},
            {std::min(face[2], face[0]), std::max(face[2], face[0])}
        };

        for (const auto& edge : edges) {
            edgeToFaceMap[edge].push_back(n);
        }
    }

    ne = edgeToFaceMap.size(); // Number of unique edges is used for computing the Euler number later elsewhere

    // disp(MSG_DEBUG,"nv = %d nf = %d ne = %d", nv, nf, ne);

    // Categorize edges based on their connectivity
    for (const auto& pair : edgeToFaceMap) {
        const auto& edge            = pair.first;
        const auto& connectedFaces  = pair.second;

        if (connectedFaces.size() == 1) {
            boundaryEdges.push_back(edge); // Edge is a boundary edge
        } else if (connectedFaces.size() > 2) {
            overconnectedEdges.push_back(edge); // Edge is an overconnected edge
        } 
        // else if (connectedFaces.size() == 2) {
        //     manifoldEdges.push_back(edge); // Edge is a manifold edge
        // }
    }

    edgesCategorized = true;

    // disp(MSG_DEBUG,"boundaryEdges:      %d",boundaryEdges.size());
    // disp(MSG_DEBUG,"overconnectedEdges: %d",overconnectedEdges.size());
}

void Surface::categorizeVertices() {

    // disp(MSG_DEBUG,"categorizeVertices");

    // printInfo();
    
    categorizeEdges();

    if (verticesCategorized==true) {
        // disp(MSG_DEBUG,"Done. vertices already categorized");
        return;
    }

    if (nv==0) {
        verticesCategorized = true;
        // disp(MSG_DEBUG,"Done. nv = 0");
        return;
    }

    // Resetting categorization vectors
    singularVertices.clear();
    boundaryVertices.clear();
    overconnectedVertices.clear();
    // manifoldVertices.clear();

    // Use sets for quick lookup
    std::unordered_set<int> boundaryVertexSet;
    std::unordered_set<int> overconnectedVertexSet;
    std::unordered_set<int> surfaceVertexSet;

    // Process boundary and overconnected edges to categorize vertices
    for (const auto& edge : boundaryEdges) {
        boundaryVertexSet.insert(edge.first);
        boundaryVertexSet.insert(edge.second);
    }
    for (const auto& edge : overconnectedEdges) {
        overconnectedVertexSet.insert(edge.first);
        overconnectedVertexSet.insert(edge.second);
    }
    
    for (int i = 0; i < nf; ++i) {
        surfaceVertexSet.insert(faces[i][0]);
        surfaceVertexSet.insert(faces[i][1]);
        surfaceVertexSet.insert(faces[i][2]);
    }

    // Assuming nv is the number of vertices and vertices is accessible
    for (int i = 0; i < nv; ++i) {

        if (surfaceVertexSet.find(i) == surfaceVertexSet.end() ) {
            singularVertices.push_back(i);
        } else if (boundaryVertexSet.find(i) != boundaryVertexSet.end()) {
            boundaryVertices.push_back(i);
        } else if (overconnectedVertexSet.find(i) != overconnectedVertexSet.end()) {
            overconnectedVertices.push_back(i);
        } 
        // else {
        //     manifoldVertices.push_back(i);
        // }
    }

    // disp(MSG_DEBUG,"singularVertices:      %d",singularVertices.size());
    // disp(MSG_DEBUG,"boundaryVertices:      %d",boundaryVertices.size());
    // disp(MSG_DEBUG,"overconnectedVertices: %d",overconnectedVertices.size());

    verticesCategorized = true;

}

void Surface::computeBoundaries() {

    // disp(MSG_DEBUG, "computeBoundaries");

    if (!boundaries.empty())
        return;
    
    categorizeVertices();

    std::unordered_map<int, std::unordered_set<int>> adjacencyList;
    for (const auto& edge : boundaryEdges) {
        adjacencyList[edge.first].insert(edge.second);
        adjacencyList[edge.second].insert(edge.first);
    }

    std::unordered_set<int> visitedVertices;
    boundaries.clear();

    for (const auto& edge : boundaryEdges) {

        if (visitedVertices.find(edge.first) == visitedVertices.end()) {

            std::vector<int> orderedVertices;
            int currentVertex = edge.first;
            orderedVertices.push_back(currentVertex);
            visitedVertices.insert(currentVertex);

            while (true) {

                int nextVertex = *adjacencyList[currentVertex].begin();
                
                if (!orderedVertices.empty() && nextVertex == orderedVertices[0]) {
                    break; // Completed a loop
                }

                orderedVertices.push_back(nextVertex);
                visitedVertices.insert(nextVertex);

                adjacencyList[currentVertex].erase(nextVertex);
                adjacencyList[nextVertex].erase(currentVertex);

                currentVertex = nextVertex;
            }

            boundaries.push_back(orderedVertices);
        }
    }

    // Find boundary lengths and areas
    for (auto& boundary : boundaries) {

        // disp(MSG_DEBUG, "Processing boundary %d with size %d", b++, boundary.size());

        double length = 0;

        for (size_t j = 0; j < boundary.size(); j++) {
            size_t nextIndex = (j + 1) % boundary.size();
            length += dist(vertices[boundary[j]], vertices[boundary[nextIndex]]);
        }
        boundaryLengths.push_back(length);

        double area   = 0;

        for (size_t j = 1; j < (boundary.size()-1); j++) {
            area += areaOfTriangle(vertices[boundary[0]], vertices[boundary[j]], vertices[boundary[j+1]]);
        }
        boundaryAreas.push_back(area);

    }

    // disp(MSG_DEBUG, "Done computeBoundaries");
}
