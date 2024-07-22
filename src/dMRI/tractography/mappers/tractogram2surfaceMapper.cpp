#include "dMRI/tractography/mappers/tractogram2surfaceMapper.h"
#include "math/gaussian.h"
#include "surface/surface_operators.h"

using namespace NIBR;

// mapping contains streamline2faceMaps for each face of the input surface
// if mapOnce is true, then a face will not have two instances from the same streamline
void NIBR::tractogram2surfaceMapper(NIBR::TractogramReader* _tractogram, NIBR::Surface* surf, std::vector<std::vector<NIBR::streamline2faceMap>>& mapping, bool mapOnce)
{

    surf->calcCentersOfFaces();

    // Make copies of tractogram for multithreader
    NIBR::TractogramReader* tractogram = new NIBR::TractogramReader[NIBR::MT::MAXNUMBEROFTHREADS()]();
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++)
        tractogram[t].copyFrom(*_tractogram);

    std::vector<std::vector<std::vector<std::vector<int>>>> surfaceGrid;

    // Index surface and create mask
    NIBR::Image<int> img;
    bool*** mask = indexSurfaceBoundary(surf,&img,&surfaceGrid,false);

    // Map tractogram on surface
    std::vector<std::vector<std::vector<streamline2faceMap>>> map2surf(NIBR::MT::MAXNUMBEROFTHREADS());
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        map2surf[t].resize(surf->nf);
    }

    auto doMapping = [&](NIBR::MT::TASK task)->void {

        int streamlineId = task.no;
        int threadNo     = task.threadId;

        // If streamline does not have a segment, then exit
        if (tractogram[threadNo].len[streamlineId]<2) return;
    
        double p0[3], p1[3], dir[3], t, length;
        
        int32_t A[3], B[3];
        
        float** streamline = tractogram[threadNo].readStreamline(streamlineId);
        
        NIBR::LineSegment seg;
        seg.id = streamlineId;

        auto addToMap=[&]()->void{
            for (auto f : surfaceGrid[A[0]][A[1]][A[2]]) {

                double pointOfIntersection[3];
                double distanceToIntersection;
                float angle = findSegmentTriangleIntersection(surf, f, seg.beg, seg.end, &pointOfIntersection[0], &distanceToIntersection);

                if (angle>0) {
                    streamline2faceMap tmp;
                    tmp.index  = seg.id;
                    tmp.dir[0] = seg.dir[0];
                    tmp.dir[1] = seg.dir[1];
                    tmp.dir[2] = seg.dir[2];
                    tmp.p[0]   = pointOfIntersection[0];
                    tmp.p[1]   = pointOfIntersection[1];
                    tmp.p[2]   = pointOfIntersection[2];
                    tmp.angle  = angle;
                    map2surf[threadNo][f].push_back(tmp);
                }

            }
        };

        // Beginning of segment in real and image space
        img.to_ijk(streamline[0],p0);
        A[0]  = std::round(p0[0]);
        A[1]  = std::round(p0[1]);
        A[2]  = std::round(p0[2]);

        // If streamline has many points and segments
        for (uint32_t i=0; i<tractogram[threadNo].len[streamlineId]-1; i++) {

            // End of segment in real and image space
            img.to_ijk(streamline[i+1],p1);
            seg.beg = streamline[i];
            seg.end = streamline[i+1];
            for (int m=0;m<3;m++) {
                seg.dir[m] = streamline[i+1][m] - streamline[i][m];
                B[m]       = std::round(p1[m]);
            }
            
            // Find segment length and direction in real space
            seg.len = norm(seg.dir);
            vec3scale(seg.dir,1.0/seg.len);
            
            // Add the first voxel in the map if it is within the mask
            if ( img.isInside(A) && mask[A[0]][A[1]][A[2]] ) {
                addToMap();
            }

            // Segment does not leave the voxel, add this segment and continue with the next one
            if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
                for (int m=0;m<3;m++) {
                    p0[m] = p1[m];
                }
                continue;
            }
            
            // Find length and direction of segment in grid space
            vec3sub(dir,p1,p0);
            length = norm(dir);
            vec3scale(dir,1.0/length);

            while (length>0.0) {

                if (rayTraceVoxel(A,p0,dir,t)) {
                    t += EPS4;
                } else {
                    t  = EPS4; // otherwise t is NAN
                }

                for (int m=0;m<3;m++) {
                    p0[m] += t*dir[m];
                    A[m]   = std::round(p0[m]);
                }

                if ( img.isInside(A) && mask[A[0]][A[1]][A[2]] ) {
                    addToMap();
                }

                length -= t;

            }

            for (int m=0;m<3;m++) {
                p0[m] = p1[m];
                A[m]  = B[m];
            }
            

        }
        
        
        for (uint32_t i=0; i<tractogram[threadNo].len[streamlineId]; i++)
            delete[] streamline[i];
        delete[] streamline;
        

    };
    NIBR::MT::MTRUN(tractogram[0].numberOfStreamlines, NIBR::MT::MAXNUMBEROFTHREADS(), "Tractogram to surface mapping", doMapping);

    // Clean up 
    for (int i = 0; i < img.imgDims[0]; i++) {
        for (int j = 0; j < img.imgDims[1]; j++) {
            delete[] mask[i][j];
        }
        delete[] mask[i];
    }
    delete[] mask;

    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        tractogram[t].destroyCopy();
    }
    delete[] tractogram;

    auto finMapping = [&](NIBR::MT::TASK task)->void {  

        int f = task.no; 

        for (int t = 1; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
            if (!map2surf[t][f].empty()) {
                map2surf[0][f].insert(map2surf[0][f].end(), map2surf[t][f].begin(), map2surf[t][f].end());
                map2surf[t][f].clear();
            }
        }
        
        if (mapOnce) {
            if (!map2surf[0][f].empty()) {
                std::sort(map2surf[0][f].begin(), map2surf[0][f].end(),[](streamline2faceMap s1,streamline2faceMap s2){return (s1.index  < s2.index) ? 1 : 0;});
                auto it = std::unique (map2surf[0][f].begin(), map2surf[0][f].end(),[](streamline2faceMap s1,streamline2faceMap s2){return (s1.index == s2.index) ? 1 : 0;});
                map2surf[0][f].erase(it, map2surf[0][f].end());
            }
        }


    };
    NIBR::MT::MTRUN(surf->nf, NIBR::MT::MAXNUMBEROFTHREADS(), finMapping);

    mapping = map2surf[0];

}
