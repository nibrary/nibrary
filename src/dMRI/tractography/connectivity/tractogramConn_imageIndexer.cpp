#include "tractogramConn_imageIndexer.h"
#include <fstream>

using namespace NIBR;

bool NIBR::SCimageIndexer::writeConn(std::string outFname) {

    std::ofstream fout;
    fout.open (outFname);

    // Write the labels
    fout << original_labels[0];
    for (size_t i=1; i<labelCnt; i++) {
        fout << ", ";
        fout << original_labels[i];
    }
    fout << std::endl;

    // Write the connectivity matrix
    for (size_t i=0; i<labelCnt; i++) {
        fout << conn[i][0].size();
        for (size_t j=1; j<labelCnt; j++) {
            fout << ", ";
            fout << conn[i][j].size();
        }
        fout << std::endl;
    }

    fout.close();

    return true;

}

// Checks if label belongs to background
bool NIBR::SCimageIndexer::isBg(int val) {
    if (bgLabels.find(val)==bgLabels.end())
        return false;
    else
        return true;
}

NIBR::SCimageIndexer::SCimageIndexer(NIBR::TractogramReader* _tractogram, NIBR::Image<int>* _img) {
    
    // Initialize tractogram
    tractogram = new NIBR::TractogramReader[NIBR::MT::MAXNUMBEROFTHREADS()]();
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++)
        tractogram[t].copyFrom(*_tractogram);

    // Initialize label image
    img = _img;
    
    // Background label is 0 for internal computation, i.e., in the input label image all values in bgLabels is set to bgVal
    // (Changing this value will break the internal computation. So let's not touch it.)
    bgVal = 0;

    // Set default parameters
    endLengthThresh = 0;
    endType         = LAST_NONE_BG_LABEL;

}

NIBR::SCimageIndexer::~SCimageIndexer() { 
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        tractogram[t].destroyCopy();
    }
    delete[] tractogram;
}


void NIBR::SCimageIndexer::run() {

    // If no bgLabel is added, add the default bgVal in the set
    if (bgLabels.empty())
        bgLabels.insert(0);

    // Find labels and create relabeled image for fast vector indexing
    std::set<int> label_set(img->data, img->data + img->voxCnt);
    labels.push_back(bgVal); // First label for internal computation is bgVal
    for (auto label: label_set) {
        if (!isBg(label)) {
            original_labels.push_back(label);
            labels.push_back(labels.size());        // Labels is then sorted as 0, 1, 2, 3,...
        }
    }

    // Convert labels in the image
    auto convertLabels = [&](NIBR::MT::TASK task)->void {
        int val = img->data[task.no];

        if (isBg(val)) {

            img->data[task.no] = bgVal;

        } else {

            for (size_t i=0; i<original_labels.size(); i++)
                if (original_labels[i] == val)
                    img->data[task.no] = labels[i+1];

        }
    };
    NIBR::MT::MTRUN(img->voxCnt,NIBR::MT::MAXNUMBEROFTHREADS(),convertLabels);

    // Create connectivity matrix
    labelCnt = original_labels.size();
    conn.resize(labelCnt, std::vector<std::set<size_t>>(labelCnt));

    NIBR::MT::MTRUN(tractogram[0].numberOfStreamlines, "Computing connectome", [&](NIBR::MT::TASK task)->void{processStreamline(task.no,task.threadId);} );
}


bool NIBR::SCimageIndexer::processStreamline(int streamlineId, uint16_t threadNo) {

    // If streamline is empty
    if (tractogram[threadNo].len[streamlineId]<2) 
        return true;

    float pi[3], pip[3], pim[3]; // p[i], p[i+1], p[i-1] - points on streamline
    double p0[3], p1[3],  dir[3], length, lengthR, lengthScale, t, endLength;

    bool stop;
    
    int A[3], B[3];

    std::vector<NIBR::Segment> end1;
    std::vector<NIBR::Segment> end2;

    auto insert2Conn = [&]()->void {
        
        int frLabel, toLabel;
        
        auto checkAndInsert = [&]()->void {
            if ( (frLabel != bgVal) && (toLabel != bgVal) ) {
                
                if (frLabel > toLabel) {
                    NIBR::MT::PROC_MX().lock();
                    conn[toLabel-1][frLabel-1].insert(streamlineId);
                    NIBR::MT::PROC_MX().unlock();
                } else {
                    NIBR::MT::PROC_MX().lock();
                    conn[frLabel-1][toLabel-1].insert(streamlineId);
                    NIBR::MT::PROC_MX().unlock();
                }
                
            }
        };


        if (endType == END_POINT_LABEL) {
            frLabel = *(int*)(end1[0].data);
            toLabel = *(int*)(end2[0].data);
            checkAndInsert();
        }

        if (endType == LONGEST_END_LABEL) {

            auto getEdgeLabel = [&](std::vector<NIBR::Segment> e)->int {

                std::sort(e.begin(), e.end(), [&](NIBR::Segment s1,NIBR::Segment s2)->bool{return (*(int*)(s1.data)) < (*(int*)(s2.data));});
            
                int   cl        = *(int*)(e[0].data);       // current label
                float cs        = 0.0;                      // current label's length
                int   edgeLabel = cl;                       // edge label, i.e., label with max label
                float ms        = 0.0;                      // max label's length

                int sl;
                for (auto s : e) {

                    sl = (*(int*)(s.data));

                    if (sl != bgVal) {

                        if ( sl == cl) {
                            cs += s.length;
                        } else {

                            if ( (ms < cs) || (edgeLabel == bgVal) ) {
                                edgeLabel  = cl;
                                ms         = cs; 
                            }

                            cl = sl;
                            cs = s.length;
                        }

                    }

                    if ( (ms < cs) || (edgeLabel == bgVal) ) {
                        edgeLabel  = cl;
                        ms         = cs; 
                    }


                }

                return edgeLabel;
            };

            frLabel = getEdgeLabel(end1);
            toLabel = getEdgeLabel(end2);
            checkAndInsert();

        }

        if (endType == FIRST_NONE_BG_LABEL) {

            auto getEdgeLabel = [&](std::vector<NIBR::Segment> e)->int {

                int edgeLabel = bgVal;

                for (auto it = e.begin(); it != e.end(); ++it) {
                    edgeLabel = (*(int*)(it->data));
                    if ( edgeLabel != bgVal)
                        return edgeLabel;
                }

                return edgeLabel;

            };

            frLabel = getEdgeLabel(end1);
            toLabel = getEdgeLabel(end2);
            checkAndInsert();

        }

        if (endType == LAST_NONE_BG_LABEL) {

            auto getEdgeLabel = [&](std::vector<NIBR::Segment> e)->int {

                int edgeLabel = bgVal;

                for (auto it = e.rbegin(); it != e.rend(); ++it) {
                    edgeLabel = (*(int*)(it->data));
                    if ( edgeLabel != bgVal)
                        return edgeLabel;
                }

                return edgeLabel;
                
            };

            frLabel = getEdgeLabel(end1);
            toLabel = getEdgeLabel(end2);
            checkAndInsert();

        }

    };

    NIBR::Segment seg;

    // Segment data here is the label computed with nearest neighbor interpolation
    auto pushToEnd = [&](int endNo)->bool {

        if ( (A[0]>-1) && (A[1]>-1) && (A[2]>-1) && A[0]<img->imgDims[0] && A[1]<img->imgDims[1] && A[2]<img->imgDims[2] )
            seg.data = img->data + img->sub2ind(A[0],A[1],A[2]);
        else
            seg.data = &bgVal;

        endLength += seg.length;
        if (endLength > endLengthThresh) {
            seg.length -= endLength - endLengthThresh;
            if (endNo==1)
                end1.push_back(seg);
            else
                end2.push_back(seg);
            stop = true;
        } else {
            if (endNo==1)
                end1.push_back(seg);
            else
                end2.push_back(seg);
            stop = false;
        }  
        return stop;
    };

    // Iterate until the whole length of segment is covered. Split the segment into subsegments if more than a single voxel is traced by the segment.
    auto splitSegment = [&](int endNo)->bool {

        // Find total segment length in image space and its direction
        vec3sub(dir,p1,p0);
        length  = norm(dir);
        vec3scale(dir,1.0f/length);
        
        // Grid lengthScale
        lengthScale = std::sqrt(dir[0]*img->pixDims[0]*dir[0]*img->pixDims[0]+dir[1]*img->pixDims[1]*dir[1]*img->pixDims[1]+dir[2]*img->pixDims[2]*dir[2]*img->pixDims[2]);

        t = 0;
        while (1) {
        
            // Segment enters the next voxel
            if (t>0) {
                
                t += EPS4;
                    
                // Calc subsegment length within the current voxel and push to end
                seg.length = t*lengthScale;
                if (lengthR < seg.length) {
                    seg.length = lengthR;
                    lengthR    = -1;
                } else {
                    lengthR   -= seg.length;
                }
                if (pushToEnd(endNo))
                    return stop;
                
                // Check if end of segment is reached
                length -= t;
                if ((length<0) || (lengthR<0)) break;
            
                // Split the segment by updating the segment beginning
                for (int m=0;m<3;m++) {
                    seg.p[m] += seg.length*seg.dir[m];
                    p0[m]    += t*dir[m];
                    A[m]      = std::round(p0[m]);
                }
                
                // If remaining part of the segment does not leave the voxel, add this part to the end and continue with the next segment
                if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
                    seg.length  = lengthR;
                    if (pushToEnd(endNo))
                        return stop;
                    break;
                }
                
            }
            
            if (rayTraceVoxel(A,p0,dir,t)) 
                continue;

            if (length<EPS4) 
                break;
            
            length -= EPS4;

            seg.length  = EPS4;
            if (pushToEnd(endNo)) return stop;
            
            for (int m=0;m<3;m++) {
                p0[m] += EPS4*dir[m];
                A[m]   = std::round(p0[m]);
            }
            
            t = 0;
        
        }

        return stop;
    };

    // PROCESS END1 : Walk from 0 -> N
    endLength   = 0;
    stop        = false;

    // Beginning of first segment and its corner in image space
    tractogram[threadNo].readPoint(streamlineId,0,pi);
    img->to_ijk(pi, p0);
    A[0] = std::round(p0[0]);
    A[1] = std::round(p0[1]);
    A[2] = std::round(p0[2]);

    if (endType == END_POINT_LABEL) {
        seg.length = 0;
        pushToEnd(1);
    } else {

        for (uint32_t i=0; i<tractogram[threadNo].len[streamlineId]-1; i++) {

            // End of segment and its corner in image space
            tractogram[threadNo].readPoint(streamlineId,i,  pi);
            tractogram[threadNo].readPoint(streamlineId,i+1,pip);
            img->to_ijk(pip, p1);
            for (int m=0;m<3;m++) {
                seg.p[m]   = pi[m];
                seg.dir[m] = pip[m] - pi[m];
                B[m]       = std::round(p1[m]);
            }

            // Find segment length and direction in real space
            lengthR = norm(seg.dir);
            vec3scale(seg.dir,1.0f/lengthR);
            
            // Segment does not leave the voxel, add this segment on the edge and continue with the next one
            if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
                seg.length  = lengthR;
                if (pushToEnd(1)) 
                    break;
                p0[0] = p1[0];
                p0[1] = p1[1];
                p0[2] = p1[2];
                continue;
            }
            
            // Segment leaves the voxel and needs to be splitted into several subsegments
            if (splitSegment(1)) 
                break;

            for (int m=0;m<3;m++) {
                p0[m] = p1[m];
                A[m]  = B[m];
            }

        }

    }

    // PROCESS END2 : Walk from N -> 0
    endLength   = 0;
    stop        = false;

    // Beginning of first segment and its corner in image space
    tractogram[threadNo].readPoint(streamlineId,tractogram[threadNo].len[streamlineId]-1, pi);
    img->to_ijk(pi, p0);
    A[0] = std::round(p0[0]);
    A[1] = std::round(p0[1]);
    A[2] = std::round(p0[2]);

    

    if (endType == END_POINT_LABEL) {
        seg.length = 0;
        pushToEnd(2);
    } else {

        for (uint32_t i=tractogram[threadNo].len[streamlineId]-1; i>0; i--) {

            // End of segment and its corner in image space
            tractogram[threadNo].readPoint(streamlineId,i,  pi);
            tractogram[threadNo].readPoint(streamlineId,i-1,pim);
            img->to_ijk(pim, p1);
            for (int m=0;m<3;m++) {
                seg.p[m]   = pi[m];
                seg.dir[m] = pim[m] - pi[m];
                B[m]       = std::round(p1[m]);
            }

            // Find segment length and direction in real space
            lengthR = norm(seg.dir);
            vec3scale(seg.dir,1.0f/lengthR);

            // Segment does not leave the voxel, add this segment on the edge and continue with the next one
            if ( (A[0]==B[0]) && (A[1]==B[1]) && (A[2]==B[2]) ) {
                seg.length  = lengthR;
                if (pushToEnd(2)) 
                    break;
                p0[0] = p1[0];
                p0[1] = p1[1];
                p0[2] = p1[2];
                continue;
            }

            // Segment leaves the voxel and needs to be splitted into several subsegments
            if (splitSegment(2)) 
                break;

            for (int m=0;m<3;m++) {
                p0[m] = p1[m];
                A[m]  = B[m];
            }
        }

    }

    insert2Conn();

    return true;

}
