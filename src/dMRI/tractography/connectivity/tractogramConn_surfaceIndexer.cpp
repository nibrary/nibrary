#include "tractogramConn_surfaceIndexer.h"
#include <fstream>

using namespace NIBR;

bool NIBR::SCsurfaceIndexer::writeConn(std::string outFname) {

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

// Lambda function to check if the label belongs to background
bool NIBR::SCsurfaceIndexer::isBg(int val) {
    if (bgLabels.find(val)==bgLabels.end())
        return false;
    else
        return true;
}

NIBR::SCsurfaceIndexer::SCsurfaceIndexer(NIBR::TractogramReader* _tractogram, NIBR::Surface* _surf, NIBR::SurfaceField *_surfLabels) {
    
    // Initialize tractogram
    tractogram = _tractogram;

    // Initialize label image
    surf            = _surf;
    surfLabels      = _surfLabels;
    
    // Background label is 0 for internal computation, i.e., in the input label image all values in bgLabels is set to bgVal
    // (Changing this value will break the internal computation. So let's not touch it.)
    bgVal = 0;

    // Set default parameters
    endLengthThresh = 0;

}

NIBR::SCsurfaceIndexer::~SCsurfaceIndexer() { 
    return;
}


void NIBR::SCsurfaceIndexer::run() {

    // If no bgLabel is added, add the default bgVal in the set
    if (bgLabels.empty())
        bgLabels.insert(0);

    // Find labels and create relabeled surface field for fast vector indexing
    std::set<int> label_set;
    faceLabels.resize(surf->nf);

    if (surfLabels->owner == NIBR::FACE) {
        for (int i=0; i<surf->nf; i++) {
            label_set.insert(surfLabels->idata[i][0]);
            faceLabels[i] = surfLabels->idata[i][0];
        }
    }
    if (surfLabels->owner == NIBR::VERTEX) {
        surf->getNeighboringFaces();
        for (int i=0; i<surf->nv; i++) {
            label_set.insert(surfLabels->idata[i][0]);
            for (auto nef : surf->neighboringFaces[i]) {
                faceLabels[nef] = surfLabels->idata[i][0];
            }
        }
    }


    labels.push_back(bgVal); // First label for internal computation is bgVal
    for (auto label: label_set) {
        if (!isBg(label)) {
            original_labels.push_back(label);
            labels.push_back(labels.size());        // Labels is then sorted as 0, 1, 2, 3,...
        }
    }

    // Create connectivity matrix
    labelCnt = original_labels.size();
    conn.resize(labelCnt, std::vector<std::set<size_t>>(labelCnt));

    // Create tractogram to surface map
    tractogram2surfaceMapper(tractogram, surf, tract2surfMap, true);

    // Compute connectome
    NIBR::MT::MTRUN(tractogram[0].numberOfStreamlines, "Computing connectome", [&](const NIBR::MT::TASK& task)->void{processStreamline(task.no,task.threadId);} );
}


bool NIBR::SCsurfaceIndexer::processStreamline(int, uint16_t) {
    
    return true;
}
