#include "tractogram_operators.h"
#include "streamline_operators.h"
#include "dMRI/tractography/io/tractogramWriter.h"
#include <cstdint>
#include <tuple>
#include <chrono>
#include <atomic>
#include "base/vectorOperations.h"

using namespace NIBR;


std::vector<std::vector<std::vector<float>>> NIBR::applyTransform(NIBR::TractogramReader* _tractogram, float M[][4]) {

    // Initialize tractogram and make copies for multithreader
    NIBR::TractogramReader* tractogram = new NIBR::TractogramReader[NIBR::MT::MAXNUMBEROFTHREADS()]();
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        tractogram[t].copyFrom(*_tractogram);
    }
    
    int N = tractogram[0].numberOfStreamlines;
    
    if (N<1) {
        return std::vector<std::vector<std::vector<float>>>();
    }

    auto out = std::vector<std::vector<std::vector<float>>>();
    out.resize(N);

    // Iterate throught the whole tractogram
    auto transform = [&](NIBR::MT::TASK task)->void{

        float** streamline = tractogram[task.threadId].readStreamline(task.no);

        int len = tractogram[task.threadId].len[task.no];

        out[task.no].reserve(len);

        for (int l=0; l<len; l++) {
            auto p = std::vector<float>(3);
            p[0] = streamline[l][0]*M[0][0] + streamline[l][1]*M[0][1] + streamline[l][2]*M[0][2] + M[0][3];
            p[1] = streamline[l][0]*M[1][0] + streamline[l][1]*M[1][1] + streamline[l][2]*M[1][2] + M[1][3];
            p[2] = streamline[l][0]*M[2][0] + streamline[l][1]*M[2][1] + streamline[l][2]*M[2][2] + M[2][3];
            out[task.no].push_back(p);
            delete[] streamline[l];
        }
        delete[] streamline;

    };

    NIBR::MT::MTRUN(N, "Applying transform", transform);

    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        tractogram[t].destroyCopy();
    }
    delete[] tractogram;

    return out;

}



std::vector<float> NIBR::getTractogramBBox(NIBR::TractogramReader* _tractogram) {

    std::vector<float> bb(6,0); 

    // Initialize tractogram and make copies for multithreader
    NIBR::TractogramReader* tractogram = new NIBR::TractogramReader[NIBR::MT::MAXNUMBEROFTHREADS()]();
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        tractogram[t].copyFrom(*_tractogram);
    }
    
    int N = tractogram[0].numberOfStreamlines;
    
    if (N<1) {
        return bb;
    }

    // Initialize bb using the first streamline
    float** firstStreamline = tractogram[0].readStreamline(0);
    bb[0] = firstStreamline[0][0];
    bb[1] = firstStreamline[0][0];
    bb[2] = firstStreamline[0][1];
    bb[3] = firstStreamline[0][1];
    bb[4] = firstStreamline[0][2];
    bb[5] = firstStreamline[0][2];

    for (uint32_t i=0; i<tractogram[0].len[0]; i++)
        delete[] firstStreamline[i];
    delete[] firstStreamline;

    // Iterate throught the whole tractogram
    auto findBB = [&](NIBR::MT::TASK task)->void{
        
        // Local bounding box for the current streamline
        std::vector<float> localBB(6);
        float** streamline = tractogram[task.threadId].readStreamline(task.no);
        localBB[0] = localBB[1] = streamline[0][0];  // x min and max
        localBB[2] = localBB[3] = streamline[0][1];  // y min and max
        localBB[4] = localBB[5] = streamline[0][2];  // z min and max

        for (uint32_t i = 1; i < tractogram[task.threadId].len[task.no]; i++) {
            localBB[0] = std::min(localBB[0], streamline[i][0]);  // Update x min
            localBB[1] = std::max(localBB[1], streamline[i][0]);  // Update x max
            localBB[2] = std::min(localBB[2], streamline[i][1]);  // Update y min
            localBB[3] = std::max(localBB[3], streamline[i][1]);  // Update y max
            localBB[4] = std::min(localBB[4], streamline[i][2]);  // Update z min
            localBB[5] = std::max(localBB[5], streamline[i][2]);  // Update z max
        }

        for (uint32_t i=0; i<tractogram[task.threadId].len[task.no]; i++)
            delete[] streamline[i];
        delete[] streamline;

        // Safely merge the local bounding box into the global bounding box
        {
            std::lock_guard<std::mutex> lock(MT::PROC_MX());

            bb[0] = std::min(bb[0], localBB[0]);
            bb[1] = std::max(bb[1], localBB[1]);
            bb[2] = std::min(bb[2], localBB[2]);
            bb[3] = std::max(bb[3], localBB[3]);
            bb[4] = std::min(bb[4], localBB[4]);
            bb[5] = std::max(bb[5], localBB[5]);
        }
        
    };
        
    NIBR::MT::MTRUN(N, NIBR::MT::MAXNUMBEROFTHREADS(), "Finding tractogram bounding box", findBB);

    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        tractogram[t].destroyCopy();
    }
    delete[] tractogram;

    return bb;

}

// tuple<diffStreamlineIdx,sameStreamlineIdx>
std::tuple<std::vector<size_t>,std::vector<size_t>> NIBR::tractogramDiff(NIBR::TractogramReader* inp_tractogram, NIBR::TractogramReader* ref_tractogram)
{

    std::vector<size_t> diffStreamlineIdx;
    std::vector<size_t> sameStreamlineIdx;

    int N = ref_tractogram->numberOfStreamlines;

    if (N==0) {
        for (size_t n=0; n<inp_tractogram->numberOfStreamlines; n++)
            diffStreamlineIdx.push_back(n);
        return std::make_tuple(diffStreamlineIdx,sameStreamlineIdx);
    }

    std::vector<std::vector<std::vector<float>>> ref = ref_tractogram->read();

    int bakMaxThreads = NIBR::MT::MAXNUMBEROFTHREADS();
	if (int(inp_tractogram->numberOfStreamlines)<NIBR::MT::MAXNUMBEROFTHREADS())
		NIBR::MT::MAXNUMBEROFTHREADS() = inp_tractogram->numberOfStreamlines;

	NIBR::TractogramReader* inp = new NIBR::TractogramReader[NIBR::MT::MAXNUMBEROFTHREADS()]();
    for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++) {
        inp[t].copyFrom(*inp_tractogram);
    }

	auto compare = [&](NIBR::MT::TASK task)->void {

		float**  streamline = NULL;
        bool     isRead     = false;
        uint32_t len        = inp_tractogram->len[task.no];

        bool     isSame     = false;

        int      n = 0;
        uint32_t l = 0;

        while ( (isSame == false) && (n<N) ) {

            if (ref_tractogram->len[n] == len) {

                if (isRead==false) {
                    streamline = inp[task.threadId].readStreamline(task.no);
                    isRead     = true;
                }
                
                for (l=0; l<len; l++) {
                    if ( (streamline[l][0]!=ref[n][l][0]) || (streamline[l][1]!=ref[n][l][1]) || (streamline[l][2]!=ref[n][l][2]) )
                        break;
                }

                isSame = (l<len) ? false : true;

            }

            n++;
        }

        if (isRead) {
            for (l=0; l<len; l++)
                delete[] streamline[l];
            delete[] streamline;
        }

        if (isSame) {
            MT::PROC_MX().lock();
            sameStreamlineIdx.push_back(task.no);
            MT::PROC_MX().unlock();
        } else {
            MT::PROC_MX().lock();
            diffStreamlineIdx.push_back(task.no);
            MT::PROC_MX().unlock();
        }


    };

	if (VERBOSE()>VERBOSE_INFO)
		NIBR::MT::MTRUN(inp_tractogram->numberOfStreamlines, "Comparing streamlines", compare);
	else
		NIBR::MT::MTRUN(inp_tractogram->numberOfStreamlines, compare);
	
	for (int t = 0; t < NIBR::MT::MAXNUMBEROFTHREADS(); t++)
        inp[t].destroyCopy();
	delete[] inp;

    NIBR::MT::MAXNUMBEROFTHREADS() = bakMaxThreads;

    return std::make_tuple(diffStreamlineIdx,sameStreamlineIdx);

}

std::vector<std::vector<std::vector<float>>> NIBR::tractogramMerge(NIBR::TractogramReader* inp1, NIBR::TractogramReader* inp2, bool checkDuplicates)
{

    std::vector<std::vector<std::vector<float>>> trk1 = inp1->read();
    std::vector<std::vector<std::vector<float>>> trk2 = inp2->read();

    if (!checkDuplicates) {
        trk1.reserve(trk1.size() + trk2.size());
        std::move(trk2.begin(), trk2.end(), std::back_inserter(trk1));
        return trk1;
    }

    std::vector<size_t> existsInInp1;

	auto compare = [&](NIBR::MT::TASK task)->void {

		const auto& streamline  = trk2[task.no];

        uint32_t len = streamline.size();
        uint32_t l   = 0;

        for (const auto& s : trk1) {

            if (s.size() == len) {
                
                for (l=0; l<len; l++) {
                    if ( (streamline[l][0]!=s[l][0]) || (streamline[l][1]!=s[l][1]) || (streamline[l][2]!=s[l][2]) )
                        break;
                }

                if (l==len) {
                    MT::PROC_MX().lock();
                    existsInInp1.push_back(task.no);
                    MT::PROC_MX().unlock();
                    return;
                }

            }

        }

    };

	if (VERBOSE()>VERBOSE_INFO)
		NIBR::MT::MTRUN(trk2.size(), "Comparing streamlines", compare);
	else
		NIBR::MT::MTRUN(trk2.size(), compare);
	
	removeIdx(trk2, existsInInp1);

    trk1.reserve(trk1.size() + trk2.size());
    std::move(trk2.begin(), trk2.end(), std::back_inserter(trk1));

    return trk1;

}


TractogramField NIBR::colorTractogram(NIBR::TractogramReader* tractogram)
{

    float*** segmentColors = new float**[tractogram->numberOfStreamlines];

    // Iterate throught the whole tractogram
    auto getColors = [&](NIBR::MT::TASK task)->void{
        auto streamline = tractogram->readStreamlineVector(task.no);
        segmentColors[task.no] = colorStreamline(streamline);
    };

    if (VERBOSE() >= VERBOSE_INFO) {
        NIBR::MT::MTRUN(tractogram->numberOfStreamlines, 1, "Computing segment colors", getColors);
    } else {
        NIBR::MT::MTRUN(tractogram->numberOfStreamlines, getColors);
    }

    TractogramField streamlineColors;
    streamlineColors.owner      = POINT_OWNER;
    streamlineColors.name       = "RGB";
    streamlineColors.datatype   = FLOAT32_DT;
    streamlineColors.dimension  = 3;
    streamlineColors.data       = reinterpret_cast<void*>(segmentColors);

    return streamlineColors;

}
