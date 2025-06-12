#include "tractogram_operators.h"
#include "streamline_operators.h"

using namespace NIBR;

StreamlineBatch NIBR::tractogramTransform(const StreamlineBatch& batch_in, float M[][4])
{
    int N = batch_in.size();
    
    if (N<1) {
        return StreamlineBatch();
    }

    StreamlineBatch out(N);

    auto transform = [&](const NIBR::MT::TASK& task)->void{

        const auto& streamline = batch_in[task.no];

        int len = streamline.size();

        Streamline streamline_out;
        streamline_out.reserve(len);

        for (int l=0; l<len; l++) {
            Point3D p;
            p[0] = streamline[l][0]*M[0][0] + streamline[l][1]*M[0][1] + streamline[l][2]*M[0][2] + M[0][3];
            p[1] = streamline[l][0]*M[1][0] + streamline[l][1]*M[1][1] + streamline[l][2]*M[1][2] + M[1][3];
            p[2] = streamline[l][0]*M[2][0] + streamline[l][1]*M[2][1] + streamline[l][2]*M[2][2] + M[2][3];
            streamline_out.emplace_back(p);
        }

        out[task.no] = std::move(streamline_out);

    };

    NIBR::MT::MTRUN(N, "Applying transform", transform);

    return out;
}



std::vector<float> NIBR::getTractogramBBox(NIBR::TractogramReader* tractogram) {

    std::vector<float> bb(6,0.0f);
    
    int N = tractogram->numberOfStreamlines;
    
    if (N<1) {
        return bb;
    }

    tractogram->reset();

    auto [success,firstStreamline,streamlineId] = tractogram->getNextStreamline();
    bb[0] = firstStreamline[0][0];
    bb[1] = firstStreamline[0][0];
    bb[2] = firstStreamline[0][1];
    bb[3] = firstStreamline[0][1];
    bb[4] = firstStreamline[0][2];
    bb[5] = firstStreamline[0][2];

    std::mutex modifier;

    // Iterate throught the whole tractogram
    auto findBB = [&]()->void{
        
        // Local bounding box for the current streamline
        std::vector<float> localBB(6);
        auto [success,streamline,streamlineId] = tractogram->getNextStreamline();

        localBB[0] = localBB[1] = streamline[0][0];  // x min and max
        localBB[2] = localBB[3] = streamline[0][1];  // y min and max
        localBB[4] = localBB[5] = streamline[0][2];  // z min and max

        for (std::size_t i = 1; i < streamline.size(); i++) {
            localBB[0] = std::min(localBB[0], streamline[i][0]);  // Update x min
            localBB[1] = std::max(localBB[1], streamline[i][0]);  // Update x max
            localBB[2] = std::min(localBB[2], streamline[i][1]);  // Update y min
            localBB[3] = std::max(localBB[3], streamline[i][1]);  // Update y max
            localBB[4] = std::min(localBB[4], streamline[i][2]);  // Update z min
            localBB[5] = std::max(localBB[5], streamline[i][2]);  // Update z max
        }

        // Safely merge the local bounding box into the global bounding box
        {
            std::lock_guard<std::mutex> lock(modifier);

            bb[0] = std::min(bb[0], localBB[0]);
            bb[1] = std::max(bb[1], localBB[1]);
            bb[2] = std::min(bb[2], localBB[2]);
            bb[3] = std::max(bb[3], localBB[3]);
            bb[4] = std::min(bb[4], localBB[4]);
            bb[5] = std::max(bb[5], localBB[5]);
        }
        
    };
        
    NIBR::MT::MTRUN(N-1, NIBR::MT::MAXNUMBEROFTHREADS(), "Finding tractogram bounding box", findBB);

    return bb;

}


// Returns bool vector marking true for streamlines which are in the inpBatch but not in the refBatch
std::vector<bool> NIBR::tractogramDiff(const StreamlineBatch& inpBatch, const StreamlineBatch& refBatch) {

    std::vector<bool> diff(inpBatch.size(),true); // Mark all as different
    std::mutex modifier;

    auto compare = [&](const NIBR::MT::TASK& task)->void {

        const auto& inpStreamline = inpBatch[task.no];

        for (const auto& refStreamline : refBatch) {
            if (areStreamlinesIdentical(inpStreamline, refStreamline)) {
                // Lock is needed because vector<bool> is not thread safe! 
                // Other vectors are. Because bool is written in bits for compactness.
                std::lock_guard lock(modifier);
                diff[task.no] = false;
                return;
            }
        }

    };
    NIBR::MT::MTRUN(inpBatch.size(), "Comparing streamlines", compare);

    return diff;
}

TractogramField NIBR::colorTractogram(NIBR::TractogramReader* tractogram)
{

    float*** segmentColors = new float**[tractogram->numberOfStreamlines];

    // Iterate throught the whole tractogram

    tractogram->reset();

    auto getColors = [&]()->void{
        auto [success,streamline,streamlineId] = tractogram->getNextStreamline();
        auto colors = colorStreamline(streamline);
        segmentColors[streamlineId] = new float*[colors.size()];
        for (std::size_t l = 0; l < colors.size(); l++) {
            segmentColors[streamlineId][l]    = new float[3];
            segmentColors[streamlineId][l][0] = colors[l][0];
            segmentColors[streamlineId][l][1] = colors[l][1];
            segmentColors[streamlineId][l][2] = colors[l][2];
        }
    };
    NIBR::MT::MTRUN(tractogram->numberOfStreamlines,"Computing streamline colors",getColors);

    TractogramField streamlineColors;
    streamlineColors.owner      = POINT_OWNER;
    streamlineColors.name       = "RGB";
    streamlineColors.datatype   = FLOAT32_DT;
    streamlineColors.dimension  = 3;
    streamlineColors.data       = reinterpret_cast<void*>(segmentColors);

    return streamlineColors;

}


struct StreamlineStats {
    float  length                  = 0.0f;
    double sumOfSquaredStepLengths = 0.0;
    size_t stepCount               = 0;
};

static StreamlineStats getStreamlineStats(const Streamline& streamline) {
    
    StreamlineStats stats;

    if (streamline.size() < 2) {
        return stats;
    }

    stats.stepCount = streamline.size() - 1;

    for (size_t i = 0; i < stats.stepCount; ++i) {
        const auto& p1 = streamline[i];
        const auto& p2 = streamline[i + 1];

        const float dx = p2[0] - p1[0];
        const float dy = p2[1] - p1[1];
        const float dz = p2[2] - p1[2];

        const float stepLength = std::sqrt(dx * dx + dy * dy + dz * dz);
        
        stats.length += stepLength;
        stats.sumOfSquaredStepLengths += static_cast<double>(stepLength * stepLength);
    }

    return stats;
}

std::tuple<bool,std::vector<float>, std::size_t, float, float> NIBR::getTractogramStats(NIBR::TractogramReader* tractogram)
{
    std::vector<float> streamlineLengths(tractogram->numberOfStreamlines);

    const float         corruptionLimit      = 100000.0f;
    std::atomic<bool>   hasCorruptedPoint    = false;
    std::mutex          modifier;

    double      totalStepLengthSum   = 0.0;
    double      totalStepLengthSumSq = 0.0;
    size_t      totalStepCount       = 0;
    size_t      totalNumberOfPoints  = 0;
    std::mutex  statUpdater;

    auto isPointCorrupted = [&](const auto& p)->bool {
        return std::abs(p[0]) > corruptionLimit || std::abs(p[1]) > corruptionLimit || std::abs(p[2]) > corruptionLimit;
    };

    tractogram->reset();

    auto getLength = [&]()->void{

        if (hasCorruptedPoint) return;

        auto [success,streamline,streamlineId]  = tractogram->getNextStreamline();
        
        for (const auto& p : streamline) {
            if (isPointCorrupted(p)) {
                std::lock_guard lock(modifier);
                disp(MSG_ERROR,"Streamline %d is corrupted.",streamlineId);
                hasCorruptedPoint = true;
                return;
            }
        }

        streamlineLengths[streamlineId] = getStreamlineLength(streamline);

        StreamlineStats stats = getStreamlineStats(streamline);

        if (stats.stepCount > 0) {
            std::lock_guard lock(statUpdater);
            totalStepLengthSum   += stats.length;
            totalStepLengthSumSq += stats.sumOfSquaredStepLengths;
            totalStepCount       += stats.stepCount;
            totalNumberOfPoints  += stats.stepCount + 1;
        }
    };
    NIBR::MT::MTRUN(tractogram->numberOfStreamlines, "Computing streamline lengths", getLength);

    float meanStepSize      = 0.0f;
    float stdDevStepSize    = 0.0f;

    if (totalStepCount > 0) {
        meanStepSize = static_cast<float>(totalStepLengthSum / totalStepCount);
        const double variance = (totalStepLengthSumSq / totalStepCount) - (static_cast<double>(meanStepSize) * static_cast<double>(meanStepSize));
        stdDevStepSize = (variance > 0) ? static_cast<float>(std::sqrt(variance)) : 0.0f;
    }

    return std::make_tuple(hasCorruptedPoint.load(), std::move(streamlineLengths), totalNumberOfPoints, meanStepSize, stdDevStepSize);
}

