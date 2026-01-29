#include <limits>
#include "base/nibr.h"      // For disp, SGNTR
#include "tractogramWriter_trx.h"
#include <trx/trx.h>

namespace NIBR {

TRXWriter::TRXWriter(std::string _filename) : filename_(std::move(_filename)) {}

TRXWriter::~TRXWriter()
{
    if (is_open_) {
        disp(MSG_WARN, "TRXWriter for %s destroyed without explicit close. Attempting to close.", filename_.c_str());
        long dummy1, dummy2;
        close(dummy1, dummy2);
    }
}

bool TRXWriter::open()
{
    streamlines_.clear();
    currentStreamlineCount_ = 0;
    currentPointCount_ = 0;
    // TRX requires knowing total sizes up front, so we buffer streamlines in memory.
    is_open_ = true;
    return true;
}

bool TRXWriter::writeBatch(const StreamlineBatch& batch)
{
    if (!is_open_) {
        disp(MSG_ERROR, "TRXWriter: File not open for writing.");
        return false;
    }

    if (batch.empty()) return true;

    streamlines_.reserve(streamlines_.size() + batch.size());

    for (const auto& streamline : batch) {
        if (streamline.empty()) {
            continue;
        }
        streamlines_.push_back(streamline);
        currentStreamlineCount_++;
        currentPointCount_ += streamline.size();
    }

    return true;
}

bool TRXWriter::close(long& finalStreamlineCount, long& finalPointCount)
{
    finalStreamlineCount = currentStreamlineCount_;
    finalPointCount      = currentPointCount_;

    if (!is_open_) {
        return true;
    }

    is_open_ = false;

    if (currentStreamlineCount_ == 0 || currentPointCount_ == 0) {
        disp(MSG_WARN, "TRXWriter: No streamlines to write for %s.", filename_.c_str());
        return true;
    }

    if (currentStreamlineCount_ > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        currentPointCount_ > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        disp(MSG_ERROR, "TRXWriter: Streamline or point count exceeds TRX limits for %s.", filename_.c_str());
        return false;
    }

    try {
        trxmmap::TrxFile<float> trx(static_cast<int>(currentPointCount_),
                                    static_cast<int>(currentStreamlineCount_),
                                    nullptr);

        if (!trx.streamlines || trx.streamlines->_offsets.size() == 0) {
            disp(MSG_ERROR, "TRXWriter: Failed to allocate TRX streamlines for %s.", filename_.c_str());
            return false;
        }

        trx.streamlines->_offsets(0, 0) = 0;

        uint64_t current_vertex = 0;
        std::size_t current_streamline = 0;

        for (const auto& streamline : streamlines_) {
            const auto tck_size = static_cast<Eigen::Index>(streamline.size());
            for (Eigen::Index i = 0; i < tck_size; ++i) {
                const auto& pos = streamline[static_cast<std::size_t>(i)];
                trx.streamlines->_data(static_cast<Eigen::Index>(current_vertex) + i, 0) = pos[0];
                trx.streamlines->_data(static_cast<Eigen::Index>(current_vertex) + i, 1) = pos[1];
                trx.streamlines->_data(static_cast<Eigen::Index>(current_vertex) + i, 2) = pos[2];
            }
            trx.streamlines->_lengths(static_cast<Eigen::Index>(current_streamline), 0) = static_cast<uint32_t>(tck_size);
            current_vertex += static_cast<uint64_t>(tck_size);
            trx.streamlines->_offsets(static_cast<Eigen::Index>(current_streamline + 1), 0) = current_vertex;
            current_streamline++;
        }

        trxmmap::save(trx, filename_, ZIP_CM_STORE);
        trx.close();
    } catch (const std::exception& e) {
        disp(MSG_ERROR, "TRXWriter: Failed to write %s. %s", filename_.c_str(), e.what());
        return false;
    }

    disp(MSG_DEBUG, "TRXWriter: Successfully closed %s. Streamlines: %ld, Points: %ld",
         filename_.c_str(), finalStreamlineCount, finalPointCount);

    return true;
}

}
