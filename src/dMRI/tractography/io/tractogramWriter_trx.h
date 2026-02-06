#pragma once

#include <string>
#include <vector>
#include "tractogramWriter.h"   // For IBatchWriter, StreamlineBatch

namespace NIBR
{
    class TRXWriter : public IBatchWriter {
    public:
        TRXWriter(std::string _filename);
        ~TRXWriter() override;

        TRXWriter(const TRXWriter&) = delete;
        TRXWriter& operator=(const TRXWriter&) = delete;

        bool open() override;
        bool writeBatch(const StreamlineBatch& batch) override;
        bool close(long& finalStreamlineCount, long& finalPointCount) override;
        const std::string& getFilename() const override { return filename_; }

    private:
        std::string             filename_;
        StreamlineBatch         streamlines_;
        std::size_t             currentStreamlineCount_ = 0;
        std::size_t             currentPointCount_      = 0;
        bool                    is_open_ = false;
    };
}
