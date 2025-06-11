#pragma once

#include <string>
#include <vector>
#include <cstdio>               // For FILE*
#include "tractogramWriter.h"   // For IBatchWriter, StreamlineBatch
#include "base/byteSwapper.h"   // For is_little_endian, swapByteOrder

namespace NIBR
{
    class TCKWriter : public IBatchWriter {
    public:
        TCKWriter(std::string _filename);
        ~TCKWriter() override;

        TCKWriter(const TCKWriter&) = delete;
        TCKWriter& operator=(const TCKWriter&) = delete;

        bool open() override;
        bool writeBatch(const StreamlineBatch& batch) override;
        bool close(long& finalStreamlineCount, long& finalPointCount) override;
        const std::string& getFilename() const override { return filename_; }

    private:
        FILE* file_ = nullptr;
        std::string filename_;
        size_t      currentStreamlineCount_ = 0;
        size_t      currentPointCount_      = 0;
        long        posCount_ = 0;          // File position for 'count' header field
        long        posFileOffset_ = 0;     // File position for 'file: . OFFSET' header field
        long        dataStartOffset_ = 0;   // Actual file offset where streamline data begins
        bool        needsByteSwap_ = false; // To handle endianness
    };
}
