#pragma once

#include <string>
#include <vector>
#include "tractogramWriter.h"

namespace NIBR
{

    class VTKBinaryWriter : public IBatchWriter {
    public:
        VTKBinaryWriter(std::string _filename);
        ~VTKBinaryWriter();
        bool open() override;
        bool writeBatch(const std::vector<std::vector<std::vector<float>>>& batch) override;
        bool close() override;
        const std::string& getFilename() const override { return filename; }
    private:
        FILE* file = nullptr;
        std::string         filename;
        std::vector<int>    lengths;
        size_t              totalPointCount = 0;
        size_t              totalStreamlineCount = 0;
        long                posPointsHeader = 0;
        long                posLinesHeader = 0;
        void writeHeaderPlaceholders();
        void writeFinalHeaders();
        void writeLinesData();
    };
    
}