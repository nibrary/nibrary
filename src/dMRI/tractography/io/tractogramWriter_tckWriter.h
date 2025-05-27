#pragma once

#include <string>
#include <vector>
#include "tractogramWriter.h"

namespace NIBR
{

    class TCKWriter : public IBatchWriter {
    public:
        TCKWriter(std::string _filename);
        ~TCKWriter();
        bool open() override;
        bool writeBatch(const std::vector<std::vector<std::vector<float>>>& batch) override;
        bool close() override;
        const std::string& getFilename() const override { return filename; }
    private:
        FILE* file = nullptr;
        std::string filename;
        size_t      totalStreamlineCount = 0;
        long        posCount = 0;
        long        posFileOffset = 0;
        long        dataStartOffset = 0;
    };
}