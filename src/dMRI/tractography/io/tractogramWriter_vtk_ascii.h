#pragma once

#include <string>
#include <vector>
#include <cstdio>               // For FILE*
#include "tractogramWriter.h"   // For IBatchWriter, StreamlineBatch, TRACTOGRAMFILEFORMAT
#include "tractogramField.h"    // For TractogramField

namespace NIBR
{
    class VTKAsciiWriter : public IBatchWriter {
    public:
        VTKAsciiWriter(std::string _filename);
        ~VTKAsciiWriter() override;

        VTKAsciiWriter(const VTKAsciiWriter&) = delete;
        VTKAsciiWriter& operator=(const VTKAsciiWriter&) = delete;

        void setVTKFields(const std::vector<TractogramField>& fields) override;

        bool open() override;
        bool writeBatch(const StreamlineBatch& batch) override;
        bool close(long& finalStreamlineCount, long& finalPointCount) override;

        const std::string& getFilename() const override { return filename_; }

    private:
        FILE* mainFile_                         = nullptr;
        std::string filename_;
        std::vector<TractogramField> fields_;
        bool fieldsSet_                         = false;

        // Temporary file handling
        FILE* tempPointsFile_                   = nullptr;
        std::string tempPointsFilename_;
        FILE* tempLinesFile_                    = nullptr;
        std::string tempLinesFilename_;


        size_t currentTotalPointCount_          = 0;
        size_t currentStreamlineCount_          = 0;
        size_t globalPointIndexOffset_          = 0;

        // Storing streamline lengths is still needed if field data is written from full arrays
        std::vector<int> streamlineLengths_;

        bool openTemporaryFiles();
        void closeTemporaryFiles(bool deleteFiles);
        bool appendAsciiFileContent(FILE* dest, const std::string& srcFilename);
    };
}
