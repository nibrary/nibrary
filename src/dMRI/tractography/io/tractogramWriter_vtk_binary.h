#pragma once

#include <string>
#include <vector>
#include <cstdio>               // For FILE*
#include "tractogramWriter.h"   // For IBatchWriter, StreamlineBatch, TRACTOGRAMFILEFORMAT
#include "tractogramField.h"    // For TractogramField

namespace NIBR
{
    class VTKBinaryWriter : public IBatchWriter {
    public:
        VTKBinaryWriter(std::string _filename);
        ~VTKBinaryWriter() override;

        VTKBinaryWriter(const VTKBinaryWriter&) = delete;
        VTKBinaryWriter& operator=(const VTKBinaryWriter&) = delete;

        void setVTKFields(const std::vector<TractogramField>& fields) override;

        bool open() override;
        bool writeBatch(const StreamlineBatch& batch) override; 
        bool close(long& finalStreamlineCount, long& finalPointCount) override;
        
        const std::string& getFilename() const override { return filename_; }

    private:
        FILE* mainFile_ = nullptr; // Main output VTK file
        std::string filename_;
        std::vector<TractogramField> fields_; 
        bool fieldsSet_ = false;

        // Temporary file handling
        FILE* tempPointsFile_ = nullptr;
        std::string tempPointsFilename_;
        FILE* tempLinesFile_ = nullptr;
        std::string tempLinesFilename_;
        
        // Placeholder for potential future implementation for sequential writing of fields
        struct TempFieldInfo {
            std::string filename;
            FILE* file_ptr = nullptr;
        };
        std::vector<TempFieldInfo> tempFieldFiles_;

        size_t currentTotalPointCount_ = 0; 
        size_t currentStreamlineCount_ = 0; 
        size_t globalPointIndexOffset_ = 0; // Tracks the starting index for points in the current batch

        // Storing streamline lengths is needed if field data is written from full arrays
        // or to correctly iterate point data fields.
        std::vector<int> streamlineLengths_; 

        
        void writeMainHeaderPlaceholders();
        bool openTemporaryFiles();
        void closeTemporaryFiles();
        bool appendFileContent(FILE* dest, const std::string& srcFilename);
        
        bool needsByteSwap_ = true; 
    };
}
