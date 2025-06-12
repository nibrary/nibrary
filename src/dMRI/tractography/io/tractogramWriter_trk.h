#pragma once

#include <string>
#include <vector>
#include <cstdio>               // For FILE*
#include "tractogramWriter.h"   // For IBatchWriter, StreamlineBatch, TRKReferenceInfo
#include "tractogramReader.h"   // For trkFileStruct

struct TRKReferenceInfo {
    short   imgDims[3];
    float   pixDims[3];
    float   ijk2xyz[4][4];
    float   xyz2ijk[4][4]; // Will be derived if not provided
    char    voxOrdr[4];

    // Information about scalars and properties to be written in the TRK header
    short   n_scalars_trk_to_write;
    char    scalar_name_trk_to_write[10][20];
    short   n_properties_trk_to_write;
    char    property_name_trk_to_write[10][20];

    // Default constructor to initialize
    TRKReferenceInfo() {
        std::memset(imgDims, 0, sizeof(imgDims));
        std::memset(pixDims, 0, sizeof(pixDims));
        for(int i=0; i<4; ++i) for(int j=0; j<4; ++j) ijk2xyz[i][j] = (i==j) ? 1.0f : 0.0f; // Identity
        std::memset(xyz2ijk, 0, sizeof(xyz2ijk)); // To be computed
        std::strcpy(voxOrdr, "LAS"); // TODO: Let's use LAS. Could be obtained from Image<T> and nibrary has functions for it. But it is not clear how to use this information in trk format.
        n_scalars_trk_to_write = 0;
        std::memset(scalar_name_trk_to_write, 0, sizeof(scalar_name_trk_to_write));
        n_properties_trk_to_write = 0;
        std::memset(property_name_trk_to_write, 0, sizeof(property_name_trk_to_write));
    }
};

namespace NIBR
{
    class TRKWriter : public IBatchWriter {
    public:
        TRKWriter(std::string _filename);
        ~TRKWriter() override;

        TRKWriter(const TRKWriter&) = delete;
        TRKWriter& operator=(const TRKWriter&) = delete;

        void setTRKReference(const Image<bool>& ref) override;

        bool open() override;
        bool writeBatch(const StreamlineBatch& batch) override;
        bool close(long& finalStreamlineCount, long& finalPointCount) override;
        const std::string& getFilename() const override { return filename_; }

    private:
        FILE* file_ = nullptr;
        std::string filename_;
        trkFileStruct header_;      // TRK file header structure
        TRKReferenceInfo refInfo_;  // Stores reference image information
        bool refInfoSet_ = false;

        size_t currentStreamlineCount_ = 0;
        size_t currentPointCount_      = 0;

        void initializeHeader(); 
    };
}
