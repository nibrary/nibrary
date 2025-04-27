#pragma once

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdint>
#include <limits>
#include <vector>
#include <stdexcept>
#include <mutex>

#include "base/nibr.h"
#include "math/core.h"
#include "image/image.h"


namespace NIBR
{

    struct trkFileStruct {

        char            id_string[6];
        short           dim[3];
        float           voxel_size[3];
        float           origin[3];
        short           n_scalars;
        char            scalar_name[10][20];
        short           n_properties;
        char            property_name[10][20];
        float           vox_to_ras[4][4];
        char            reserved[444];
        char            voxel_order[4];
        char            pad2[4];
        float           image_orientation_patient[6];
        char            pad1[2];
        unsigned char   invert_x;
        unsigned char   invert_y;
        unsigned char   invert_z;
        unsigned char   swap_xy;
        unsigned char   swap_yz;
        unsigned char   swap_zx;
        int             n_count;
        int             version;
        int             hdr_size;

    };

    struct Segment {
        int    streamlineNo;
        float  p[3];
        float  dir[3];
        float  length;
        void*  data;
    };

    typedef enum {
        VTK_ASCII,
        VTK_BINARY,
        TCK,
        TRK
    } TRACTOGRAMFILEFORMAT;

    struct TractogramField;

    class TractogramReader {
        
    public:
        
        TractogramReader(const std::string& _fileName, std::size_t _preloadCount = 100000);
        ~TractogramReader();

        // Copying not allowed because file handle, cache, mutex cannot be triviallycopied 
        TractogramReader(const TractogramReader& obj) = delete;
        TractogramReader& operator=(const TractogramReader&) = delete;

        bool isOpen() {return file != NULL;}
        void printInfo();

        template <typename T>
        void setReferenceImage(const Image<T>* ref);
        
        std::vector<std::vector<std::vector<float>>>    read(); // Reads all the streamlines        
        float**                                         readStreamline(std::size_t n);
        std::vector<std::vector<float>>                 readStreamlineVector(std::size_t n);
        std::vector<Point>                              readStreamlinePoints(std::size_t n);

        FILE                   *file                = NULL;
        std::string             fileName            = "";
        std::string             fileDescription     = "";
        TRACTOGRAMFILEFORMAT    fileFormat          = VTK_BINARY;
        std::size_t             numberOfPoints      = 0;
        std::size_t             numberOfStreamlines = 0;
        uint32_t*               len                 = NULL;
        long                    endPosOfStreamlines = 0;

        // Nifti assumes voxel center to be 0,0,0.
        // TRK assumes 0,0,0 to be one of the corners. But which corner? 
        // That depends on the image orientation, which makes it difficult to compute.
        // The ijk2xyz and xyz2ijk belong to the reference image.
        // TODO: For now, we will assume that the reference image was in LAS orientation. This needs to be properly computed in the future.
        // For example using the approach in Nipy: https://github.com/nipy/nibabel/blob/d23a8336582d07d683d3af73cf097c3be2de9af8/nibabel/orientations.py#L356
        //
        // With LAS assumption: 
        // - when reading a trk file, we will subtract 0.5 from all values
        // - when writing a trk file, we will add 0.5 to all values

        short                   imgDims[3];
        float                   pixDims[3];
        char                    voxOrdr[4];         // We assume this to be LAS for now
        float                   ijk2xyz[4][4];
        float                   xyz2ijk[4][4];

        

    private:


        long*                   streamlinePos       = NULL;      // file positions for first points of streamlines

        mutable std::mutex      reader_mutex;       // Mutex for thread safety

        // TRK specific
        short                   n_scalars_trk    = 0; // TRK file format extension
        short                   n_properties_trk = 0; // TRK file format extension

        std::vector<TractogramField>        findTractogramFields(); // used for printing info only

        // Preloading cache state
        std::size_t                         preloadCount        = 0;
        std::size_t                         preloadedStartIndex = std::numeric_limits<std::size_t>::max();
        std::vector<std::vector<Point>>     preloadedStreamlines;

        // Helper functions

        // Not locked
        bool initReader();
        void cleanupMemory();
        void resetCache();

        // Locked
        bool  loadBatchContaining(std::size_t streamlineIndex);
        bool  readStreamlinePointsFromFile(FILE* bf, std::size_t n, std::vector<Point>& points);

        int   missedCacheCounter = 0;

        // Parallel batch loading
        std::vector<FILE*> batchFile;

    };

    // TODO: Compute this properly and not assume LAS.
    template<typename T>
    void TractogramReader::setReferenceImage(const Image<T>* ref) 
    {
        std::lock_guard<std::mutex> lock(reader_mutex);

        for(int i=0;i<3;++i) {
            for(int j=0;j<4;++j) {
                ijk2xyz[i][j] = ref->ijk2xyz[i][j];
                xyz2ijk[i][j] = ref->xyz2ijk[i][j];
            }
        }

        for(int j=0;j<3;++j) {
            ijk2xyz[3][j] = 0;
            xyz2ijk[3][j] = 0;
            imgDims[j]    = ref->imgDims[j];
            pixDims[j]    = ref->pixDims[j];
        }

        ijk2xyz[3][3] = 1;
        xyz2ijk[3][3] = 1;

        std::strcpy(voxOrdr,"LAS");             

    }

}


