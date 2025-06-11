#pragma once

#include <fstream>
#include <deque>
#include <atomic>
#include <mutex>
#include <future>
#include "base/nibr.h"
#include "math/core.h"
#include "image/image.h"
#include "dMRI/tractography/tractogram.h"

namespace NIBR
{

    #pragma pack(push, 1)
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
    #pragma pack(pop)

    typedef enum {
        UNKNOWN_TRACTOGRAM_FORMAT,
        VTK_ASCII_3,
        VTK_BINARY_3,
        VTK_ASCII_4,
        VTK_BINARY_4,
        VTK_ASCII_5,
        VTK_BINARY_5,
        VTK_ASCII,
        VTK_BINARY,
        TCK,
        TRK
    } TRACTOGRAMFILEFORMAT;

    class TractogramReader {
        
        public:
            
            TractogramReader(std::string _fileName, bool _preload = false);
            ~TractogramReader();

            TractogramReader(const TractogramReader& obj) = delete;             // Disable the copy constructor
            TractogramReader& operator=(const TractogramReader& obj) = delete;  // Disable copy assignment operator

            void printInfo();
            bool isReady()      const {return isInitialized;}
            bool isPreloaded()  const {return isPreloadMode;}
            void reset();

            std::tuple<bool, Streamline, size_t>    getNextStreamline();                        // Returns a tuple: {success, streamline, streamline_index}
            StreamlineBatch                         getNextStreamlineBatch(size_t batchSize);
            const Streamline&                       getStreamline(std::size_t n);               // Preloaded mode specific
            Tractogram                              getTractogram();                            // Complete tractogram reader
            size_t                                  getNumberOfStreamlines() const { return numberOfStreamlines; }
            const std::vector<uint64_t>&            getNumberOfPoints();                        // Size of streamline n is out[n] - out[n-1]. The last element, out[numberOfStreamlines+1] is the total number of points

            std::string             fileName;
            std::string             fileDescription;
            TRACTOGRAMFILEFORMAT    fileFormat;
            std::size_t             numberOfStreamlines;

            // TRK specific
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

            short   imgDims[3];
            float   pixDims[3];
            char    voxOrdr[4];         // We assume this to be LAS for now
            float   ijk2xyz[4][4];
            float   xyz2ijk[4][4];
            
            template<typename T>
            void    setReferenceImage(Image<T>* ref);            

        private:

            FILE* file;
            bool  initReader(std::string _fileName, bool _preload);

            // Core I/O logic, reads a batch from the file. Not thread-safe by itself.
            StreamlineBatch readBatchFromFile(size_t batchSize);
        
            // Producer-consumer members
            void producerLoop();
            std::thread             producerThread;
            std::deque<Streamline>  streamline_buffer;
            std::mutex              buffer_mutex;
            std::condition_variable buffer_cv;
            std::atomic<bool>       stop_producer{false};
            std::atomic<bool>       producer_finished{false};
            size_t                  buffer_capacity;
            size_t                  buffer_low_water_mark;
            bool                    isPreloadMode;
            
            // State variables
            bool  isInitialized = false;
            char* readerBuffer  = nullptr;
            
            // File position tracking
            std::atomic<size_t>     streamlines_read_from_file{0};
            std::atomic<size_t>     consumed_streamline_count{0};
            long                    firstStreamlinePos   = 0;
            long                    currentStreamlinePos = 0;

            std::vector<uint64_t>   numberOfPoints; // Number of points

            // VTK specific
            bool            vtk_needsByteSwap   = true;
            long            vtk_currOffsetPos   = 0;
            int64_t         vtk_prevOffset      = 0;
            long            vtk_firstOffsetPos  = 0;
            int64_t         vtk_firstOffset     = 0;
            bool            vtk_isV5            = false;
            bool            vtk5_64bitOffsets   = true;

            // TCK specific
            std::vector<char>   tck_read_buffer;
            size_t              tck_buffer_offset = 0;

            // TRK specific
            short           n_scalars_trk;          // TRK file format extension
            short           n_properties_trk;       // TRK file format extension

    };

    template<typename T>
    void TractogramReader::setReferenceImage(Image<T>* ref) 
    {
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

        std::strcpy(voxOrdr,"LAS");             // TODO: Compute this properly and not assume LAS.

    }

}

