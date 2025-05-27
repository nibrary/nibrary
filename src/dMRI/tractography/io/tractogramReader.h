#pragma once

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <float.h>
#include <cstdint>
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
        TRK,
        VTP
    } TRACTOGRAMFILEFORMAT;

    class TractogramReader {
        
        public:
            
            TractogramReader();
            TractogramReader(std::string _fileName);
            TractogramReader(std::string _fileName, bool _usePreload);

            ~TractogramReader();
            TractogramReader(const TractogramReader& obj);
            void copyFrom(const TractogramReader& obj);
            void destroyCopy();
            void printInfo();
            bool isPreloaded() {return usePreload;}
            
            bool     initReader(std::string _fileName);
            bool     initReader(std::string _fileName, bool _usePreload);

            template<typename T>
            void     setReferenceImage(Image<T>* ref);
            
            float**  readStreamline(std::size_t n);          // Core function to read streamlines
            void     deleteStreamline(float**, std::size_t);

            bool     readNextBatch(size_t batchSize, std::vector<std::vector<std::vector<float>>>& batch_out);
            void     resetBatchReader() {currentStreamlineIdx = 0;}
            
            std::vector<std::vector<float>> readStreamlineVector(std::size_t n);
            std::vector<Point>              readStreamlinePoints(std::size_t n);

            FILE                   *file;
            std::string             fileName;
            std::string             fileDescription;
            TRACTOGRAMFILEFORMAT    fileFormat;
            
            std::size_t             numberOfPoints;
            std::size_t             numberOfStreamlines;
            uint32_t*               len;

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
            
            long*                   streamlinePos;      // file positions for first points of streamlines


            // --- VTP Specific ---
            long                    points_xml_offset_       = -1;
            long                    connectivity_xml_offset_ = -1;
            long                    offsets_xml_offset_      = -1;
            bool                    vtp_is_little_endian_    = true;
            std::string             vtp_current_tag_; 
        
        private:

            bool                    readToMemory();

            // TRK specific
            short                   n_scalars_trk;      // TRK file format extension
            short                   n_properties_trk;   // TRK file format extension

            // loading the whole tractogram to memory
            float*                  fileBuffer       = NULL;
            std::vector<float**>*   streamlineBuffer = NULL;
            bool                    usePreload       = false;    // if we load the whole thing in memory
            bool                    finishedLoading  = false;    //this is set to true after reading the whole thing in memmory

            // batch reading
            size_t                  currentStreamlineIdx = 0; // Tracks the next streamline to read

            // --- VTP Specific ---
            bool                    parseVTPHeaderWithExpat(); // <-- New name
            template<typename T>
            bool                    readVTPAppendedData(long offset, std::vector<T>& data);
            bool                    readVTPPoints(int64_t index, float* point);

            long                    appended_data_start_pos_ = 0;
            
            bool                    vtp_needs_swap_          = false;
            long                    vtp_points_file_offset_  = 0;
            std::vector<int64_t>    vtp_offsets_;
            std::vector<int64_t>    vtp_connectivity_;
            
            
                      

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

    template<typename T>
    bool NIBR::TractogramReader::readVTPAppendedData(long offset, std::vector<T>& data) {
        if (appended_data_start_pos_ == 0 || offset < 0) return false; // Use member var

        fseek(file, appended_data_start_pos_ + offset, SEEK_SET); // Use member var

        uint64_t data_size_bytes = 0;
        if (fread(&data_size_bytes, sizeof(uint64_t), 1, file) != 1) return false;
        if (vtp_needs_swap_) swapByteOrder(data_size_bytes); // Use member var

        size_t num_elements = data_size_bytes / sizeof(T);
        data.resize(num_elements);

        if (fread(data.data(), sizeof(T), num_elements, file) != num_elements) return false;

        if (vtp_needs_swap_) { // Use member var
            for (T& val : data) {
                swapByteOrder(val); 
            }
        }
        return true;
    }

}


