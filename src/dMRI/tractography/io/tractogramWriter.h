#pragma once

#include <string>
#include <vector>
#include "base/nibr.h"
#include "tractogramReader.h"
#include "tractogramField.h"

namespace NIBR
{

    class VTKBinaryWriter {
    public:
        VTKBinaryWriter(std::string _filename);
        ~VTKBinaryWriter();

        bool open();
        bool writeBatch(const std::vector<std::vector<std::vector<float>>>& batch);
        bool close();

    private:
        FILE* file = nullptr;
        std::string         filename;
        std::vector<int>    lengths;              // Store lengths of all streamlines
        size_t              totalPointCount = 0;
        size_t              totalStreamlineCount = 0;
        long                posPointsHeader = 0;  // Position to write POINTS count
        long                posLinesHeader = 0;   // Position to write LINES count

        void writeHeaderPlaceholders();
        void writeFinalHeaders();
        void writeLinesData();
    };

    
    // Other writer functions
    bool writeTractogram            (std::string out_fname,NIBR::TractogramReader* tractogram);
    bool writeTractogram            (std::string out_fname,std::vector<std::vector<std::vector<float>>>& tractogram);
    template<typename T>
    bool writeTractogram            (std::string out_fname,std::vector<std::vector<std::vector<float>>>& tractogram,Image<T>* refImg);
    bool writeTractogram            (std::string out_fname,NIBR::TractogramReader* tractogram,std::vector<size_t>& idx);
    bool writeTractogram            (std::string out_fname,std::string inp_fname,std::vector<size_t>& idx);
    bool writeTractogram            (std::string out_kept_fname,std::string out_rmvd_fname,std::string inp_fname,std::vector<size_t>& idx);

    bool writeTractogram_VTK_binary (std::string out_fname,NIBR::TractogramReader* tractogram);
    bool writeTractogram_VTK_binary (std::string out_fname,NIBR::TractogramReader* tractogram,std::vector<size_t>& idx);
    bool writeTractogram_VTK_binary (std::string out_fname,std::vector<std::vector<std::vector<float>>>& tractogram);

    bool writeTractogram_VTK_ascii  (std::string out_fname,NIBR::TractogramReader* tractogram);
    bool writeTractogram_VTK_ascii  (std::string out_fname,NIBR::TractogramReader* tractogram,std::vector<size_t>& idx);
    bool writeTractogram_VTK_ascii  (std::string out_fname,std::vector<std::vector<std::vector<float>>>& tractogram);

    bool writeTractogram_TRK        (std::string out_fname,NIBR::TractogramReader* tractogram);
    bool writeTractogram_TRK        (std::string out_fname,NIBR::TractogramReader* tractogram,std::vector<size_t>& idx);
    template<typename T>
    bool writeTractogram_TRK        (std::string out_fname,std::vector<std::vector<std::vector<float>>>& tractogram,Image<T>* refImg);

    bool writeTractogram_TCK        (std::string out_fname,NIBR::TractogramReader* tractogram);
    bool writeTractogram_TCK        (std::string out_fname,NIBR::TractogramReader* tractogram,std::vector<size_t>& idx);
    bool writeTractogram_TCK        (std::string out_fname,std::vector<std::vector<std::vector<float>>>& tractogram);

    // The following only works for VTK files
    bool writeTractogram(std::string fname,std::vector<std::vector<std::vector<float>>>& tractogram,std::vector<TractogramField>& fields);

}
