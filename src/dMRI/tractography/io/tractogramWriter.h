#pragma once

#include <string>
#include <vector>
#include "base/nibr.h"
#include "tractogramReader.h"
#include "tractogramField.h"

namespace NIBR
{    

    class IBatchWriter {
    public:
        virtual ~IBatchWriter() = default;
        virtual bool open() = 0;
        virtual bool writeBatch(const std::vector<std::vector<std::vector<float>>>& batch) = 0;
        virtual bool close() = 0;
        virtual const std::string& getFilename() const = 0;
    };

    class TractogramWriter {
    public:
        TractogramWriter(std::string _filename);
        ~TractogramWriter();

        TractogramWriter(const TractogramWriter&) = delete;
        TractogramWriter& operator=(const TractogramWriter&) = delete;

        bool open();
        bool writeBatch(const std::vector<std::vector<std::vector<float>>>& batch);
        bool close();
        bool isOpen() const  { return is_open_ && !is_closed_; }
        bool isValid() const { return pImpl_ != nullptr; }


    private:
        std::string filename_;
        std::unique_ptr<IBatchWriter> pImpl_;
        bool is_open_   = false;
        bool is_closed_ = false;
    };

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
