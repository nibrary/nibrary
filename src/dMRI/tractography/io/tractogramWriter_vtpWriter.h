#pragma once

#include <string>
#include <vector>
#include "tractogramWriter.h"

namespace NIBR
{

    class VTPAppendedWriter : public IBatchWriter {
    public:
        VTPAppendedWriter(std::string _filename);
        ~VTPAppendedWriter();

        bool open() override;
        bool writeBatch(const std::vector<std::vector<std::vector<float>>>& batch) override;
        bool close() override;
        const std::string& getFilename() const override { return final_filename_; }

    private:
        std::string final_filename_;
        std::string temp_points_filename_;
        FILE* temp_points_file_ = nullptr;
        
        std::vector<int64_t> lengths_; // Store streamline lengths
        int64_t              total_point_count_ = 0;
        int64_t              total_streamline_count_ = 0;

        bool write_xml_header(FILE* f, int64_t offset_conn, int64_t offset_offs);
        bool write_appended_data(FILE* f);
        bool write_connectivity(FILE* f);
        bool write_offsets(FILE* f);
    };


}