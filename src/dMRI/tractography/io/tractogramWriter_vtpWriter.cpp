#include <cstdio>
#include <vector>
#include <cstdint>
#include <numeric>

#include "tractogramWriter_vtpWriter.h"

using namespace NIBR;


VTPAppendedWriter::VTPAppendedWriter(std::string _filename) 
    : final_filename_(_filename) 
{
    // Create a unique-ish temp filename
    temp_points_filename_ = final_filename_ + ".points.tmp";
}

VTPAppendedWriter::~VTPAppendedWriter() {
    // Ensure temp files are closed and removed if something went wrong
    if (temp_points_file_ != nullptr) {
        fclose(temp_points_file_);
        std::remove(temp_points_filename_.c_str());
    }
}

bool VTPAppendedWriter::open() {
    // Open the temporary file for points
    temp_points_file_ = fopen(temp_points_filename_.c_str(), "wb");
    if (temp_points_file_ == nullptr) {
        disp(MSG_ERROR, "Failed to open temporary points file: %s", temp_points_filename_.c_str());
        return false;
    }
    lengths_.clear();
    total_point_count_ = 0;
    total_streamline_count_ = 0;
    disp(MSG_DEBUG, "Opened temp file %s for VTP writing.", temp_points_filename_.c_str());
    return true;
}

bool VTPAppendedWriter::writeBatch(const std::vector<std::vector<std::vector<float>>>& batch) {
    if (temp_points_file_ == nullptr) {
        disp(MSG_ERROR, "Cannot write VTP batch, temp file not open.");
        return false;
    }

    for (const auto& streamline : batch) {
        int64_t len = streamline.size();
        if (len > 0) {
            lengths_.push_back(len);
            total_streamline_count_++;
            total_point_count_ += len;

            // Write points to temp file
            for (const auto& point : streamline) {
                fwrite(point.data(), sizeof(float), 3, temp_points_file_);
            }
        }
    }

    return !ferror(temp_points_file_);
}


bool VTPAppendedWriter::write_xml_header(FILE* f, int64_t offset_conn, int64_t offset_offs) {
    fprintf(f, "<?xml version=\"1.0\"?>\n");
    fprintf(f, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(f, "  <PolyData>\n");
    fprintf(f, "    <Piece NumberOfPoints=\"%lld\" NumberOfVerts=\"0\" NumberOfLines=\"%lld\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n", 
            (long long)total_point_count_, (long long)total_streamline_count_);
    
    // Points
    fprintf(f, "      <Points>\n");
    fprintf(f, "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\" />\n");
    fprintf(f, "      </Points>\n");

    // Lines (Connectivity & Offsets)
    fprintf(f, "      <Lines>\n");
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\" />\n", (long long)offset_conn);
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%lld\" />\n", (long long)offset_offs);
    fprintf(f, "      </Lines>\n");

    // TODO: Add PointData and CellData sections if needed in the future

    fprintf(f, "    </Piece>\n");
    fprintf(f, "  </PolyData>\n");
    fprintf(f, "  <AppendedData encoding=\"raw\">\n");
    fprintf(f, "   _"); // Start of appended data marker (note space before for alignment)

    return !ferror(f);
}

bool VTPAppendedWriter::write_connectivity(FILE* f) {
    int64_t current_point = 0;
    const int64_t buffer_size = 10000; // Write in chunks
    std::vector<int64_t> buffer;
    buffer.reserve(buffer_size);

    while(current_point < total_point_count_) {
        buffer.clear();
        int64_t points_to_write = std::min(buffer_size, total_point_count_ - current_point);
        for(int64_t i = 0; i < points_to_write; ++i) {
            buffer.push_back(current_point + i);
        }
        fwrite(buffer.data(), sizeof(int64_t), buffer.size(), f);
        current_point += points_to_write;
    }
    return !ferror(f);
}

bool VTPAppendedWriter::write_offsets(FILE* f) {
    std::vector<int64_t> offsets;
    offsets.reserve(total_streamline_count_);
    int64_t current_offset = 0;
    for (int64_t len : lengths_) {
        current_offset += len;
        offsets.push_back(current_offset);
    }
    fwrite(offsets.data(), sizeof(int64_t), offsets.size(), f);
    return !ferror(f);
}


bool VTPAppendedWriter::write_appended_data(FILE* f) {
    
    // --- 1. Write Points ---
    FILE* temp_in = fopen(temp_points_filename_.c_str(), "rb");
    if (!temp_in) return false;

    fseek(temp_in, 0, SEEK_END);
    uint64_t size_points = ftell(temp_in);
    fseek(temp_in, 0, SEEK_SET);

    fwrite(&size_points, sizeof(uint64_t), 1, f);
    
    char copy_buffer[8192];
    size_t bytes_read;
    while ((bytes_read = fread(copy_buffer, 1, sizeof(copy_buffer), temp_in)) > 0) {
        fwrite(copy_buffer, 1, bytes_read, f);
    }
    fclose(temp_in);
    
    // --- 2. Write Connectivity ---
    uint64_t size_conn = total_point_count_ * sizeof(int64_t);
    fwrite(&size_conn, sizeof(uint64_t), 1, f);
    if (!write_connectivity(f)) return false;

    // --- 3. Write Offsets ---
    uint64_t size_offs = total_streamline_count_ * sizeof(int64_t);
    fwrite(&size_offs, sizeof(uint64_t), 1, f);
    if (!write_offsets(f)) return false;

    return !ferror(f);
}


bool VTPAppendedWriter::close() {
    if (temp_points_file_ == nullptr) {
        disp(MSG_WARN, "VTP writer already closed or never opened.");
        return true; 
    }

    // Close the temp file first
    fclose(temp_points_file_);
    temp_points_file_ = nullptr;

    // Calculate sizes and offsets
    FILE* temp_check = fopen(temp_points_filename_.c_str(), "rb");
    if (!temp_check) { 
        disp(MSG_ERROR,"Could not reopen temp file to check size."); 
        std::remove(temp_points_filename_.c_str()); 
        return false; 
    }
    fseek(temp_check, 0, SEEK_END);
    uint64_t size_points = ftell(temp_check);
    fclose(temp_check);

    uint64_t size_conn = total_point_count_ * sizeof(int64_t);
    // uint64_t size_offs = total_streamline_count_ * sizeof(int64_t);
    uint64_t header_size = sizeof(uint64_t);

    uint64_t offset_conn = header_size + size_points;
    uint64_t offset_offs = offset_conn + header_size + size_conn;

    // Open the final VTP file
    FILE* out_file = fopen(final_filename_.c_str(), "w"); // Use "w" - text for XML header
    if (!out_file) {
        disp(MSG_ERROR, "Failed to open final VTP file: %s", final_filename_.c_str());
        std::remove(temp_points_filename_.c_str());
        return false;
    }

    // Write XML header
    if (!write_xml_header(out_file, offset_conn, offset_offs)) {
        disp(MSG_ERROR, "Failed to write VTP XML header.");
        fclose(out_file);
        std::remove(temp_points_filename_.c_str());
        return false;
    }

    // IMPORTANT: Reopen in binary mode to append data, or use fseek and 'wb' from start.
    // Simpler: Finish text write, then reopen in 'ab' (append binary).
    fflush(out_file);
    fclose(out_file);
    out_file = fopen(final_filename_.c_str(), "ab"); // Append Binary
    if (!out_file) {
        disp(MSG_ERROR, "Failed to reopen VTP file for binary append: %s", final_filename_.c_str());
        std::remove(temp_points_filename_.c_str());
        return false;
    }

    // Write Appended Data
    if (!write_appended_data(out_file)) {
        disp(MSG_ERROR, "Failed to write VTP appended data.");
        fclose(out_file);
        std::remove(temp_points_filename_.c_str());
        return false;
    }

    // Write XML Footer - Reopen in 'a' (append text)
    fclose(out_file);
    out_file = fopen(final_filename_.c_str(), "a"); // Append Text
    fprintf(out_file, "\n  </AppendedData>\n");
    fprintf(out_file, "</VTKFile>\n");
    fclose(out_file);

    // Clean up temp file
    std::remove(temp_points_filename_.c_str());

    disp(MSG_DEBUG, "Successfully closed VTP file %s", final_filename_.c_str());
    return true;
}