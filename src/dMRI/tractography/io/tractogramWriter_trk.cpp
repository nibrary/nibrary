#include <cstring>           // For std::memcpy, std::strcpy, std::memset
#include "base/nibr.h"       // For disp, SGNTR
#include "math/core.h"       // For applyTransform, inverseAffine
#include "tractogramWriter_trk.h"

namespace NIBR {

TRKWriter::TRKWriter(std::string _filename) : filename_(std::move(_filename)), refInfoSet_(false) 
{
    std::memset(&header_, 0, sizeof(trkFileStruct));
    std::strcpy(header_.id_string, "TRACK");
    header_.version = 2;                        // TrackVis version 2
    header_.hdr_size = sizeof(trkFileStruct);   // Should be 1000
}

TRKWriter::~TRKWriter() 
{
    if (file_ != nullptr) {
        disp(MSG_WARN, "TRKWriter for %s destroyed without explicit close. Attempting to close.", filename_.c_str());
        long dummy1, dummy2;
        close(dummy1, dummy2);
    }
}

void TRKWriter::setTRKReference(const Image<bool>& ref) 
{
    
    for(int i=0; i<3; ++i) {
        refInfo_.imgDims[i] = ref.imgDims[i];
        refInfo_.pixDims[i] = ref.pixDims[i];
    }
    for(int i=0; i<3; ++i) {
        for(int j=0; j<4; ++j) {
            refInfo_.ijk2xyz[i][j] = ref.ijk2xyz[i][j];
        }
    }
    refInfo_.ijk2xyz[3][0] = 0.0f;
    refInfo_.ijk2xyz[3][1] = 0.0f;
    refInfo_.ijk2xyz[3][2] = 0.0f;
    refInfo_.ijk2xyz[3][3] = 1.0f;

    if (!NIBR::inverseAffine(refInfo_.ijk2xyz, refInfo_.xyz2ijk)) {
        disp(MSG_WARN, "Failed to compute xyz2ijk transform for TRK writer. Using identity.");
        for(int i=0; i<4; ++i) for(int j=0; j<4; ++j) refInfo_.xyz2ijk[i][j] = (i==j) ? 1.0f : 0.0f;
    }

    // TODO: Let's use LAS. Could be obtained from Image<T> and nibrary has functions for it. 
    // But it is not clear how to use this information in trk format.
    std::strcpy(refInfo_.voxOrdr, "LAS");

    if (NIBR::isAffineZero(refInfo_.xyz2ijk)) {
        if (!NIBR::inverseAffine(refInfo_.ijk2xyz, refInfo_.xyz2ijk)) {
            disp(MSG_WARN, "TRKWriter: Failed to compute xyz2ijk from ijk2xyz for reference info. Using identity for xyz2ijk.");
            for(int i=0; i<4; ++i) for(int j=0; j<4; ++j) refInfo_.xyz2ijk[i][j] = (i==j) ? 1.0f : 0.0f;
        }
    }
    
    refInfoSet_ = true;
    initializeHeader();
}

void TRKWriter::initializeHeader() 
{

    if (!refInfoSet_) {
        disp(MSG_WARN, "TRKWriter: Reference info not set. Header will use defaults.");
        for(int i=0; i<4; ++i) for(int j=0; j<4; ++j) header_.vox_to_ras[i][j] = (i==j) ? 1.0f : 0.0f;
        std::strcpy(header_.voxel_order, "LAS");
        return;
    }

    for(int i=0; i<3; ++i) {
        header_.dim[i] = refInfo_.imgDims[i];
        header_.voxel_size[i] = refInfo_.pixDims[i];
    }

    for(int i=0; i<4; ++i) {
        for(int j=0; j<4; ++j) {
            header_.vox_to_ras[i][j] = refInfo_.ijk2xyz[i][j];
        }
    }
    std::strcpy(header_.voxel_order, refInfo_.voxOrdr);

    // For now, forcefully disable scalars and properties to prevent corrupt files.
    // Data writing part is not yet implemented.
    header_.n_scalars    = 0;
    header_.n_properties = 0;

    if (refInfo_.n_scalars_trk_to_write > 0 || refInfo_.n_properties_trk_to_write > 0) {
        disp(MSG_WARN, "TRKWriter: Writing scalar/property data is not yet implemented. n_scalars and n_properties will be set to 0 in the header.");
    }

}


bool TRKWriter::open() 
{

    if (!refInfoSet_) {
        disp(MSG_FATAL, "TRKWriter: Opening file %s without TRK reference information. Header might be incomplete/default.", filename_.c_str());
        return false;
    }

    file_ = fopen(filename_.c_str(), "wb+");

    if (file_ == nullptr) {
        disp(MSG_ERROR, "TRKWriter: Cannot open output file: %s", filename_.c_str());
        return false;
    }

    // Write the header (with n_count = 0 initially)
    header_.n_count = 0; // Will be updated in close()
    if (fwrite(&header_, sizeof(trkFileStruct), 1, file_) != 1) {
        disp(MSG_ERROR, "TRKWriter: Failed to write header to %s.", filename_.c_str());
        fclose(file_);
        file_ = nullptr;
        return false;
    }

    currentStreamlineCount_ = 0;
    currentPointCount_ = 0;
    return true;
}

bool TRKWriter::writeBatch(const StreamlineBatch& batch) 
{

    if (file_ == nullptr) {
        disp(MSG_ERROR, "TRKWriter: File not open for writing.");
        return false;
    }

    for (const auto& streamline_world : batch) {

        if (streamline_world.empty()) {
            continue;
        }

        int32_t num_points = static_cast<int32_t>(streamline_world.size());

        // Write number of points for this streamline
        if (fwrite(&num_points, sizeof(int32_t), 1, file_) != 1) {
            disp(MSG_ERROR, "TRKWriter: Failed to write streamline length to %s.", filename_.c_str());
            return false;
        }

        // For each point, transform to voxel space, add 0.5, and write
        Point3D p_vox;
        for (const auto& p_world : streamline_world) {

            applyTransform(p_vox, p_world, refInfo_.xyz2ijk);

            p_vox[0] += 0.5f;
            p_vox[1] += 0.5f;
            p_vox[2] += 0.5f;

            if (fwrite(p_vox.data(), sizeof(float), 3, file_) != 3) {
                disp(MSG_ERROR, "TRKWriter: Failed to write point data to %s.", filename_.c_str());
                return false;
            }
            currentPointCount_++;
        }
        
        // TODO: Write per-point scalars if header_.n_scalars > 0
        // This requires field data to be passed and processed.
        // For now, assuming no per-point scalars are written beyond coordinates.
        if (header_.n_scalars > 0) {
            // Example: float zero_scalar[header_.n_scalars] = {0};
            // for (int i=0; i<num_points; ++i) fwrite(zero_scalar, sizeof(float), header_.n_scalars, file_);
            // This part needs actual scalar data.
            // For now, we'll skip writing scalars if data isn't available.
            // A robust implementation would require TractogramField-like data.
        }
        currentStreamlineCount_++;
    }
    // TODO: Write per-streamline properties if header_.n_properties > 0
    // This would happen after all points of a streamline (or after the loop for all streamlines in batch).

    return true;
}

bool TRKWriter::close(long& finalStreamlineCount, long& finalPointCount) 
{

    if (file_ == nullptr) {
        finalStreamlineCount = currentStreamlineCount_;
        finalPointCount      = currentPointCount_;
        return true;
    }

    // Update the streamline count in the header
    header_.n_count = static_cast<int32_t>(currentStreamlineCount_);

    // Seek to the beginning of the file to rewrite the header
    if (fseek(file_, 0, SEEK_SET) != 0) {
        disp(MSG_ERROR, "TRKWriter: Failed to seek to beginning of file %s to update header.", filename_.c_str());
        // Continue to try to close
    } else {
        if (fwrite(&header_, sizeof(trkFileStruct), 1, file_) != 1) {
            disp(MSG_ERROR, "TRKWriter: Failed to write updated header to %s.", filename_.c_str());
            // Continue to try to close
        }
    }

    if (fclose(file_) != 0) {
        disp(MSG_ERROR, "TRKWriter: Error closing file %s.", filename_.c_str());
        file_ = nullptr;
        finalStreamlineCount = currentStreamlineCount_;
        finalPointCount      = currentPointCount_;
        return false;
    }

    file_ = nullptr;
    finalStreamlineCount = currentStreamlineCount_;
    finalPointCount      = currentPointCount_;
    disp(MSG_DEBUG, "TRKWriter: Successfully closed %s. Streamlines: %ld, Points: %ld", filename_.c_str(), finalStreamlineCount, finalPointCount);
    return true;
}

}
