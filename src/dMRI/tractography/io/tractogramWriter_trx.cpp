#include "base/nibr.h"      // For disp, SGNTR
#include "zip.h"
#include "tractogramWriter_trx.h"

namespace NIBR {

TRXWriter::TRXWriter(std::string _filename) : filename_(std::move(_filename)) {}

TRXWriter::~TRXWriter()
{
    if (is_open_) {
        disp(MSG_WARN, "TRXWriter for %s destroyed without explicit close.", filename_.c_str());
    }
}

void TRXWriter::setTRXReference(const Image<bool>& ref)
{
    ref_affine_ = Eigen::Matrix4f::Identity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            ref_affine_(i, j) = ref.ijk2xyz[i][j];
    ref_dims_ = { static_cast<uint16_t>(ref.imgDims[0]),
                  static_cast<uint16_t>(ref.imgDims[1]),
                  static_cast<uint16_t>(ref.imgDims[2]) };
    has_reference_ = true;
}

void TRXWriter::setTRXFields(const std::vector<TractogramField>& fields)
{
    fields_ = fields;
}

bool TRXWriter::open()
{
    stream_ = trx::TrxStream(dtype_);
    if (has_reference_) {
        stream_.set_voxel_to_rasmm(ref_affine_);
        stream_.set_dimensions(ref_dims_);
    }
    stream_.set_metadata_mode(trx::TrxStream::MetadataMode::OnDisk);
    is_open_ = true;
    return true;
}

bool TRXWriter::writeBatch(const StreamlineBatch& batch)
{
    if (!is_open_) {
        disp(MSG_ERROR, "TRXWriter: File not open for writing.");
        return false;
    }

    for (const auto& streamline : batch) {
        if (streamline.empty()) continue;
        stream_.push_streamline(streamline);
        lengths_.push_back(streamline.size());
    }

    return true;
}

bool TRXWriter::close(long& finalStreamlineCount, long& finalPointCount)
{
    if (!is_open_) {
        finalStreamlineCount = 0;
        finalPointCount      = 0;
        return true;
    }

    is_open_ = false;

    finalStreamlineCount = static_cast<long>(stream_.num_streamlines());
    finalPointCount      = static_cast<long>(stream_.num_vertices());

    if (finalStreamlineCount == 0) {
        disp(MSG_WARN, "TRXWriter: No streamlines to write for %s.", filename_.c_str());
        return true;
    }

    // Push DPS (STREAMLINE_OWNER) fields
    for (const auto& field : fields_) {
        if (field.owner != STREAMLINE_OWNER) continue;
        if (field.datatype != FLOAT32_DT)    continue;
        if (field.data == nullptr)            continue;

        float** data = reinterpret_cast<float**>(field.data);

        std::vector<float> flat;
        flat.reserve(static_cast<size_t>(finalStreamlineCount) * field.dimension);
        for (size_t s = 0; s < static_cast<size_t>(finalStreamlineCount); ++s)
            for (int d = 0; d < field.dimension; ++d)
                flat.push_back(data[s][d]);

        const std::string dps_name = (field.dimension > 1)
            ? field.name + "." + std::to_string(field.dimension)
            : field.name;

        try {
            stream_.push_dps_from_vector(dps_name, "float32", flat);
        } catch (const std::exception& e) {
            disp(MSG_ERROR, "TRXWriter: Failed to push DPS field '%s': %s",
                 field.name.c_str(), e.what());
        }
    }

    // Push POINT_OWNER fields as DPV before finalizing
    for (const auto& field : fields_) {
        if (field.owner != POINT_OWNER)    continue;
        if (field.datatype != FLOAT32_DT)  continue;
        if (field.data == nullptr)          continue;

        float*** data = reinterpret_cast<float***>(field.data);

        // Flatten to row-major interleaved: [p0_d0, p0_d1, ..., p1_d0, ...]
        std::vector<float> flat;
        flat.reserve(static_cast<size_t>(finalPointCount) * field.dimension);
        for (size_t s = 0; s < lengths_.size(); ++s) {
            if (data[s] == nullptr) {
                // Streamline color was not computed (e.g. allocation failure) — fill with zeros
                for (size_t p = 0; p < lengths_[s]; ++p)
                    for (int d = 0; d < field.dimension; ++d)
                        flat.push_back(0.0f);
                continue;
            }
            for (size_t p = 0; p < lengths_[s]; ++p) {
                for (int d = 0; d < field.dimension; ++d) {
                    flat.push_back(data[s][p][d]);
                }
            }
        }

        // Embed n_cols in the name for multi-column DPV (e.g. "RGB.3")
        // so the zip entry becomes dpv/RGB.3.float32, which TRX readers parse correctly.
        const std::string dpv_name = (field.dimension > 1)
            ? field.name + "." + std::to_string(field.dimension)
            : field.name;

        try {
            stream_.push_dpv_from_vector(dpv_name, "float32", flat);
        } catch (const std::exception& e) {
            disp(MSG_ERROR, "TRXWriter: Failed to push DPV field '%s': %s",
                 field.name.c_str(), e.what());
        }
    }

    try {
        const trx::TrxScalarType out_dtype = (dtype_ == "float16")
            ? trx::TrxScalarType::Float16
            : trx::TrxScalarType::Float32;
        stream_.finalize(filename_, out_dtype, ZIP_CM_STORE);
    } catch (const std::exception& e) {
        disp(MSG_ERROR, "TRXWriter: Failed to write %s. %s", filename_.c_str(), e.what());
        return false;
    }

    // Append groups after finalize (append_groups_to_zip opens the zip without truncating)
    if (!groups_.empty()) {
        try {
            trx::append_groups_to_zip(filename_, groups_);
        } catch (const std::exception& e) {
            disp(MSG_ERROR, "TRXWriter: Failed to write groups to %s: %s",
                 filename_.c_str(), e.what());
        }
    }

    disp(MSG_DEBUG, "TRXWriter: Successfully closed %s. Streamlines: %ld, Points: %ld",
         filename_.c_str(), finalStreamlineCount, finalPointCount);

    return true;
}

}
