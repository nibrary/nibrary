#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>
#include <cstdint>
#include <Eigen/Dense>
#include <trx/trx.h>
#include "tractogramWriter.h"   // For IBatchWriter, StreamlineBatch
#include "tractogramField.h"

namespace NIBR
{
    class TRXWriter : public IBatchWriter {
    public:
        TRXWriter(std::string _filename);
        ~TRXWriter() override;

        TRXWriter(const TRXWriter&) = delete;
        TRXWriter& operator=(const TRXWriter&) = delete;

        bool open() override;
        bool writeBatch(const StreamlineBatch& batch) override;
        bool close(long& finalStreamlineCount, long& finalPointCount) override;
        const std::string& getFilename() const override { return filename_; }

        void setTRXDtype(const std::string& dtype) override { dtype_ = dtype; }
        void setTRXReference(const Image<bool>& ref) override;
        void setTRXFields(const std::vector<TractogramField>& fields) override;
        void setGroups(std::map<std::string, std::vector<uint32_t>> groups) { groups_ = std::move(groups); }

    private:
        std::string                                      filename_;
        std::string                                      dtype_            = "float16";
        trx::TrxStream                                   stream_;
        bool                                             is_open_          = false;
        bool                                             has_reference_    = false;
        Eigen::Matrix4f                                  ref_affine_       = Eigen::Matrix4f::Identity();
        std::array<uint16_t,3>                           ref_dims_         = {0, 0, 0};
        std::vector<TractogramField>                     fields_;
        std::vector<size_t>                              lengths_;   // point count per streamline, in write order
        std::map<std::string, std::vector<uint32_t>>     groups_;
    };
}
