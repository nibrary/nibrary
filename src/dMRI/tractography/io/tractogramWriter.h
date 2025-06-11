#pragma once

#include <string>
#include <vector>
#include <memory>
#include "base/nibr.h"
#include "dMRI/tractography/tractogram.h"
#include "image/image.h"
#include "tractogramReader.h"
#include "tractogramField.h"

namespace NIBR
{

    // Interface for batch writing operations
    class IBatchWriter {
    public:
        virtual ~IBatchWriter() = default;

        virtual bool open() = 0;
        virtual bool writeBatch(const StreamlineBatch& batch) = 0;
        virtual bool close(long& finalStreamlineCount, long& finalPointCount) = 0;

        virtual const std::string& getFilename() const = 0;

        // Optional methods for setting context, with default no-op implementations
        virtual void setTRKReference(const Image<bool>& /*reference*/) {}
        virtual void setVTKFields(const std::vector<TractogramField>& /*fields*/) {}
    };

    // Main TractogramWriter class using Pimpl idiom
    class TractogramWriter {
    public:
        TractogramWriter(std::string _filename);
        ~TractogramWriter();

        // Disable copy and assignment
        TractogramWriter(const TractogramWriter&) = delete;
        TractogramWriter& operator=(const TractogramWriter&) = delete;

        // Methods to set context before opening
        template<typename T>
        void setTRKReference(const Image<T>& refImg) {if (pImpl_) { Image<bool> tmp; tmp.createFromTemplate(refImg,false); pImpl_->setTRKReference(tmp);} else {disp(MSG_WARN, "Cannot set TRK reference info: writer implementation is null.");}}
        void setVTKIsAscii(bool isAscii = true); // Has to be used immediately after the constructor
        void setVTKFields(const std::vector<TractogramField>& fields);

        bool open();
        bool writeStreamline(const Streamline& streamline);
        bool writeBatch(const StreamlineBatch& batch);
        bool writeBatch(StreamlineBatch&& batch);
        bool close(); // Returns true on success

        bool isOpen() const  { return is_open_ && !is_closed_; }
        bool isValid() const { return pImpl_ != nullptr; }

        long getFinalStreamlineCount() const { return finalStreamlineCount_; }
        long getFinalPointCount() const { return finalPointCount_; }


    private:
        std::string filename_;
        std::unique_ptr<IBatchWriter> pImpl_;
        bool is_open_                       = false;
        bool is_closed_                     = false;
        long finalStreamlineCount_          = 0;
        long finalPointCount_               = 0;
        TRACTOGRAMFILEFORMAT fileFormat_    = UNKNOWN_TRACTOGRAM_FORMAT;

        // For ascynchronous double buffering
        void asyncWriteLoop();
        
        StreamlineBatch asyncBufferA_;
        StreamlineBatch asyncBufferB_;
        
        StreamlineBatch* fillAsyncBuffer_{&asyncBufferA_};
        StreamlineBatch* writeAsyncBuffer_{&asyncBufferB_};

        std::thread             writerThread_;
        std::mutex              asyncMutex_;
        std::condition_variable asyncCond_;
        std::atomic<bool>       write_pending_{false};
        std::atomic<bool>       stop_requested_{false};
        std::atomic<bool>       async_write_failed_{false};

    };

    // Helper free functions for common writing tasks

    // Write a whole tractogram from a TractogramReader
    bool writeTractogram(std::string out_fname, NIBR::TractogramReader* reader);

    // Write a whole tractogram from an in-memory Tractogram object
    bool writeTractogram(std::string out_fname, const Tractogram& tractogram);

    // Write a tractogram with associated field data (primarily for VTK)
    bool writeTractogram(std::string out_fname, const Tractogram& tractogram, const std::vector<TractogramField>& fields);

    // Write a subset of streamlines from a TractogramReader using indices
    bool writeTractogram(std::string out_fname, NIBR::TractogramReader* reader, const std::vector<size_t>& idx_to_keep);

    // Write a subset of streamlines from one file to another using indices
    bool writeTractogram(std::string out_fname, std::string inp_fname, const std::vector<size_t>& idx_to_keep);

    // Write kept and removed streamlines to two separate files
    bool writeTractogram(std::string out_kept_fname, std::string out_rmvd_fname, std::string inp_fname, const std::vector<size_t>& idx_to_keep);

    // For TRK, reference image information must be provided separately if not using a reader
    template<typename T>
    bool writeTractogram(std::string out_fname, const Tractogram& tractogram, const Image<T>& refImg);

}
