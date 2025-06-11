#include <cstring>                 // For strcpy, memset
#include <algorithm>               // For std::sort, std::unique for idx processing
#include "base/vectorOperations.h" // For removeIdx
#include "math/core.h"             // For inverseAffine if needed
#include "tractogramWriter.h"
#include "tractogramWriter_tck.h"
#include "tractogramWriter_trk.h"
#include "tractogramWriter_vtk_binary.h"
#include "tractogramWriter_vtk_ascii.h"

#define WRITE_BUFFER_SIZE   20000

namespace NIBR {

TractogramWriter::TractogramWriter(std::string _filename) : filename_(std::move(_filename)) 
{
    
    std::string ext = getFileExtension(filename_);

    if (ext == "tck") {
        pImpl_      = std::make_unique<TCKWriter>(filename_);
        fileFormat_ = TCK;
        disp(MSG_DEBUG, "Selected TCKWriter for %s", filename_.c_str());
    } else if (ext == "trk") {
        pImpl_      = std::make_unique<TRKWriter>(filename_);
        fileFormat_ = TRK;
        disp(MSG_DEBUG, "Selected TRKWriter for %s", filename_.c_str());
    } else if (ext == "vtk") {
        pImpl_      = std::make_unique<VTKBinaryWriter>(filename_);
        fileFormat_ = VTK_BINARY_3;
        disp(MSG_DEBUG, "Selected VTKBinaryWriter for %s", filename_.c_str());
    } else {
        pImpl_      = nullptr;
        fileFormat_ = UNKNOWN_TRACTOGRAM_FORMAT;
        disp(MSG_ERROR, "Unsupported output file extension: .%s for file: %s. Cannot create writer.", ext.c_str(), filename_.c_str());
    }

    asyncBufferA_.reserve(WRITE_BUFFER_SIZE);
    asyncBufferB_.reserve(WRITE_BUFFER_SIZE);
}

TractogramWriter::~TractogramWriter() 
{
    if (is_open_ && !is_closed_) {
        disp(MSG_WARN, "TractogramWriter for %s destroyed without explicit close(). Flushing and closing.", filename_.c_str());
        close();
    }
}

void TractogramWriter::setVTKIsAscii(bool isAscii) 
{
    if ((isAscii == true) && (fileFormat_ == VTK_BINARY_3)) {
        pImpl_ = std::make_unique<VTKAsciiWriter>(filename_);
        disp(MSG_DEBUG, "Selected VTKAsciiWriter for %s", filename_.c_str());
    }

    if ((isAscii == true) && (fileFormat_ != VTK_BINARY_3)) {
        disp(MSG_ERROR, "Unsupported output file extension. VTKAsciiWriter only support .vtk. %s", filename_.c_str());
        return;
    }
}

void TractogramWriter::setVTKFields(const std::vector<TractogramField>& fields) 
{
    if (pImpl_) {
        pImpl_->setVTKFields(fields);
    } else {
        disp(MSG_WARN, "Cannot set VTK fields: writer implementation is null.");
    }
}

bool TractogramWriter::open() 
{
    if (!pImpl_) {
        disp(MSG_ERROR, "Cannot open %s: Writer implementation is null (unsupported format or creation failed).", filename_.c_str());
        return false;
    }
    if (is_open_) {
        disp(MSG_WARN, "File %s is already open.", filename_.c_str());
        return true;
    }
    if (is_closed_) {
        disp(MSG_ERROR, "File %s was already closed and cannot be reopened.", filename_.c_str());
        return false;
    }

    is_open_ = pImpl_->open();
    if (!is_open_) {
        disp(MSG_ERROR, "Failed to open file: %s", filename_.c_str());
    }

    // Reset state and start the background writer thread
    stop_requested_ = false;
    write_pending_  = false;
    writerThread_   = std::thread(&TractogramWriter::asyncWriteLoop, this);

    return is_open_;
}

void TractogramWriter::asyncWriteLoop()
{
    while (true) {
        std::unique_lock<std::mutex> lock(asyncMutex_);
        asyncCond_.wait(lock, [this] { 
            return write_pending_.load() || stop_requested_.load(); 
        });

        if (stop_requested_.load() && !write_pending_.load()) {
            break; // Correct way to exit the loop
        }

        lock.unlock();

        if (!writeAsyncBuffer_->empty()) {
            if (!pImpl_->writeBatch(*writeAsyncBuffer_)) {
                async_write_failed_ = true;
            }
        }
        
        lock.lock();
        writeAsyncBuffer_->clear(); // Clear the buffer after writing
        write_pending_ = false;
        asyncCond_.notify_one(); 
    }
}

bool TractogramWriter::writeStreamline(const NIBR::Streamline& streamline) 
{
    if (!is_open_ || is_closed_) {
        disp(MSG_ERROR, "Cannot write streamline: File not open or already closed.", filename_.c_str());
        return false;
    }

    std::unique_lock<std::mutex> lock(asyncMutex_);

    if (async_write_failed_.load()) {
        disp(MSG_ERROR, "Async writer failed. Aborting writeStreamline.");
        return false;
    }

    fillAsyncBuffer_->push_back(streamline);

    // If the fill buffer is full, swap and notify the writer
    if (fillAsyncBuffer_->size() >= WRITE_BUFFER_SIZE) {
        
        // Wait for the background writer to become free
        asyncCond_.wait(lock, [this]{ return !write_pending_.load(); });
        
        // Swap buffers
        std::swap(fillAsyncBuffer_, writeAsyncBuffer_);
        
        // Signal writer thread
        write_pending_ = true;
        asyncCond_.notify_one();
    }

    return true;
}

bool TractogramWriter::writeBatch(const StreamlineBatch& batch) 
{
    if (!is_open_ || is_closed_) {
        disp(MSG_ERROR, "Cannot write batch: File not open or already closed.", filename_.c_str());
        return false;
    }
    
    if (batch.empty()) return true;

    std::unique_lock<std::mutex> lock(asyncMutex_);

    if (async_write_failed_.load()) {
        disp(MSG_ERROR, "Async writer failed. Aborting writeBatch.");
        return false;
    }

    asyncCond_.wait(lock, [this]{ return !write_pending_.load(); });

    *fillAsyncBuffer_ = batch; 

    std::swap(fillAsyncBuffer_, writeAsyncBuffer_);
    
    write_pending_ = true;
    lock.unlock();
    asyncCond_.notify_one();

    return true;
}

bool TractogramWriter::writeBatch(StreamlineBatch&& batch)
{
    if (!is_open_ || is_closed_) {
        disp(MSG_ERROR, "Cannot write batch: File not open or already closed.", filename_.c_str());
        return false;
    }
    
    if (batch.empty()) return true;

    std::unique_lock<std::mutex> lock(asyncMutex_);

    if (async_write_failed_.load()) {
        disp(MSG_ERROR, "Async writer failed. Aborting writeBatch.");
        return false;
    }

    asyncCond_.wait(lock, [this]{ return !write_pending_.load(); });

    *fillAsyncBuffer_ = std::move(batch);

    std::swap(fillAsyncBuffer_, writeAsyncBuffer_);

    write_pending_ = true;
    lock.unlock();
    asyncCond_.notify_one();

    return true;
}

bool TractogramWriter::close()
{
    if (!is_open_ && is_closed_) return true;

    std::unique_lock<std::mutex> lock(asyncMutex_);

    // Immediately prevent any new writes from other threads
    is_open_ = false;

    // Check for prior failure before proceeding
    if (async_write_failed_.load()) {        
        disp(MSG_ERROR, "Cannot close file cleanly, an async write operation failed previously.");
    } else {

        // Before waiting, check if there's a final partial batch to write.
        if (!fillAsyncBuffer_->empty()) {
            
            // Wait for the writer to be free
            asyncCond_.wait(lock, [this]{ return !write_pending_.load(); });
            
            // Swap the final partial batch into the write buffer
            std::swap(fillAsyncBuffer_, writeAsyncBuffer_);
            
            // Signal this final write
            write_pending_ = true;
            asyncCond_.notify_one();
        }
    }

    // Now, wait for any pending write
    asyncCond_.wait(lock, [this]{ return !write_pending_.load(); });

    // Signal the thread to stop.
    stop_requested_ = true;
    lock.unlock();
    asyncCond_.notify_one();

    // Wait for the writer thread to terminate.
    if (writerThread_.joinable()) {
        writerThread_.join();
    }

    // Now it's safe to call the underlying implementation's close
    is_closed_ = pImpl_->close(finalStreamlineCount_, finalPointCount_);
    is_open_   = false; // Ensure state is consistent

    if (!is_closed_) {
        disp(MSG_ERROR, "Failed to properly close file: %s", filename_.c_str());
    } else {
        disp(MSG_DEBUG, "Successfully closed %s. Final Streamlines: %ld, Final Points: %ld", filename_.c_str(), finalStreamlineCount_, finalPointCount_);
    }

    // Report failure if either the underlying close failed OR an async write failed.
    return is_closed_ && !async_write_failed_.load();
}


// Helper free functions implementations
bool writeTractogram(std::string out_fname, NIBR::TractogramReader* reader) 
{

    if (!reader || !reader->isReady()) {
        disp(MSG_FATAL, "Input TractogramReader is not valid or not initialized.");
        return false;
    }

    TractogramWriter writer(out_fname);
    if (!writer.isValid())  return false;
    if (!writer.open())     return false;

    // Set TRK reference info if the output is TRK and reader has it
    if (getFileExtension(out_fname) == "trk") {
        disp(MSG_FATAL, "Missing reference image for trk output.");
        return false;
    }

    int batch_count = 0;
    reader->reset();

    while (true) {
        
        StreamlineBatch batch = reader->getNextStreamlineBatch(WRITE_BUFFER_SIZE);

        if (batch.empty()) break; 

        if (!writer.writeBatch(std::move(batch))) {
            disp(MSG_ERROR, "Failed to write batch %d to %s.", batch_count, out_fname.c_str());
            writer.close();
            return false;
        }

        batch_count++;
        if (VERBOSE_LEVEL() == VERBOSE_DEBUG) {
            if (batch_count % 10 == 0) {
                disp(MSG_DEBUG, "Written %d batches to %s...", batch_count, out_fname.c_str());
            }
        }
    }

    return writer.close();
}

bool writeTractogram(std::string out_fname, const Tractogram& tractogram) 
{
    
    if (getFileExtension(out_fname) == "trk") {
        disp(MSG_FATAL, "Missing reference image for trk output.");
        return false;
    }
    
    TractogramWriter writer(out_fname);
    if (!writer.isValid())  return false;
    if (!writer.open())     return false;

    const size_t single_batch_threshold = WRITE_BUFFER_SIZE * 5;

    if (tractogram.size() > single_batch_threshold) {
        disp(MSG_DEBUG, "Input tractogram is large (%zu streamlines), writing in chunks.", tractogram.size());
        StreamlineBatch current_batch;
        current_batch.reserve(WRITE_BUFFER_SIZE);
        for(size_t i = 0; i < tractogram.size(); ++i) {
            current_batch.push_back(tractogram[i]);
            if (current_batch.size() >= WRITE_BUFFER_SIZE) {
                if (!writer.writeBatch(std::move(current_batch))) {
                    writer.close(); return false;
                }
                current_batch.clear();
            }
        }
        if (!current_batch.empty()) {
            if (!writer.writeBatch(std::move(current_batch))) {
                 writer.close(); return false;
            }
        }
    } else {
         if (!tractogram.empty()) {
            if (!writer.writeBatch(tractogram)) {
                writer.close(); return false;
            }
        }
    }


    return writer.close();
}

template<typename T>
bool writeTractogram(std::string out_fname, const Tractogram& tractogram, const Image<T>& refImg) 
{
    
    if (getFileExtension(out_fname) != "trk") {
        disp(MSG_WARN, "Reference image provided for non-TRK file type (%s). It will be ignored.", out_fname.c_str());
        return writeTractogram(out_fname, tractogram);
    }

    TractogramWriter writer(out_fname);
    if (!writer.isValid()) return false;

    writer.setTRKReference(refImg);
    if (!writer.open())    return false;
    
    const size_t single_batch_threshold = WRITE_BUFFER_SIZE * 5;

    if (tractogram.size() > single_batch_threshold) {
        disp(MSG_DEBUG, "Input tractogram is large (%zu streamlines), writing in chunks.", tractogram.size());
        StreamlineBatch current_batch;
        current_batch.reserve(WRITE_BUFFER_SIZE);
        for(size_t i = 0; i < tractogram.size(); ++i) {
            current_batch.push_back(tractogram[i]);
            if (current_batch.size() >= WRITE_BUFFER_SIZE) {
                if (!writer.writeBatch(std::move(current_batch))) {
                    writer.close(); return false;
                }
                current_batch.clear();
            }
        }
        if (!current_batch.empty()) {
            if (!writer.writeBatch(std::move(current_batch))) {
                 writer.close(); return false;
            }
        }
    } else {
         if (!tractogram.empty()) {
            if (!writer.writeBatch(tractogram)) {
                writer.close(); return false;
            }
        }
    }

    return writer.close();
}

template bool writeTractogram<bool>         (std::string out_fname, const Tractogram& tractogram, const Image<bool>&          refImg);
template bool writeTractogram<int8_t>       (std::string out_fname, const Tractogram& tractogram, const Image<int8_t>&        refImg);
template bool writeTractogram<int16_t>      (std::string out_fname, const Tractogram& tractogram, const Image<int16_t>&       refImg);
template bool writeTractogram<int32_t>      (std::string out_fname, const Tractogram& tractogram, const Image<int32_t>&       refImg);
template bool writeTractogram<int64_t>      (std::string out_fname, const Tractogram& tractogram, const Image<int64_t>&       refImg);
template bool writeTractogram<uint8_t>      (std::string out_fname, const Tractogram& tractogram, const Image<uint8_t>&       refImg);
template bool writeTractogram<uint16_t>     (std::string out_fname, const Tractogram& tractogram, const Image<uint16_t>&      refImg);
template bool writeTractogram<uint32_t>     (std::string out_fname, const Tractogram& tractogram, const Image<uint32_t>&      refImg);
template bool writeTractogram<uint64_t>     (std::string out_fname, const Tractogram& tractogram, const Image<uint64_t>&      refImg);
template bool writeTractogram<float>        (std::string out_fname, const Tractogram& tractogram, const Image<float>&         refImg);
template bool writeTractogram<double>       (std::string out_fname, const Tractogram& tractogram, const Image<double>&        refImg);
template bool writeTractogram<long double>  (std::string out_fname, const Tractogram& tractogram, const Image<long double>&   refImg);


bool writeTractogram(std::string out_fname, const Tractogram& tractogram, const std::vector<TractogramField>& fields) 
{
    
    if (getFileExtension(out_fname) != "vtk" ){
        disp(MSG_ERROR, "Field data can be written only for vtk output.");
        return false;
    }

    TractogramWriter writer(out_fname);
    if (!writer.isValid()) return false;

    writer.setVTKFields(fields);
    if (!writer.open())    return false;

    const size_t single_batch_threshold = WRITE_BUFFER_SIZE * 5;

    if (tractogram.size() > single_batch_threshold) {
        disp(MSG_DEBUG, "Input tractogram is large (%zu streamlines), writing in chunks.", tractogram.size());
        StreamlineBatch current_batch;
        current_batch.reserve(WRITE_BUFFER_SIZE);
        for(size_t i = 0; i < tractogram.size(); ++i) {
            current_batch.push_back(tractogram[i]);
            if (current_batch.size() >= WRITE_BUFFER_SIZE) {
                if (!writer.writeBatch(std::move(current_batch))) {
                    writer.close(); return false;
                }
                current_batch.clear();
            }
        }
        if (!current_batch.empty()) {
            if (!writer.writeBatch(std::move(current_batch))) {
                 writer.close(); return false;
            }
        }
    } else {
         if (!tractogram.empty()) {
            if (!writer.writeBatch(tractogram)) {
                writer.close(); return false;
            }
        }
    }
    
    return writer.close();
}


bool writeTractogram(std::string out_fname, NIBR::TractogramReader* reader, const std::vector<size_t>& idx_to_keep_const) 
{
    
    if (!reader || !reader->isReady()) {
        disp(MSG_FATAL, "Input TractogramReader is not valid or not initialized.");
        return false;
    }

    if (idx_to_keep_const.empty()) {
        disp(MSG_WARN, "Empty index list provided to writeTractogram. Output file %s will be empty.", out_fname.c_str());
        // Create an empty file by opening and closing the writer
        TractogramWriter writer(out_fname);
        if (!writer.isValid()) return false;
        
        if (getFileExtension(out_fname) == "trk") {
            disp(MSG_FATAL, "Can't output selected streamlines with trk format.");
            return false;
        }
        if (!writer.open()) return false;
        return writer.close();
    }

    // Sort and unique indices for efficient processing if reader is not preloaded
    std::vector<size_t> idx_to_keep = idx_to_keep_const;
    std::sort(idx_to_keep.begin(), idx_to_keep.end());
    idx_to_keep.erase(std::unique(idx_to_keep.begin(), idx_to_keep.end()), idx_to_keep.end());

    TractogramWriter writer(out_fname);
    if (!writer.isValid()) return false;

    
    if (getFileExtension(out_fname) == "trk") {
        disp(MSG_FATAL, "Can't output selected streamlines with trk format.");
        return false;
    }
    // TODO: Handle VTK fields if subsetting. This is complex as field data also needs subsetting.
    // Current implementation does not transfer fields for indexed writing.
    if (!writer.open())   return false;

    reader->reset();

    if (reader->isPreloaded()) {
        StreamlineBatch out_batch;
        out_batch.reserve(idx_to_keep.size());
        for (size_t target_idx : idx_to_keep) {
            if (target_idx < reader->numberOfStreamlines) {
                out_batch.push_back(reader->getStreamline(target_idx));
            } else {
                disp(MSG_WARN, "Index %zu out of bounds (total streamlines: %zu). Skipping.", target_idx, reader->numberOfStreamlines);
            }
        }
        if (!out_batch.empty()) {
            if (!writer.writeBatch(std::move(out_batch))) {
                writer.close(); return false;
            }
        }
    } else {
        // Streaming mode: read batches and pick streamlines
        StreamlineBatch read_batch_buffer;
        StreamlineBatch write_batch_buffer;
        write_batch_buffer.reserve(WRITE_BUFFER_SIZE);

        size_t current_reader_streamline_idx = 0;
        size_t current_target_idx_pos = 0; // Position in sorted idx_to_keep

        while (current_target_idx_pos < idx_to_keep.size()) {
            size_t target_streamline_global_idx = idx_to_keep[current_target_idx_pos];

            // Read batches until we reach or pass the target_streamline_global_idx
            bool batch_read_ok = false;
            while (current_reader_streamline_idx + read_batch_buffer.size() <= target_streamline_global_idx) {
                current_reader_streamline_idx += read_batch_buffer.size();
                read_batch_buffer = reader->getNextStreamlineBatch(WRITE_BUFFER_SIZE);
                if (read_batch_buffer.size() <= 0) {
                    goto end_streaming_loop; // End of input or error
                }
                batch_read_ok = true;
            }
            if (!batch_read_ok && read_batch_buffer.empty()) { // If we didn't read a new batch and buffer is empty
                read_batch_buffer = reader->getNextStreamlineBatch(WRITE_BUFFER_SIZE);
                if (read_batch_buffer.size() <= 0) {
                    goto end_streaming_loop; // End of input or error
                }
            }


            // Process the current read_batch_buffer for any matching indices
            for (size_t i = 0; i < read_batch_buffer.size() && current_target_idx_pos < idx_to_keep.size(); ++i) {
                size_t streamline_in_batch_global_idx = current_reader_streamline_idx + i;
                if (streamline_in_batch_global_idx == idx_to_keep[current_target_idx_pos]) {
                    write_batch_buffer.push_back(read_batch_buffer[i]);
                    current_target_idx_pos++;
                    if (write_batch_buffer.size() >= WRITE_BUFFER_SIZE) {
                        if(!writer.writeBatch(std::move(write_batch_buffer))) { writer.close(); return false; }
                        write_batch_buffer.clear();
                    }
                } else if (streamline_in_batch_global_idx > idx_to_keep[current_target_idx_pos]) {
                    // This can happen if indices are skipped or out of order, but we sorted.
                    // It means the target index was missed or doesn't exist.
                    // This shouldn't happen if idx_to_keep is sorted and reader->numberOfStreamlines is respected.
                    disp(MSG_WARN, "Missed target index %zu while streaming. Current global index: %zu", idx_to_keep[current_target_idx_pos], streamline_in_batch_global_idx);
                    current_target_idx_pos++; // Move to next target to avoid infinite loop
                }
            }
        }
        end_streaming_loop:;
        if (!write_batch_buffer.empty()) {
            if(!writer.writeBatch(std::move(write_batch_buffer))) { writer.close(); return false; }
        }
        if (current_target_idx_pos < idx_to_keep.size()) {
            disp(MSG_WARN, "Not all requested streamlines were found/written from %s. Found %zu of %zu.",
                reader->fileName.c_str(), current_target_idx_pos, idx_to_keep.size());
        }
    }

    return writer.close();
}


bool writeTractogram(std::string out_fname, std::string inp_fname, const std::vector<size_t>& idx_to_keep) 
{
    NIBR::TractogramReader reader(inp_fname.c_str());
    if (!reader.isReady()) {
        disp(MSG_FATAL,"Failed to initialize reader for input file: %s", inp_fname.c_str());
        return false;
    }
    return writeTractogram(out_fname, &reader, idx_to_keep);
}

bool writeTractogram(std::string out_kept_fname, std::string out_rmvd_fname, std::string inp_fname, const std::vector<size_t>& idx_to_keep_const) 
{
    
    NIBR::TractogramReader reader(inp_fname.c_str()); // Open reader once
    if (!reader.isReady()) {
        disp(MSG_FATAL,"Failed to initialize reader for input file: %s", inp_fname.c_str());
        return false;
    }

    // Write kept streamlines
    disp(MSG_INFO, "Writing kept streamlines to: %s", out_kept_fname.c_str());
    if (!writeTractogram(out_kept_fname, &reader, idx_to_keep_const)) {
        disp(MSG_ERROR, "Failed to write kept streamlines.");
        return false;
    }
    disp(MSG_INFO, "Finished writing kept streamlines.");

    // Prepare indices for removed streamlines
    std::vector<bool> is_kept(reader.numberOfStreamlines, false);
    for (size_t original_idx : idx_to_keep_const) {
        if (original_idx < reader.numberOfStreamlines) {
            is_kept[original_idx] = true;
        }
    }

    std::vector<size_t> idx_to_remove;
    idx_to_remove.reserve(reader.numberOfStreamlines); // Max possible size
    for (size_t i = 0; i < reader.numberOfStreamlines; ++i) {
        if (!is_kept[i]) {
            idx_to_remove.push_back(i);
        }
    }
    
    // Reset reader to read again for the removed streamlines
    reader.reset(); 

    disp(MSG_INFO, "Writing removed streamlines to: %s", out_rmvd_fname.c_str());
    if (!writeTractogram(out_rmvd_fname, &reader, idx_to_remove)) {
        disp(MSG_ERROR, "Failed to write removed streamlines.");
        return false;
    }
    disp(MSG_INFO, "Finished writing removed streamlines.");

    return true;
}


}
