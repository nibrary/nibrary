#include "dMRI/tractography/io/tractogramReader.h"
#include "tractogramField.h"
#include <cstdint>
#include <cstring>
#include <iterator>

// Default buffer size for streaming mode
#define LOW_BUFFER_CAPACITY     200000
#define DISC_IO_BATCH_SIZE      20000
#define READER_BUFFER_SIZE       1 * 1024 * 1024
#define TCK_CHUNK_SIZE          16 * 1024 * 1024

using namespace NIBR;
			
NIBR::TractogramReader::TractogramReader(std::string _fileName, bool _preload) 
{
    initReader(_fileName, _preload);
}

NIBR::TractogramReader::~TractogramReader() 
{
    stop_producer = true;
    buffer_cv.notify_all();
    if (producerThread.joinable()) {
        producerThread.join();
    }
    if (file != nullptr) fclose(file);
    if (readerBuffer != nullptr) delete[] readerBuffer;
}

bool NIBR::TractogramReader::initReader(std::string _fileName, bool _preload) 
{

    if (isInitialized) return true;

    // On POSIX systems system, try to optimize file reading
    std::string extension = getFileExtension(_fileName);

    if ((extension == "tck") || (extension == "trk")){
        prepSequentialRead(_fileName);
    } else if (extension == "vtk") {
        prepRandomRead(_fileName);
    }


	fileName 	= _fileName;
	file 		= fopen(fileName.c_str(), "rb");

	if (file == nullptr) {
        disp(MSG_ERROR, "Failed to open file: %s", fileName.c_str());
        return false;
    }

    // Set up the file I/O buffer
    readerBuffer = new char[READER_BUFFER_SIZE];
    if (setvbuf(file, readerBuffer, _IOFBF, READER_BUFFER_SIZE) != 0) {
        delete[] readerBuffer;
        readerBuffer = nullptr;
    }

    
    // Parse header
	const std::size_t strLength = 256;
	char dummy[strLength];

	disp(MSG_DEBUG,"File extension: %s", extension.c_str());

	numberOfStreamlines = 0;

	if (extension == "tck") {

		fileFormat 		= TCK;

        fgets(dummy, strLength, file);
        fileDescription = std::string(dummy);
        fileDescription.pop_back();

		std::string tmps;
		long pos = 0;

		do {
			fgets(dummy, strLength, file);
			tmps = std::string(dummy);

			std::size_t column = tmps.find_last_of(":");

			if (column != std::string::npos) {
				std::string key = tmps.substr(0, column);

				if (key == "count") {	
					numberOfStreamlines = std::atoi((tmps.substr(column + 1)).c_str());
				} 
				if (key == "file") {
					pos = std::atoi((tmps.substr(column + 3)).c_str());
				} 
			}

		} while (tmps != "END\n");

		currentStreamlinePos = pos;
		firstStreamlinePos   = currentStreamlinePos;
		// List of points until NaN NaN NaN

	} else 
    
    if (extension == "trk") {
		
		fileFormat 		= TRK;
		fileDescription = "trk file";

		trkFileStruct trkFileHeader;

		std::fread(&trkFileHeader, sizeof(trkFileStruct), 1, file);

		for(int i=0;i<4;++i)
			for(int j=0;j<4;++j)
				ijk2xyz[i][j] = trkFileHeader.vox_to_ras[i][j];

		for(int i=0;i<3;++i) {
			imgDims[i] = trkFileHeader.dim[i];
			pixDims[i] = trkFileHeader.voxel_size[i];
		}
		std::strcpy(voxOrdr,"LAS");

		bool allZero = true;

		for(int i=0;i<3;++i)
			for(int j=0;j<3;++j)
				if (ijk2xyz[i][j] != 0)
					allZero = false;

		if ( allZero || (trkFileHeader.version != 2) ) {
			disp(MSG_FATAL,"Trk file version must be 2. Old trk files are not accepted.");
			return false;
		}

		inverseAffine(ijk2xyz, xyz2ijk);

		n_scalars_trk 		= trkFileHeader.n_scalars;
		n_properties_trk	= trkFileHeader.n_properties;
		numberOfStreamlines	= trkFileHeader.n_count;

		currentStreamlinePos = ftell(file);
		firstStreamlinePos   = currentStreamlinePos;
		// First value is int, which is the length of streamlines
		// Then comes points

	} else

	if (extension == "vtk") {

		std::string vtkFormat;
		
		int majorVersion = 0, minorVersion = 0;
		std::fscanf(file, "# vtk DataFile Version %d.%d\n", &majorVersion, &minorVersion);
		disp(MSG_DEBUG,"Vtk version is: %d.%d", majorVersion,minorVersion);

		vtk_isV5 = (majorVersion >= 5);

		fgets(dummy, strLength, file);                        // file description
		fileDescription = std::string(dummy);
		fileDescription = fileDescription.substr(0, fileDescription.size() - 1);
		std::fscanf(file, "%s ", dummy);                      // ascii or binary
		vtkFormat = std::string(dummy);
		fgets(dummy, strLength, file);                        // always DATASET POLYDATA, we skip to check this for now
		
		std::size_t totalNumberOfPoints;
		std::fscanf(file, "%*s %zu %*s ", &totalNumberOfPoints);   // number of points and datatype, we assume datatype is float and skip checking this
		disp(MSG_DEBUG,"totalNumberOfPoints: %d", totalNumberOfPoints);

		currentStreamlinePos = ftell(file);
		firstStreamlinePos   = currentStreamlinePos;


		if ((vtkFormat == "ascii") || (vtkFormat == "ASCII"))        {if(majorVersion == 3) fileFormat = VTK_ASCII_3;  if(majorVersion == 4) fileFormat = VTK_ASCII_4;  if(majorVersion == 5) fileFormat = VTK_ASCII_5; }
		else if ((vtkFormat == "binary") || (vtkFormat == "BINARY")) {if(majorVersion == 3) fileFormat = VTK_BINARY_3; if(majorVersion == 4) fileFormat = VTK_BINARY_4; if(majorVersion == 5) fileFormat = VTK_BINARY_5;}
		else return false;

		// Skip points
		if ((fileFormat == VTK_BINARY_3) || (fileFormat == VTK_BINARY_4) || (fileFormat == VTK_BINARY_5)){
			std::fseek(file, sizeof(float) * totalNumberOfPoints * 3, SEEK_CUR);
			int tmp = std::fgetc(file); if (tmp != '\n') std::ungetc(tmp, file); // Make sure to go end of the line
		}

		if ((fileFormat == VTK_ASCII_3) || (fileFormat == VTK_ASCII_4) || (fileFormat == VTK_ASCII_5)){
			for (std::size_t i = 0; i < totalNumberOfPoints; i++)
				std::fscanf(file, "%*f %*f %*f ");
		}

		vtk_needsByteSwap = is_little_endian();

		if (vtk_needsByteSwap) disp(MSG_DEBUG, "Swapping bytes");
			
		if (totalNumberOfPoints > 0) {

			// --- Seek to LINES and parse ---
			while (fgets(dummy, strLength, file)) {
				if (strncmp(dummy, "LINES", 5) == 0) break;
				if (feof(file)) { disp(MSG_ERROR, "EOF reached before finding LINES."); return false; }
			}

			if (sscanf(dummy, "LINES %zu %*u", &numberOfStreamlines) != 1) {
				disp(MSG_ERROR, "Failed to parse LINES header: %s", dummy);
				return false;
			}

			if (vtk_isV5) {

                disp(MSG_DEBUG,"Parsing v5");

                numberOfStreamlines -= 1;

				if (numberOfStreamlines < 1) {
					disp(MSG_DEBUG,"No lines in file");
					numberOfStreamlines = 0;
				} else {

					std::size_t numberOfOffsets = numberOfStreamlines + 1;

					disp(MSG_DEBUG,"Streamlines: %zu, Connectivity: %zu, Offsets: %zu", numberOfStreamlines, totalNumberOfPoints, numberOfOffsets);

					// --- Handle OFFSETS ---
					vtk5_64bitOffsets = false;
					while (fgets(dummy, strLength, file)) {
						if (strncmp(dummy, "OFFSETS", 7) == 0) break;
						if (feof(file)) { disp(MSG_ERROR, "EOF reached before finding OFFSETS."); return false; }
					}
					if (strstr(dummy, "vtktypeint64")) {
						vtk5_64bitOffsets = true;
						disp(MSG_DEBUG, "Offsets are 64-bit.");
					} else if (strstr(dummy, "vtktypeint32")) {
						vtk5_64bitOffsets = false;
						disp(MSG_DEBUG, "Offsets are 32-bit.");
					} else {
						disp(MSG_ERROR, "Could not determine OFFSETS data type in V5 file: %s", dummy);
						return false;
					}

					// Dummy read for the first offset, which is always 0.
					if ((fileFormat == VTK_BINARY_3) || (fileFormat == VTK_BINARY_4) || (fileFormat == VTK_BINARY_5)){
						if (vtk5_64bitOffsets) {
							int64_t offset;
							std::fread(&offset, sizeof(int64_t), 1, file);
							if (vtk_needsByteSwap) swapByteOrder(offset);
							vtk_prevOffset = offset;
						} else {
							int32_t offset;
							std::fread(&offset, sizeof(int32_t), 1, file);
							if (vtk_needsByteSwap) swapByteOrder(offset);
							vtk_prevOffset  = offset;
						}
						// Consume potential newline after binary read
						int tmp_char = std::fgetc(file); if (tmp_char != '\n' && tmp_char != EOF) std::ungetc(tmp_char, file);
					} else {
						int64_t offset;
						std::fscanf(file, "%lld", (long long*)&offset);
						vtk_prevOffset = offset;
					}

					vtk_currOffsetPos  = ftell(file);
					vtk_firstOffsetPos = vtk_currOffsetPos;
					vtk_firstOffset = vtk_prevOffset;
				}

            } else {

				if (numberOfStreamlines < 1) {
					disp(MSG_DEBUG,"No lines in file");
					numberOfStreamlines = 0;
				} else {
					disp(MSG_DEBUG,"numberOfStreamlines: %d", numberOfStreamlines);
					vtk_currOffsetPos  = ftell(file);
					vtk_firstOffsetPos = vtk_currOffsetPos;
				}

			}

		} else {
			numberOfStreamlines = 0;
		}


	} else {
        disp(MSG_ERROR,"Unknown file extension. Only .vtk, .tck and .trk extensions are supported.");
        isInitialized = false;
        return false;
    }

    currentStreamlinePos = firstStreamlinePos;

    // Configure the buffer and producer thread
    isPreloadMode           = _preload;


    if (isPreloadMode) {
        buffer_capacity       = numberOfStreamlines;
        buffer_low_water_mark = numberOfStreamlines;
    } else {
        buffer_low_water_mark = LOW_BUFFER_CAPACITY;
        buffer_capacity       = buffer_low_water_mark + DISC_IO_BATCH_SIZE + 1;
    }

    stop_producer               = false;
    streamlines_read_from_file  = 0;
    consumed_streamline_count   = 0;
    producer_finished           = false;
    producerThread              = std::thread(&TractogramReader::producerLoop, this);
    isInitialized               = true;

	return true;
}

void NIBR::TractogramReader::reset()
{
    disp(MSG_DEBUG, "Reseting tractogramReader.");

    if (!isInitialized) {
        disp(MSG_ERROR, "Cannot reset a reader that was not successfully initialized.");
        return;
    }

    if (isPreloadMode) {
        
        disp(MSG_DEBUG, "Performing fast in-memory reset for preload mode.");
        std::lock_guard<std::mutex> lock(buffer_mutex);
        consumed_streamline_count = 0;

    } else {

        // Signal the producer to stop its loop
        disp(MSG_DEBUG, "Signaled the producer to stop its loop");
        stop_producer = true;
        
        // Notify the producer in case it is sleeping in a wait() call
        disp(MSG_DEBUG, "Notified the producer in case it is sleeping in a wait() call");
        buffer_cv.notify_all();

        // Wait for the producer thread to actually finish and exit
        disp(MSG_DEBUG, "Wait for the producer thread to actually finish and exit...");
        if (producerThread.joinable()) {
            producerThread.join();
        }
        disp(MSG_DEBUG, "Done");

        // Now that the producer is stopped, we can safely reset all state variables.
        {
            // Lock the mutex to ensure no consumer thread can access the buffer while we clear it.
            disp(MSG_DEBUG, "Reset all state variables");
            std::lock_guard<std::mutex> lock(buffer_mutex);
            streamline_buffer.clear();
        }

        // Reset all counters and position trackers
        streamlines_read_from_file = 0;
        consumed_streamline_count  = 0;
        producer_finished          = false;
        currentStreamlinePos       = firstStreamlinePos; // Go back to the start of streamline data

        // Reset VTK-specific variables
        if ((fileFormat >= VTK_ASCII_3) && (fileFormat <= VTK_BINARY_5)) {
            vtk_currOffsetPos = vtk_firstOffsetPos;
            vtk_prevOffset    = vtk_firstOffset;
        }

        // Reset TCK-specific variables
        tck_read_buffer.clear();
        tck_buffer_offset = 0;

        // Physically seek the file pointer back to the start of the streamline data
        disp(MSG_DEBUG, "Seek the file pointer back to the start of the streamline data");
        fseek(file, firstStreamlinePos, SEEK_SET);

        // With a clean state, we can now launch a new producer thread. Allow the new producer to run
        disp(MSG_DEBUG, "Launch a new producer thread");
        stop_producer = false;

        // Launch the new producer, which will start pre-buffering from the beginning
        producerThread = std::thread(&TractogramReader::producerLoop, this);
    }

    disp(MSG_DEBUG, "TractogramReader has been reset.");
}

void NIBR::TractogramReader::producerLoop() 
{

    disp(MSG_DEBUG, "Started producer.");

    while (!stop_producer) {

        disp(MSG_DEBUG, "Buffering...");
        
        if (!isPreloadMode) {
            std::unique_lock<std::mutex> lock(buffer_mutex);
            buffer_cv.wait(lock, [this] {
                return (streamline_buffer.size() <= buffer_low_water_mark) || stop_producer;
            });
        }

        if (stop_producer) break;

        if (streamlines_read_from_file >= numberOfStreamlines) break;

        size_t batchReadSize = 0;

        {
            std::lock_guard<std::mutex> lock(buffer_mutex);
            batchReadSize = buffer_capacity - streamline_buffer.size();
        }
        
        // If batchReadSize is 0, it means the buffer is full. 
        // In preload mode, this could cause a busy-loop. We yield to prevent it.
        if (batchReadSize == 0) {
             // Give other threads a chance to run.
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            continue; // Go back to the top of the while loop to re-check everything
        }

        StreamlineBatch batch = readBatchFromFile(std::min(batchReadSize, size_t(DISC_IO_BATCH_SIZE)));
        
        if (batch.empty()) break;

        {
            std::lock_guard<std::mutex> lock(buffer_mutex);
            for (auto& s : batch) {
                streamline_buffer.push_back(std::move(s));
            }
            disp(MSG_DEBUG, "%d streamlines buffered", batch.size());
        }
        buffer_cv.notify_all();
    }

    disp(MSG_DEBUG, "Producer finished.");

    producer_finished = true;
    buffer_cv.notify_all();
}


std::tuple<bool, Streamline, size_t> NIBR::TractogramReader::getNextStreamline() 
{
    std::unique_lock<std::mutex> lock(buffer_mutex);

    // Wait for buffer has data, or producer is finished.
    buffer_cv.wait(lock, [this] {
        return !streamline_buffer.empty() || producer_finished;
    });

    // If we woke up and the buffer is still empty, it must mean the producer is finished.
    if (streamline_buffer.empty()) {
        return {false, Streamline(), 0};
    }

    size_t streamline_idx = consumed_streamline_count++;

    if (isPreloadMode) {
        return {true, streamline_buffer.at(streamline_idx), streamline_idx};
    } else {
        Streamline s = std::move(streamline_buffer.front());
        streamline_buffer.pop_front();
        
        lock.unlock();
        buffer_cv.notify_one(); // Notify producer that space is available

        return {true, std::move(s), streamline_idx};
    }
}


StreamlineBatch NIBR::TractogramReader::getNextStreamlineBatch(size_t batchSize)
{
    StreamlineBatch batch_out;

    if (batchSize == 0) return batch_out;

    batch_out.reserve(batchSize);

    for (size_t i = 0; i < batchSize; ++i) {
        
        auto [success, streamline, streamline_idx] = this->getNextStreamline();

        if (success) {
            batch_out.push_back(std::move(streamline));
        } else {
            break; 
        }
    }

    return batch_out;
}


const Streamline& NIBR::TractogramReader::getStreamline(size_t n) 
{
    if (!isPreloadMode) {
        disp(MSG_FATAL, "getStreamline(n) is only available in preload mode.");
        static const Streamline empty; return empty;
    }
    if (n >= numberOfStreamlines) {
        disp(MSG_FATAL, "Streamline index %zu is out of bounds.", n);
        static const Streamline empty; return empty;
    }
    
    std::unique_lock<std::mutex> lock(buffer_mutex);
    buffer_cv.wait(lock, [this, n]{ return n < streamline_buffer.size() || producer_finished; });

    if (n >= streamline_buffer.size()) {
        disp(MSG_FATAL, "Producer finished but streamline index %zu not found.", n);
        static const Streamline empty; return empty;
    }
    
    return streamline_buffer.at(n);
}


// It is NOT thread-safe by itself and must be called from the producerLoop.
StreamlineBatch NIBR::TractogramReader::readBatchFromFile(size_t batchSize) 
{
    disp(MSG_DEBUG,"Reading batch from file...");

    StreamlineBatch batch_out;

    if (batchSize == 0 || streamlines_read_from_file >= numberOfStreamlines) {
        return batch_out;
    }

    size_t actualBatchSize = std::min(batchSize, numberOfStreamlines - streamlines_read_from_file.load());
    if (actualBatchSize == 0) return batch_out;
    batch_out.reserve(actualBatchSize);

    switch (fileFormat) {

        case TCK:
        {

            // This is the main loop to build the requested batch of streamlines.
            for (std::size_t i = 0; i < actualBatchSize; ++i) {

                Streamline streamline;
                bool endOfFile = false;

                // This loop builds a single streamline point by point.
                while (true) {

                    // Check if our current in-memory buffer has enough data for one point.
                    if ((tck_read_buffer.size() - tck_buffer_offset) < sizeof(Point3D)) {

                        // --- Buffer is exhausted, time to read a new chunk from the disk ---

                        // 1. Preserve any leftover bytes from the previous read.
                        size_t leftover_bytes = tck_read_buffer.size() - tck_buffer_offset;
                        if (leftover_bytes > 0) {
                            std::memmove(tck_read_buffer.data(), tck_read_buffer.data() + tck_buffer_offset, leftover_bytes);
                        }
                        tck_read_buffer.resize(leftover_bytes);
                        tck_buffer_offset = 0;

                        // 2. Read the next large chunk from the file.
                        // We seek first to ensure our file position is correct.
                        fseek(file, currentStreamlinePos, SEEK_SET);
                        size_t current_buffer_size = tck_read_buffer.size();
                        tck_read_buffer.resize(current_buffer_size + TCK_CHUNK_SIZE);
                        
                        size_t bytes_read = fread(tck_read_buffer.data() + current_buffer_size, 1, TCK_CHUNK_SIZE, file);
                        
                        // Update our master file position and resize buffer to actual bytes read.
                        currentStreamlinePos += bytes_read;
                        tck_read_buffer.resize(current_buffer_size + bytes_read);

                        // 3. If after reading, we STILL don't have enough for a point, we've hit the end.
                        if (tck_read_buffer.size() < sizeof(Point3D)) {
                            endOfFile = true;
                            break;
                        }
                    }

                    // --- We have enough data in the buffer, so we can parse a point ---

                    Point3D p;
                    // Copy a point directly from our fast in-memory buffer.
                    memcpy(p.data(), tck_read_buffer.data() + tck_buffer_offset, sizeof(Point3D));
                    
                    // Advance our position within the buffer.
                    tck_buffer_offset += sizeof(Point3D);

                    // TCK uses NaN to separate streamlines and Inf to terminate the file.
                    if (std::isnan(p[0]) || std::isinf(p[0])) {
                        break; // End of this streamline.
                    }
                    streamline.emplace_back(p);
                }

                if (!streamline.empty()) {
                    batch_out.emplace_back(std::move(streamline));
                }

                if (endOfFile) {
                    disp(MSG_DEBUG,"Reached end of TCK file.");
                    break; // Stop trying to read more streamlines for this batch.
                }
            }
            break;
        }

        case TRK:
        {
            // Read All at Once: Perform a pre-pass to find the total byte size of the batch,
            // read the entire block into memory, then parse the memory buffer.
            // This is much faster for I/O than many small reads.

            long batchStartPosition = ftell(file);
            std::vector<int32_t> batch_lengths;
            batch_lengths.reserve(actualBatchSize);

            // 1. Pre-pass to get lengths and calculate total byte size of the batch
            for (std::size_t i = 0; i < actualBatchSize; ++i) {
                int32_t len = 0;
                if (fread(&len, sizeof(int32_t), 1, file) < 1) break;
                batch_lengths.push_back(len);
                // Seek past the data for this streamline to get to the next length field
                fseek(file, len * (sizeof(Point3D) + n_scalars_trk * sizeof(float)), SEEK_CUR);
            }
            long batchEndPosition = ftell(file);
            uint64_t totalBytesToRead = batchEndPosition - batchStartPosition;

            // 2. Single bulk read of the entire batch data
            std::vector<char> rawBuffer(totalBytesToRead);
            fseek(file, batchStartPosition, SEEK_SET);
            fread(rawBuffer.data(), 1, totalBytesToRead, file);

            // 3. Parse the in-memory buffer
            char* bufferPtr = rawBuffer.data();
            for (int32_t len : batch_lengths) {
                
                // The length is in the buffer, but we already have it. Skip it.
                bufferPtr += sizeof(int32_t);

                Streamline streamline;
                streamline.reserve(len);
                for (int j = 0; j < len; ++j) {
                    Point3D p_vox, p_world;

                    // Copy point from buffer instead of reading from file
                    std::memcpy(p_vox.data(), bufferPtr, sizeof(Point3D));
                    bufferPtr += sizeof(Point3D);

                    p_vox[0] -= 0.5f;
                    p_vox[1] -= 0.5f;
                    p_vox[2] -= 0.5f;
                    applyTransform(p_world, p_vox, ijk2xyz);
                    streamline.emplace_back(p_world);

                    // Skip scalars in the buffer
                    if (n_scalars_trk > 0) {
                        bufferPtr += n_scalars_trk * sizeof(float);
                    }
                }
                batch_out.push_back(std::move(streamline));
            }

            currentStreamlinePos = ftell(file);
            break;
        }

        case VTK_BINARY_3:
        case VTK_BINARY_4:
        case VTK_BINARY_5:
        {
            // Step 1: Read streamline lengths for the batch.
            std::vector<uint32_t> batch_lengths;
            batch_lengths.reserve(actualBatchSize);
            fseek(file, vtk_currOffsetPos, SEEK_SET);

            if (vtk_isV5) {
                std::vector<int64_t> offsets(actualBatchSize);
                if (vtk5_64bitOffsets) {
                    fread(offsets.data(), sizeof(int64_t), actualBatchSize, file);
                    if (vtk_needsByteSwap) for(auto& v : offsets) swapByteOrder(v);
                } else {
                    std::vector<int32_t> offsets32(actualBatchSize);
                    fread(offsets32.data(), sizeof(int32_t), actualBatchSize, file);
                    if (vtk_needsByteSwap) for(auto& v : offsets32) swapByteOrder(v);
                    std::copy(offsets32.begin(), offsets32.end(), offsets.begin());
                }

                int64_t last_offset = vtk_prevOffset;
                for (int64_t off : offsets) {
                    batch_lengths.push_back(static_cast<uint32_t>(off - last_offset));
                    last_offset = off;
                }
                vtk_prevOffset = last_offset;
            } else {
                for (std::size_t i = 0; i < actualBatchSize; ++i) {
                    uint32_t len;
                    fread(&len, sizeof(uint32_t), 1, file);
                    if (vtk_needsByteSwap) swapByteOrder(len);
                    batch_lengths.push_back(len);
                    fseek(file, len * sizeof(uint32_t), SEEK_CUR); // Skip connectivity indices
                }
            }

            vtk_currOffsetPos = ftell(file);

            // Step 2: Read all points for the batch in one go.
            uint32_t totalPointsInBatch = 0;
            for (uint32_t len : batch_lengths) totalPointsInBatch += len;

            std::vector<Point3D> pointBuffer;
            if (totalPointsInBatch > 0) {
                pointBuffer.resize(totalPointsInBatch);
                fseek(file, currentStreamlinePos, SEEK_SET);
                fread(pointBuffer.data(), sizeof(Point3D), totalPointsInBatch, file);
                currentStreamlinePos = ftell(file);

                if (vtk_needsByteSwap) {
                    for (auto& p : pointBuffer) {
                        swapByteOrder(p[0]); swapByteOrder(p[1]); swapByteOrder(p[2]);
                    }
                }
            }

            // Step 3: Distribute points into the output batch.
            int pointBufferOffset = 0;
            for (uint32_t len : batch_lengths) {
                Streamline streamline(pointBuffer.begin() + pointBufferOffset, pointBuffer.begin() + pointBufferOffset + len);
                batch_out.push_back(std::move(streamline));
                pointBufferOffset += len;
            }

            currentStreamlinePos = ftell(file);
            break;
        }

        case VTK_ASCII_3:
        case VTK_ASCII_4:
        case VTK_ASCII_5:
        {

            // Step 1: Read lengths
            std::vector<uint32_t> batch_lengths;
            batch_lengths.reserve(actualBatchSize);
            fseek(file, vtk_currOffsetPos, SEEK_SET);

            if (vtk_isV5) {
                int64_t last_offset = vtk_prevOffset;
                for (std::size_t i = 0; i < actualBatchSize; ++i) {
                    int64_t off = 0;
                    fscanf(file, "%lld", (long long*)&off);
                    batch_lengths.push_back(static_cast<uint32_t>(off - last_offset));
                    last_offset = off;
                }
                vtk_prevOffset = last_offset;

            } else {
                for (std::size_t i = 0; i < actualBatchSize; ++i) {
                    uint32_t len = 0;
                    fscanf(file, "%u ", &len);
                    batch_lengths.push_back(len);
                    for (uint32_t j = 0; j < len; ++j) fscanf(file, "%*d ");
                }
            }
            vtk_currOffsetPos = ftell(file);

            // Step 2: Read points
            uint32_t totalPointsInBatch = 0;
            for (uint32_t len : batch_lengths) totalPointsInBatch += len;
            
            std::vector<Point3D> pointBuffer;
            if (totalPointsInBatch > 0) {
                pointBuffer.resize(totalPointsInBatch);
                fseek(file, currentStreamlinePos, SEEK_SET);
                
                // Read points for this batch
                for (uint32_t i = 0; i < totalPointsInBatch; ++i) {
                    fscanf(file, "%f %f %f ", &pointBuffer[i][0], &pointBuffer[i][1], &pointBuffer[i][2]);
                }
                currentStreamlinePos = ftell(file);
                
            }
            
            // Step 3: Distribute
            int pointBufferOffset = 0;
            for (uint32_t len : batch_lengths) {
                Streamline streamline(pointBuffer.begin() + pointBufferOffset, pointBuffer.begin() + pointBufferOffset + len);
                batch_out.push_back(std::move(streamline));
                pointBufferOffset += len;
            }
            currentStreamlinePos = ftell(file);
            break;
        }

        default:
            // Unknown format, do nothing.
            break;
    }

    // Update the count of streamlines read from the file
    streamlines_read_from_file += batch_out.size();

    disp(MSG_DEBUG,"Done");

    return batch_out;
}


Tractogram NIBR::TractogramReader::getTractogram() {
    Tractogram all_streamlines;
    all_streamlines.reserve(numberOfStreamlines);

    reset();

    if (isPreloadMode) {
        std::unique_lock<std::mutex> lock(buffer_mutex);
        buffer_cv.wait(lock, [this] { return producer_finished && streamline_buffer.size() == numberOfStreamlines; });
        all_streamlines.assign(streamline_buffer.begin(), streamline_buffer.end());
        return all_streamlines;
    }

    while (true) {
        auto [success, s, idx] = getNextStreamline();
        if (!success) break;
        all_streamlines.push_back(std::move(s));
    }
    return all_streamlines;
}

const std::vector<uint64_t>& NIBR::TractogramReader::getNumberOfPoints()
{
    if (!numberOfPoints.empty()) return numberOfPoints;

    numberOfPoints.resize(numberOfStreamlines + 1, 0);
    
    // In preload mode, we need the full buffer to calculate this from memory
    if (isPreloadMode) {
        std::unique_lock<std::mutex> lock(buffer_mutex);
        buffer_cv.wait(lock, [this] { return producer_finished && streamline_buffer.size() == numberOfStreamlines; });
        
        uint64_t cumulativeSum = 0;
        for (size_t i = 0; i < numberOfStreamlines; ++i) {
            cumulativeSum           += streamline_buffer.at(i).size();
            numberOfPoints[i + 1]    = cumulativeSum;
        }
        return numberOfPoints;
    }

    // We will open a new file for this purpose
    auto numFile = fopen(fileName.c_str(), "rb");
    if (numFile == nullptr) {
        disp(MSG_ERROR, "Failed to open file %s.", fileName.c_str());
        return numberOfPoints;
    }

    // Use a 64-bit integer for the cumulative sum to prevent potential overflow
    // with very large tractograms before casting to uint64_t for each element.
    uint64_t cumulativeSum = 0;

    switch (fileFormat) {

        case TRK:
        {
            fseek(numFile, firstStreamlinePos, SEEK_SET);
            for (std::size_t i = 0; i < numberOfStreamlines; ++i) {
                int32_t streamlineLength = 0;
                if (fread(&streamlineLength, sizeof(int32_t), 1, numFile) != 1) {
                    disp(MSG_ERROR, "Failed to read streamline length from .trk file.");
                    numberOfPoints.clear(); // Invalidate results on error
                    break;
                }
                cumulativeSum           += streamlineLength;
                numberOfPoints[i + 1]    = cumulativeSum;
                
                // Seek past the point data and associated scalars to the next length header.
                fseek(numFile, (long)streamlineLength * (sizeof(Point3D) + n_scalars_trk * sizeof(float)), SEEK_CUR);
            }
            break;
        }

        case TCK:
        {
            // TCK format requires reading through the file to count points, as lengths are not stored explicitly.
            fseek(numFile, firstStreamlinePos, SEEK_SET);
            for (std::size_t i = 0; i < numberOfStreamlines; ++i) {
                uint64_t pointsInStreamline = 0;
                while (true) {
                    float p[3];
                    if (fread(p, sizeof(float), 3, numFile) < 3) { // End of file
                        break;
                    }
                    if (std::isnan(p[0])) { // End of streamline marker
                        break;
                    }
                    pointsInStreamline++;
                }
                cumulativeSum           += pointsInStreamline;
                numberOfPoints[i + 1]    = cumulativeSum;
            }
            break;
        }

        case VTK_BINARY_5:
        case VTK_ASCII_5:
        {
            // Modern VTK (v5+) uses an OFFSETS array, which is already a cumulative sum.
            fseek(numFile, vtk_firstOffsetPos, SEEK_SET);
            
            // The first offset (0) is handled by vector initialization. We read the remaining N offsets.
            if (fileFormat == VTK_BINARY_5) {
                if (vtk5_64bitOffsets) {
                    std::vector<int64_t> offsets(numberOfStreamlines);
                    if (fread(offsets.data(), sizeof(int64_t), numberOfStreamlines, numFile) == numberOfStreamlines) {
                        if (vtk_needsByteSwap) for(auto& v : offsets) swapByteOrder(v);
                        for(size_t i=0; i < numberOfStreamlines; ++i) numberOfPoints[i+1] = offsets[i];
                    }
                } else { // 32-bit offsets
                    std::vector<int32_t> offsets(numberOfStreamlines);
                    if (fread(offsets.data(), sizeof(int32_t), numberOfStreamlines, numFile) == numberOfStreamlines) {
                        if (vtk_needsByteSwap) for(auto& v : offsets) swapByteOrder(v);
                        for(size_t i=0; i < numberOfStreamlines; ++i) numberOfPoints[i+1] = offsets[i];
                    }
                }
            } else { // ASCII
                for (std::size_t i = 0; i < numberOfStreamlines; ++i) {
                    int64_t offset_val = 0;
                    if (fscanf(numFile, "%lld", (long long*)&offset_val) != 1) {
                         disp(MSG_ERROR, "Failed to read offset from ASCII VTK file.");
                         numberOfPoints.clear(); break;
                    }
                    numberOfPoints[i + 1] = offset_val;
                }
            }
            break;
        }

        case VTK_BINARY_3:
        case VTK_BINARY_4:
        case VTK_ASCII_3:
        case VTK_ASCII_4:
        {
            // Legacy VTK stores the length of each line individually.
            fseek(numFile, vtk_firstOffsetPos, SEEK_SET);
            bool is_binary = (fileFormat == VTK_BINARY_3 || fileFormat == VTK_BINARY_4);

            for (std::size_t i = 0; i < numberOfStreamlines; ++i) {
                uint32_t streamlineLength = 0;
                if (is_binary) {
                    if (fread(&streamlineLength, sizeof(uint32_t), 1, numFile) != 1) {
                        disp(MSG_ERROR, "Failed to read streamline length from binary VTK file.");
                        numberOfPoints.clear(); break;
                    }
                    if (vtk_needsByteSwap) swapByteOrder(streamlineLength);
                    // Skip connectivity indices
                    fseek(numFile, (long)streamlineLength * sizeof(uint32_t), SEEK_CUR);
                } else { // ASCII
                    if (fscanf(numFile, "%u", &streamlineLength) != 1) {
                        disp(MSG_ERROR, "Failed to read streamline length from ASCII VTK file.");
                        numberOfPoints.clear(); break;
                    }
                    // Skip connectivity indices
                    for (uint32_t j = 0; j < streamlineLength; ++j) fscanf(numFile, "%*d");
                }
                cumulativeSum           += streamlineLength;
                numberOfPoints[i + 1]    = cumulativeSum;
            }
            break;
        }
        
        default:
            disp(MSG_ERROR,"Unknown file format.");
            numberOfPoints.clear();
            break;
    }

    // Restore the file pointer to its original position.
    fclose(numFile);

    return numberOfPoints;
}

void NIBR::TractogramReader::printInfo() {

	disp(MSG_INFO,"Tractogram info");
	std::cout << "\033[32m";

	std::cout << "File name: " << fileName << std::endl;

	std::cout << "Format:                   ";
	if (fileFormat==NIBR::TCK)               std::cout << "tck"             << std::endl << std::flush;
	else if (fileFormat==NIBR::TRK)          std::cout << "trk"             << std::endl << std::flush;
	else if (fileFormat==NIBR::VTK_ASCII_3)  std::cout << "vtk v3 (ascii)"  << std::endl << std::flush;
	else if (fileFormat==NIBR::VTK_BINARY_3) std::cout << "vtk v3 (binary)" << std::endl << std::flush;
    else if (fileFormat==NIBR::VTK_ASCII_4)  std::cout << "vtk v4 (ascii)"  << std::endl << std::flush;
	else if (fileFormat==NIBR::VTK_BINARY_4) std::cout << "vtk v4 (binary)" << std::endl << std::flush;
    else if (fileFormat==NIBR::VTK_ASCII_5)  std::cout << "vtk v5 (ascii)"  << std::endl << std::flush;
	else if (fileFormat==NIBR::VTK_BINARY_5) std::cout << "vtk v5 (binary)" << std::endl << std::flush;
	else {
		std::cout << "unknown"      << std::endl << std::flush;
		return;
	}

	std::cout << "Description:              "  << fileDescription      << std::endl << std::flush;
	std::cout << "Streamline count:         "  << numberOfStreamlines  << std::endl << std::flush;

	if (fileFormat==NIBR::TRK) {

		std::cout << "ijk2xyz: " << std::endl << std::flush;
		
		std::cout << "   " << ijk2xyz[0][0] << " " << ijk2xyz[0][1] << " " << ijk2xyz[0][2] << " " << ijk2xyz[0][3] << std::endl << std::flush;
		std::cout << "   " << ijk2xyz[1][0] << " " << ijk2xyz[1][1] << " " << ijk2xyz[1][2] << " " << ijk2xyz[1][3] << std::endl << std::flush;
		std::cout << "   " << ijk2xyz[2][0] << " " << ijk2xyz[2][1] << " " << ijk2xyz[2][2] << " " << ijk2xyz[2][3] << std::endl << std::flush;
		std::cout << "   " << ijk2xyz[3][0] << " " << ijk2xyz[3][1] << " " << ijk2xyz[3][2] << " " << ijk2xyz[3][3] << std::endl << std::flush;

	}

	
	if ( (fileFormat==NIBR::VTK_BINARY_3) || (fileFormat==NIBR::VTK_BINARY_4) ){
		std::vector<NIBR::TractogramField> fields = findTractogramFields(*this);
		std::cout << "Field count: " << fields.size() << std::endl << std::flush;
		int i = 1;
		for (auto f : fields) {
			std::cout << "   " << i++ << ". " << f.name.c_str() << ": " << ((f.owner == POINT_OWNER) ? "Point" : "Streamline") << " field with dim " << f.dimension << std::endl << std::flush;
		}
	}

	

	std::cout << "\033[0m";

	return;


}
