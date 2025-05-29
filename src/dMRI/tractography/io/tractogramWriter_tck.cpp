#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <cmath>

#include "base/verbose.h"
#include "base/dataTypeHandler.h"
#include "tractogramWriter_tck.h"

using namespace NIBR;

TCKWriter::TCKWriter(std::string _filename) : filename(_filename) {}

TCKWriter::~TCKWriter() {
    if (file != nullptr) {
        fclose(file);
    }
}

bool TCKWriter::open() {
    file = fopen(filename.c_str(), "wb+"); // Use "wb+" for binary write + read/seek
    if (file == NULL) {
        disp(MSG_ERROR, "Cannot open output file: %s", filename.c_str());
        return false;
    }

    int written;

    written = fprintf(file, "mrtrix tracks\n");
    written = fprintf(file, "datatype: Float32LE\n");

    posCount = ftell(file);
    written = fprintf(file, "count: %20lu\n", (unsigned long)0); // Padded placeholder

    posFileOffset = ftell(file);
    written = fprintf(file, "file: . %20ld\n", (long)0);        // Padded placeholder

    written = fprintf(file, "END\n");

    // Ensure END is followed by a newline, some readers might need it.
    // Check if the last char was newline, if not, add it.
    // However, fprintf adds it, so ftell should be correct.
    dataStartOffset = ftell(file);

    if (written < 0) {
        disp(MSG_ERROR, "Failed to write TCK header placeholders.");
        fclose(file);
        file = nullptr;
        return false;
    }

    // Rewrite the file offset placeholder with the correct data offset now.
    // We will rewrite count later.
    fseek(file, posFileOffset, SEEK_SET);
    fprintf(file, "file: . %20ld\n", dataStartOffset);
    fseek(file, dataStartOffset, SEEK_SET); // Go back to where data should start

    return true;
}

bool TCKWriter::writeBatch(const std::vector<std::vector<std::vector<float>>>& batch) {
    if (file == NULL) return false;

    const float NAN_val = NAN; // Get NAN value
    const float NAN_arr[3] = {NAN_val, NAN_val, NAN_val};

    for (const auto& streamline : batch) {
        int len = streamline.size();
        if (len > 0) {
            // Write points for the streamline
            for (const auto& point : streamline) {
                fwrite(point.data(), sizeof(float), 3, file);
            }
            // Write NAN separator
            fwrite(NAN_arr, sizeof(float), 3, file);
            totalStreamlineCount++;
        }
    }

    // Check for write errors (optional, checking ferror(file) periodically)
    if (ferror(file)) {
        disp(MSG_ERROR, "Error writing TCK batch data.");
        return false;
    }

    return true;
}

bool TCKWriter::close() {
    if (file == NULL) return false;

    const float INF_val = INFINITY;
    const float INF_arr[3] = {INF_val, INF_val, INF_val};

    // Write the final INF marker (following your original writer's convention)
    // Note: Many TCK readers just expect EOF or the last NAN.
    // If INF causes issues, remove this line and ensure the last write was NAN.
    fwrite(INF_arr, sizeof(float), 3, file);

    // Seek back and write the final count
    fseek(file, posCount, SEEK_SET);
    fprintf(file, "count: %20lu\n", totalStreamlineCount);

    // Close the file
    fclose(file);
    file = nullptr;
    return true;
}
