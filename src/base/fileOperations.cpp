#include <fstream>
#include <random>
#include <regex>
#include <fcntl.h>
#include <thread>
#include <system_error>
#include "config.h"
#include "fileOperations.h"
#include "verbose.h"
#include "multithreader.h"

using namespace NIBR;

// Returns true if file exists
bool NIBR::existsFile(const std::string& name) 
{
    std::filesystem::path p(name);

    // Checks if the path exists and is a regular file (not a directory)
    return (std::filesystem::exists(p) && std::filesystem::is_regular_file(p));  
}

// Returns true if file exists
bool NIBR::existsFolder(const std::string& name) 
{
    std::filesystem::path p(name);

    // Checks if the path exists and is a directory
    return (std::filesystem::exists(p) && std::filesystem::is_directory(p));  
}

bool NIBR::deleteFile(const std::string& fileName) 
{

    std::filesystem::path filePath(fileName);

    // Checks if the path exists and is a regular file (not a directory)
    if (std::filesystem::exists(filePath) && std::filesystem::is_regular_file(filePath)) {
        try {
            std::filesystem::remove(filePath);
            return true;
        } catch (const std::filesystem::filesystem_error& e) {
            disp(MSG_ERROR, ("Error deleting the file: " + std::string(e.what())).c_str() );
            return false;
        }
    } else {
        disp(MSG_WARN, ("File does not exist: " + fileName).c_str() );
    }

    return true;

}

bool NIBR::copyFile(const std::string& sourceFileName, const std::string& destinationFileName) 
{
    std::filesystem::path sourceFilePath(sourceFileName);
    std::filesystem::path destinationFilePath(destinationFileName);

    disp(MSG_DEBUG, ("Source file path: " + sourceFilePath.string()).c_str());
    disp(MSG_DEBUG, ("Destination file path: " + destinationFilePath.string()).c_str());

    // Checks if the source path exists and is a regular file (not a directory)
    if (std::filesystem::exists(sourceFilePath) && std::filesystem::is_regular_file(sourceFilePath)) {
        try {
            // Create directories in the destination path if they do not exist
            disp(MSG_DEBUG, ("Creating directories for destination path: " + destinationFilePath.parent_path().string()).c_str());
            std::filesystem::create_directories(destinationFilePath.parent_path());
            
            disp(MSG_DEBUG, "Copying file...");
            std::filesystem::copy(sourceFilePath, destinationFilePath, std::filesystem::copy_options::overwrite_existing);
            disp(MSG_DEBUG, "File copied successfully.");
            return true;
        } catch (const std::filesystem::filesystem_error& e) {
            disp(MSG_ERROR, ("Error copying the file: " + std::string(e.what())).c_str() );
            return false;
        }
    } else {
        disp(MSG_WARN, ("Source file does not exist: " + sourceFileName).c_str() );
    }

    return false;
}


std::string NIBR::getFileExtension(std::string filePath) 
{

    std::filesystem::path p(filePath);
    std::string str = p.filename().string();

    // Find the last occurrence of '.' in the given string
    std::size_t pos = str.find_last_of(".");

    // If a period is found 
    if (pos != std::string::npos) 
    {
        // Extract extension after '.'
        std::string extension = str.substr(pos+1);

        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
        
        // Check if extension is "gz"
        if (extension == "gz")
        {
            // Find the second to last '.'
            std::size_t second_to_last_pos = str.substr(0, pos).find_last_of(".");
            if (second_to_last_pos != std::string::npos)
            {
                // Return the extension from before '.gz' to end of string 
                return str.substr(second_to_last_pos + 1);
            }
        }
        return extension;
    }

    // If no period is found
    return "";
    
}

std::string NIBR::removeFileExtension(std::string filePath)
{

    // Convert string to path
    std::filesystem::path p(filePath);
    std::string str = p.filename().string();

    // Find the last occurrence of '.' in the given string
    std::size_t pos = str.find_last_of(".");

    // If a period is found 
    if (pos != std::string::npos)
    {
        // Extract the extension after '.'
        std::string extension = str.substr(pos+1);
        
        // If the extension is 'gz', we need to check if there's another . preceding it
        if (extension == "gz")
        {
            // Find the second to last '.'
            std::size_t second_to_last_pos = str.substr(0, pos).find_last_of(".");
            if (second_to_last_pos != std::string::npos)
            {
                // The extension is till second '.' for '.gz' files;
                // so return the filename till the second last '.'
                return str.substr(0, second_to_last_pos);
            }
        }
        // Regular case - return the part of the string before the last '.'
        return str.substr(0, pos);
    }

    // If no period is found, return the entire filename
    return str;

}


std::string NIBR::getFolderPath(std::string filePath) {
    
    std::filesystem::path p(filePath);

    // Check if the input path is a directory or a file
    if (std::filesystem::is_directory(p)) {
        return p.string();
    } else {
        p = p.parent_path();
        return p.string();
    }
}

bool NIBR::makeFolder(std::string dirPath) {
    
    std::filesystem::path p(dirPath);
    
    if (std::filesystem::exists(p)) {
        if (!std::filesystem::is_directory(p)) {
            disp(MSG_FATAL,"Path points to a file: %s", dirPath.c_str());
            return false;
        }
    } else {
        if (!std::filesystem::create_directories(p)) {
            disp(MSG_FATAL,"Failed to create directory: %s", dirPath.c_str());
            return false;
        }
    }

    return true;
}

std::string NIBR::generateRandomFileName() {

    // Generate a file name of length 16
    std::random_device rd;
    std::default_random_engine rng(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    std::stringstream ss;
    ss << std::hex << std::setfill('0');
    for (int i = 0; i < 8; ++i) {
        ss << std::setw(2) << dist(rng);
    }

    return ss.str();

}

std::string NIBR::generateRandomFileNameInThisPath(std::string filePath) {
    std::filesystem::path p(filePath);

    // Check if the input path is a directory or a file
    if (std::filesystem::is_directory(p)) {
        p /= generateRandomFileName();
    } else {
        p = p.parent_path();
        p /= generateRandomFileName();
    }

    return p.string();
}

std::string NIBR::generateRandomFileNameInTempPath() {    
    std::filesystem::path tempPath = std::filesystem::temp_directory_path();
    tempPath /= generateRandomFileName();
    return tempPath.string();
}


// Returns false if values can't be read
std::tuple<bool,std::vector<Point3D>> NIBR::readTripletsFromTextFile(std::string fname)
{

    std::vector<Point3D> out;

    std::string   values;
    std::ifstream inpFile(fname);

    if (!inpFile.good()) {
        disp(MSG_ERROR,"Can't read values from %s", fname.c_str());
        return std::make_tuple(false,out);
    }

    int lineNo = 1;
    bool readError = false;
    while(std::getline(inpFile,values)) {
        if (!values.empty()) {
            std::stringstream xyz(values);
            float x,y,z;
            if (xyz.good()) xyz >> x; else { readError = true; break; }
            if (xyz.good()) xyz >> y; else { readError = true; break; }
            if (xyz.good()) xyz >> z; else { readError = true; break; }
            out.push_back({x,y,z});
        }
        lineNo++;
    }
    inpFile.close();
    if (readError) {
        disp(MSG_ERROR,"Can't read values from %s line %d", fname.c_str(),lineNo);
        return std::make_tuple(false,out);
    }

    return std::make_tuple(true,out);

}

std::string wildcardToRegex(const std::string& wildcard) {
    std::string regexPattern;
    for(char c : wildcard) {
        switch(c) {
            case '*': regexPattern += ".*"; break;
            case '?': regexPattern += "."; break;
            case '.': regexPattern += "\\."; break;
            default: regexPattern += c; break;
        }
    }
    return regexPattern;
}

std::vector<std::string> NIBR::getMatchingFiles(const std::string& pattern) {

    std::vector<std::string> result;

    std::filesystem::path p(pattern);
    std::filesystem::path directory;

    if (p.has_parent_path()) {
        directory = p.parent_path();
    } else {
        directory = std::filesystem::current_path();
    }
    disp(MSG_DEBUG,"directory: %s", directory.c_str());
    
    const auto filenameWildcard = p.filename().string();
    const auto filenameRegexPattern = wildcardToRegex(filenameWildcard);
    disp(MSG_DEBUG,"filenameRegexPattern: %s", filenameRegexPattern.c_str());
    std::regex filenameRegex(filenameRegexPattern, std::regex::ECMAScript | std::regex::icase);

    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory)) {
        disp(MSG_ERROR,"Invalid directory: %s", directory.c_str());
        return result;
    }

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (std::filesystem::is_regular_file(entry) && std::regex_match(entry.path().filename().string(), filenameRegex)) {
            result.push_back(entry.path().string());
        }
    }

    return result;
}

bool NIBR::areFilesIdentical(const std::string& fileName1, const std::string& fileName2) 
{
    
    std::ifstream file1(fileName1, std::ifstream::binary | std::ifstream::ate);
    std::ifstream file2(fileName2, std::ifstream::binary | std::ifstream::ate);

    // Check if either file failed to open
    if (!file1.is_open() || !file2.is_open()) {
        if (!file1.is_open()) {
            std::cerr << "Error opening file: " << fileName1 << std::endl;
        }
        if (!file2.is_open()) {
            std::cerr << "Error opening file: " << fileName2 << std::endl;
        }
        return false;
    }

    // Get the size of each file
    std::ifstream::pos_type fileSize1 = file1.tellg();
    std::ifstream::pos_type fileSize2 = file2.tellg();

    // If the file sizes are different, they are not identical
    if (fileSize1 != fileSize2) {
        return false;
    }

    // Move the file pointers back to the beginning of the files
    file1.seekg(0, std::ifstream::beg);
    file2.seekg(0, std::ifstream::beg);

    // Define a buffer size for reading chunks of the files
    const size_t bufferSize = 1024 * 1024 * 10; // 10 MB
    std::vector<char> buffer1(bufferSize);
    std::vector<char> buffer2(bufferSize);

    // Read and compare the files in chunks
    while (file1.read(buffer1.data(), bufferSize) && file2.read(buffer2.data(), bufferSize)) {
        if (std::memcmp(buffer1.data(), buffer2.data(), bufferSize) != 0) {
            return false;
        }
    }

    // After the loop, check the number of bytes read in the last read operation
    if (file1.gcount() != file2.gcount() || std::memcmp(buffer1.data(), buffer2.data(), file1.gcount()) != 0) {
        return false;
    }

    // If all chunks are identical, the files are identical
    return true;
}

void NIBR::prepSequentialRead(const std::string& filename) 
{

#ifndef _WIN32

    disp(MSG_DETAIL, "Hinting kernel to preload %s using posix_fadvise...", filename.c_str());

    int fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
        disp(MSG_FATAL, "Could not open file %s", filename.c_str());
        return;
    }

    // Advise the kernel that file access will be sequential
    #if defined(__APPLE__)
        fcntl(fd, F_RDAHEAD, 1);
    #else
        posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
    #endif

    close(fd);

#endif

    return;
}

void NIBR::prepRandomRead(const std::string& filename)
{

#ifndef _WIN32

    disp(MSG_DETAIL, "Hinting kernel for random access on %s using posix_fadvise...", filename.c_str());

    int fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
        disp(MSG_FATAL, "Could not open file %s", filename.c_str());
        return;
    }

    // Advise the kernel that file access will be random
    #if defined(__APPLE__)
        fcntl(fd, F_NOCACHE, 1);
    #else
        posix_fadvise(fd, 0, 0, POSIX_FADV_RANDOM);
    #endif

    close(fd);
    
#endif

    return;
}
