#pragma once

#include <filesystem>
#include "base/nibr.h"
#include "math/core.h"
#include "base/stringOperations.h"


namespace NIBR 
{

    bool existsFile(const std::string& name);
    bool existsFolder(const std::string& name);
    
    bool deleteFile(const std::string& fileName);
    bool copyFile(const std::string& sourceFileName, const std::string& destinationFileName); // Overwrites existing files

    std::string getFileExtension(std::string filePath);
    std::string removeFileExtension(std::string filePath);
    std::string getFolderPath(std::string filePath);
    bool        makeFolder(std::string filePath);

    std::string generateRandomFileName();
    std::string generateRandomFileNameInTempPath();
    std::string generateRandomFileNameInThisPath(std::string filePath);

    std::tuple<bool,std::vector<Point3D>> readTripletsFromTextFile(std::string fname);

    std::vector<std::string> getMatchingFiles(const std::string& pattern);

    bool areFilesIdentical(const std::string& fileName1, const std::string& fileName2);
    
    void prepSequentialRead(const std::string& filename);
    void prepRandomRead(const std::string& filename);

    template<class T>
    std::tuple<bool, std::vector<T>> readValuesFromFile(std::string fname) {

        CNumericLocaleGuard lockScopeForNumericReading;

        auto f = fopen(fname.c_str(), "rb");
        if (f == NULL) {
            return std::make_tuple(false, std::vector<T>());
        }

        std::vector<T> out;
        T value;
        size_t read_size = sizeof(T);
        while (fread(&value, read_size, 1, f) == 1) {
            out.push_back(value);
        }

        fclose(f);
        return std::make_tuple(true, out);
    }

    // Explicit template instantiations for int, float, and double
    template std::tuple<bool, std::vector<int>>    readValuesFromFile<int>(std::string fname);
    template std::tuple<bool, std::vector<float>>  readValuesFromFile<float>(std::string fname);
    template std::tuple<bool, std::vector<double>> readValuesFromFile<double>(std::string fname);


}