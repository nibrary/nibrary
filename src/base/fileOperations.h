#pragma once

#include <filesystem>
#include "math/core.h"


namespace NIBR 
{

    bool existsFile(const std::string& name);
    bool existsFolder(const std::string& name);
    
    bool deleteFile(const std::string& fileName);

    std::string getFileExtension(std::string filePath);
    std::string removeFileExtension(std::string filePath);
    std::string getFolderPath(std::string filePath);
    bool        makeFolder(std::string filePath);

    std::string generateRandomFileName();
    std::string generateRandomFileNameInTempPath();
    std::string generateRandomFileNameInThisPath(std::string filePath);

    std::tuple<bool,std::vector<Point>> readTripletsFromTextFile(std::string fname);
    std::tuple<bool,std::vector<std::vector<float>>> readTripletsFromTextFileTo2DVector(std::string fname);

    std::vector<std::string> getMatchingFiles(const std::string& pattern);


}