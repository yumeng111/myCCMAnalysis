#include <string>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>

#include "CCMAnalysis/io/IOUtils.h"

CCMFileType DetermineFileType(std::string fname) {
    if(IsRootFile(fname)) {
        return CCMFileType::ROOT;
    } else {
        return CCMFileType::RawBinary;
    }
}

bool IsRootFile(std::string const & fname) {
    // Check if the first 4 bytes are "root"
    std::ifstream input(fname.c_str(), std::ios::binary);
    constexpr unsigned int n_bytes = 4;
    char bytes[n_bytes+1];
    bytes[n_bytes] = 0;
    input.read(bytes, n_bytes);

    return std::string(bytes) == std::string("root");
}

bool FileExists(std::string const & name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
