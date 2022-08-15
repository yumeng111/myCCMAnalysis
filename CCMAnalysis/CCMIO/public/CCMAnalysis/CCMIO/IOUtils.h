#ifndef CCMAnalysis_IOUtils_H
#define CCMAnalysis_IOUtils_H

#include <string>
#include <fstream>

enum class CCMFileType {
    ROOT,
    RawBinary
};

CCMFileType DetermineFileType(std::string fname);
bool IsRootFile(std::string const & fname);
bool FileExists(std::string const & name);

#endif // CCMAnalysis_IOUtils_H
