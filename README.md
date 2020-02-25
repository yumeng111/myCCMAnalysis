Requires C++11 and to read the processed files using ROOT6.
Also need CMAKE 3


# Compile
This code is designed to be compilable on a linux or mac

To compile, go into build and run


`cmake ../`

if it finishes successfully then run

`make`

If for some reason cmake is not successful, you might need to change the following lines in CMakeLists.txt
```cmake
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/root/cmake)
```


# Output and running with ROOT
Add the following lines to your ~/.bash_profile, ~/.bashrc (if you are running tcsh then syntax needs to change and you put the code into ~/.tcshrc)

```bash
export CCMINSTALL="/Users/rtthorn/Documents/CCM/code" # change so it matches where the files are located on your computer

export LD_LIBRARY_PATH="$CCMINSTALL/build:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="$CCMINSTALL/build:$DYLD_LIBRARY_PATH"
export PATH="$PATH:$CCMINSTALL/build"
```
Put the following lines in your ~/.rootlogon.C
```c++
{
  //Load libraries needed for CCM analysis
  TString temp(gSystem->GetIncludePath());
  temp += "-I${CCMINSTALL} -I${CCMINSTALL}/build";
  gSystem->SetIncludePath(temp.Data());
  TString temp2(".include");
  temp2 += gSystem->Getenv("CCMINSTALL");
  temp2 += "/build";
  gROOT->ProcessLine(temp2.Data());
  TString temp4(".include");
  temp4 += gSystem->Getenv("CCMINSTALL");
  gROOT->ProcessLine(temp4.Data());
  TString temp3(".include ");
  gROOT->ProcessLine(temp3.Data());

  // use the following for a mac
  gSystem->Load("libEvent.dylib");

  // use the following for linux
  //gSystem->Load("libEvent.so");

}
```

If this does not work for you, please let me know.


# Mapping information
The class PMTMapInfo handles a map of PMTInformation. PMTInformation is a class that contains all the information
about a given PMT that you would want in your analysis.
The process files create a branch in the tree with the PMTInformation that was used to generate the file.
You can use the tree to fill the PMT map with `PMTMapInfo::FillPMTMap(TTree * tree)`, or you can use
`PMTMapInfo::DefaultFillPMTMap()` that will use the mapping_master*.csv files in your $CCMINSTALL directory
to populate the map. You can always pass your own .csv file with `PMTMapInfo::FillPMTMap(std::istream& infile)`
as long as it has the same format as the default .csv files.



# Example Code
- The `EnergyCalibration` executable (`energyCalibration.cc`) shows how to read in the raw binary data and populate the Pulses class.
- The `FindEvents` executable (`findEvents.cc`) shows how to build events from the pulses
- The `ApplyPreCuts` executable (`applyPreCuts.cc`) shows how make cuts on the events that were found



Let me know if you have any questions.
Tyler Thornton
P-25 LANL
rtthorn@lanl.gov
