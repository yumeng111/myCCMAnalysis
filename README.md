Requires C++11 and to read the processed files using ROOT6.
Also need CMAKE3 and xerces-c


# Compile
This code is designed to be compilable on a linux or mac

To compile:
create a build directory where the source code resides
go into build and run


`cmake ../`
or
`cmake3 ../`

if it finishes successfully then run

`make`

# Output and running with ROOT
Add the following lines to your ~/.bash_profile, ~/.bashrc (if you are running tcsh then syntax needs to change and you put the code into ~/.tcshrc)

```bash
export CCMINSTALL="/Users/rtthorn/Documents/CCM/code" # change so it matches where the files are located on your computer
export CCMINSTALLBuild=$CCMINSTALL/build" # change so it matches where the files are located on your computer

export LD_LIBRARY_PATH="$CCMINSTALLBuild:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="$CCMINSTALLBuild:$DYLD_LIBRARY_PATH"
export PATH="$PATH:$CCMINSTALLBuild"
```
Put the following lines in your ~/.rootlogon.C
```c++
{
  //Load libraries needed for CCM analysis
  TString temp(gSystem->GetIncludePath());
  temp += "-I${CCMINSTALL} -I${CCMINSTALLBuild}";
  gSystem->SetIncludePath(temp.Data());
  TString temp2(".include");
  temp2 += gSystem->Getenv("CCMINSTALLBuild");
  gROOT->ProcessLine(temp2.Data());
  TString temp4(".include");
  temp4 += gSystem->Getenv("CCMINSTALLBuild");
  gROOT->ProcessLine(temp4.Data());
  TString temp3(".include ");
  gROOT->ProcessLine(temp3.Data());

  // use the following for linux
  gSystem->Load("${CCMINSTALLBuild}/src/utils/libCCMUtils.so");
  gSystem->Load("${CCMINSTALLBuild}/src/ds/libCCMDS.so");
  gSystem->Load("${CCMINSTALLBuild}/src/spe/libCCMSPE.so");
  gSystem->Load("${CCMINSTALLBuild}/src/modules/framework/libCCMXML.so");
  gSystem->Load("${CCMINSTALLBuild}/src/modules/reco/libCCMReco.so");

}
```

If this does not work for you, please let me know.


# Mapping information
The class PMTMapInfo handles a map of PMTInformation. PMTInformation is a class that contains all the information
about a given PMT that you would want in your analysis.
The process files create a branch in the tree with the PMTInformation that was used to generate the file.
You can use the tree to fill the PMT map with `PMTMapInfo::FillPMTMap(TTree * tree)`, or you can use
`PMTMapInfo::DefaultFillPMTMap()` that will use the mappings/mapping_master*.csv files in your $CCMINSTALL directory
to populate the map. You can always pass your own .csv file with `PMTMapInfo::FillPMTMap(std::istream& infile)`
as long as it has the same format as the default .csv files.



# Module Code
- The `CCMAnalysis` executable reads in the xml config file (examples are in configFiles) and runs the module you want to run
- Current Modules
  - CCMConvertBinary2ROOT (convert binary to ROOT files and creates pulses)
  - CCMNearlineDiag (calculate SPE rates and calculate SPE calibrations)
  - CCMFindEvents (find events and information about them)
  - CCMApplyPreCuts (apply cuts on the events that were found)



Let me know if you have any questions.
Tyler Thornton
P-25 LANL
rtthorn@lanl.gov

