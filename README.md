Requires C++11 and to read the processed files using ROOT6.
Also need CMAKE3 and xerces-c

This README is for the main analysis code. For the README on the simulation see the README.md
in the simulationCCM directory.


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
Add the following lines to your ~/.bash_profile, ~/.bashrc (if you are running tcsh then syntax needs to change and you put the code into ~/.tcshrc). Also you must compile for gcc version 4.8.5 when running on LANL HPC computers, thanks to the version used when compiling ROOT.

```bash
#must have the XERCES environmental variables when running on LANL HPC Computers
export XERCESINSTALL=$PROJECT/Software/xerces-c-3.2.2_install
export XERCESINCLUDE=$XERCESINSTALL/include
export XERCESLIB=$XERCESINSTALL/lib

export CCMINSTALL=$HOME/CCMCode/CCM_analysis_sw_dev_xmlConfig # change for your directory
export CCMINSTALLBuild=$CCMINSTALL/build # change for your directory
export LD_LIBRARY_PATH=$XERCESLIB:$CCMINSTALLBuild:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$XERCESINSTALL/lib:$CCMINSTALLBuild:$DYLD_LIBRARY_PATH
export PATH=$XERCESINSTALL/bin:$CCMINSTALLBuild:$PATH

# the following is for HPC SLRUM job submission only.
# I had created directories to copy rawFiles to and processFile to be stored
# that made it easier for me in creating job submission scripts, but  is not needed
export LUSTRETMP=/lustre/scratch3/turquoise/rtthorn
export CCMRAW=$LUSTRETMP/rawFiles
export CCMPROCESSED=$LUSTRETMP/processedFiles

# need the following to compile the code
module load cmake/3.14.0
```

If you are running on the LANL HPC and you want to use the default mappings and calibration files then you should also add
```bash
export CCMPROJECT=/usr/projects/w20_ccm_lanl
```
to your $HOME/.bash_profile or $HOME/.bashrc file otherwise you will not be able to use the default options when speficifing
the location of the mappings, and calibration files.

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
  gSystem->Load("${CCMINSTALLBuild}/src/io/libCCMIO.so");
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
  - CCMFindPulses(convert binary to ROOT files and creates pulses)
  - CCMSPECalc (calculate SPE calibrations)
  - CCMSRateCalc (calculate SPE rates)
  - CCMFindEvents (find events and information about them)
  - CCMNa22Cuts (apply cuts on the events that were found)
- <b>If you want to develope your own module or suggest changes to an existing, please let me know.</b>

# Analysis Directory
Each user should create a personal directory under $CCMINSTALL/analysis to store their root and python scripts.
Please do not directly modify someone elses scripts, copy to your directory and change from there.



Let me know if you have any questions.<br>
Tyler Thornton<br>
P-25 LANL<br>
rtthorn@lanl.gov<br>

