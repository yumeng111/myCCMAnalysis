This is the README for the CCM simulation built by Edward. Please email any questions to e.dunton@columbia.edu.



#To make the simulation:

Ensure geant4 is loaded correctly. <b>On the HPC run 'module list' and make sure the version of gcc loaded is 8.4.5, plus cmake 3.14.0.</b> Then run (or have in your .bashrc) the line:

```bash
#Geant4
source /usr/projects/w20_ccm_lanl/Software/geant4.10.06_install/bin/geant4.sh
```

Next, create a build directory. I typically make mine in the simulation directory, so that the path to the build is:

`$CCMINSTALL/simulationCCM/build/`

From there you need to run an extensive cmake command, pointed at the simulationCCM directory. Mine would be:

`cmake -DGeant4_DIR=/usr/projects/w20_ccm_lanl/Software/geant4.10.06_install/lib64/Geant4-10.6.0 -DXercesC_LIBRARY=/usr/projects/w20_ccm_lanl/Software/xerces-c-3.2.2_install/lib/libxerces-c-3.2.so /usr/projects/w20_ccm_lanl/Software/geant4.10.06/ -DXercesC_INCLUDE_DIR=/usr/projects/w20_ccm_lanl/Software/xerces-c-3.2.2_install/include -DCMAKE_PREFIX_PATH=/usr/projects/w20_ccm_lanl/Software/geant4.10.06_install/ ../`

For convenience I added the following to my .bashrc:

```bash
alias g4cmake="cmake -DGeant4_DIR=/usr/projects/w20_ccm_lanl/Software/geant4.10.06_install/lib64/Geant4-10.6.0 -DXercesC_LIBRARY=/usr/projects/w20_ccm_lanl/Software/xerces-c-3.2.2_install/lib/libxerces-c-3.2.so /usr/projects/w20_ccm_lanl/Software/geant4.10.06/ -DXercesC_INCLUDE_DIR=/usr/projects/w20_ccm_lanl/Software/xerces-c-3.2.2_install/include -DCMAKE_PREFIX_PATH=/usr/projects/w20_ccm_lanl/Software/geant4.10.06_install/"
```

So my new cmake command is just:

`g4cmake ../`

Which is much more convenient. After cmake has run, assuming it ran without any issues, stay within the /build/ directory and run
"make" (no options should be needed). This will compile the various .cc files in the /src/ directory, one by one, so if it fails you see where the issue came from. Once that has completed it should produce the executable "simulationCCM". 

Once the cmake command has been run once, it should not be necessary to run cmake again even if you modify the simulation source code files. Just make will usually more than suffice for alterations to the simulation.
 
To run the simulation, just choose a macro file and run:

`./simulationCCM runlaser.mac`

You could also add the following to your .bashrc
```bash
export LD_LIBRARY_PATH=$CCMINSTALL/simulationCCM/build:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$CCMINSTALL/simulationCCM/build:$DYLD_LIBRARY_PATH
export PATH=$CCMINSTALL/simulationCCM/build:$PATH
```
so you can run the the executable from any directory

Where runlaser is replaced by whatever macro file you want. Macro files use the available commands to run a variety of simulations without modifying the source codes or rebuilding. Currently the simulation is capable of running the Laser, Sodium, or other calibration source of arbitrary energy (assuming it produces only a single particle) without remaking. 

OVERLAP WARNINGS:

One note about running: with any version of the simulation there will be a large number of 'Overlap Warnings' that show up in the output stream. If you inspect them more closely, they will all be relatively small (a few mm at most) and there will be a note this is just a warning. The reason this comes up is that in placing the PMTs along the outside of the detector, we are placing a flat surface (the bottom of the PMT semisphere) against a curved one (the edge of the Fiducial volume). Thus there is inevitably either a bit of overlap, or a small hole between the PMT and the edge of the Fiducial cylinder. Neither situation causes a major problem, but in making the simulation I chose to prefer a procedural issue (overlap warning) to a possible physics issue (gaps behind the PMTs). This does not affect the running of the simulation at all, and can be safely ignored. The only issue it causes is cluttering the verbose pre-run outputs, which is not a significant issue for mosts cases.



#Running the Simulation: macro files:

When running geant4 simulations, a large number of ui commands are normally required. As such it is advisible to prepare a macro file that contains these commands and thus allows you to run longer simulations rather than trying to run one by one through user input. The current simulation has 4 major macro files included, each of which demonstrates some feature of the simulation.

     runlaser.mac is the macro file to run the laser at both wavelengths (532 and 213) and all three positions, producing a number of events for each of these simulations and placing them all in different output files. 

     runcalib.mac is a macro file that will run the Sodium and Cobalt calibrations in sequence, placing the outputs in two seperate macro files.

     run1.mac is a simple macro file that keeps everything default, so is least likely to cause segfaults if that comes up as a problem. It can be used as a test to see if the issue is in the running of the code or the actual construction; if the latter, changing the defaults in the deteector construction sometimes works. If the former, there is a more complex issue.

     init_vis.mac and vis.mac are essentially a pair of linked codes, where one runs a few setups before invoking the other. These are the files that control the visualization (requires X11 forwarding; make sure that is turned on or you're running on a local machine) and will be automatically called if you run just the executable.

You can add as many macro files as you want, to do any number of things. These four just provide an example and describe a few of the commands.
      


#Analyzing the Simulation output:

The simulation currently outputs text files that record the start of an event (when the calibration source is generated/laser is created/particle gun is fired) and the photons that enter a pmt. It also records non-optical photons which deposit energy in the Fiducial volume, which can be used to calculate the total energy seen by the Fiducial volume during an event (useful for calibration sources). 

The start of an event is recorded simply as 'Start Event' so searching for that in the output file will let you know where each event is started. As geant simulates events one by one, this works as an effective means of determining where one event ends and the next begins. The photons that impact PMTs are a bit more complicated. Each photon gets its own line which records the PMT hit, the energy of the photon, the time it hits, and the cosine of the incident angle. These values can be used to calculate a probability that the photon gets turned into a photo-electron, by using the energy dependance of the PMT response function and the angular dependance, which is essentially sqrt(cos(theta)) on top of the energy dependance. 

Currently included in the build directory are a set of python files, calibfunctions.py and fastAnalysis.py. the first contains a large number of defined functions for reading an output file, turning its output into waveforms, analyzing those waveforms to find events (analysis events, seperate from geant4 events), then plotting sets of events to get spectrums. It also contains functions for use with the laser, including calculating row ratios and coated/uncoated ratios. fastAnalysis.py, on the other hand, is an example code that uses the functions defined in order to analyze a set of laser result files and return a bunch of useful information about the result. 



#To modify the simulation:

Modifying the simulation requires going into the source code and making modifications. This is necessary if you want to change the physics of the liquid Argon, reduce the TPB efficiency, remove specific PMTs, alter the reflectivity of the foils, or a number of other things. More specifics on what each part of the code can do are included in the comments of each individual code. 

After making an alteration, one should go to the build directory and run 'make' to remake the executable with the changes. If that works, in theory running with a macro will work as well. One change that can be used to test this is exchanging the current 'detectorConstruction.cc' with one of the '*Detector.file' files which contain examples of different detector designs that we have been considering. Do this by running 'cp idealDetector.file detectorConstruction.cc' (exchange 'idealDetector' for the detector of your choice) so as to not lose or overwrite the .file files. 

The example files are as follows:

    idealDetector.file, which provides the geometry and physics for an ideal, best case scenario detector with no lAr or TPB issues. As such it can be used for optimal expectations, showing what the detector could do in the best possible case.

    paintedFoilsDetector.file, which is exactly the same as idealDetector but with two kinds of TPB to represent the different efficiencies between painted and evaporative coated tpb. This file can be used to see what the best case scenario was for our current detector. It is also an example of how to manipulate the TPB properties to produce lower efficiency in two different ways.

    LBOCdetector.file replaces the tpb on the foils with a Light guide Bars Optical Cylinder (LBOC). This file also contains a number of possible liquid Argon absorption sets to represent different levels of contamination. They are commented, with descriptors explaining how much contamination each set represents. The LBOCdetector also has options for single vs triple cylinder configuration through the same flag as was used to turn the tpbfoils on and off with the foil-based detectors.

    bestFitDetector.file is the current best fit detector to the various calibration values obtained from data (laser row ratios, sodium and cobalt energy scale). It is significantly more complex than the other files, including such features as unsmooth geometries, layers of liquid Argon, three types of liquid Argon, multiple sections of tpb, and a large number of options in comments for controlling all of those things. Unfortunately it is not a very good fit yet, so it is not advised to used this detector for anything other than obtaining approximations of the energy scale of various energy calibration sources. It is, however, still the best fit to the current data and thus included for that reason. 



