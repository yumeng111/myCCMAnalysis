/* 
CCM simulation primary code

This code calls the header and code files from the /include and /src directories respectively. 
It then initializes the overall simulation and manages the threading and visualization. 

It should not be necessary to manipulate this code for most simulation functions. 
The only except is in modifying some scintillation parameters or physics. 
For example, changing from nuclear to electromagnetic scintillation requires a modification around line 70.
Otherwise this is just a container code that calls the more specific headers and specific codes. 
*/

#include "G4Types.hh"
#include "CCMAnalysis/CCMDetectorSimulation/detectorConstruction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/actionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "TROOT.h"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
//the #include's call the various headers necessary for running the code. 

//main method called when running the simulation.
int main(int argc, char** argv)
{
  ROOT::EnableThreadSafety();

  //initialization of the user interface
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  
  //defines the randomization engine.
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  //initialize the threading, finding whether or not multiple threads are available and using 2 if there are.
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif
  
  //create a new instance of the detectorConstruction. 
  //The commented line is a debugging tool if compilation or running fails with a segfault
  detectorConstruction* det = new detectorConstruction();
  runManager->SetUserInitialization(det);
  //G4cout << "maincodenote: created detector" << G4endl;

  //Initialize the physics. The chosen physics lists are specific for scintillation. 
  //The FTFP_BERT physics list used as the base has most nuclear processes for normal particles 
  //(photons, neutrons, electrons, etc) but has some approximations that one should be cautious of. 
  //As such this simulation is best used for simulating scintillation, not nucleonic interactions. 
  G4VModularPhysicsList* physics = new FTFP_BERT;
  physics->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  //G4cout << "maincodenote: installed main physics" << G4endl;
  
  //The Following lines set more specific parameters for scintillation. 
  //Specifically, lines 71 and 72 hold the parameters for nucleonic (no triple) and electromagnetic scintillation.
  opticalPhysics->SetWLSTimeProfile("delta");
  opticalPhysics->SetScintillationYieldFactor(1.0);//for e/m
  //opticalPhysics->SetScintillationYieldFactor(.25);//for nucleon
  opticalPhysics->SetScintillationExcitationRatio(0.0);
  
  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);
  
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);
  //G4cout << "maincodenote: installed optical physics" << G4endl;

  physics->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physics);

  //Initialize the physics action using the above physics lists. 
  runManager->SetUserInitialization(new actionInitialization());
  //G4cout << "maincodenote: initialized runmanager" << G4endl;
  
  //Define the visualization managers. It should not be necessary to modify anything below this point. 
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  //G4cout << "maincodenote: started visualization" << G4endl;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  //Defines the first argument after the executable as a macro file to run. 
  if ( !  ui ) {
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UImanager->ApplyCommand(command+filename);
  }
  //calls the visualization macro if none is provided, starts a user interface session. 
  else {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  //G4cout << "maincodenote: done" << G4endl;
  
  //deletes the managers to prevent memory leaks. 
  delete visManager;
  delete runManager;
  return 0;
}
