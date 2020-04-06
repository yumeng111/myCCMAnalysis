/*
  Primary Generator for CCM simulation

This code create the starting particles for each event. It currently has three major options:
Laser: setting the laser to true produces a large number of optical photons for each event, 
       skipping the scintillation step.
Sodium: setting Sodium to true creates events with the proper Sodium-22 decay chain,
        specifically a 1.2 MeV gamma, plus two back to back 511 keV gammas 90% of the time.
Neither: setting both to false allows the creation of other kinds of initial events, 
         such as 120 keV gammas for cobalt, or neutrons/Ar-39 decay events. 
 */
#include "primaryGenerator.hh"
#include "detectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//constructor: sets the base parameters. 
//Anything set here can be changed by a macro file, and will be overwritten by one.
primaryGenerator::primaryGenerator()
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(0)
{
  //create a particle gun that creates 1 particle at a time (so each particle can be controlled individually).
  G4int nparticles = 1;
  fParticleGun = new G4ParticleGun(nparticles);

  //set default particle kinematics
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(1*eV);
}

//Deconstructor
primaryGenerator::~primaryGenerator()
{
  delete fParticleGun;
}

//Generate primaries function. 
//Anything defined here will NOT be overwritten by a macro files.
void primaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
  //This function is called at the beginning of an event

  //Booleans to be flagged if running either a Sodium or Laser simulation.
  G4bool Sodium = false;
  G4bool Laser = true;

  //generates several randoms for position, momentum, and probabilities. 
  G4double phi = G4RandFlat::shoot(-0.28, 6.002);
  G4double theta = G4RandFlat::shoot(-0.0005, 3.142);
  G4double check = G4RandFlat::shoot(0.1,1.1);
  G4double radi = G4RandFlat::shoot(0.01, 0.1)*cm;

  //Define the position based on the randoms generated above. current settings: within the laser plates.
  G4double ypos = radi*sin(phi);
  G4double xpos = radi*cos(phi);
  G4double zpos = 0.1*cm; 

  //Generate the initial momentum based on the random variables above. 
  G4double momx = sin(theta)*cos(phi);
  G4double momy = sin(theta)*sin(phi);
  G4double momz = cos(theta);

  //defines a variable for rodHeight, necessary for the laser.
  G4double rodHeight = 0;

  //Obtain the rod Height from the detector to produce the photons at the right location.
  const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  rodHeight = detector->GetRodHeight();
  
  //determine if the Laser or Sodium is on from the detector. Note: Laser overwrites Sodium, so only one runs if both are on.
  Laser = detector->GetfLaser();
  if (!Laser) {
    Sodium = detector->GetfSodium();
  }
  
  //begin the primary generation for laser runs.
  if (Laser) {
    zpos = rodHeight -.95*cm;
    
    //Make the initial particles opticalphotons.
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="opticalphoton");
    fParticleGun->SetParticleDefinition(particle);
    
    //get the particle energy entered either through the above lines or from a macro.
    //then set the number of photons to be generated to 1000 for 532 and 10000 for 213
    //note: more needed for 213 due to non-ideal detector killing most of them.
    G4double partEneg = fParticleGun->GetParticleEnergy();
    G4int totalops = 1000;
    if (partEneg > 5*eV) {
      totalops = 10000;
    }

    //run through the photons, generating them one at a time.
    for (G4int i=0;i<totalops;++i){
      //randomize the angle for each photon based on the approx. gaussian after the diffuser plates.
      phi = G4RandFlat::shoot(-0.28, 6.002);
      theta = G4RandGauss::shoot(0,12.8);
      if (check < 0.2) {
	//add a flat component to the generated angle to better match the diffuser distribution.
	theta = G4RandFlat::shoot(0.01,90);
      }
      //redfine the angle in radians and use it to generate initial photon angles.
      theta = CLHEP::pi/180.*theta;
      radi = G4RandFlat::shoot(0.001,1.33)*cm;
      momx = sin(theta)*cos(phi);
      momy = sin(theta)*sin(phi);
      momz = -cos(theta);

      //define the x and y position based on random angles earlier.
      xpos = radi*cos(phi);
      ypos = radi*sin(phi);
      
      //generate the photon through the particle gun.
      fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }
  }

  //if Sodium, overwrite any input from the macro file and create the 1.2 MeV gamma at the position (0,0,0).
  if (Sodium) { 
    fParticleGun->SetParticleEnergy(1.2*MeV);
    xpos = 0*cm;
    ypos = 0*cm;
    Laser = false;
  }
    
  //if not Laser, generate a primary for a single particle (Sodium 1.2 MeV or other from macro input).
  if (!Laser) {
    fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  //if Sodium, determine if the 2 511 keVs are going to be generated and do so.
  if (Sodium && check > 0.2) {
    //create a new momentum angle.
    phi = G4RandFlat::shoot(-0.28,6);
    theta = G4RandFlat::shoot(-0.14,3);
    momx = sin(theta)*cos(phi);
    momy = sin(theta)*sin(phi);
    momz = cos(theta);

    G4double momx2 = -momx;
    G4double momy2 = -momy;
    G4double momz2 = -momz;

    //generate the 511's back to back.
    fParticleGun->SetParticleEnergy(511*keV);

    fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  
    fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx2, momy2, momz2));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }//*/
  
  //output at the start of an event.
  G4cout << "primarygeneratornote: Start Event" << G4endl;

  //more complete output for starting the event, showing momentum and position initiated.
  //  G4cout << "S0" << "  particleinit  " << momx << "  " << momy << "  " << momz << "  " <<
  // xpos/cm << "  " << ypos/cm << "  " << zpos/cm <<G4endl;

}
