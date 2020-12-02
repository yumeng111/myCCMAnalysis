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
#include "G4IonTable.hh"
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
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  momx = 0.0;
  momy = 0.0;
  momz = 1.0;
  partEneg = 1.0*keV;
  nameString = "undefined name";

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
  G4bool Ar39 = false;
  G4bool darkMatter = false;
  nameString = "Cobalt";

  //generates several randoms for position, momentum, and probabilities. 
  G4double phi = G4RandFlat::shoot(-0.28, 6.002);
  G4double theta = G4RandFlat::shoot(-0.0005, 3.142);
  G4double check = G4RandFlat::shoot(0.1,1.1);
  G4double radi = G4RandFlat::shoot(0.01, 0.1)*cm;

  //Define the position based on the randoms generated above. current settings: within the laser plates.
  ypos = radi*sin(phi);
  xpos = radi*cos(phi);
  zpos = 0.1*cm; 

  //Generate the initial momentum based on the random variables above. 
  momx = sin(theta)*cos(phi);
  momy = sin(theta)*sin(phi);
  momz = cos(theta);
  partEneg = fParticleGun->GetParticleEnergy()/keV;
    
  //defines a variable for rodHeight, necessary for the laser.
  G4double rodHeight = 0;

  //defines the Ar39 energy and probability spectra
  G4double Ar39eneg[28];
  Ar39eneg[0] = 0.01;
  for (int i=1;i<28;++i){ Ar39eneg[i] = 20.0*i; }
  G4double Ar39prob[] = { .195, .395, .605,   .825, 1.05, 1.28,  1.515, 1.755, 2.0, 2.25, 
			  2.50, 2.75, 2.995, 3.235, 3.47, 3.695, 3.905, 4.1, 4.28, 4.445, 
			  4.59, 4.715, 4.815, 4.89, 4.94, 4.97,  4.99,  5.0 };

  //Obtain the rod Height from the detector to produce the photons at the right location.
  const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //detectorConstruction* detector = new detectorConstruction;
  
  /*  if (detector->IsRandom()) {
    G4String randoms = detector->GetRandoms();
    G4cout << randoms << G4endl;
    }*/

  rodHeight = detector->GetRodHeight();
  zpos = zpos + rodHeight;
  
  //determine if the Laser or Sodium is on from the detector. Note: Laser overwrites Sodium, so only one runs if both are on.
  Laser = detector->GetfLaser();
  if (!Laser) {
    Sodium = detector->GetfSodium();
    Ar39 = detector->GetfAr39();
    darkMatter = detector->GetDarkMatter();
    //    G4cout << "Primary found sodium as " << Sodium << G4endl;
  }
  
  //begin the primary generation for laser runs.
  if (Laser) {
    nameString = "Laser";
    zpos = rodHeight -.95*cm;
    
    //Make the initial particles opticalphotons.
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="opticalphoton");
    fParticleGun->SetParticleDefinition(particle);
    
    //get the particle energy entered either through the above lines or from a macro.
    //then set the number of photons to be generated to 1000 for 532 and 10000 for 213
    //note: more needed for 213 due to non-ideal detector killing most of them.
    G4int totalops = 1000;
    if (partEneg*keV > 5*eV) {
      totalops = 4000;
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

  //output at the start of an event.
  //G4cout << "primarygeneratornote: Start Event" << G4endl;

  if (Ar39) {
    nameString = "Argon39";
    G4double test = G4RandFlat::shoot(0.00001,5.0);
    G4double energy = 200.0;
    for (int i=0;i<28;++i) {
      if (test < Ar39prob[i]) {
	energy = Ar39eneg[i];
	break;
      }
    }
    
    radi = G4RandFlat::shoot(0.001,75.0)*cm;
    ypos = radi*sin(phi);
    xpos = radi*cos(phi);
    zpos = G4RandFlat::shoot(-50.0, 50.0)*cm;
    partEneg = G4RandFlat::shoot(energy,(energy+20.0));
    fParticleGun->SetParticleEnergy(partEneg*keV);
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
    fParticleGun->SetParticleDefinition(particle);
    
    //G4cout << "InitialConditions:" << '\t' << xpos  << '\t' << ypos << '\t' << zpos << '\t' << partEneg << G4endl;
  }

  if (darkMatter) {
    nameString = "DarkMatter";
    G4double DMeneg[137];
    for (int i=0;i<96;++i){ DMeneg[i] = 10.0*i+50.0; }
    for (int i=96;i<137;++i){ DMeneg[i] = 100.0*(i-95)+1000.0; }
    G4double DMprob[] = { 0.185888,   0.303587,   0.390679,   0.4586973,  0.5136096,
			  0.5592977,  0.5980322,  0.6311708,  0.6597705,  0.6851431,
			  0.7072998,  0.7272171,  0.7450814,  0.7611439,  0.7758223,
			  0.7889518,  0.8011255,  0.8122191,  0.8225298,  0.83205342,
			  0.84091939, 0.84906535, 0.8566666,  0.86383961, 0.87039379,
			  0.87663503, 0.88247568, 0.88795573, 0.89303931, 0.89790348,
			  0.90246058, 0.90669474, 0.91077831, 0.91451835, 0.91815368,
			  0.92155666, 0.92474904, 0.92780378, 0.93075793, 0.93357443,
			  0.9361674,  0.93862801, 0.94099745, 0.9432563,  0.94542279,
			  0.94743516, 0.94938989, 0.95125109, 0.9530017,  0.95466348,
			  0.95625291, 0.95779645, 0.95921411, 0.96066412, 0.96200766,
			  0.96333061, 0.96455474, 0.96574475, 0.96687182, 0.96793712,
			  0.96895125, 0.969940672,0.970870093,0.971797749,0.972715993,
			  0.97357306, 0.974384833,0.9751919,  0.975900731,0.976659562,
			  0.977403099,0.978078988,0.978746642,0.979394884,0.979983125,
			  0.980562543,0.981143137,0.981649613,0.982176089,0.9826808,
			  0.983166687,0.983610221,0.984048461,0.984499054,0.984921411,
			  0.985346709,0.985718477,0.98608201, 0.986430249,0.986777311,
			  0.987133197,0.987475553,0.987802027,0.988111442,0.988402621,
			  0.990951471,0.992793842,0.994185033,0.995288574,0.996145642,
			  0.996804472,0.997335065,0.997780952,0.998137426,0.998451547,
			  0.998705667,0.998928611,0.999086848,0.999245085,0.999375086,
			  0.999480381,0.999575676,0.999636853,0.999703325,0.99975509,
			  0.999799796,0.999836267,0.999864503,0.999890974,0.999909798,
			  0.999930386,0.999943327,0.999957445,0.999972151,0.999976857,
			  0.999982739,0.999987445,0.999989798,0.99999274, 0.999994504,
			  0.999996269,0.999997445,0.999998622,0.99999921, 0.999999798, 1.0};
   
    G4double test = G4RandFlat::shoot(0.1,1.1);
    test = test - 0.1;
    G4double energy = 200.0;
    G4double energy2 = 200.0;
    for (int i=0;i<137;++i) {
      if (test < DMprob[i]) {
	energy = DMeneg[i];
	energy2 = DMeneg[i+1];
	break;
      }
    }

    G4int Z = 18, A = 40;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;

    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
 
    radi = G4RandFlat::shoot(0.001,80.0)*cm;
    ypos = radi*sin(phi);
    xpos = radi*cos(phi);
    zpos = G4RandFlat::shoot(-65.0, 65.0)*cm;
    partEneg = G4RandFlat::shoot(energy,energy2);
    fParticleGun->SetParticleEnergy(partEneg*keV);
    
    //G4cout << "InitialConditions:" << '\t' << xpos  << '\t' << ypos << '\t' << zpos << '\t' << partEneg << G4endl;
  }

  //if Sodium, overwrite any input from the macro file and create the 1.2 MeV gamma at the position (0,0,0).
  if (Sodium) {
    nameString = "Sodium";
    partEneg=1200;
    fParticleGun->SetParticleEnergy(partEneg*keV);
    xpos = 0*cm;
    ypos = 0*cm;
    Laser = false;
    //G4cout << "First Sodium Loop" << G4endl;
  }
    
  //if not Laser, generate a primary for a single particle (Sodium 1.2 MeV or other from macro input).
  if (!Laser) {
    fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  //if Sodium, determine if the 2 511 keVs are going to be generated and do so.
  if (Sodium && check > 0.2) {
    //G4cout << "Start Second Sodium Loop" << G4endl;
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
    //G4cout << "End Second Sodium Loop" << G4endl;
  }//*/
  
  //more complete output for starting the event, showing momentum and position initiated.
  //  G4cout << "S0" << "  particleinit  " << momx << "  " << momy << "  " << momz << "  " <<
  // xpos/cm << "  " << ypos/cm << "  " << zpos/cm <<G4endl;

}
