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
#include "globals.hh"
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
  G4bool Cosmic = false;
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

  //Obtain the rod Height from the detector to produce the photons at the right location.
  const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  
  rodHeight = detector->GetRodHeight();
  zpos = zpos + rodHeight;
  
  //determine if the Laser or Sodium is on from the detector. Note: Laser overwrites Sodium, so only one runs if both are on.
  Laser = detector->GetfLaser();
  if (!Laser) {
    Sodium = detector->GetfSodium();
    Ar39 = detector->GetfAr39();
    darkMatter = detector->GetDarkMatter();
    Cosmic = detector->GetfCosmic();
    //    G4cout << "Primary found sodium as " << Sodium << G4endl;
  }

  //call methods to generate events of non-laser types. the order of priority is as below, so only one kind of event can be generated at a time. See below methods for details.
  if (Ar39) {
    shootArgon(anEvent);
  } else if (darkMatter) {
    shootDarkMatter(anEvent);
  } else if (Sodium) {
    shootSodium(anEvent);
  } else if (Cosmic) {
    shootCosmic(anEvent);
  }
  
  //begin the primary generation for laser runs.
  if (Laser) {
    nameString = "Laser";
    zpos = rodHeight -.95*cm; //the Laser generatio point is a bit under a cm above the rodHeight position (in the diffuser glass).
    
    //Make the initial particles opticalphotons.
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="opticalphoton");
    fParticleGun->SetParticleDefinition(particle);
    
    //get the particle energy entered either through the above lines or from a macro.
    //then set the number of photons to be generated to 2000 for 532 and 10000 for 213
    //note: more needed for 213 due to non-ideal detector killing most of them.
    G4int totalops = 2000;
    if (partEneg*keV > 5*eV) {
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
}

//method to generate Argon39 Events
void primaryGenerator::shootArgon(G4Event* anEvent) {

  //defines the Ar39 energy and probability spectra
  G4double Ar39eneg[28];
  Ar39eneg[0] = 0.01;
  for (int i=1;i<28;++i){ Ar39eneg[i] = 20.0*i; }
  G4double Ar39prob[] = { .195, .395, .605,   .825, 1.05, 1.28,  1.515, 1.755, 2.0, 2.25, 
			  2.50, 2.75, 2.995, 3.235, 3.47, 3.695, 3.905, 4.1, 4.28, 4.445, 
			  4.59, 4.715, 4.815, 4.89, 4.94, 4.97,  4.99,  5.0 };
  
  nameString = "Argon39";
  
  //pulls an energy value from the Ar39 energy distribution
  G4double test = G4RandFlat::shoot(0.00001,5.0);
  G4double energy = 200.0;
  for (int i=0;i<28;++i) {
    if (test < Ar39prob[i]) {
      energy = Ar39eneg[i];
      break;
    }
  }
  partEneg = G4RandFlat::shoot(energy,(energy+20.0));
  fParticleGun->SetParticleEnergy(partEneg*keV);
  
  //defines a position for the generation of the argon39 recoil electron
  G4double phi = G4RandFlat::shoot(-0.28, 6.002);
  G4double radi2 = G4RandFlat::shoot(0.001,12735.1225);
  radi2 = sqrt(radi2)*cm;
  ypos = radi2*sin(phi);
  xpos = radi2*cos(phi);
  zpos = G4RandFlat::shoot(-70.79, 70.79)*cm;

  //defines the Ar39 electron
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
    
  //shoots the particle
  fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//defines a method for generating a dark matter event
void primaryGenerator::shootDarkMatter(G4Event* anEvent) {

  //define the Dark Matter energy probability distribution
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
  
  //pulls an upper and lower energy from the Dark Matter energy table
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
  
  //defines the recoiling argon nucleus 
  G4int Z = 18, A = 40;
  G4double ionCharge   = 0.*eplus;
  G4double excitEnergy = 0.*keV;
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
  fParticleGun->SetParticleDefinition(ion);
  fParticleGun->SetParticleCharge(ionCharge);
  
  //defines the initial position and energy for the recoiling nucleus
  G4double radi2 = G4RandFlat::shoot(0.001,12735.1225);
  radi2 = sqrt(radi2)*cm;
  G4double phi = G4RandFlat::shoot(-0.28, 6.002);
  ypos = radi2*sin(phi);
  xpos = radi2*cos(phi);
  zpos = G4RandFlat::shoot(-70.79, 70.79)*cm;
  partEneg = G4RandFlat::shoot(energy,energy2);
  fParticleGun->SetParticleEnergy(partEneg*keV);
  
  //shoots the particle
  fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//defines a method for generating a Sodium22 event
void primaryGenerator::shootSodium(G4Event* anEvent) {
  //rewrites the particle energy for the particle creation
  nameString = "Sodium";
  partEneg=1200;
  fParticleGun->SetParticleEnergy(partEneg*keV);
  G4double check = G4RandFlat::shoot(0.1,1.1);


  //shoots the particle
  fParticleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //if Sodium, determine if the 2 511 keVs are going to be generated and do so.
  if (check > 0.2) {
    //G4cout << "Start Second Sodium Loop" << G4endl;
    //create a new momentum angle.
    G4double phi = G4RandFlat::shoot(-0.28,6);
    G4double theta = G4RandFlat::shoot(-0.14,3);
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
  }
}

//method to generate a cosmic muon event
void primaryGenerator::shootCosmic(G4Event* anEvent) {
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
  //  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  
  fParticleGun->SetParticleDefinition(particle);

  // Sphere model

  // This model was used in B. Olmos Yáñez and A. A. Aguilar-Arevalo, "A method to measure the integral vertical intensity and angular distribution of atmospheric
  // muons with a stationary plastic scintillator bar detector", Nuclear Inst. and Methods in Physics Research, A 987 (2021) 164870, doi.org/10.1016/j.nima.2020.164870.

  double R = 200.0;// Sphere radius
  double px = 100.0;// Plane dimensions
  double py = 100.0;//

  double theta;// Theta (zenithal) angle

  bool a = true;
  while (a) {// Loop to sample theta by acceptance-rejection method
    double p = G4RandFlat::shoot(0.0, 1.0);//
    theta = G4RandFlat::shoot(0.0, 3.1416/2);
    double q = pow(cos(theta), 2)*sin(theta);
    if (p<=q) a = false;
  }

  double phi = G4RandFlat::shoot(0.0, 2*3.1416);// Phi angle uniformly sampled
  double X = R*sin(theta)*cos(phi);// Coordinates (X,Y,Z) over the sphere
  double Y = R*sin(theta)*sin(phi);//
  double Z = R*cos(theta);

  double u = G4RandFlat::shoot(-px/2, px/2);// Point over the plane uniformly sampled
  double v = G4RandFlat::shoot(-py/2, py/2);//
  xpos = (X+u*cos(theta)*cos(phi)-v*sin(phi))*cm;// Muon initial position (x0,y0,z0)
  ypos = (Y+u*cos(theta)*sin(phi)+v*cos(phi))*cm;//
  zpos = (Z-u*sin(theta))*cm;


  momx = -sin(theta)*cos(phi);// Muon unitary direction vector
  momy = -sin(theta)*sin(phi);//
  momz = -cos(theta);


  fParticleGun->SetParticlePosition(G4ThreeVector(xpos,ypos,zpos));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));

  // *****   Muons   *****

  // Smith & Duller model
  // This model gives the muon differential intensity at ground level

  // Model taken from J. A. Smith, N. M. Duller, "Effects of pi meson decay-absorption phenomena on the high-energy mu meson zenithal variation near sea level",
  // Journal of Geophysical Research 64 (12) (1959) 2297-2305. doi:10.1029/JZ064i012p02297.

  double Emin = 1e0;// Min kinetic energy
  double Emax = 1e5;// Max kinetic energy
  const int ee = 10000;// Energy-array size
  double ES[ee];// Energy array from Emin to Emax
  double dE = (Emax-Emin)/ee;// Steps

  double Eu;// Kinetic energy variable
  double Au = 2e9;                      // Parameter values taken from S. Chatzidakis, S. Chrysikopoulou, L. Tsoukalas, "Developing a cosmic ray muon sampling capability 
  double gu = 2.645;                    // for muon tomography and monitoring applications", Nuclear Instruments and Methods in Physics Research Section A: Accelerators, 
  double ru = 0.76;                     // Spectrometers, Detectors and Associated Equipment 804 (2015) 33-42. dx.doi.org/10.1016/j.nima.2015.09.033
  double au = 2.5;//
  double y0u = 1000.0;
  double bmu = 0.80;
  double cu = 299792458.0e2;
  double mmu = 105.7/pow(cu,2);
  double t0mu = 2.2e-6;
  double r0u = 0.00129;
  double Epu;
  double Bmu = bmu*mmu*y0u*cu/(t0mu*r0u);// Equation (15) in Chatzidakis
  double Pmu;
  double lpu = 120.0;
  double bu = 0.771;
  double mpu = 139.6/pow(cu,2);
  double t0pu = 2.6e-8;
  double jpu = mpu*y0u*cu/(t0pu*r0u);// Equation (11) in Chatzidakis
  
  for (int j=0; j<ee; j++) {    // This loop creates the energy array from Emin to Emax in steps of size dE
    Eu = Emin+j*dE;// Kinetic energy
    Epu = (Eu+au*y0u*(1.0/cos(theta)-0.100))/ru;// Definitions (Equation (8) in Smith)
    Pmu = pow(0.100*cos(theta)*(1-(au*(y0u/cos(theta)-100)/(ru*Epu))),(Bmu/((ru*Epu+100*au)*cos(theta))));//             (Equation (9) in Smith)
    ES[j] = Au*(pow(Epu,-gu))*Pmu*lpu*bu*jpu/(Epu*cos(theta)+bu*jpu);// Energy-array elements (Equation (17) in Chatzidakis)
  }

  int nbins = ee;
  G4RandGeneral GenDist(ES,nbins);          // Energy distribution for sampling
  //  double E = Emin + (GenDist.shoot())*(Emax-Emin);   // Energy sampling

  double add = GenDist.shoot();
  double E = Emin + (add)*(Emax-Emin);   // Energy sampling
  
  fParticleGun->SetParticleEnergy(E*MeV);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  //  G4cout << "Generated Muon: " << E << "\t" << G4ThreeVector(momx,momy,momz) << '\t' << G4ThreeVector(xpos/cm,ypos/cm,zpos/cm) << G4endl; 

}
