/*
  stepping action for CCM simulation

This code describes the functions to be applied at every step of the simulation. 
It also contains most of the output generation, currently sent to a text file through the G4cout.
Adjustment to the output, perhaps to better fit the .root data style, should be done through this code.

Current active functions:
 */
#include "steppingAction.hh"
#include "eventAction.hh"
#include "detectorConstruction.hh"
#include "primaryGenerator.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>
#include <string>

//Constructor
steppingAction::steppingAction(eventAction* eventaction)
  : G4UserSteppingAction(),
    fEventAction(eventaction)
{}

//Deconstructor
steppingAction::~steppingAction()
{}

//Stepping action: describes the actions to be taken at every stepof the simulation.
//G4Step: an event is made up of ultiple steps, each describing every possible physics process that a particle can take.
//Some example steps: OpAbsorption, elastic scattering, reflection, , Traveling (change logical volume)
//Note: as the stepping action is called for every step, it is advisible to put simple escapes as early as possible to reduce run time.
void steppingAction::UserSteppingAction(const G4Step* step)
{
  G4String particlen = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  //  G4cout << "got particle type " << particlen << G4endl;

  //debugging output lines
  //  G4cout << "volumes defined"<< G4endl;
  //G4cout << "maincodenote: start stepping Action" << G4endl;
  
  //get the volume at the end of the step. If none, return.
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() == NULL){
    //G4cout << "particle volume pointer = null" <<G4endl;
    return;
  }

  //Defines variables for the current volume, volume before the step was made, and origin volume of the track.
  G4LogicalVolume* volume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* volumep = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* volumei = step->GetTrack()->GetOriginTouchableHandle()->GetVolume()->GetLogicalVolume();

  //  G4cout <<"got particle volume" << G4endl;

  //Get the names of the three volumes defined above. 
  G4String volname = volume->GetName();
  G4String volpname = volumep->GetName();
  G4String voliname = volumei->GetName();

  //kill all particles that go through the frame. Needs to be turned off before trying any veto simulations
  if (volname=="Frame"){
    //G4cout << volname << '\t' << kinEn/eV << G4endl;
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }  

  //Get the first three characters of the post and pre step volumes (for checking if PMT)
  G4String testn = volname.substr(0,3);
  G4String testnp = volpname.substr(0,3);

  //  G4cout <<"got particle volume name" << testn << G4endl;

  //Get the process and particle name.
  G4String process  = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //define variables for particle momentum, position, pmt position, energy desposited, kinetic energy, and more.
  G4double px = 0;
  G4double py = 0;
  G4double pz = 0;
  G4double xx = 0;
  G4double xy = 0;
  G4double xz = 0;
  G4double pmtx = 0;
  G4double pmty = 0;
  G4double pmtz = 0;
  G4double kinEn = 0;
  G4double kinEnpre = 0;
  G4double kinEnpost = 0;
  G4ThreeVector mom = G4ThreeVector(0,0,0);
  G4ThreeVector pos = G4ThreeVector(0,0,0);
  G4double time = 0;
  G4double angle = 0;
  G4double magpos = 0;
  G4double dot = 0;

  //obtain the kinetic energy of the particle.
  kinEn = step->GetTrack()->GetKineticEnergy();

  //if the particle is not an optical photon and is in the Fiducial volume, 
  // calculate the energy deposited during the step and output it. 
  if (particlen != "opticalphoton" && volpname=="Fiducial" && voliname != "Fiducial") {
    kinEnpre = step->GetPreStepPoint()->GetKineticEnergy();
    kinEnpost = step->GetPostStepPoint()->GetKineticEnergy();
    G4double Endep = kinEnpre-kinEnpost;

    //G4cout << particlen << '\t' << Endep/keV << G4endl;
    return;
  }

  //If the particle is a optical photon, and has just transferred from the reflector foil on the top or bottom,
  //alter the direction to induce a slight randomness due to the unsmoothness of the foils there.
  if ( (particlen == "opticalphoton" && volname=="TPBfoilb" && volpname == "ptfefoil") || (particlen == "opticalphoton" && volname=="tpbbcone" && volpname == "ptfebcone") ) {
    const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    mom = step->GetPostStepPoint()->GetMomentumDirection();
    px = mom.x();
    py = mom.y();
    pz = mom.z();
    //G4cout << "test pz>0: " << px << '\t' << py  << '\t' << pz << G4endl;
    
    //creates three random variables, two angles and a test.
    //the gaussian angles are for a slight rotation biased towards maintaing the same direction, 
    //the flats commented after them for complete randomness
    G4double randwide = detector->GetUnsmooth();
    G4double thmax = randwide/2.0;
    G4double randtest = G4RandFlat::shoot(0.1,1.1);
    G4double phi = G4RandGauss::shoot(0,randwide);//G4RandFlat::shoot(-.002,6.28);//
    G4double theta = G4RandGauss::shoot(0,thmax);//G4RandFlat::shoot(.01,1.57);//
    
    //Use the two angles to calculate a new outgoing angle based on a rotation from the original
    G4double c1 = cos(phi);
    G4double c2 = cos(theta);
    G4double s1 = sin(phi);
    G4double s2 = sin(theta);
    G4double px1 = px;
    G4double py1 = py;
    G4double pz1 = pz;
    
    px = c2*px1-s2*py1;
    py = c1*s2*px1+c1*c2*py1-s1*pz;
    pz = s1*s2*px1+s1*c2*py1+c1*pz;

    //make sure the z direction of the output matches the z direction of the input.
    if ((pz < 0 && pz1 > 0) || (pz > 0 && pz1 < 0)) {
      pz = -pz;
    }
    //define the proportion of photons deflected thus. Currently, 100% of photons are deflected.
    if (randtest > 0) {
      /*
	//this section: for complete randomness (flat angles above).
      px = sin(theta)*cos(phi);
      py = sin(theta)*sin(phi);
      if (pz > 0) {
	pz = cos(theta);
      } else {
	pz = -cos(theta);
	}//*/
      mom = G4ThreeVector(px,py,pz);
      step->GetPostStepPoint()->SetMomentumDirection(mom);
      //G4cout << "set px,y,z: " << px << '\t' << py  << '\t' << pz << G4endl;
    }    
  }//*/
  
  //Check if the volume is a PMT, then get the energy and ingoing angle, as well as the full name of the pmt.
  if (testn == "PMT" || testnp == "PMT") {
    //skip any photons above 5 eV (wavelength < 250 nm) as the PMTs can't detect them.
    if (kinEn >= 5*eV) {
      //px = 0;
      return;
    }

    //get the pmt row and column. rows 0 and 6 are for top and bottom of CCM200
    G4int col = 0;
    G4int row = 0;
    G4int row2 = 0;
    G4bool coated = false;
    if (testn == "PMT") {
      std::istringstream iss (volname.substr(volname.find('C')+1,volname.find('R')));
      iss >> col; 
      std::istringstream isr (volname.substr(volname.find('R')+1,volname.find('R')+2));
      isr >> row;
      if (row == 0 || row == 6) {
	std::istringstream isr2 (volname.substr(volname.find('R')+2,volname.find('R')+3));
	isr2 >> row2;
      }	
      if (volname.find("coat") != std::string::npos) {
	coated = true;
      }
    } else if (testnp == "PMT") {
      std::istringstream iss (volpname.substr(volpname.find('C')+1,volpname.find('R')));
      iss >> col; 
      std::istringstream isr (volpname.substr(volpname.find('R')+1,volpname.find('R')+2));
      isr >> row;
      if (row == 0 || row == 6) {
	std::istringstream isr2 (volname.substr(volname.find('R')+2,volname.find('R')+3));
	isr2 >> row2;
      }	
      if (volname.find("coat") != std::string::npos) {
	coated = true;
      }
    }

    //Get the location of the PMT based on row and column. 
    //Caution: the simulation has row numbering backwards from the mapping, with row 1 on bottom and row 5 on top.
    angle = (2*CLHEP::pi/24)*(col-1);
    pmtx = 96*std::cos(angle);
    pmty = 96*std::sin(angle);
    pmtz = (row-3)*23.11;

    //get pmt location of rows 0 and 6 (bottom and top, ccm200)
    if (row == 0 || row == 6) {
      pmtz = -68.0;
      if (row == 6) {
	pmtz = 68.0;
	row = 60+row2;
      }else {
	row = 10+row2;
      }
      if (row2 == 4) {
	angle = (col-1)*2*CLHEP::pi/20;
	pmtx = 85.5*std::cos(angle);
	pmty = 85.5*std::sin(angle);
      } else if (row2 == 3) {
	angle = (col-1)*2*CLHEP::pi/15;
	pmtx = 64.2*std::cos(angle);
	pmty = 64.2*std::sin(angle);
      } else if (row2 == 2) {
	angle = (col-1)*2*CLHEP::pi/10;
	pmtx = 42.0*std::cos(angle);
	pmty = 42.0*std::sin(angle);
      } else if (row2 == 1) {
	angle = (col-1)*2*CLHEP::pi/5;
	pmtx = 20.0*std::cos(angle);
	pmty = 20.0*std::sin(angle);
      }
    }
    
    //get the momentum direction and position of the photon when it enters the PMT.
    mom = step->GetTrack()->GetMomentumDirection();
    px = mom.x();
    py = mom.y();
    pz = mom.z();
    pos = step->GetTrack()->GetPosition();
    xx = pos.x();
    xy = pos.y();
    xz = pos.z();
    time = step->GetTrack()->GetGlobalTime();

    //define the radial vector of the PMT matched to the location of photon impact.
    pmtx = xx-pmtx;
    pmty = xy-pmty;
    pmtz = xz-pmtz;
    
    //Calculate the ingoing angle based on radial vector and momentum vector.
    magpos = std::sqrt(pmtx*pmtx+pmty*pmty+pmtz*pmtz);
    //dot = cos(ingoing angle).
    dot = (pmtx*px+pmty*py+pmtz*pz)/magpos;

    //longer potential output with all variables (position, momentum, time, angle)
    /*G4cout << volname << "\t" << particlen << "\t" << kinEn/eV << "\t" << 
      px << "\t" << py << "\t" << pz << "\t" << xx/cm << "\t" << xy/cm << "\t" <<
      xz/cm << "\t" << time/ns << "\t" << dot << G4endl;
    */

    //terminate any photon absorbed by the PMT.
    step->GetTrack()->SetTrackStatus(fStopAndKill);
  
    //abbreviated output with only the relevant information (pmt, photon energy, time, and angle)
    fEventAction->AddHit(row,col,coated,kinEn/eV,time/ns,dot);
    //mctruth->AddHitInformation(row,col,coated,kinEn/eV,time/ns,dot);
    //G4cout << volname << '\t' << kinEn/eV << '\t' << time/ns << '\t' << dot << G4endl;
    return;
  }
  //  G4cout << "maincodenote: end stepping Action" << G4endl;
}
