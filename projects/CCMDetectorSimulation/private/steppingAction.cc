/*
  stepping action for CCM simulation

This code describes the functions to be applied at every step of the simulation. 
It also contains most of the output generation, currently sent to a text file through the G4cout.
Adjustment to the output, perhaps to better fit the .root data style, should be done through this code.

Current active functions:
 */
#include "CCMAnalysis/CCMDetectorSimulation/steppingAction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/eventAction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/detectorConstruction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/primaryGenerator.hh"

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
  //  G4cout << "Step Made" << G4endl;
  //get the volume at the end of the step. If none, return (indicates particle was just created and so should be ignored).
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() == NULL){
    return;
  }

  //Defines variables for the current volume, volume before the step was made, and origin volume of the track.
  G4String particlen = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4LogicalVolume* volume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* volumep = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* volumei = step->GetTrack()->GetOriginTouchableHandle()->GetVolume()->GetLogicalVolume();

  //Get the names of the three volumes defined above. 
  G4String volname = volume->GetName();
  G4String volpname = volumep->GetName();
  G4String voliname = volumei->GetName();

  if (volname == "expHall" || volpname == "expHall") {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }

  //Get the first three characters of the post and pre step volumes (for checking if PMT)
  G4String testn = volname.substr(0,3);
  G4String testnp = volpname.substr(0,3);

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
  G4double topx = 0;
  G4double topy = 0;
  G4double topz = 0;
  G4double kinEn = 0;
  G4ThreeVector mom = G4ThreeVector(0,0,0);
  G4ThreeVector pos = G4ThreeVector(0,0,0);
  G4double time = 0;
  G4double angle = 0;
  G4double magpos = 0;
  G4double dot = 0;
  G4int col = 0;
  G4int row = 0;
  G4int row2 = 0;
  G4bool coated = false;

  //obtain the kinetic energy of the particle.
  kinEn = step->GetTrack()->GetKineticEnergy();
  time = step->GetTrack()->GetGlobalTime();
  if (time >= 20000*ns) {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }

  if (process == "nCapture") {
    kinEn = step->GetPreStepPoint()->GetKineticEnergy();
    G4cout << process << '\t' << kinEn << '\n';
  }

  G4int sNum = step->GetTrack()->GetCurrentStepNumber();
  if ((particlen == "gamma" || particlen == "e-") && sNum==1) {
    //get the original process.
    G4String creatorProcess = step->GetTrack()->GetCreatorProcess()->GetProcessName();

    if (creatorProcess == "nCapture") {
      G4cout << particlen << '\t' << creatorProcess << '\t' << kinEn << '\n';
    }
  }
  //if the particle is not an optical photon and is in the Fiducial volume, 
  // calculate the energy deposited during the step and output it. 
  //currently off for root based output, lacking the mechanics to currently use it.
  /*if (particlen != "opticalphoton" && volpname=="Fiducial" && voliname != "Fiducial") {
    G4double kinEnpre = 0;
    G4double kinEnpost = 0;
    kinEnpre = step->GetPreStepPoint()->GetKineticEnergy();
    kinEnpost = step->GetPostStepPoint()->GetKineticEnergy();
    G4double Endep = kinEnpre-kinEnpost;

    //G4cout << particlen << '\t' << Endep/keV << G4endl;
    return;
    }//*/
  
  //Check if the volume is a PMT, then get the energy and ingoing angle, as well as the full name of the pmt.
  if (particlen == "opticalphoton" && testn == "PMT") {
    //skip any photons above 5 eV (wavelength < 250 nm) as the PMTs can't detect them.
    if (kinEn >= 5*eV) {
      step->GetTrack()->SetTrackStatus(fStopAndKill);
      return;
    }

    //get the pmt row and column. rows 0 and 6 are for top and bottom of CCM200
    if (testn == "PMT") {
      std::istringstream iss (volname.substr(volname.find('C')+1,volname.find('R')));
      iss >> col; 
      std::istringstream isr (volname.substr(volname.find('R')+1,volname.find('R')+2));
      isr >> row;
      if (volname.find("coat") != std::string::npos) {
	coated = true;
      }
    } else if (testnp == "PMT") {
      std::istringstream iss (volpname.substr(volpname.find('C')+1,volpname.find('R')));
      iss >> col; 
      std::istringstream isr (volpname.substr(volpname.find('R')+1,volpname.find('R')+2));
      isr >> row;
      if (volname.find("coat") != std::string::npos) {
	coated = true;
      }
    }

    //G4cout << row << '\t' << col << G4endl;
    //Get the location of the PMT based on row and column. 
    //Caution: the simulation has row numbering backwards from the mapping, with row 1 on bottom and row 5 on top.//No longer applies, old comment.
    angle = (2*CLHEP::pi/24)*(col-1);
    pmtx = 96*std::cos(angle);
    pmty = 96*std::sin(angle);
    pmtz = (3-row)*22.86;

    //Get the position of the top of the pmt
    topx = 90*std::cos(angle);
    topy = 90*std::sin(angle);
    topz = (3-row)*22.86;

    //get pmt location of rows 0 and 6 (top and bottom, ccm200)
    if (row == 0 || row == 6) {
      pmtz = 61.6;
      if (row == 6) {
	pmtz = -61.6;
      }
      int col2 = col%100;
      row2 = col/100;
      //G4cout << col << '\t' << col2 << '\t' << row  << '\t' << row2 << G4endl;
      if (row2 == 4) {
	angle = (col2-1)*2*CLHEP::pi/20;
	pmtx = 85.5*std::cos(angle);
	pmty = 85.5*std::sin(angle);
      } else if (row2 == 3) {
	angle = (col2-1)*2*CLHEP::pi/15;
	pmtx = 64.2*std::cos(angle);
	pmty = 64.2*std::sin(angle);
      } else if (row2 == 2) {
	angle = (col2-1)*2*CLHEP::pi/10;
	pmtx = 42.0*std::cos(angle);
	pmty = 42.0*std::sin(angle);
      } else if (row2 == 1) {
	angle = (col2-1)*2*CLHEP::pi/5;
	pmtx = 20.0*std::cos(angle);
	pmty = 20.0*std::sin(angle);
      }
      topx = pmtx;
      topy = pmty;
      topz = pmtz*.9;
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

    //define the radial vector of the PMT matched to the location of photon impact.
    G4double vectx = xx-pmtx;
    G4double vecty = xy-pmty;
    G4double vectz = xz-pmtz;
    
    //Calculate the ingoing angle based on radial vector and momentum vector.
    magpos = std::sqrt(vectx*vectx+vecty*vecty+vectz*vectz);
    //dot = cos(ingoing angle).
    dot = (vectx*px+vecty*py+vectz*pz)/magpos;

    //define the pointing vector of the PMT
    vectx = pmtx-topx;
    vecty = pmty-topy;
    vectz = pmtz-topz;

    //Calculate the ingoing angle based on pointing vector and momentum vector.
    magpos = std::sqrt(vectx*vectx+vecty*vecty+vectz*vectz);
    G4double dot2 = (vectx*px+vecty*py+vectz*pz)/magpos;
    
    //combine the angles of radial and popinting vector to get the overall cos2 angle
    dot = dot*std::sqrt(std::sqrt(dot2));

    //terminate any photon absorbed by the PMT.
    step->GetTrack()->SetTrackStatus(fStopAndKill);
  
    //get the original process.
    G4String creatorProcess = step->GetTrack()->GetCreatorProcess()->GetProcessName();

    //abbreviated output with only the relevant information (pmt, photon energy, time, and angle)
    fEventAction->AddHit(row,col,coated,kinEn/eV,time/ns,dot,creatorProcess);
    //G4cout << row << '\t' << col << '\t' << coated << '\t' << kinEn/eV << '\t' << time/ns << '\t' << dot << G4endl;
    //G4cout << volname << '\t' << kinEn/eV << '\t' << time/ns << '\t' << dot << G4endl;
    return;
  }

  //If the particle is a optical photon, and has just transferred from the reflector foil on the top or bottom,
  //alter the direction to induce a slight randomness due to the unsmoothness of the foils there.
  //if ( (particlen == "opticalphoton" && volname=="TPBfoilb" && volpname == "ptfefoil") || (particlen == "opticalphoton" && volname=="tpbbcone" && volpname == "ptfebcone") ) {
  
  //CCM200 alter path of visible photons in tpb volume of any sort. 
  /*if ( (particlen == "opticalphoton" && volname.find("TPB") != std::string::npos) || (particlen == "opticalphoton" && volname.find("tpb") != std::string::npos) ) {
    const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    //G4cout << "Stepping calls to primaries: " << detector->GetfSodium() << '\t' << detector->GetfLaser() << '\t' << detector->GetDarkMatter() << '\t' << detector->GetALP() << '\n';

    if (kinEn >= 4*eV) {
      return;
    }

    mom = step->GetPostStepPoint()->GetMomentumDirection();
    px = mom.x();
    py = mom.y();
    pz = mom.z();
    
    //creates three random variables, two angles and a test.
    //the gaussian angles are for a slight rotation biased towards maintaing the same direction, 
    G4double randwide = detector->GetUnsmooth();
    G4double test = G4RandFlat::shoot(0.1,1.1);
    if (test > randwide+0.1) {
      //G4double thmax = randwide/2.0;
      //G4double phi = G4RandGauss::shoot(0,randwide);
      //G4double theta = G4RandGauss::shoot(0,thmax);
      
      //CCM200 randmization of photons in TPB: total random as in white paint. 
      G4double pi = CLHEP::pi;
      G4double phi = G4RandFlat::shoot(0.0,pi);
      G4double theta = G4RandFlat::shoot(0.0,2*pi);
    
      //Use the two angles to calculate a new outgoing angle based on a rotation from the original
      G4double px1 = px;
      G4double py1 = py;
      G4double pz1 = pz;
      
      px = cos(theta)*px1-sin(theta)*py1;
      py = cos(phi)*sin(theta)*px1+cos(phi)*cos(theta)*py1-sin(phi)*pz;
      pz = sin(phi)*sin(theta)*px1+sin(phi)*cos(theta)*py1+cos(phi)*pz;
      
      //make sure the z direction of the output matches the z direction of the input.
      if ((pz < 0 && pz1 > 0) || (pz > 0 && pz1 < 0)) {
	pz = -pz;
      }
      mom = G4ThreeVector(px,py,pz);
      step->GetPostStepPoint()->SetMomentumDirection(mom);
    }
    }//*/

  //  G4cout << "maincodenote: end stepping Action" << G4endl;
}
