/*
   stepping action for CCM simulation

   This code describes the functions to be applied at every step of the simulation. 
   It also contains most of the output generation, currently sent to a text file through the G4cout.
   Adjustment to the output, perhaps to better fit the .root data style, should be done through this code.

   Current active functions:
*/
#include "g4-larsim/g4classes/G4CCMUserSteppingAction.h"
#include "g4-larsim/g4classes/G4CCMUserEventAction.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"
#include "g4-larsim/g4classes/G4CCMPrimaryGenerator.h"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>
#include <string>

//Constructor
G4CCMUserSteppingAction::G4CCMUserSteppingAction(G4CCMUserEventAction* eventaction)
    : G4UserSteppingAction(),
    fEventAction(eventaction)
{}

//Deconstructor
G4CCMUserSteppingAction::~G4CCMUserSteppingAction()
{}

// Stepping action: describes the actions to be taken at every stepof the simulation.
// G4Step: an event is made up of ultiple steps, each describing every possible physics process that a particle can take.
// Some example steps: OpAbsorption, elastic scattering, reflection, , Traveling (change logical volume)
// Note: as the stepping action is called for every step, it is advisible to put simple escapes as early as possible to reduce run time.
void G4CCMUserSteppingAction::UserSteppingAction(const G4Step* step)
{  
    // get the volume at the end of the step. If none, return (indicates particle was just created and so should be ignored).
    if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() == NULL){
        return;
    }

    // Defines variables for the current volume, volume before the step was made, and origin volume of the track.
    G4String particlen = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4LogicalVolume* volume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4LogicalVolume* volumep = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4LogicalVolume* volumei = step->GetTrack()->GetOriginTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Get the names of the three volumes defined above. 
    G4String volname = volume->GetName();
    G4String volpname = volumep->GetName();
    G4String voliname = volumei->GetName();

    if (volname == "expHall" || volpname == "expHall") {
        step->GetTrack()->SetTrackStatus(fStopAndKill);
        return;
    }

    // Get the first three characters of the post and pre step volumes (for checking if PMT)
    G4String testn = volname.substr(0,3);
    G4String testnp = volpname.substr(0,3);

    // Get the process and particle name.
    G4String process  = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    // define variables for particle momentum, position, pmt position, energy desposited, kinetic energy, and more.
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

    // obtain the kinetic energy of the particle.
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

    // Check if the volume is a PMT, then get the energy and ingoing angle, as well as the full name of the pmt.
    if (particlen == "opticalphoton" && testn == "PMT") {
        // skip any photons above 5 eV (wavelength < 250 nm) as the PMTs can't detect them.
        if (kinEn >= 5*eV) {
            step->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
        }

        // get the pmt row and column. rows 0 and 6 are for top and bottom of CCM200
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

        // Get the location of the PMT based on row and column. 
        angle = (2*CLHEP::pi/24)*(col-1);
        pmtx = 96*std::cos(angle);
        pmty = 96*std::sin(angle);
        pmtz = (3-row)*22.86;

        // Get the position of the top of the pmt
        topx = 90*std::cos(angle);
        topy = 90*std::sin(angle);
        topz = (3-row)*22.86;

        // get pmt location of rows 0 and 6 (top and bottom, ccm200)
        if (row == 0 || row == 6) {
            pmtz = 61.6;
            if (row == 6) {
                pmtz = -61.6;
            }
            int col2 = col%100;
            row2 = col/100;
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

        // get the momentum direction and position of the photon when it enters the PMT.
        mom = step->GetTrack()->GetMomentumDirection();
        px = mom.x();
        py = mom.y();
        pz = mom.z();
        pos = step->GetTrack()->GetPosition();
        xx = pos.x();
        xy = pos.y();
        xz = pos.z();

        // define the radial vector of the PMT matched to the location of photon impact.
        G4double vectx = xx-pmtx;
        G4double vecty = xy-pmty;
        G4double vectz = xz-pmtz;

        // Calculate the ingoing angle based on radial vector and momentum vector.
        magpos = std::sqrt(vectx*vectx+vecty*vecty+vectz*vectz);
        dot = (vectx*px+vecty*py+vectz*pz)/magpos;

        // define the pointing vector of the PMT
        vectx = pmtx-topx;
        vecty = pmty-topy;
        vectz = pmtz-topz;

        // Calculate the ingoing angle based on pointing vector and momentum vector.
        magpos = std::sqrt(vectx*vectx+vecty*vecty+vectz*vectz);
        G4double dot2 = (vectx*px+vecty*py+vectz*pz)/magpos;

        // combine the angles of radial and popinting vector to get the overall cos2 angle
        dot = dot*std::sqrt(std::sqrt(dot2));

        // terminate any photon absorbed by the PMT.
        step->GetTrack()->SetTrackStatus(fStopAndKill);

        // get the original process.
        G4String creatorProcess = step->GetTrack()->GetCreatorProcess()->GetProcessName();

        // abbreviated output with only the relevant information (pmt, photon energy, time, and angle)
        fEventAction->AddHit(row,col,coated,kinEn/eV,time/ns,dot,creatorProcess);
        return;
    }
}
