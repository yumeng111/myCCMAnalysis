// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...

#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4CCMScintSD.h"
#include "g4-larsim/g4classes/G4CCMRunManager.h"

#include <G4Run.hh>
#include <G4Timer.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4ParticleGun.hh>
#include <G4EventManager.hh>

G4CCMRunManager::G4CCMRunManager(): G4MTRunManager() {}

void G4CCMRunManager::SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) {
    BeamOn(1);
}

void G4CCMRunManager::SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries) {
    BeamOn(primaries.size());
}
