#ifndef G4CCMRUNMANAGER_H
#define G4CCMRUNMANAGER_H

#include <dataclasses/I3Map.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>

#include <g4-larsim/g4classes/G4CCMPMTHit.h>
#include <g4-larsim/g4classes/G4CCMScintHit.h>

#include <icetray/CCMPMTKey.h>

#include <simclasses/CCMMCPE.h>

#include <vector>

#include <G4TaskRunManager.hh>
#include <G4DataVector.hh>
#include <G4SystemOfUnits.hh>

/**
 * Implementation of G4RunManager
 */
class G4CCMRunManager: public G4TaskRunManager {
public:
    G4CCMRunManager();

    static G4CCMRunManager* GetCCMRunManager() {return (G4CCMRunManager*)GetRunManager();}

    void SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries);
    void SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries);
};

#endif
