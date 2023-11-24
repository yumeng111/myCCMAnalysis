#ifndef _TOPSIM_G4INTERFACE_H_
#define _TOPSIM_G4INTERFACE_H_

#include <g4-larsim/g4classes/G4CCMRunManager.h>
#include <icetray/I3Logging.h>

#ifdef G4VIS_USE
class G4VisManager;
#endif

class I3Particle;
class G4IceTopTank;
class G4IceScintStation;
class G4CCMDetectorConstruction;

/**
 * Top-level class to handle Geant4. All global things are initialized here (run manager, visualization manager, detector construction, physics list and user actions).
 */

class G4Interface {
public:
    G4Interface(const std::string& visMacro="");
    ~G4Interface();

    //  Static method which returns the singleton pointer to this class
    static G4Interface* GetInstance() {return g4Interface_;}

    /// Add a tank to the geometry. Should not be called after initialized.
    void InstallTank(G4IceTopTank* tank);

    /// Add a scintillator to the geometry. Should not be called after initialized.
    void InstallScintillator(G4IceScintStation* scint);

    /// Initialize event. Most Geant4 global things are initialized the first time this is called.
    void InitializeEvent();
    /// To be called after simulating each IceTray event.
    void TerminateEvent();
    /// Simulate a single particle (InitializeEvent must be called first)
    void InjectParticle(const I3Particle& particle);

private:
    void Initialize();

    static G4Interface* g4Interface_;

    G4CCMRunManager runManager_;

#ifdef G4VIS_USE
    G4VisManager* visManager_;
#endif

    G4CCMDetectorConstruction* detector_;
    bool initialized_;
    bool eventInitialized_;
    std::string visMacro_;
    bool triangulatedSnow_;
    bool useScintillator_;

    SET_LOGGER("G4Interface");
};

#endif
