#ifndef G4CCMActionInitialization_h
#define G4CCMActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class G4CCMDetectorConstruction;

class G4CCMActionInitialization : public G4VUserActionInitialization {
public:
    G4CCMActionInitialization(G4CCMDetectorConstruction const * det);
    ~G4CCMActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    const G4CCMDetectorConstruction* fDetector = nullptr;
};

#endif

