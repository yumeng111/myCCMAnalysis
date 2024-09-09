#ifndef G4CCMVisActionInitialization_h
#define G4CCMVisActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class G4CCMDetectorConstruction;

class G4CCMVisActionInitialization : public G4VUserActionInitialization {
public:
    G4CCMVisActionInitialization(G4CCMDetectorConstruction const * det);
    ~G4CCMVisActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    const G4CCMDetectorConstruction* fDetector = nullptr;
};

#endif

