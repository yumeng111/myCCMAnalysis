#ifndef G4CCMCerenkov_h
#define G4CCMCerenkov_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4VProcess.hh"

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

struct PhotonNumInfo {
    G4double avgPhotons;
    G4double omegaMin;
    G4double omegaMax;
};

class G4Material;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4Step;
class G4Track;
class G4VParticleChange;

class G4CCMCerenkov : public G4VProcess {
public:
    explicit G4CCMCerenkov(const G4String& processName = "G4CCMCerenkov",
            G4ProcessType type          = fElectromagnetic);
    ~G4CCMCerenkov();

    explicit G4CCMCerenkov(const G4CCMCerenkov& right);

    G4CCMCerenkov& operator=(const G4CCMCerenkov& right) = delete;

    G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
    // Returns true -> 'is applicable', for all charged particles
    // except short-lived particles.

    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;
    // Build table at a right time

    void PreparePhysicsTable(const G4ParticleDefinition& part) override;
    void Initialise();

    G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);
    // Returns the discrete step limit and sets the 'StronglyForced'
    // condition for the DoIt to be invoked at every step.

    G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double,
            G4ForceCondition*) override;
    // Returns the discrete step limit and sets the 'StronglyForced'
    // condition for the DoIt to be invoked at every step.

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
            const G4Step& aStep) override;
    // This is the method implementing the Cerenkov process.

    //  no operation in  AtRestDoIt and  AlongStepDoIt
    virtual G4double AlongStepGetPhysicalInteractionLength(
            const G4Track&, G4double, G4double, G4double&, G4GPILSelection*) override
    {
        return -1.0;
    };

    virtual G4double AtRestGetPhysicalInteractionLength(
            const G4Track&, G4ForceCondition*) override
    {
        return -1.0;
    };

    //  no operation in  AtRestDoIt and  AlongStepDoIt
    virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override
    {
        return nullptr;
    };

    virtual G4VParticleChange* AlongStepDoIt(const G4Track&,
            const G4Step&) override
    {
        return nullptr;
    };

    void SetTrackSecondariesFirst(const G4bool state);
    // If set, the primary particle tracking is interrupted and any
    // produced Cerenkov photons are tracked next. When all have
    // been tracked, the tracking of the primary resumes.

    G4bool GetTrackSecondariesFirst() const;
    // Returns the boolean flag for tracking secondaries first.

    void SetMaxBetaChangePerStep(const G4double d);
    // Set the maximum allowed change in beta = v/c in % (perCent) per step.

    G4double GetMaxBetaChangePerStep() const;
    // Returns the maximum allowed change in beta = v/c in % (perCent)

    void SetMaxNumPhotonsPerStep(const G4int NumPhotons);
    // Set the maximum number of Cerenkov photons allowed to be generated during
    // a tracking step. This is an average ONLY; the actual number will vary
    // around this average. If invoked, the maximum photon stack will roughly be
    // of the size set. If not called, the step is not limited by the number of
    // photons generated.

    G4int GetMaxNumPhotonsPerStep() const;
    // Returns the maximum number of Cerenkov photons allowed to be
    // generated during a tracking step.

    void SetStackPhotons(const G4bool);
    // Call by the user to set the flag for stacking the scint. photons

    G4bool GetStackPhotons() const;
    // Return the boolean for whether or not the scint. photons are stacked

    G4int GetNumPhotons() const;
    // Returns the current number of scint. photons (after PostStepDoIt)

    G4PhysicsTable* GetPhysicsTable() const;
    // Returns the address of the physics table.

    void DumpPhysicsTable() const;
    // Prints the physics table.

    PhotonNumInfo GetAverageNumberOfPhotons(const G4double charge, const G4double beta, const G4Material* aMaterial, G4MaterialPropertyVector* Rindex) const;

    void DumpInfo() const override {ProcessDescription(G4cout);};
    void ProcessDescription(std::ostream& out) const override;

    void SetVerboseLevel(G4int);
    // sets verbosity

    void SetPhotonSamplingFactor(G4double factor);
    // sets the scaling factor for the number of photons produced

    G4double GetPhotonSamplingFactor() const;
    // returns the scaling factor for the number of photons produced

protected:
    G4PhysicsTable* thePhysicsTable;
    boost::shared_ptr<std::vector<std::vector<G4double>>> RefractiveIndexValsVectors = boost::make_shared<std::vector<std::vector<G4double>>>();
    boost::shared_ptr<std::vector<std::vector<G4double>>> RefractiveIndexEnergyVectors = boost::make_shared<std::vector<std::vector<G4double>>>();

private:
    G4double fMaxBetaChange;

    G4double fPhotonSamplingFactor;

    G4int fMaxPhotons;
    G4int fNumPhotons;

    G4bool fStackingFlag;
    G4bool fTrackSecondariesFirst;

    G4int secID = -1;  // creator modelID

};

inline G4bool G4CCMCerenkov::GetTrackSecondariesFirst() const
{
    return fTrackSecondariesFirst;
}

inline G4double G4CCMCerenkov::GetMaxBetaChangePerStep() const
{
    return fMaxBetaChange;
}

inline G4int G4CCMCerenkov::GetMaxNumPhotonsPerStep() const { return fMaxPhotons; }

inline G4bool G4CCMCerenkov::GetStackPhotons() const { return fStackingFlag; }

inline G4int G4CCMCerenkov::GetNumPhotons() const { return fNumPhotons; }

inline G4PhysicsTable* G4CCMCerenkov::GetPhysicsTable() const
{
    return thePhysicsTable;
}

inline G4double G4CCMCerenkov::GetPhotonSamplingFactor() const { return fPhotonSamplingFactor; }

#endif /* G4CCMCerenkov_h */

