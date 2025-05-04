#ifndef CCM200RESPONSE_H
#define CCM200RESPONSE_H
// standard library stuff
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/g4classes/G4Interface.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3ServiceBase.h"

#include "phys-services/I3RandomService.h"

#include "simclasses/CCMMCPE.h"
#include "simclasses/PhotonSummary.h"
#include "simclasses/DetectorResponseConfig.h"

#include <vector>
#include <string>
#include <algorithm>

class CCM200Response : public CCMDetectorResponse {
private:
    CCM200Response operator= (const CCM200Response& rhs);

    std::string visMacroFile_;

    bool SaveAllEnergyLossesTree_;

    bool VetoSDSaveEnergyLossesVector_;
    bool VetoSDSaveEnergyLossesTree_;
    bool VetoSDPruneTree_;

    bool InteriorSDSaveEnergyLossesVector_;
    bool InteriorSDSaveEnergyLossesTree_;
    bool InteriorSDPruneTree_;

    bool KillNeutrinos_ = false;
    bool KillPhotons_ = false;
    bool KillScintillation_ = false;
    bool KillCherenkov_ = false;

    bool TimeCut_ = false;
    bool DetailedPhotonTracking_ = false;
    bool TrackParticles_ = false;
    bool TrackEnergyLosses_ = false;

    double SimulateNuclearRecoils_;
    double G4RangeCut_;
    double G4EDepMin_;
    double G4ETrackingMin_;

    bool RecordHits_;
    bool SourceRodIn_; // place source rod into detector geometry
    double SourceRodLocation_; // location of end of source rod
    bool CobaltSourceRun_; // true for cobalt57
    bool SodiumSourceRun_; // true for sodium22
    bool TrainingSource_; // true for training source data
    double DecayX_; // if training source data, provide X position of decay
    double DecayY_; // if training source data, provide Y position of decay
    double DecayZ_; // if training source data, provide Z position of decay
    double Rayleigh128_; // set rayl scattering length for 128nm light
    double UVAbsA_;
    double UVAbsB_;
    double UVAbsD_;
    double UVAbsScaling_;
    double WLSNPhotonsEndCapFoil_;
    double WLSNPhotonsSideFoil_;
    double WLSNPhotonsPMT_;
    double EndCapFoilTPBThickness_;
    double SideFoilTPBThickness_;
    double PMTTPBThickness_;
    double TPBAbsTau_;
    double TPBAbsNorm_;
    double TPBAbsScale_;
    double MieGG_;
    double MieRatio_;
    double Normalization_;
    double PhotonSampling_;
    long RandomSeed_; // random seed for geant4

    std::shared_ptr<G4Interface> g4Interface_ = nullptr;

    SET_LOGGER("CCM200Response");

public:

    CCM200Response(const I3Context& context);
    virtual ~CCM200Response() override;

    virtual void Configure() override;
    virtual I3FrameObjectPtr GetSimulationConfiguration() override;

    // initialize geant4 detector
    virtual void Initialize() override;

    virtual void SimulateEvent(I3Particle const & primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries, I3MCTreePtr veto_tree, I3MCTreePtr inner_tree, I3VectorI3ParticlePtr veto_vector, I3VectorI3ParticlePtr inner_vector) override;
    virtual void SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> veto_trees, std::vector<I3MCTreePtr> inner_trees, std::vector<I3VectorI3ParticlePtr> veto_vectors, std::vector<I3VectorI3ParticlePtr> inner_vectors) override;

    virtual void DestroyInterface() override;
};

#endif // CCM200RESPONSE_H


