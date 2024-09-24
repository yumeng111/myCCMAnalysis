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

#include <vector>
#include <string>
#include <algorithm>

class CCM200Response : public CCMDetectorResponse {
private:
    CCM200Response operator= (const CCM200Response& rhs);

    std::string visMacroFile_;

    bool PMTSDStatus_; // turn PMT SD on/off
    bool LArSDStatus_; // turn fiducial LAr SD on/off
    bool SourceRodIn_; // place source rod into detector geometry
    double SourceRodLocation_; // location of end of source rod
    bool CobaltSourceRun_; // true for cobalt57
    bool SodiumSourceRun_; // true for sodium22
    double SingletTau_; // set LAr singlet time constant
    double TripletTau_; // set LAr triplet time constant
    double Rayleigh128_; // set rayl scattering length for 128nm light
    double UVAbsLength_; // set uv abs length at 128nm
    double WLSNPhotonsEndCapFoil_;
    double WLSNPhotonsSideFoil_;
    double WLSNPhotonsPMT_;
    double EndCapFoilTPBThickness_;
    double SideFoilTPBThickness_;
    double PMTTPBThickness_;
    bool TimeCut_; // true ends all events after 200 nsec
    bool KillCherenkov_; // true turns off cherenkov light
    long RandomSeed_; // random seed for geant4

    std::shared_ptr<G4Interface> g4Interface_ = nullptr;

    SET_LOGGER("CCM200Response");

public:

    CCM200Response(const I3Context& context);
    virtual ~CCM200Response() override;

    virtual void Configure() override;

    // initialize geant4 detector
    virtual void Initialize() override;

    virtual void SimulateEvent(I3Particle const & primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) override;
    virtual void SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries) override;

    virtual void DestroyInterface() override;
};

#endif // CCM200RESPONSE_H


