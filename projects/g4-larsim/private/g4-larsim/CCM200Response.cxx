// standard library stuff

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCM200Response.h"
#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/g4classes/G4Interface.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3SingleServiceFactory.h"

#include "phys-services/I3RandomService.h"
#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>
#include <algorithm>

CCM200Response::CCM200Response(const I3Context& context) :
    CCMDetectorResponse(context), PMTSDStatus_(true), LArSDStatus_(true), SourceRodIn_(false), SourceRodLocation_(0.0 * I3Units::cm),
    CobaltSourceRun_(false), SodiumSourceRun_(false), SingletTau_(8.2 * I3Units::nanosecond), TripletTau_(743.0 * I3Units::nanosecond),
    Rayleigh128_(95.0 * I3Units::cm), UVAbsLength_(55.0 * I3Units::cm), WLSNPhotonsEndCapFoil_(0.605), WLSNPhotonsSideFoil_(0.605), WLSNPhotonsPMT_(0.605),
    EndCapFoilTPBThickness_(0.00278035 * I3Units::mm), SideFoilTPBThickness_(0.00278035 * I3Units::mm), PMTTPBThickness_(0.00203892 * I3Units::mm),
    TimeCut_(true), KillCherenkov_(false), RandomSeed_(0){
    AddParameter("PMTSDStatus", "true if tracking photon hits on PMTs", PMTSDStatus_);
    AddParameter("LArSDStatus", "true if tracking scintillation depositions in fiducial LAr", LArSDStatus_);
    AddParameter("SourceRodIn", "true if we want to simulate the sodium source rod", SourceRodIn_);
    AddParameter("SourceRodLocation", "z location of the end of the sodium source rod", SourceRodLocation_);
    AddParameter("CobaltSourceRun", "true if we want to simulate cobalt source pellet", CobaltSourceRun_);
    AddParameter("SodiumSourceRun", "true if we want to simulate sodium source pellet", SodiumSourceRun_);
    AddParameter("SingletTimeConstant", "LAr singlet tau", SingletTau_);
    AddParameter("TripletTimeConstant", "LAr triplet tau", TripletTau_);
    AddParameter("Rayleigh128Length", "Rayleigh scattering length for 128nm light", Rayleigh128_);
    AddParameter("UVAbsLength", "set UV absorption length at 128nm", UVAbsLength_);
    AddParameter("WLSNPhotonsEndCapFoil", "mean number of photons produced per WLS for TPB foils on the end caps of the detector", WLSNPhotonsEndCapFoil_);
    AddParameter("WLSNPhotonsSideFoil", "mean number of photons produced per WLS for TPB foils on the sides of the detector", WLSNPhotonsSideFoil_);
    AddParameter("WLSNPhotonsPMT", "mean number of photons produced per WLS for TPB on PMTs", WLSNPhotonsPMT_);
    AddParameter("EndCapFoilTPBThickness", "thickness of TPB on the endcap foil", EndCapFoilTPBThickness_);
    AddParameter("SideFoilTPBThickness", "thickness of TPB on the side foil", SideFoilTPBThickness_);
    AddParameter("PMTTPBThickness", "thickness of TPB on the PMTs", PMTTPBThickness_);
    AddParameter("TimeCut", "only track events up to 200nsec", TimeCut_);
    AddParameter("KillCherenkov", "turn cherenkov light on/off", KillCherenkov_);
    AddParameter("RandomSeed", "seed for geant4 random generator", RandomSeed_);
}

void CCM200Response::Configure() {
    GetParameter("PMTSDStatus", PMTSDStatus_);
    GetParameter("LArSDStatus", LArSDStatus_);
    GetParameter("SourceRodIn", SourceRodIn_);
    GetParameter("SourceRodLocation", SourceRodLocation_);
    GetParameter("CobaltSourceRun", CobaltSourceRun_);
    GetParameter("SodiumSourceRun", SodiumSourceRun_);
    GetParameter("SingletTimeConstant", SingletTau_);
    GetParameter("TripletTimeConstant", TripletTau_);
    GetParameter("Rayleigh128Length", Rayleigh128_);
    GetParameter("UVAbsLength", UVAbsLength_);
    GetParameter("WLSNPhotonsEndCapFoil", WLSNPhotonsEndCapFoil_);
    GetParameter("WLSNPhotonsSideFoil", WLSNPhotonsSideFoil_);
    GetParameter("WLSNPhotonsPMT", WLSNPhotonsPMT_);
    GetParameter("EndCapFoilTPBThickness", EndCapFoilTPBThickness_);
    GetParameter("SideFoilTPBThickness", SideFoilTPBThickness_);
    GetParameter("PMTTPBThickness", PMTTPBThickness_);
    GetParameter("TimeCut", TimeCut_);
    GetParameter("KillCherenkov", KillCherenkov_);
    GetParameter("RandomSeed", RandomSeed_);
}

CCM200Response::~CCM200Response() {}

void CCM200Response::Initialize() {

    if (g4Interface_ == nullptr){
        g4Interface_ = G4Interface::GetInstance();
    }

    if (!PMTSDStatus_ and !LArSDStatus_){
        log_warn("Oops! Both sensitive detectors are turned off!");
    }

    // let's let's construct the detector
    g4Interface_->InstallDetector(PMTSDStatus_, LArSDStatus_, SourceRodIn_, SourceRodLocation_, CobaltSourceRun_, SodiumSourceRun_, 
                                  SingletTau_, TripletTau_, Rayleigh128_, UVAbsLength_, WLSNPhotonsEndCapFoil_, WLSNPhotonsSideFoil_, WLSNPhotonsPMT_,
                                  EndCapFoilTPBThickness_, SideFoilTPBThickness_, PMTTPBThickness_, TimeCut_, KillCherenkov_, RandomSeed_);
}

I3FrameObjectPtr CCM200Response::GetSimulationConfiguration() {
    DetectorResponseConfigPtr config = boost::make_shared<DetectorResponseConfig>();
    config->rayleigh_scattering_length_ = Rayleigh128_;
    config->pmt_tpb_qe_ = WLSNPhotonsPMT_;
    config->endcap_tpb_qe_ = WLSNPhotonsEndCapFoil_;
    config->side_tpb_qe_ = WLSNPhotonsSideFoil_;
    config->pmt_tpb_thickness_ = PMTTPBThickness_;
    config->endcap_tpb_thickness_ = EndCapFoilTPBThickness_;
    config->side_tpb_thickness_ = SideFoilTPBThickness_;
    return config;
}

void CCM200Response::SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) {
    g4Interface_->SimulateEvent(primary, tree, mcpeseries);
}

void CCM200Response::SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries) {
    g4Interface_->SimulateEvents(primaries, trees, mcpeseries);
}

void CCM200Response::DestroyInterface() {

    g4Interface_->DestroyInstance();
}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);


