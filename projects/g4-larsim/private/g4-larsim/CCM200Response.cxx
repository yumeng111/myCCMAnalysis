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
#include "icetray/I3ServiceBase.h"
#include "icetray/I3SingleServiceFactory.h"

#include "phys-services/I3RandomService.h"
#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>
#include <algorithm>


CCM200Response::CCM200Response(const I3Context& context) :
    CCMDetectorResponse(context), PMTSDStatus_(true), LArSDStatus_(true), SodiumSourceRun_(false), SodiumSourceLocation_(0.0 * I3Units::cm),
    SingletTau_(8.2 * I3Units::nanosecond), TripletTau_(743.0 * I3Units::nanosecond), UVAbsStatus_(true), TimeCut_(true), CerenkovControl_(false){
    AddParameter("PMTSDStatus", "true if tracking photon hits on PMTs", PMTSDStatus_);
    AddParameter("LArSDStatus", "true if tracking scintillation depositions in fiducial LAr", LArSDStatus_);
    AddParameter("SodiumSourceRun", "true if we want to simulate the sodium source rod + pellet", SodiumSourceRun_);
    AddParameter("SodiumSourceLocation", "z location of the end of the sodium source rod", SodiumSourceLocation_);
    AddParameter("SingletTimeConstant", "LAr singlet tau", SingletTau_);
    AddParameter("TripletTimeConstant", "LAr triplet tau", TripletTau_);
    AddParameter("UVAbsLenStatus", "turn uv abs on/off", UVAbsStatus_);
    AddParameter("TimeCut", "only track events up to 200nsec", TimeCut_);
    AddParameter("CerenkovControl", "turn cerenkov light on/off", CerenkovControl_);
}

void CCM200Response::Configure() {
    GetParameter("PMTSDStatus", PMTSDStatus_);
    GetParameter("LArSDStatus", LArSDStatus_);
    GetParameter("SodiumSourceRun", SodiumSourceRun_);
    GetParameter("SodiumSourceLocation", SodiumSourceLocation_);
    GetParameter("SingletTimeConstant", SingletTau_);
    GetParameter("TripletTimeConstant", TripletTau_);
    GetParameter("UVAbsLenStatus", UVAbsStatus_);
    GetParameter("TimeCut", TimeCut_);
    GetParameter("CerenkovControl", CerenkovControl_);
}

CCM200Response::~CCM200Response() {
    if (G4Interface::GetInstance()) {
        delete g4Interface_;
    }
}

void CCM200Response::Initialize() {

    g4Interface_ = G4Interface::GetInstance();
    if (!g4Interface_) {
        g4Interface_ = new G4Interface(visMacroFile_);
    }

    if (!PMTSDStatus_ and !LArSDStatus_){
        log_warn("Oops! Both sensitive detectors are turned off!");
    }

    // let's let's construct the detector
    g4Interface_->InstallDetector(PMTSDStatus_, LArSDStatus_, SodiumSourceRun_, SodiumSourceLocation_, SingletTau_, TripletTau_, UVAbsStatus_,
                                  TimeCut_, CerenkovControl_);
    g4Interface_->InitializeRun();

}

void CCM200Response::BeginEvent(const I3Particle& primary) {
    // inject our particle
    g4Interface_->InjectParticle(primary);

}

void CCM200Response::EndEvent(I3MCTreePtr & LArEnergyDep, boost::shared_ptr<CCMMCPESeriesMap> & CCMMCPEMap, PhotonSummarySeriesPtr & photon_summary_series,
                              boost::shared_ptr<I3Map<int, size_t>> & photon_summary_series_map ) {

    g4Interface_->TerminateEvent(); // this ends event and grabs salient information from geant4

    // now let's call function in G4Interface to grab the map between CCMPMTKey and std::vector<CCMMCPE>
    LArEnergyDep = g4Interface_->GetLArEnergyDep(); 
    CCMMCPEMap = g4Interface_->GetCCMMCPEMap();
    photon_summary_series = g4Interface_->GetPhotonSummarySeries();
    photon_summary_series_map = g4Interface_->GetPhotonSummaryMap();

}

void CCM200Response::TerminateRun() {

    g4Interface_->TerminateRun(); 
}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);


