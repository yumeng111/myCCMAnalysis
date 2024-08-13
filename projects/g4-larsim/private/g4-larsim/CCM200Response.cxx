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
    CCMDetectorResponse(context), PMTSDStatus_(true), LArSDStatus_(true), SourceRodIn_(false), SourceRodLocation_(0.0 * I3Units::cm),
    CobaltSourceRun_(false), SodiumSourceRun_(false), SingletTau_(8.2 * I3Units::nanosecond), TripletTau_(743.0 * I3Units::nanosecond),
    Rayleigh128_(95.0 * I3Units::cm), UVAbsLength_(55.0 * I3Units::cm), TimeCut_(true), KillCherenkov_(false), RandomSeed_(0){
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
                                  SingletTau_, TripletTau_, Rayleigh128_,
                                  UVAbsLength_, TimeCut_, KillCherenkov_, RandomSeed_);
    g4Interface_->InitializeRun();

}

void CCM200Response::BeginEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) {
    // inject our particle
    g4Interface_->InjectParticle(primary, tree, mcpeseries);
}

void CCM200Response::EndEvent() {
    g4Interface_->TerminateEvent(); // this ends event and grabs salient information from geant4
}

void CCM200Response::TerminateRun() {

    g4Interface_->TerminateRun(); 
}

void CCM200Response::DestroyInterface() {

    g4Interface_->DestroyInstance(); 
}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);


