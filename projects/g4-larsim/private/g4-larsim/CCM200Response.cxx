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
    CCMDetectorResponse(context), PMTSDStatus_(true), LArSDStatus_(true)   {
    AddParameter("PMTSDStatus", "true if tracking photon hits on PMTs", PMTSDStatus_);
    AddParameter("LArSDStatus", "true if tracking scintillation depositions in fiducial LAr", LArSDStatus_);
}

void CCM200Response::Configure() {
    GetParameter("PMTSDStatus", PMTSDStatus_);
    GetParameter("LArSDStatus", LArSDStatus_);
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
    g4Interface_->InstallDetector(PMTSDStatus_, LArSDStatus_);

}

void CCM200Response::BeginEvent(const I3Particle& primary) {
    g4Interface_->InitializeEvent();

    // inject our particle
    g4Interface_->InjectParticle(primary);

}

void CCM200Response::EndEvent() {

    g4Interface_->TerminateEvent(); // this ends event and grabs salient information from geant4

    // now let's call function in G4Interface to grab the map between CCMPMTKey and std::vector<CCMMCPE>
    CCMMCPEMap = g4Interface_->GetCCMMCPEMap();  
    LArEnergyDep = g4Interface_->GetLArEnergyDep();

}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);


