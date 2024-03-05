// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/CCM200Response.h"
#include <g4-larsim/g4classes/G4Interface.h>
#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <simclasses/CCMMCPE.h>
#include <icetray/I3SingleServiceFactory.h>


CCM200Response::CCM200Response(const I3Context& context) :
    CCMDetectorResponse(context), geometry_name_("CCMGeometry"), mc_tree_name_("I3MCTree"), output_reco_pulse_name_("SimulationRecoPulses") {
    AddParameter("InputMCTreeName", "Name of the MCTree in the frame.", mc_tree_name_);
    AddParameter("OutputRecoPulseMapName", "Name of the CCMRecoPulse series to save to the frame.", output_reco_pulse_name_);
}

void CCM200Response::Configure() {
    GetParameter("InputMCTreeName", mc_tree_name_);
    GetParameter("OutputRecoPulseMapName", output_reco_pulse_name_);
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

    // let's let's construct the detector
    g4Interface_->InstallDetector();

}

void CCM200Response::BeginEvent(const I3Particle& primary)
{
    g4Interface_->InitializeEvent();
    
    // inject our particle
    g4Interface_->InjectParticle(particle);
}


void CCM200Response::EndEvent() {
    // let's grab the photon hits from our event and then save

    g4Interface_->TerminateEvent();

    std::map<OMKey, double>::const_iterator pePerVEM_iter;
    for (pePerVEM_iter = pePerVEM_.begin(); pePerVEM_iter != pePerVEM_.end(); ++pePerVEM_iter) {
        hitHC.GetHitHisto(pePerVEM_iter->first).Scale(scalingFactor);
        cherHitCollection.GetHitHisto(pePerVEM_iter->first).Scale(scalingFactor);
    }

    log_trace_stream(tankKey_ << " reached VEM-threshold, particles before: "
                              << particlesBeforeThreshold_ << ", after: "
                              << particlesAfterThreshold_ << ", scaling-factor: "
                              << scalingFactor);
    }
}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);
