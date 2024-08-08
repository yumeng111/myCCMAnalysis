
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3EventHeader.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/physics/I3ScintRecoPulseSeriesMap.h"

#include "g4-larsim/CCMSimulator.h"

#include "icetray/I3Bool.h"
#include "icetray/I3Frame.h"
#include "icetray/I3Context.h"

#include <stdexcept>

I3_MODULE(CCMSimulator);

// Name of the IncEventID boolean in the frame (see I3MCTimeGeneratorService)
const std::string CCMSimulator::INC_ID_NAME = "MCTimeIncEventID";

CCMSimulator::CCMSimulator(const I3Context& context): I3Module(context) {
    injector_ = CCMParticleInjectorPtr();
    response_ = CCMDetectorResponsePtr();

    responseServiceName_ = "CCM200Response";
    AddParameter("ResponseServiceName", "Name of the detector response service.", responseServiceName_);

    injectorServiceName_ = "CCMSimpleInjector";
    AddParameter("InjectorServiceName", "Name of the injector service.", injectorServiceName_);

    mcPrimaryName_ = "MCPrimary";
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);

    PMTHitSeriesName_ = "PMTMCHitsMap";
    AddParameter("PMTHitSeriesName", "Name of the resulting PMT hit map in the frame.", PMTHitSeriesName_);

    LArMCTreeName_ = "LArMCTree";
    AddParameter("LArMCTreeName", "Name of the MC tree containing energy depositions in LAr", LArMCTreeName_);

    PhotonSummarySeriesName_ = "PhotonSummarySeries";
    AddParameter("PhotonSummarySeriesName", "Name of the photon summary series containing optical photon hits in LAr", PhotonSummarySeriesName_);

    AddOutBox("OutBox");
}

CCMSimulator::~CCMSimulator() { }

void CCMSimulator::Configure() {
    log_info("Configuring the CCMSimulator:");

    GetParameter("InjectorServiceName", injectorServiceName_);
    log_info("+ Injector Service: %s", injectorServiceName_.c_str());
    injector_ = GetContext().Get<CCMParticleInjectorPtr>(injectorServiceName_);
    if(!injector_) log_fatal("No injector service \"%s\" in context", injectorServiceName_.c_str());

    GetParameter("ResponseServiceName", responseServiceName_);
    log_info("+ Response Service: %s", responseServiceName_.c_str());
    response_ = GetContext().Get<CCMDetectorResponsePtr>(responseServiceName_);
    if(!response_) log_fatal("No response service \"%s\" in context", responseServiceName_.c_str());

    GetParameter("PrimaryName", mcPrimaryName_);
    log_info("+ MC primary: %s", mcPrimaryName_.c_str());

    GetParameter("PMTHitSeriesName", PMTHitSeriesName_);
    log_info("+ PMT hit series : %s", PMTHitSeriesName_.c_str());

    GetParameter("LArMCTreeName", LArMCTreeName_);
    log_info("+ LAr MC Tree : %s", LArMCTreeName_.c_str());

    GetParameter("PhotonSummarySeriesName", PhotonSummarySeriesName_);
    log_info("+ PhotonSummarySeries : %s", PhotonSummarySeriesName_.c_str());

    // initialize injector and response services
    //injector_->Configure();
    //response_->Configure();
    response_->Initialize();

    // set our controls over SD
    PMTSDStatus_ = response_->GetPMTSDStatus();
    LArSDStatus_ = response_->GetLArSDStatus();
}


// we don't have DetectorStatus frames at this moment...
void CCMSimulator::DetectorStatus(I3FramePtr frame) {
    // Configure and initialize response and injector services
    //response_->Configure();
    //response_->Initialize();

    //injector_->SetResponseService(response_);
    //injector_->Configure();

    PushFrame(frame);
}

void CCMSimulator::Simulation(I3FramePtr frame) {
    seen_s_frame_ = true;
    FillSimulationFrame(frame);
    PushFrame(frame);
}

void CCMSimulator::FillSimulationFrame(I3FramePtr frame) {
    // put mcTree into frame
    I3FrameObjectPtr obj = injector_->GetSimulationConfiguration();
    frame->Put("InjectorConfiguration", obj);
}

void CCMSimulator::DAQ(I3FramePtr frame) {
    if(not seen_s_frame_) {
        I3FramePtr sim_frame = boost::make_shared<I3Frame>(I3Frame::Simulation);
        FillSimulationFrame(sim_frame);
        PushFrame(sim_frame);
        seen_s_frame_ = true;
    }

    log_debug("   Simulating CCM");

    // let's grab the mcPrimary from the injector
    I3MCTreePtr injection_tree = injector_->GetMCTree();
    frame->Put("I3MCTree", injection_tree);

    I3MCTreePtr edep_tree = boost::make_shared<I3MCTree>(*injection_tree);
    CCMMCPESeriesMapPtr mcpeseries_map = boost::make_shared<CCMMCPESeriesMap>();

    // Iterate over all particles in the MCTree
    typename I3MCTree::fast_const_iterator tree_iter(*injection_tree), tree_end=injection_tree->cend_fast();
    for(;tree_iter != tree_end; tree_iter++) {
        I3Particle const & particle = *tree_iter;

        // Tell the response service of a new event
        response_->BeginEvent(particle, edep_tree, mcpeseries_map);

        // Tell the response service that the event has ended
        // this will also populate the map between CCMPMTKey and std::vector<CCMMCPE> to save to frame
        // also grab hits in fiducial argon for voxelization
        // note -- if SD not enabled, these will just be empty objects
        response_->EndEvent();
    }

    // sort mcpeseries_map by time
    for (CCMMCPESeriesMap::iterator it = mcpeseries_map->begin(); it != mcpeseries_map->end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.g4_time < b.g4_time; });
    }

    frame->Put(PMTHitSeriesName_, mcpeseries_map);
    frame->Put(LArMCTreeName_, edep_tree);

    PushFrame(frame);
}

void CCMSimulator::Finish() {
    // terminate geant4
    response_->TerminateRun();
    // destruct g4 interface
    response_->DestroyInterface();
}

