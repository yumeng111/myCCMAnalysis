
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

CCMSimulator::CCMSimulator(const I3Context& context): I3Module(context) {
    response_ = CCMDetectorResponsePtr();

    responseServiceName_ = "CCM200Response";
    AddParameter("ResponseServiceName", "Name of the detector response service.", responseServiceName_);

    input_mc_tree_name_ = "I3MCTree";
    AddParameter("InputMCTreeName", "Name of the input MC tree in the frame.", input_mc_tree_name_);

    PMTHitSeriesName_ = "PMTMCHitsMap";
    AddParameter("PMTHitSeriesName", "Name of the resulting PMT hit map in the frame.", PMTHitSeriesName_);

    LArMCTreeName_ = "LArMCTree";
    AddParameter("LArMCTreeName", "Name of the MC tree containing energy depositions in LAr", LArMCTreeName_);

    PhotonSummarySeriesName_ = "PhotonSummarySeries";
    AddParameter("PhotonSummarySeriesName", "Name of the photon summary series containing optical photon hits in LAr", PhotonSummarySeriesName_);

    AddOutBox("OutBox");
}

void CCMSimulator::Configure() {
    log_info("Configuring the CCMSimulator:");

    GetParameter("ResponseServiceName", responseServiceName_);
    log_info("+ Response Service: %s", responseServiceName_.c_str());
    response_ = GetContext().Get<CCMDetectorResponsePtr>(responseServiceName_);
    if(!response_) log_fatal("No response service \"%s\" in context", responseServiceName_.c_str());

    GetParameter("InputMCTreeName", input_mc_tree_name_);
    log_info("+ Input MC Tree : %s", input_mc_tree_name_.c_str());

    GetParameter("PMTHitSeriesName", PMTHitSeriesName_);
    log_info("+ PMT hit series : %s", PMTHitSeriesName_.c_str());

    GetParameter("LArMCTreeName", LArMCTreeName_);
    log_info("+ LAr MC Tree : %s", LArMCTreeName_.c_str());

    GetParameter("PhotonSummarySeriesName", PhotonSummarySeriesName_);
    log_info("+ PhotonSummarySeries : %s", PhotonSummarySeriesName_.c_str());

    // initialize the response services
    response_->Initialize();
}

void CCMSimulator::Simulation(I3FramePtr frame) {
    seen_s_frame_ = true;
    // Do nothing with the simulation frame so far
    // Could later be used to load systematic information
    PushFrame(frame);
}

void CCMSimulator::DAQ(I3FramePtr frame) {
    if(!seen_s_frame_) {
        log_fatal("No simulation frame seen before DAQ frame");
    }

    // let's grab the mcPrimary from the injector
    I3MCTreeConstPtr injection_tree = frame->Get<I3MCTreeConstPtr>(input_mc_tree_name_);
    if(!injection_tree) {
        log_fatal("No MCTree found in frame with name %s", input_mc_tree_name_.c_str());
    }

    I3MCTreePtr edep_tree = boost::make_shared<I3MCTree>(*injection_tree);
    CCMMCPESeriesMapPtr mcpeseries_map = boost::make_shared<CCMMCPESeriesMap>();


    log_debug("Simulating CCM");
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

