
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
    PhotonSummarySeriesPtr AllPhotonSummarySeries;
    boost::shared_ptr<I3Map<int, size_t>> AllPhotonSummaryMap;
    size_t max_id = 0;

    // Iterate over all particles in the MCTree
    typename I3MCTree::fast_const_iterator tree_iter(*injection_tree), tree_end=injection_tree->cend_fast();
    for(;tree_iter != tree_end; tree_iter++) {
        I3Particle const & particle = *tree_iter;

        // Tell the response service of a new event
        response_->BeginEvent(particle, edep_tree);

        // Tell the response service that the event has ended
        // this will also populate the map between CCMPMTKey and std::vector<CCMMCPE> to save to frame
        // also grab hits in fiducial argon for voxelization
        // note -- if SD not enabled, these will just be empty objects

        // The following are uninitialized because they get passed by reference to the response service
        // and are overwritten there
        boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap;
        PhotonSummarySeriesPtr photon_summary_series;
        boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map;

        response_->EndEvent(CCMMCPEMap, photon_summary_series, photon_summary_series_map);

        // Order of these matters
        // Recursively add energy depositions from LArEnergyDep to tree
        //AppendSubTree(edep_tree, particle.GetID(), LArEnergyDep);
        // Update photon summary map with new parent/track ids and offset indices
        AppendPhotonSummaryMap(AllPhotonSummaryMap, photon_summary_series_map, max_id, AllPhotonSummarySeries->size());
        // Append photon summary series
        AppendPhotonSummarySeries(AllPhotonSummarySeries, particle, photon_summary_series);
        // Append MCPE map, updating time, parent_id, and track_id
        AppendMCPEMap(mcpeseries_map, particle, CCMMCPEMap, max_id);
    }

    // sort mcpeseries_map by time
    for (CCMMCPESeriesMap::iterator it = mcpeseries_map->begin(); it != mcpeseries_map->end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.g4_time < b.g4_time; });
    }

    frame->Put(PMTHitSeriesName_, mcpeseries_map);
    frame->Put(LArMCTreeName_, edep_tree);
    frame->Put(PhotonSummarySeriesName_, AllPhotonSummarySeries);
    frame->Put("PhotonSummaryMap", AllPhotonSummaryMap);

    PushFrame(frame);
}

void CCMSimulator::AppendSubTree(I3MCTreePtr tree, I3ParticleID primary_id, I3MCTreeConstPtr subtree) {
    // Append the subtree to the tree
    std::vector<const I3Particle*> subtree_primaries = I3MCTreeUtils::GetPrimariesPtr(subtree);
    std::vector<const I3Particle*> subtree_daughters;
    for (const I3Particle* primary : subtree_primaries) {
        I3MCTreeUtils::AppendChild(*tree, primary_id, *primary);
        std::vector<const I3Particle*> daughters = I3MCTreeUtils::GetDaughtersPtr(subtree, primary->GetID());
        subtree_daughters.insert(subtree_daughters.end(), daughters.begin(), daughters.end());
    }
    subtree_primaries = subtree_daughters;
    subtree_daughters.clear();
    while(subtree_primaries.size() > 0) {
        for (const I3Particle* primary : subtree_primaries) {
            I3Particle const * parent = I3MCTreeUtils::GetParentPtr(subtree, primary->GetID());
            I3MCTreeUtils::AppendChild(*tree, parent->GetID(), *primary);
            std::vector<const I3Particle*> daughters = I3MCTreeUtils::GetDaughtersPtr(subtree, primary->GetID());
            subtree_daughters.insert(subtree_daughters.end(), daughters.begin(), daughters.end());
        }
        subtree_primaries = subtree_daughters;
        subtree_daughters.clear();
    }
}

void CCMSimulator::AppendPhotonSummaryMap(boost::shared_ptr<I3Map<int, size_t>> parent_map, boost::shared_ptr<I3Map<int, size_t>> sub_map, size_t max_id, size_t offset) {
    for(I3Map<int, size_t>::const_iterator it = sub_map->begin(); it != sub_map->end(); ++it) {
        parent_map->insert(std::make_pair(it->first + max_id, it->second + offset));
    }
}

void CCMSimulator::AppendPhotonSummarySeries(PhotonSummarySeriesPtr parent_series, I3Particle const & parent_particle, PhotonSummarySeriesConstPtr sub_series) {
    double parent_time = parent_particle.GetTime();

    parent_series->reserve(parent_series->size() + sub_series->size());
    for(PhotonSummary const & ps : *sub_series) {
        parent_series->push_back(ps);
        parent_series->back().g4_time += parent_time;
    }
}

void CCMSimulator::AppendMCPEMap(CCMMCPESeriesMapPtr parent_map, I3Particle const & parent_particle, CCMMCPESeriesMapConstPtr sub_map, size_t & max_id) {

    size_t sub_map_max_id = 0;
    for(CCMMCPESeriesMap::const_iterator it = sub_map->begin(); it != sub_map->end(); ++it) {
        for(CCMMCPE const & mpe : it->second) {
            sub_map_max_id = std::max(sub_map_max_id, mpe.parent_id);
            sub_map_max_id = std::max(sub_map_max_id, mpe.track_id);
        }
    }

    double parent_time = parent_particle.GetTime();
    for (CCMMCPESeriesMap::const_iterator it = sub_map->begin(); it != sub_map->end(); ++it) {
        CCMMCPESeriesMap::iterator parent_it = parent_map->find(it->first);
        if (parent_it == parent_map->end()) {
            parent_map->insert(*it);
        } else {
            parent_it->second.reserve(parent_it->second.size() + it->second.size());
            for(CCMMCPE const & mpe : it->second) {
                parent_it->second.push_back(mpe);
                CCMMCPE & p = parent_it->second.back();
                p.g4_time += parent_time;
                p.parent_id += max_id;
                p.track_id += max_id;
            }
        }
    }

    max_id += sub_map_max_id + 1;
}

void CCMSimulator::Finish() {
    // terminate geant4
    response_->TerminateRun();
    // destruct g4 interface
    response_->DestroyInterface();
}



