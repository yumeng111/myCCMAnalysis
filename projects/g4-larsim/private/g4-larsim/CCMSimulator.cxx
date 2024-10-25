
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
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"

#include <stdexcept>

#define i3_log(format, ...) log_trace("%s: " format, this->GetName().c_str(), ##__VA_ARGS__)

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

    multithreaded_ = false;
    AddParameter("Multithreaded", "Run the simulation in multithreaded mode.", multithreaded_);

    batch_size_ = 1;
    AddParameter("BatchSize", "Number of events to simulate in one batch.", batch_size_);

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

    GetParameter("Multithreaded", multithreaded_);
    log_info("+ Multithreaded : %s", multithreaded_ ? "true" : "false");

    GetParameter("BatchSize", batch_size_);
    log_info("+ BatchSize : %zu", batch_size_);

    // initialize the response services
    response_->Initialize();
}


void CCMSimulator::Process() {

    i3_log("%s", __PRETTY_FUNCTION__);

    if (inbox_)
        log_trace("%zu frames in inbox", inbox_->size());

    I3FramePtr frame = PopFrame();

    if (!frame) {
        return;
    }

    methods_t::iterator miter = methods_.find(frame->GetStop());

    if (miter != methods_.end()) {
        miter->second(frame);
        return;
    }

    if(not multithreaded_) {
        ProcessNormally(frame);
        return;
    }

    if(frame->GetStop() == I3Frame::DAQ) {
        frame_queue_.push_back(frame);
        ++n_daq_frames_;
        if(n_daq_frames_ == batch_size_) {
            ndaqcall_ += n_daq_frames_;
            n_daq_frames_ = 0;
            // Reached the batch size, process all frames in the queue
            DAQMultiThreaded();
        }
    } else {
        if(frame->GetStop() == I3Frame::Geometry or frame->GetStop() == I3Frame::Calibration or frame->GetStop() == I3Frame::DetectorStatus or frame->GetStop() == I3Frame::Simulation) {
            if(frame_queue_.size() == 0) {
                ProcessNormally(frame);
            } else {
                // Got a frame that should be processed before the DAQ frame
                // Process all frames in the queue
                DAQMultiThreaded();
                ProcessNormally(frame);
            }
        } else {
            frame_queue_.push_back(frame);
        }
    }
}

void CCMSimulator::ProcessNormally(I3FramePtr frame) {
    if(frame->GetStop() == I3Frame::Physics && ShouldDoPhysics(frame)) {
        ++nphyscall_;
        Physics(frame);
    } else if(frame->GetStop() == I3Frame::Geometry && ShouldDoGeometry(frame))
        Geometry(frame);
    else if(frame->GetStop() == I3Frame::Calibration && ShouldDoCalibration(frame))
        Calibration(frame);
    else if(frame->GetStop() == I3Frame::DetectorStatus && ShouldDoDetectorStatus(frame))
        DetectorStatus(frame);
    else if(frame->GetStop() == I3Frame::Simulation && ShouldDoSimulation(frame))
        Simulation(frame);
    else if(frame->GetStop() == I3Frame::DAQ && ShouldDoDAQ(frame)) {
        ++ndaqcall_;
        DAQ(frame);
    } else if(ShouldDoOtherStops(frame))
        OtherStops(frame);
}

void CCMSimulator::FillSimulationFrame(I3FramePtr frame) {
    I3FrameObjectPtr obj = response_->GetSimulationConfiguration();
    frame->Put("DetectorConfiguration", obj);
}

void CCMSimulator::Simulation(I3FramePtr frame) {
    seen_s_frame_ = true;
    FillSimulationFrame(frame);
    PushFrame(frame);
}

void CCMSimulator::DAQ(I3FramePtr frame) {
    if(multithreaded_) {
        // Process all frames in the queue
        // Most recent DAQ frame is already in the queue
        DAQMultiThreaded();
    } else {
        DAQSingleThreaded(frame);
    }
}

void CCMSimulator::DAQSingleThreaded(I3FramePtr frame) {
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
        // Simulate the event with the response
        response_->SimulateEvent(particle, edep_tree, mcpeseries_map);
    }

    // sort mcpeseries_map by time
    for (CCMMCPESeriesMap::iterator it = mcpeseries_map->begin(); it != mcpeseries_map->end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.g4_time < b.g4_time; });
    }

    frame->Put(PMTHitSeriesName_, mcpeseries_map);
    frame->Put(LArMCTreeName_, edep_tree);

    PushFrame(frame);
}

void CCMSimulator::DAQMultiThreaded() {
    // Pop everything until the first DAQ frame, processing / passing other frame types as appropriate
    while(frame_queue_.size() > 0 && frame_queue_.front()->GetStop() != I3Frame::DAQ) {
        ProcessNormally(frame_queue_.front());
        frame_queue_.pop_front();
    }

    size_t max_queue_index = 0;
    std::deque<I3FramePtr> daq_frames;
    // Get all the DAQ frames until the next frame that supersedes them
    for(size_t i=0; i<frame_queue_.size(); ++i) {
        I3FramePtr frame = frame_queue_.at(i);
        if(frame->GetStop() == I3Frame::Geometry or frame->GetStop() == I3Frame::Calibration or frame->GetStop() == I3Frame::DetectorStatus or frame->GetStop() == I3Frame::Simulation) {
            break;
        } else if(frame->GetStop() == I3Frame::DAQ) {
            daq_frames.push_back(frame);
            max_queue_index = i;
        }
    }

    if(!seen_s_frame_) {
        log_fatal("No simulation frame seen before DAQ frame");
    }

    std::vector<size_t> particles_per_event;
    std::vector<I3Particle> particles;
    std::vector<I3MCTreePtr> edep_trees;
    std::vector<CCMMCPESeriesMapPtr> mcpeseries_maps;

    // Create input/ouput objects for each event
    for(size_t i=0; i<daq_frames.size(); ++i) {
        I3FramePtr frame = daq_frames.at(i);
        // let's grab the mcPrimary from the injector
        I3MCTreeConstPtr injection_tree = frame->Get<I3MCTreeConstPtr>(input_mc_tree_name_);
        if(!injection_tree) {
            log_fatal("No MCTree found in frame with name %s", input_mc_tree_name_.c_str());
        }

        size_t n_particles = 0;

        typename I3MCTree::fast_const_iterator tree_iter(*injection_tree), tree_end=injection_tree->cend_fast();
        for(;tree_iter != tree_end; ++tree_iter) {
            ++n_particles;
            particles.push_back(*tree_iter);
            edep_trees.push_back(boost::make_shared<I3MCTree>(*injection_tree));
            mcpeseries_maps.push_back(boost::make_shared<CCMMCPESeriesMap>());
        }
        particles_per_event.push_back(n_particles);
    }

    log_debug("Simulating CCM");

    response_->SimulateEvents(particles, edep_trees, mcpeseries_maps);

    std::vector<I3MCTreePtr> final_edep_trees;
    std::vector<CCMMCPESeriesMapPtr> final_mcpeseries_maps;

    // Now need to merge everything into individual data structures for each event
    size_t particle_idx = 0;
    for(size_t i=0; i<daq_frames.size(); ++i) {
        I3MCTreePtr edep_tree = edep_trees.at(particle_idx);
        CCMMCPESeriesMapPtr mcpeseries_map = mcpeseries_maps.at(particle_idx);

        final_edep_trees.push_back(edep_tree);
        final_mcpeseries_maps.push_back(mcpeseries_map);
        for(size_t j=1; j<particles_per_event.at(i); ++j) {
            MergeEDepTrees(edep_tree, edep_trees.at(particle_idx+j), particles.at(particle_idx+j));
            MergeMCPESeries(mcpeseries_map, mcpeseries_maps.at(particle_idx+j));
        }

        // sort mcpeseries_map by time
        for (CCMMCPESeriesMap::iterator it = mcpeseries_map->begin(); it != mcpeseries_map->end(); ++it) {
            std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.g4_time < b.g4_time; });
        }

        particle_idx += particles_per_event.at(i);
    }

    // Put it all back into the DAQ frames
    for(size_t i=0; i<daq_frames.size(); ++i) {
        I3FramePtr frame = daq_frames.at(i);
        frame->Put(PMTHitSeriesName_, final_mcpeseries_maps.at(i));
        frame->Put(LArMCTreeName_, final_edep_trees.at(i));
    }

    size_t original_frame_queue_size = frame_queue_.size();

    size_t daq_frames_pushed = 0;
    // Push the frames
    for(size_t i=0; i<original_frame_queue_size and daq_frames_pushed<daq_frames.size(); ++i) {
        I3FramePtr frame = frame_queue_.front();
        if(frame->GetStop() == I3Frame::Geometry or frame->GetStop() == I3Frame::Calibration or frame->GetStop() == I3Frame::DetectorStatus or frame->GetStop() == I3Frame::Simulation) {
            break;
        }
        frame_queue_.pop_front();
        PushFrame(frame);
        if(frame->GetStop() == I3Frame::DAQ) {
            ++daq_frames_pushed;
        }
    }

    // Pop everything until the first DAQ frame, processing / passing other frame types as appropriate
    while(frame_queue_.size() > 0 && frame_queue_.front()->GetStop() != I3Frame::DAQ) {
        ProcessNormally(frame_queue_.front());
        frame_queue_.pop_front();
    }
}

void CCMSimulator::Finish() {
    if(frame_queue_.size() > 0) {
        // Process all frames in the queue
        DAQMultiThreaded();
    }
    // destruct g4 interface
    response_->DestroyInterface();
}

void CCMSimulator::MergeMCPESeries(CCMMCPESeriesMapPtr mcpeseries_dest, CCMMCPESeriesMapPtr mcpeseries_source) {
    // Iterate over PMTs in source map
    for (CCMMCPESeriesMap::iterator it = mcpeseries_source->begin(); it != mcpeseries_source->end(); ++it) {
        // Find the corresponding PMT in the destination map
        CCMMCPESeriesMap::iterator it_dest = mcpeseries_dest->find(it->first);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == mcpeseries_dest->end()) {
            mcpeseries_dest->insert(std::make_pair(it->first, CCMMCPESeries()));
            // Update the iterator so it points to our new entry
            it_dest = mcpeseries_dest->find(it->first);
        }

        // Reference to the destination
        CCMMCPESeries & dest_series = it_dest->second;
        dest_series.insert(dest_series.end(), it->second.begin(), it->second.end());
    }
}

void CCMSimulator::MergeEDepTrees(I3MCTreePtr dest, I3MCTreePtr source, I3Particle primary) {
    I3Particle * source_particle = I3MCTreeUtils::GetParticlePtr(source, primary.GetID());
    I3Particle * dest_particle = I3MCTreeUtils::GetParticlePtr(dest, primary.GetID());

    if(source_particle == NULL) {
        log_fatal("Source particle not found in source tree");
    }
    if(dest_particle == NULL) {
        log_fatal("Source particle not found in destination tree");
    }

    std::vector<I3Particle *> daughters = I3MCTreeUtils::GetDaughtersPtr(source, source_particle->GetID());
    std::deque<std::tuple<I3Particle *, I3Particle *>> source_children(daughters.size());
    for(size_t i=0; i<daughters.size(); ++i) {
        source_children[i] = std::make_tuple(source_particle, daughters.at(i));
    }

    while(source_children.size() > 0) {
        I3Particle * source_parent = std::get<0>(source_children.front());
        I3Particle * source_child = std::get<1>(source_children.front());
        source_children.pop_front();

        I3MCTreeUtils::AppendChild(*dest, source_parent->GetID(), *source_child);
        daughters = I3MCTreeUtils::GetDaughtersPtr(source, source_child->GetID());

        for(size_t i=0; i<daughters.size(); ++i) {
            source_children.push_back(std::make_tuple(source_child, daughters.at(i)));
        }
    }
}
