
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

    // initialize the response services
    response_->Initialize();
}

void CCMSimulator::Process() {
    i3_log("%s", __PRETTY_FUNCTION__);

    if (inbox_)
        log_trace("%zu frames in inbox", inbox_->size());

    I3FramePtr frame = PopFrame();

    if (!frame)
        return;

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
        ++n_daq_frames_;
        if(n_daq_frames_ == batch_size_) {
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
        else {
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

void CCMSimulator::Simulation(I3FramePtr frame) {
    seen_s_frame_ = true;
    // Do nothing with the simulation frame so far
    // Could later be used to load systematic information
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

void CCMSimulator::DAQMultiThreaded(I3FramePtr frame) {
    // Pop everything until the first DAQ frame, processing / passing other frame types as appropriate
    while(frame_queue_.size() > 0 && frame_queue_.front()->GetStop() != I3Frame::DAQ) {
        ProcessNormally(frame_queue_.front());
        frame_queue_.pop_front();
    }

    std::deque<I3FramePtr> daq_frames;
    // Get all the DAQ frames until the next frame that supersedes them
    for(size_t i=0; i<frame_queue_.size(); ++i) {
        I3FramePtr frame = frame_queue_.at(i);
        if(frame->GetStop() == I3Frame::Geometry or frame->GetStop() == I3Frame::Calibration or frame->GetStop() == I3Frame::DetectorStatus or frame->GetStop() == I3Frame::Simulation) {
            break;
        } else if(frame->GetStop() == I3Frame::DAQ) {
            daq_frames.push_back(frame);
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

    // TODO
    // Now need to merge everything into individual data structures for each event
    size_t particle_idx = 0;
    for(size_t i=0; i<daq_frames.size(); ++i) {
        CCMMCPESeriesMapPtr> mcpeseries_map = mcpeseries_maps.at(particle_idx);
        for(size_t j=1; i<particles_per_event.at(i); ++j) {

            mcpeseries_map->insert(mcpeseries_map->end(), mcpeseries_maps.at(particle_idx+i)->begin(), mcpeseries_maps.at(particle_idx+i)->end());
        }

        ++particle_idx;
    }

    // Put it all back into the DAQ frames
    for(size_t i=0; i<daq_frames.size(); ++i) {
        I3FramePtr frame = daq_frames.at(i);
        frame->Put(PMTHitSeriesName_, final_mcpeseries_maps.at(i));
        frame->Put(LArMCTreeName_, final_edep_trees.at(i));
    }

    // Push the frames
    for(size_t i=0; i<frame_queue_.size(); ++i) {
        I3FramePtr frame = frame_queue_.at(i);
        if(frame->GetStop() == I3Frame::Geometry or frame->GetStop() == I3Frame::Calibration or frame->GetStop() == I3Frame::DetectorStatus or frame->GetStop() == I3Frame::Simulation) {
            break;
        }
        frame_queue_.pop_front();
        PushFrame(frame);
    }

    // Pop everything until the first DAQ frame, processing / passing other frame types as appropriate
    while(frame_queue_.size() > 0 && frame_queue_.front()->GetStop() != I3Frame::DAQ) {
        ProcessNormally(frame_queue_.front());
        frame_queue_.pop_front();
    }
}

void CCMSimulator::Finish() {
    // terminate geant4
    response_->TerminateRun();
    // destruct g4 interface
    response_->DestroyInterface();
}

