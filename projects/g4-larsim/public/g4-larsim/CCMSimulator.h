#ifndef _CCMSIMULATOR_CCMSIMULATOR_H_
#define _CCMSIMULATOR_CCMSIMULATOR_H_

#include "dataclasses/physics/I3Particle.h"

#include "g4-larsim/CCMDetectorResponse.h"

#include "icetray/I3Module.h"

#include "simclasses/CCMMCPE.h"

/**
 * \brief The CCMSimulator module handles the whole simulation and
 * stores the results in the frame.
 */

class CCMSimulator : public I3Module  {
public:
    CCMSimulator(const I3Context& context);
    ~CCMSimulator() = default;

    void Configure() override;

    void Process() override;
    void ProcessNormally(I3FramePtr frame);

    void DAQSingleThreaded(I3FramePtr frame);
    void DAQMultiThreaded();
    void DAQ(I3FramePtr frame) override;

    void Simulation(I3FramePtr frame) override;
    void FillSimulationFrame(I3FramePtr frame);

    void Finish() override;

    static void MergeMCPESeries(CCMMCPESeriesMap & mcpeseries_dest, CCMMCPESeriesMap const & mcpeseries_source);
    static void MergeEDepTrees(I3MCTreePtr dest, I3MCTreePtr source, I3Particle primary);
private:
    std::string responseServiceName_ = "CCM200Response";
    std::string input_mc_tree_name_ = "I3MCTree";
    std::string configuration_name_ = "DetectorResponseConfig";
    std::string PMTHitSeriesName_ = "PMTMCHitsMap";
    std::string LArMCTreeName_ = "LArMCTree";

    bool multithreaded_ = false;
    size_t batch_size_ = 1;
    size_t n_threads_ = 1;
    bool multi_particle_output_ = false;
    bool per_event_output_ = true;

    std::deque<I3FramePtr> frame_queue_;
    bool seen_s_frame_ = false;
    unsigned int n_daq_frames_ = 0;

    CCMDetectorResponsePtr response_;

    SET_LOGGER("CCMSimulator");
};

#endif
