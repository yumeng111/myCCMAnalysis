#ifndef _CCMSIMULATOR_CCMSIMULATOR_H_
#define _CCMSIMULATOR_CCMSIMULATOR_H_

#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/physics/I3Particle.h"

#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/I3Module.h"

#include "simclasses/CCMMCPE.h"
#include "simclasses/PhotonSummary.h"

/**
 * \brief The CCMSimulator module handles the whole simulation and
 * stores the results in the frame.
 */

class CCMSimulator : public I3Module  {
public:
    CCMSimulator(const I3Context& context);
    ~CCMSimulator() = default;

    void Configure();

    void Process();
    void ProcessNormally(I3FramePtr frame);

    void DAQSingleThreaded(I3FramePtr frame);
    void DAQMultiThreaded();
    void DAQ(I3FramePtr frame);

    void Simulation(I3FramePtr frame);

    void Finish();
private:
    std::string responseServiceName_;
    std::string input_mc_tree_name_;

    std::string PMTHitSeriesName_;
    std::string LArMCTreeName_;
    std::string PhotonSummarySeriesName_;

    bool multithreaded_;
    int batch_size_;
    std::deque<I3FramePtr> frame_queue_;

    bool seen_s_frame_ = false;
    unsigned int n_daq_frames_ = 0;

    CCMDetectorResponsePtr response_;

    SET_LOGGER("CCMSimulator");
};

#endif
