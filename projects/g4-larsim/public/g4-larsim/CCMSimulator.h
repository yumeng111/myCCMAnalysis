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
        ~CCMSimulator();

        void Configure();
        void DetectorStatus(I3FramePtr frame);
        void DAQ(I3FramePtr frame);
        void Simulation(I3FramePtr frame);
        void FillSimulationFrame(I3FramePtr frame);
        void Finish();
    private:
  
        std::string injectorServiceName_;
        std::string responseServiceName_;
        std::string mcPrimaryName_;
        std::string PMTHitSeriesName_;
        std::string LArMCTreeName_; 
        std::string PhotonSummarySeriesName_; 

        bool seen_s_frame_ = false;

        CCMParticleInjectorPtr injector_;
        CCMDetectorResponsePtr response_;

        static void AppendSubTree(I3MCTreePtr tree, I3ParticleID primary_id, I3MCTreeConstPtr subtree);
        static void AppendPhotonSummaryMap(boost::shared_ptr<I3Map<int, size_t>> parent_map, boost::shared_ptr<I3Map<int, size_t>> sub_map, size_t max_id, size_t offset);
        static void AppendPhotonSummarySeries(PhotonSummarySeriesPtr parent_series, I3Particle const & parent_particle, PhotonSummarySeriesConstPtr sub_series);
        static void AppendMCPEMap(CCMMCPESeriesMapPtr parent_map, I3Particle const & parent_particle, CCMMCPESeriesMapConstPtr sub_map, size_t & max_id);

        bool PMTSDStatus_; // turn PMT SD on/off
        bool LArSDStatus_; // turn fiducial LAr SD on/off
  
        static const std::string INC_ID_NAME;
  
        SET_LOGGER("CCMSimulator");
};

#endif
