#ifndef _CCMSIMULATOR_CCMSIMULATOR_H_
#define _CCMSIMULATOR_CCMSIMULATOR_H_

#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/physics/I3Particle.h"

#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/I3Module.h"

#include "simclasses/CCMMCPE.h"

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

    private:
  
        std::string injectorServiceName_;
        std::string responseServiceName_;
        std::string mcPrimaryName_;
        std::string PMTHitSeriesName_;
        std::string LArMCTreeName_; 
        I3MCTreePtr mcTree_;

        CCMParticleInjectorPtr injector_;
        CCMDetectorResponsePtr response_;
        boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap = boost::make_shared<CCMMCPESeriesMap> ();
        I3MCTreePtr LArEnergyDep = boost::make_shared<I3MCTree>();
        
        bool PMTSDStatus_; // turn PMT SD on/off
        bool LArSDStatus_; // turn fiducial LAr SD on/off
  
        static const std::string INC_ID_NAME;
  
        SET_LOGGER("CCMSimulator");
};

#endif
