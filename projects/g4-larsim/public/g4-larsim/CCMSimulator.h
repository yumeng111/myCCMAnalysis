#ifndef _CCMSIMULATOR_CCMSIMULATOR_H_
#define _CCMSIMULATOR_CCMSIMULATOR_H_


#include <icetray/I3Module.h>
#include <g4-larsim/CCMParticleInjector.h>
#include <g4-larsim/CCMDetectorResponse.h>
#include <dataclasses/physics/I3Particle.h>
#include <simclasses/CCMMCPE.h>
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
        std::string hitSeriesName_;

        CCMParticleInjectorPtr injector_;
        CCMDetectorResponsePtr response_;
        boost::shared_ptr<I3Map<CCMPMTKey, std::vector<CCMMCPE>>> CCMMCPEMap = boost::make_shared<I3Map<CCMPMTKey, std::vector<CCMMCPE>>> ();
  
        static const std::string INC_ID_NAME;
  
        SET_LOGGER("CCMSimulator");
};

#endif
