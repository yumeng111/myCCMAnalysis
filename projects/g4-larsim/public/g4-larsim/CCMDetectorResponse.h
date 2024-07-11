#ifndef CCMDETECTORRESPONSE_H
#define CCMDETECTORRESPONSE_H

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "icetray/I3Context.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3Configuration.h"
#include "icetray/I3PointerTypedefs.h"
#include "icetray/I3SingleServiceFactory.h"

#include "simclasses/CCMMCPE.h"

class CCMDetectorResponse : public I3ServiceBase {
    public:
        // construct self & declare configuration parameters
        CCMDetectorResponse(const I3Context &context): I3ServiceBase(context) {}

        // cleanup
        virtual ~CCMDetectorResponse() = default;

        // get configuration parameters
        virtual void Configure() = 0;

        // initialize geant4 detector
        virtual void Initialize() = 0;

        // begin events
        virtual void BeginEvent(const I3Particle& primary) = 0;
        
        // end events
        virtual void EndEvent(I3MCTreePtr & LArEnergyDep, boost::shared_ptr<CCMMCPESeriesMap> & CCMMCPEMap, I3MCTreePtr &  photon_summary) = 0;

        // terminate run
        virtual void TerminateRun() = 0;

        // get SD status
        virtual bool GetPMTSDStatus() = 0;
        virtual bool GetLArSDStatus() = 0;

    protected:
        template <class ParamType>
        void AddParameter(const std::string& name,
                        const std::string& description,
                        const ParamType& defaultValue)
        { I3ServiceBase::AddParameter<ParamType>(name, description, defaultValue); }

        template <class ParamType>
        void GetParameter(const std::string& name, ParamType& value) const
        { I3ServiceBase::GetParameter<ParamType>(name, value); }

        const I3Context& GetContext()
        { return I3ServiceBase::GetContext(); }
        
        SET_LOGGER("CCMDetectorResponse");
};

I3_POINTER_TYPEDEFS(CCMDetectorResponse);

#endif
