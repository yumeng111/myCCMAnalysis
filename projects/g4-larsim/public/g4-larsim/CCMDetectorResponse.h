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
#include "simclasses/PhotonSummary.h"

class CCMDetectorResponse : public I3ServiceBase {
    public:
        // construct self & declare configuration parameters
        CCMDetectorResponse(const I3Context &context): I3ServiceBase(context) {}

        // cleanup
        virtual ~CCMDetectorResponse() = default;

        // get configuration parameters
        virtual void Configure() = 0;
        
        // get config information
        virtual I3FrameObjectPtr GetSimulationConfiguration() = 0;

        // initialize geant4 detector
        virtual void Initialize() = 0;

        virtual void SimulateEvent(I3Particle const & primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) = 0;
        virtual void SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries) = 0;

        virtual void DestroyInterface() = 0;

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
