#ifndef CCMPARTICLEINJECTOR_H
#define CCMPARTICLEINJECTOR_H

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "icetray/I3Context.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3Configuration.h"
#include "icetray/I3PointerTypedefs.h"
#include "icetray/I3SingleServiceFactory.h"

class CCMParticleInjector : public I3ServiceBase {
    public:
        // construct self & declare configuration parameters
        CCMParticleInjector(const I3Context &context): I3ServiceBase(context) {}

        // cleanup
        virtual ~CCMParticleInjector() = default;

        // get configuration parameters
        virtual void Configure() = 0;

        // get I3MCTree and Configuration for simulation
        virtual I3MCTreePtr GetMCTree() = 0;
        virtual I3FrameObjectPtr GetSimulationConfiguration() = 0;

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

        SET_LOGGER("CCMParticleInjector");
};

I3_POINTER_TYPEDEFS(CCMParticleInjector);

#endif
