#ifndef CCMParticleInjector_H
#define CCMParticleInjector_H
// standard library stuff

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"

#include "phys-services/I3RandomService.h"

#include <vector>
#include <string>
#include <algorithm>

class CCMParticleInjector : public I3ConditionalModule {
public:
    CCMParticleInjector(const I3Context& context);
    ~CCMParticleInjector() = default;

    void Configure();
    void Simulation(I3FramePtr frame);
    void DAQ(I3FramePtr frame);

    virtual I3MCTreePtr GetMCTree() = 0;
    virtual I3FrameObjectPtr GetSimulationConfiguration() = 0;

    void FillSimulationFrame(I3FramePtr frame);

private:
    bool seen_s_frame_ = false;

    std::string mcPrimaryName_;
    std::string output_mc_tree_name_;
};

#endif // CCMParticleInjector_H

