#ifndef SODIUMSOURCEINJECTOR_H
#define SODIUMSOURCEINJECTOR_H
// standard library stuff

#include "simclasses/SodiumInjectorConfig.h"

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3ServiceBase.h"

#include "phys-services/I3RandomService.h"

#include <vector>
#include <string>
#include <algorithm>

class SodiumSourceInjector : public CCMParticleInjector {
public:
    SodiumSourceInjector(const I3Context& context);
    ~SodiumSourceInjector() = default;

    virtual void Configure() override;
    virtual I3MCTreePtr GetMCTree() override;
    virtual I3FrameObjectPtr GetSimulationConfiguration() override;
private:
    double z_position_;
    double inset_ = 0.25 * I3Units::cm;
    double pellet_radius_ = 0.4 * I3Units::cm;
    double pellet_height_ = 0.3 * I3Units::cm;

    std::string randomServiceName_;
    I3RandomServicePtr randomService_;
};

#endif // SODIUMSOURCEINJECTOR_H

