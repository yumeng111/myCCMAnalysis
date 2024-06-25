#ifndef SODIUMSOURCEINJECTOR_H
#define SODIUMSOURCEINJECTOR_H
// standard library stuff

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
private:
    SodiumSourceInjector operator= (const SodiumSourceInjector& rhs);

    double z_position_;
    std::string mcPrimaryName_;
    std::string output_mc_tree_name_;
    I3Particle::ParticleType particleType_;
    std::string randomServiceName_;
    I3RandomServicePtr randomService_;

    SET_LOGGER("SodiumSourceInjector");

public:

    SodiumSourceInjector(const I3Context& context);
    virtual ~SodiumSourceInjector() override = default;

    virtual void Configure() override;
    virtual void FillMCTree(I3FramePtr frame);
    virtual I3MCTreePtr GetMCTree() override;
};

#endif // SODIUMSOURCEINJECTOR_H


