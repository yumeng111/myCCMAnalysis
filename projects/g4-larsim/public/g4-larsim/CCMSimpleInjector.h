#ifndef CCMSIMPLEINJECTOR_H
#define CCMSIMPLEINJECTOR_H
// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

#include <icetray/I3Frame.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3ServiceBase.h>
#include "g4-larsim/CCMParticleInjector.h"

class CCMSimpleInjector : public CCMParticleInjector {
private:
    CCMSimpleInjector operator= (const CCMSimpleInjector& rhs);

    double energy_;
    I3Vector<double> location_;
    I3Particle::ParticleType particle_type_;
    std::string mcPrimaryName_;
    std::string output_mc_tree_name_;

    SET_LOGGER("CCMSimpleInjector");

public:

    CCMSimpleInjector(const I3Context& context);
    virtual ~CCMSimpleInjector() override = default;

    void Configure();

    virtual void FillMCTree(I3FramePtr frame);
    virtual I3MCTreePtr GetMCTree() override;
    virtual I3FrameObjectPtr GetSimulationConfiguration() override;
};

#endif // CCMSIMPLEINJECTOR_H


