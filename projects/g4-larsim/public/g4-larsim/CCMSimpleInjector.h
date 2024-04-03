#ifndef CCMSIMPLEINJECTOR_H
#define CCMSIMPLEINJECTOR_H
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

class CCMSimpleInjector : public CCMParticleInjector {
private:
    CCMSimpleInjector operator= (const CCMSimpleInjector& rhs);

    double energy_;
    I3Vector<double> location_;
    I3Vector<double> direction_;
    std::string typeName_;
    std::string mcPrimaryName_;
    std::string output_mc_tree_name_;
    I3Particle::ParticleType particleType_;
    I3Particle::ParticleType GetParticleType(const std::string& typeName);

    SET_LOGGER("CCMSimpleInjector");

public:

    CCMSimpleInjector(const I3Context& context);
    virtual ~CCMSimpleInjector() override = default;

    virtual void Configure() override;
    virtual void FillMCTree(I3FramePtr frame);
    virtual I3MCTreePtr GetMCTree() override;
    //virtual I3FrameObjectPtr GetSimulationConfiguration() override;
};

#endif // CCMSIMPLEINJECTOR_H


