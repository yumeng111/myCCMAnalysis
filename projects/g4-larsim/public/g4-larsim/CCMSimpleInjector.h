#ifndef CCMSIMPLEINJECTOR_H
#define CCMSIMPLEINJECTOR_H
// standard library stuff

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"

#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"

#include "phys-services/I3RandomService.h"

#include <string>

class CCMSimpleInjector : public CCMParticleInjector {
public:
    CCMSimpleInjector(const I3Context& context);
    ~CCMSimpleInjector() = default;

    virtual void Configure() override;

    virtual I3MCTreePtr GetMCTree() override;
    virtual std::map<std::string, I3FrameObjectPtr> GetSimulationConfiguration() override;
private:
    bool seen_s_frame_ = false;

    std::string mcPrimaryName_ = "MCPrimary";
    double energy_;
    double time_;
    I3Vector<double> location_;
    I3Vector<double> direction_;
    std::string typeName_;
    I3Particle::ParticleType particleType_;
public:
    I3Particle::ParticleType GetParticleType(const std::string& typeName);
};

#endif // CCMSIMPLEINJECTOR_H


