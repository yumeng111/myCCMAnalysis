#ifndef LASERINJECTOR_H
#define LASERINJECTOR_H
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

class LaserInjector : public CCMParticleInjector {
private:
    LaserInjector operator= (const LaserInjector& rhs);

    size_t fiber_;
    double wavelength_;
    size_t nphotons_; 
    I3Vector<double> profile_;
    I3Vector<double> fiber_location_;
    double mu1;
    double sig1;
    double mu2;
    double sig2;
    double sample1;
    double sample2;
    double combined_sample;
    std::string mcPrimaryName_;
    std::string output_mc_tree_name_;
    I3Particle::ParticleType particleType_;
    std::string randomServiceName_;
    I3RandomServicePtr randomService_;

    static const std::unordered_map<size_t, I3Vector<double>> fibertoLocation;

    SET_LOGGER("LaserInjector");

public:

    LaserInjector(const I3Context& context);
    virtual ~LaserInjector() override = default;

    virtual void Configure() override;
    virtual void FillMCTree(I3FramePtr frame);
    virtual I3MCTreePtr GetMCTree() override;
};

#endif // LASERINJECTOR_H


