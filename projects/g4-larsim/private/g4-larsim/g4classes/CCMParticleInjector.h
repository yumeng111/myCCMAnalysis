#ifndef CCMPARTICLEINJECTOR_H
#define CCMPARTICLEINJECTOR_H

#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3Context.h>
#include <icetray/I3Configuration.h>
#include <icetray/I3ServiceBase.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/I3SingleServiceFactory.h>

class CCMParticleInjector : public I3ServiceBase
{
 public:

  CCMParticleInjector(I3Configuration& config, const I3Context& context):
    context_(context),
    configuration_(config)
  {}

  virtual ~CCMParticleInjector() {}

  virtual void Configure() {}
  virtual I3MCTreePtr GetI3MCTree() {}
  virtual I3FramePtr GetI3ConfigurationFrame() {}

 private:

  const I3Context& context_;
  I3Configuration& configuration_;

  SET_LOGGER("CCMParticleInjector");
};


I3_POINTER_TYPEDEFS(CCMParticleInjector);

#endif
