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

class CCMParticleInjector : public I3ServiceBase {
 private: 
 public:


  // construct self & declare configuration parameters
  CCMParticleInjector(const I3Context &context): I3ServiceBase(context) {}

  // cleanup
  virtual ~CCMParticleInjector() {}

  // get configuration parameters
  virtual void Configure() {}

  // get IMCTree and Configuration for simulation  
  I3MCTreePtr GetI3MCTree() {return 0;}
  I3FramePtr GetI3ConfigurationFrame() {return 0;}

  SET_LOGGER("CCMParticleInjector");
};


//I3_POINTER_TYPEDEFS(CCMParticleInjector);

#endif
