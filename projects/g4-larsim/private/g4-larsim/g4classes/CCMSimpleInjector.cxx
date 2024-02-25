// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

#include "g4-larsim/g4classes/CCMSimpleInjector.h"
#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3SingleServiceFactory.h>

CCMSimpleInjector::CCMSimpleInjector(const I3Context& context): I3ServiceBase(context)
{
    AddParameter("ParticleEnergy", "energy of particle to inject into Geant4 in MeV", default_energy_);
    AddParameter("ParticleLocation", "location of particle to inject into Geant4 in m", default_location_);
    AddParameter("ParticleType", "type of particle to inject into Geant4", default_particle_type_);
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", default_mcPrimaryName_);
}

CCMSimpleInjector::~CCMSimpleInjector(){}

void CCMSimpleInjector::Configure(){
    GetParameter("ParticleEnergy", energy_);
    GetParameter("ParticleLocation", location_);
    GetParameter("ParticleType", particle_type_);
    GetParameter("PrimaryName", mcPrimaryName_);
}


void CCMSimpleInjector::FillMCTree(){
  
    // first let's create our MC tree
    I3ParticlePtr mcPrimary(new I3Particle(I3Particle::Null, I3Particle::unknown));
    I3MCTreePtr mcTree = I3MCTreePtr();
    I3MCTreeUtils::AddPrimary(*mcTree, *mcPrimary);

    // let's create and fill our I3Particle
    I3Particle particle;
    particle.SetType(particle_type_);
    particle.SetEnergy(energy_);
    I3MCTreeUtils::AppendChild(*mcTree, *mcPrimary, particle);

    // now save to a frame
    I3FramePtr frame;
    frame->Put(mcPrimaryName_, mcPrimary);
}

I3MCTreePtr CCMSimpleInjector::GetI3MCTree(){
    return mcTree;
}

I3FramePtr CCMSimpleInjector::GetI3ConfigurationFrame(){
    return frame; 
}

typedef
I3SingleServiceFactory<CCMSimpleInjector,CCMParticleInjector>
CCMSimpleInjectorFactory;
I3_SERVICE_FACTORY( CCMSimpleInjectorFactory );
