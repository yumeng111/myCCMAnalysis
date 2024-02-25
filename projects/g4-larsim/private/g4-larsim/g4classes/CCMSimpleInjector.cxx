// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

#include "CCMSimpleInjector.h"
#include "CCMParticleInjector.h"
#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3SingleServiceFactory.h>

CCMSimpleInjector::CCMSimpleInjector(const I3Context& context) :
    CCMParticleInjector(context), energy_(1.0 * I3Units::MeV), location_({0.0, 0.0, 0.0}),
    particle_type_(I3Particle::ParticleType::EMinus), mcPrimaryName_("CCMMCPrimary"), output_mc_tree_name_("I3MCTree") {
    AddParameter("ParticleEnergy", "energy of particle to inject into Geant4 in MeV", energy_);
    AddParameter("ParticleLocation", "location of particle to inject into Geant4 in m", location_);
    AddParameter("ParticleType", "type of particle to inject into Geant4", particle_type_);
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
}

void CCMSimpleInjector::Configure() {
    GetParameter("ParticleEnergy", energy_);
    GetParameter("ParticleLocation", location_);
    GetParameter("ParticleType", particle_type_);
    GetParameter("PrimaryName", mcPrimaryName_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
}

void CCMSimpleInjector::FillMCTree(I3FramePtr frame) {

    I3MCTreePtr mcTree = GetMCTree();

    std::vector<I3Particle> primaries = I3MCTreeUtils::GetPrimaries(*mcTree);
    I3VectorI3ParticlePtr primary(new I3Vector<I3Particle>(primaries.begin(), primaries.end()));

    // now save to a frame
    frame->Put(mcPrimaryName_, primary);
    frame->Put(output_mc_tree_name_, mcTree);
}

I3MCTreePtr CCMSimpleInjector::GetMCTree() {
    // first let's create our MC tree
    I3ParticlePtr mcPrimary(new I3Particle(I3Particle::Null, I3Particle::unknown));
    I3MCTreePtr mcTree = I3MCTreePtr();
    I3MCTreeUtils::AddPrimary(*mcTree, *mcPrimary);

    // let's create and fill our I3Particle
    I3Particle particle;
    particle.SetType(particle_type_);
    particle.SetEnergy(energy_);
    I3MCTreeUtils::AppendChild(*mcTree, *mcPrimary, particle);

    return mcTree;
}

I3FrameObjectPtr CCMSimpleInjector::GetSimulationConfiguration() {
    return nullptr;
}

typedef I3SingleServiceFactory<CCMSimpleInjector,CCMParticleInjector> CCMSimpleInjectorFactory;

I3_SERVICE_FACTORY(CCMSimpleInjectorFactory);
