// standard library stuff

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCMSimpleInjector.h"
#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3SingleServiceFactory.h"

#include "phys-services/I3RandomService.h"

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#define PARTICLE_TYPE_MAP_APPEND(r, data, t) types[BOOST_PP_STRINGIZE(t)] = I3Particle::t;
namespace {
  const std::map<std::string, I3Particle::ParticleType>& AvailableTypes()
  {
    static std::map<std::string, I3Particle::ParticleType> types;
    if (types.size() == 0) {
      // default types for backward compatibility
      types["mu-"] = I3Particle::MuMinus;
      types["mu+"] = I3Particle::MuPlus;
      types["e-"] = I3Particle::EMinus;
      types["e+"] = I3Particle::EPlus;
      types["gamma"] = I3Particle::Gamma;
      types["p+"] = I3Particle::PPlus;
      types["p-"] = I3Particle::PMinus;
      types["n"] = I3Particle::Neutron;
      types["nbar"] = I3Particle::NeutronBar;
      types["pi0"] = I3Particle::Pi0;
      types["pi+"] = I3Particle::PiPlus;
      types["pi-"] = I3Particle::PiMinus;
      BOOST_PP_SEQ_FOR_EACH(PARTICLE_TYPE_MAP_APPEND, ~, I3PARTICLE_H_I3Particle_ParticleType);
    }
    return types;
  }
}

CCMSimpleInjector::CCMSimpleInjector(const I3Context& context) :
    CCMParticleInjector(context), energy_(1.0 * I3Units::MeV), location_({0.0, 0.0, 0.0}), direction_({1.0, 0., 0.}),
    typeName_("e-"), mcPrimaryName_("CCMMCPrimary"), output_mc_tree_name_("CCMMCTree") {
    AddParameter("ParticleEnergy", "energy of particle to inject into Geant4 in MeV", energy_);
    AddParameter("ParticleLocation", "location of particle to inject into Geant4 in m", location_);
    AddParameter("ParticleDirection", "direction of particle to inject into Geant4", direction_);
    AddParameter("ParticleType", "type of particle to inject into Geant4", typeName_);
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
}

void CCMSimpleInjector::Configure() {
    GetParameter("ParticleEnergy", energy_);
    GetParameter("ParticleLocation", location_);
    GetParameter("ParticleDirection", direction_);
    //GetParameter("ParticleType", particle_type_);
   
    std::string typeName;
    GetParameter("ParticleType", typeName);
    particleType_ = GetParticleType(typeName);
    if(particleType_ == I3Particle::unknown) {
        std::ostringstream msg;
        msg << "Invalid particle type " << typeName.c_str() << ". Available types are:\n";
        for (std::map<std::string, I3Particle::ParticleType>::const_iterator it = AvailableTypes().begin(); it != AvailableTypes().end(); ++it) {
            msg << "   " << it->first << "\n";
        }
        log_fatal("%s",msg.str().c_str());
    }
    log_info("+ Particle type: %s", typeName.c_str());
    
    GetParameter("PrimaryName", mcPrimaryName_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
}

I3Particle::ParticleType CCMSimpleInjector::GetParticleType(const std::string& typeName) {
    std::map<std::string, I3Particle::ParticleType>::const_iterator p = AvailableTypes().find(typeName);
    if (p != AvailableTypes().end())
        return p->second;
    return I3Particle::unknown;
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
    I3MCTreePtr mcTree = boost::make_shared<I3MCTree>();

    // let's create and fill our I3Particle
    // this will be the primary in our mcTree
    I3Particle primary(particleType_);
    primary.SetEnergy(energy_);
    primary.SetPos(location_[0], location_[1], location_[2]);
    primary.SetDir(direction_[0], direction_[1], direction_[2]);
    
    I3MCTreeUtils::AddPrimary(*mcTree, primary);
    return mcTree;
}

//I3FrameObjectPtr CCMSimpleInjector::GetSimulationConfiguration() {
//    return nullptr;
//}

typedef I3SingleServiceFactory<CCMSimpleInjector,CCMParticleInjector> CCMSimpleInjectorFactory;

I3_SERVICE_FACTORY(CCMSimpleInjectorFactory);


