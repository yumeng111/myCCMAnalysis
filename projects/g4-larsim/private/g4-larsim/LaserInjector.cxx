// standard library stuff

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/LaserInjector.h"
#include "g4-larsim/CCMParticleInjector.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3SingleServiceFactory.h"

#include "phys-services/I3GSLRandomService.h"

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

const std::unordered_map<size_t, I3Vector<double>> LaserInjector::fibertoLocation = {{0, I3VectorDouble([10.16 * I3Units::cm, 2.54 * I3Units::cm, 58.0 * I3Units::cm])},
                                                                                     {1, I3VectorDouble([-13.261461129557242 * I3Units::cm, 100.73079791557913 * I3Units::cm, 0.0 * I3Units::cm])},
                                                                                     {2, I3VectorDouble([93.86616050314673 * I3Units::cm, -38.88063672829312 * I3Units::cm, 0.0 * I3Units::cm])},
                                                                                     {3, I3VectorDouble([-80.6046993735895 * I3Units::cm, -61.85016118728599 * I3Units::cm, 0.0 * I3Units::cm])},
                                                                                     {5, I3VectorDouble([10.16 * I3Units::cm, -2.54 * I3Units::cm, -58.0 * I3Units::cm])}};
LaserInjector::LaserInjector(const I3Context& context) :
    CCMParticleInjector(context), fiber_(0), wavelength_(532 * I3Units::nanometer), nphotons_(1e3), profile_({0.0, 1.0, 0.0, 12.0}), 
    mcPrimaryName_("CCMMCPrimary"), output_mc_tree_name_("CCMMCTree"), randomServiceName_("") {
    AddParameter("LaserFiber", "fiber you to simulate", fiber_);
    AddParameter("LaserWavelength", "wavlength of laser to simulate", wavelength_);
    AddParameter("NPhotonsProduced", "number of photons produced by laser pulse", nphotons_);
    AddParamter("DiffuserProfile", "parameters of two gaussian diffuser profile (ex: [sig1, mu1, sig2, mu2])", profile_);
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
    randomService_ = I3RandomServicePtr();
    AddParameter("RandomServiceName", "Name of the random service in the context. If empty default random service will be used.", randomServiceName_);
}

void LaserInjector::Configure() {
    GetParameter("LaserFiber", fiber_);
    GetParamter("LaserWavelength", wavelength_);
    GetParameter("NPhotonsProduced", nphotons_);
    GetParameter("DiffuserProfile", profile_);
    GetParameter("ParticleLocation", location_);
    GetParameter("ParticleDirection", direction_);
  
    // set our diffuser profile
    mu1 = profile_[0];
    sig1 = profile_[1];
    mu2 = profile_[2];
    sig2 = profile_[3];

    particleType_ = I3Particle::Gamma; // co-opting I3Particle::Gamma to refer to optical photons

    // now let's check fibertoLocation to grab location of fiber
    for (auto it = fibertoLocation.begin(); it != fibertoLocation.end(); ++it) {
        size_t this_fiber = it->first;
        if (this_fiber == fiber_){
            fiber_location_ = fibertoLocation.at(this_fiber);
        }
    }
    GetParameter("PrimaryName", mcPrimaryName_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
    GetParameter("RandomServiceName", randomServiceName_);
    if(randomServiceName_.empty()) {
        randomService_ = I3RandomServicePtr(new I3GSLRandomService(0));
        log_info("+ Random service: I3GSLRandomService  (default)");
    }
    else {
        randomService_ = GetContext().Get<I3RandomServicePtr>(randomServiceName_);
        if(randomService_) log_info("+ Random service: %s  (EXTERNAL)",  randomServiceName_.c_str());
        else log_fatal("No random service \"%s\" in context!", randomServiceName_.c_str());
    }
}

void LaserInjector::FillMCTree(I3FramePtr frame) { return nullptr; }

I3MCTreePtr LaserInjector::GetMCTree() {
    // first let's create our MC tree
    I3MCTreePtr mcTree = boost::make_shared<I3MCTree>();

    // let's create and fill our I3Particle
    // we want to fill out mcTree with a primary for the number of photons we are simulating
    // this will be the primary in our mcTree
    for (size_t photon_it = 0; photon_it < nphotons_; photon_it++){
        I3Particle primary(particleType_);
        primary.SetEnergy(energy_);
        primary.SetPos(fiber_location_[0], fiber_location_[1], fiber_location_[2]);
        // let's sample from our distributions to get particle direction
        sample1 = randomService_->Gaus(mu1, sig1);
        sample2 = randomService_->Gaus(mu2, sig2);
        combined_sample = sample1 + sample2; // this is our angle!!! now convert into direction vector

        //primary.SetDir(direction_[0], direction_[1], direction_[2]);
        
        I3MCTreeUtils::AddPrimary(*mcTree, primary);
    }
    
    return mcTree;
}


typedef I3SingleServiceFactory<LaserInjector,CCMParticleInjector> LaserInjectorFactory;

I3_SERVICE_FACTORY(LaserInjectorFactory);


