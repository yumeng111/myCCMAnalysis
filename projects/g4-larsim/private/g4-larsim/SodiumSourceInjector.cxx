// standard library stuff

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/SodiumSourceInjector.h"
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
#include <random>

SodiumSourceInjector::SodiumSourceInjector(const I3Context& context) :
    CCMParticleInjector(context), z_position_(0.0 * I3Units::cm), mcPrimaryName_("CCMMCPrimary"), output_mc_tree_name_("CCMMCTree"), randomServiceName_("") {
    AddParameter("SourceZPosition", "z location of source pellet", z_position_);
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
    randomService_ = I3RandomServicePtr();
    AddParameter("RandomServiceName", "Name of the random service in the context. If empty default random service will be used.", randomServiceName_);
}

void SodiumSourceInjector::Configure() {
    GetParameter("SourceZPosition", z_position_);
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

void SodiumSourceInjector::FillMCTree(I3FramePtr frame) {}

I3MCTreePtr SodiumSourceInjector::GetMCTree() {
    // first let's create our MC tree
    I3MCTreePtr mcTree = boost::make_shared<I3MCTree>();

    // let's make the location of our event randomly within our sodium pellet
    // set up geometry of sodium source pellet
    double inset = 0.25 * I3Units::cm;
    double pellet_radius = 0.4 * I3Units::cm;
    double pellet_height = 0.3 * I3Units::cm;

    // now throw events
    // note -- using c++ random distributions, but should probably change to use I3GSLRandomService at some point,,.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis_angle(0.0, 2.0*M_PI); // uniform distribution 0 - 2pi
    std::uniform_real_distribution<double> dis_z(z_position_ + inset, z_position_ + inset + pellet_height); // uniform distribution across z position of sodium pellet 
    std::uniform_real_distribution<double> dis_radius_squared(0.0, pellet_radius * pellet_radius);

    for (size_t p = 0; p < 2; p++){
        // let's create and fill our I3Particle
        I3Particle primary(I3Particle::Na22Nucleus);

        double theta_pos = dis_angle(gen);
        double r = std::sqrt(dis_radius_squared(gen));
        double x = r * std::cos(theta_pos);
        double y = r * std::sin(theta_pos);
        double z = dis_z(gen);

        std::cout << "event location = " << x << ", " << y << ", " << z << std::endl;
        primary.SetPos(x, y, z);
        primary.SetEnergy(0.0 * I3Units::MeV); 
        primary.SetDir(0.0, 0.0, 0.0); // doesnt really matter
        
        I3MCTreeUtils::AddPrimary(*mcTree, primary);
    }

    return mcTree;
}


typedef I3SingleServiceFactory<SodiumSourceInjector,CCMParticleInjector> SodiumSourceInjectorFactory;

I3_SERVICE_FACTORY(SodiumSourceInjectorFactory);


