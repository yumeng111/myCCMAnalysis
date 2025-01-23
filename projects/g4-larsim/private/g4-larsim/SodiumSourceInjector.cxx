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
    CCMParticleInjector(context), z_position_(0.0 * I3Units::cm), randomServiceName_(""), training_(false),
    decay_x_(0.0 * I3Units::cm), decay_y_(0.0 * I3Units::cm), decay_z_(0.0 * I3Units::cm){
    AddParameter("SourceZPosition", "z location of source pellet", z_position_);
    AddParameter("Inset", "Inset of the source pellet", inset_);
    AddParameter("PelletRadius", "Radius of the source pellet", pellet_radius_);
    AddParameter("PelletHeight", "Height of the source pellet", pellet_height_);
    randomService_ = I3RandomServicePtr();
    AddParameter("RandomServiceName", "Name of the random service in the context. If empty default random service will be used.", randomServiceName_);
    AddParameter("Training", "If true, allow sodium pellet to be anywhere in detector training purposes", training_);
    AddParameter("DecayX", "if simulating training data, will use decay x for x pos of soidum", decay_x_);
    AddParameter("DecayY", "if simulating training data, will use decay y for y pos of soidum", decay_y_);
    AddParameter("DecayZ", "if simulating training data, will use decay z for z pos of soidum", decay_z_);
}

void SodiumSourceInjector::Configure() {
    CCMParticleInjector::Configure();
    GetParameter("SourceZPosition", z_position_);
    GetParameter("Inset", inset_);
    GetParameter("PelletRadius", pellet_radius_);
    GetParameter("PelletHeight", pellet_height_);
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
    GetParameter("Training", training_);
    GetParameter("DecayX", decay_x_);
    GetParameter("DecayY", decay_y_);
    GetParameter("DecayZ", decay_z_);
}

I3MCTreePtr SodiumSourceInjector::GetMCTree() {
    // first let's create our MC tree
    I3MCTreePtr mcTree = boost::make_shared<I3MCTree>();

    // let's make the location of our event randomly within our sodium pellet

    // let's create and fill our I3Particle
    I3Particle primary(I3Particle::Na22Nucleus);

    double theta_pos;
    double r;
    double x;
    double y;
    double z;

    if (training_){
        // ok we want a sodium position for training purposes!
        x = decay_x_;
        y = decay_y_;
        z = decay_z_;
    } else {
        // want realistic sodium position (ie in the sodium pellet)
        theta_pos = randomService_->Uniform(0.0, 2.0*M_PI);
        r = std::sqrt(randomService_->Uniform(0.0, pellet_radius_ * pellet_radius_));
        x = r * std::cos(theta_pos);
        y = r * std::sin(theta_pos);
        z = randomService_->Uniform(z_position_ + inset_, z_position_ + inset_ + pellet_height_);
    }

    primary.SetTime(0.0);
    primary.SetPos(x, y, z);
    primary.SetEnergy(0.0 * I3Units::MeV);
    primary.SetDir(0.0, 0.0, 0.0); // doesnt really matter

    I3MCTreeUtils::AddPrimary(*mcTree, primary);

    return mcTree;
}

I3FrameObjectPtr SodiumSourceInjector::GetSimulationConfiguration() {
    SodiumInjectorConfigPtr config = boost::make_shared<SodiumInjectorConfig>();
    config->z_position_ = z_position_;
    config->inset_ = inset_;
    config->pellet_radius_ = pellet_radius_;
    config->pellet_height_ = pellet_height_;
    return config;
}

I3_MODULE(SodiumSourceInjector);

