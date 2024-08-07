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

std::ostream & SodiumInjectorConfig::Print(std::ostream& oss) const {
    oss << "SodiumInjectorConfig: " << std::endl;
    oss << "  z_position: " << z_position_ << std::endl;
    oss << "  inset: " << inset_ << std::endl;
    oss << "  pellet_radius: " << pellet_radius_ << std::endl;
    oss << "  pellet_height: " << pellet_height_ << std::endl;
    return oss;
}

bool SodiumInjectorConfig::operator==(const SodiumInjectorConfig& rhs) const {
    return (z_position_ == rhs.z_position_ &&
            inset_ == rhs.inset_ &&
            pellet_radius_ == rhs.pellet_radius_ &&
            pellet_height_ == rhs.pellet_height_);
}

std::ostream& operator<<(std::ostream& oss, SodiumInjectorConfig const & bcm) {
    return bcm.Print(oss);
}

std::ostream& operator<<(std::ostream& oss, SodiumInjectorConfig & bcm) {
    return bcm.Print(oss);
}

SodiumSourceInjector::SodiumSourceInjector(const I3Context& context) :
    CCMParticleInjector(context), z_position_(0.0 * I3Units::cm), mcPrimaryName_("CCMMCPrimary"), output_mc_tree_name_("CCMMCTree"), randomServiceName_("") {
    AddParameter("SourceZPosition", "z location of source pellet", z_position_);
    AddParameter("Inset", "Inset of the source pellet", inset_);
    AddParameter("PelletRadius", "Radius of the source pellet", pellet_radius_);
    AddParameter("PelletHeight", "Height of the source pellet", pellet_height_);
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
    randomService_ = I3RandomServicePtr();
    AddParameter("RandomServiceName", "Name of the random service in the context. If empty default random service will be used.", randomServiceName_);
}

void SodiumSourceInjector::Configure() {
    GetParameter("SourceZPosition", z_position_);
    GetParameter("Inset", inset_);
    GetParameter("PelletRadius", pellet_radius_);
    GetParameter("PelletHeight", pellet_height_);
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

    // let's create and fill our I3Particle
    I3Particle primary(I3Particle::Na22Nucleus);

    double theta_pos = randomService_->Uniform(0.0, 2.0*M_PI);
    double r = std::sqrt(randomService_->Uniform(0.0, pellet_radius_ * pellet_radius_));
    double x = r * std::cos(theta_pos);
    double y = r * std::sin(theta_pos);
    double z = randomService_->Uniform(z_position_ + inset_, z_position_ + inset_ + pellet_height_);

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

typedef I3SingleServiceFactory<SodiumSourceInjector,CCMParticleInjector> SodiumSourceInjectorFactory;

I3_SERVICE_FACTORY(SodiumSourceInjectorFactory);


