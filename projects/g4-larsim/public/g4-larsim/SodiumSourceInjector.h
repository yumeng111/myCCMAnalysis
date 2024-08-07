#ifndef SODIUMSOURCEINJECTOR_H
#define SODIUMSOURCEINJECTOR_H
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

static const unsigned sodiuminjectorconfig_version_ = 1;
class SodiumInjectorConfig : public I3FrameObject {
public:
    double z_position_;
    double inset_;
    double pellet_radius_;
    double pellet_height_;

    SodiumInjectorConfig() = default;
    virtual ~SodiumInjectorConfig() override = default;

    std::ostream& Print(std::ostream&) const;

    bool operator==(const SodiumInjectorConfig& rhs) const;
private:
    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

std::ostream& operator<<(std::ostream& oss, SodiumInjectorConfig const & bcm);
std::ostream& operator<<(std::ostream& oss, SodiumInjectorConfig & bcm);

template <class Archive>
void SodiumInjectorConfig::save(Archive& ar, unsigned version) const {
    if (version>sodiuminjectorconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of SodiumInjectorConfig class.",version,sodiuminjectorconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("z_position", z_position_);
    ar & make_nvp("inset", inset_);
    ar & make_nvp("pellet_radius", pellet_radius_);
    ar & make_nvp("pellet_height", pellet_height_);
}

template <class Archive>
void SodiumInjectorConfig::load(Archive& ar, unsigned version) {
    if (version>sodiuminjectorconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of SodiumInjectorConfig class.",version,sodiuminjectorconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("z_position", z_position_);
    ar & make_nvp("inset", inset_);
    ar & make_nvp("pellet_radius", pellet_radius_);
    ar & make_nvp("pellet_height", pellet_height_);
}

I3_CLASS_VERSION(SodiumInjectorConfig, sodiuminjectorconfig_version_);
I3_POINTER_TYPEDEFS(SodiumInjectorConfig);
I3_DEFAULT_NAME(SodiumInjectorConfig);

class SodiumSourceInjector : public CCMParticleInjector {
private:
    SodiumSourceInjector operator= (const SodiumSourceInjector& rhs);

    double z_position_;
    double inset_ = 0.25 * I3Units::cm;
    double pellet_radius_ = 0.4 * I3Units::cm;
    double pellet_height_ = 0.3 * I3Units::cm;

    std::string mcPrimaryName_;
    std::string output_mc_tree_name_;
    std::string randomServiceName_;
    I3RandomServicePtr randomService_;

    SET_LOGGER("SodiumSourceInjector");

public:

    SodiumSourceInjector(const I3Context& context);
    virtual ~SodiumSourceInjector() override = default;

    virtual void Configure() override;
    virtual void FillMCTree(I3FramePtr frame);
    virtual I3MCTreePtr GetMCTree() override;
    virtual I3FrameObjectPtr GetSimulationConfiguration() override;
};

#endif // SODIUMSOURCEINJECTOR_H


