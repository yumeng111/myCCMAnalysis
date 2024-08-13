#ifndef SodiumInjectorConfig_H_INCLUDED
#define SodiumInjectorConfig_H_INCLUDED

#include <icetray/I3DefaultName.h>
#include "icetray/I3FrameObject.h"
#include "dataclasses/Utility.h"

#include <string>
#include <iostream>
#include <sstream>

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

#endif // SodiumInjectorConfig_H_INCLUDED

