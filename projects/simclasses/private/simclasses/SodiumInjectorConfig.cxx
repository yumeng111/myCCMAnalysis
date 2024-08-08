#include "simclasses/SodiumInjectorConfig.h"

#include "icetray/I3FrameObject.h"

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

I3_SERIALIZABLE(SodiumInjectorConfig);
