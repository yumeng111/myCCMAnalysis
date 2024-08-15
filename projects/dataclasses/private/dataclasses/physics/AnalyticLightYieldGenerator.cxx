#include <dataclasses/physics/AnalyticLightYieldGenerator.h>
#include <ostream>


const std::unordered_map<AnalyticLightYieldGenerator::LArLightProfileType, std::string> AnalyticLightYieldGenerator::LArLightProfileTypeToName = 
                                                                                                     {{AnalyticLightYieldGenerator::LArLightProfileType::Unknown, "Unknown"},
                                                                                                      {AnalyticLightYieldGenerator::LArLightProfileType::Full, "Full"},
                                                                                                      {AnalyticLightYieldGenerator::LArLightProfileType::Simplified, "Simplified"}};
std::ostream& operator<<(std::ostream& os, const AnalyticLightYieldGenerator& pe) {
    return(pe.Print(os));
}

std::ostream& AnalyticLightYieldGenerator::Print(std::ostream& os) const{
    os << "[ AnalyticLightYieldGenerator::"
        << "\n  Ratio Singlet :" << Rs
        << "\n  Ratio Triplet :" << Rt
        << "\n  Singlet Time Constant :" << tau_s
        << "\n  Triplet Time Constant :" << tau_t
        << "\n  Recombination Time Constant :" << tau_rec
        << "\n  TPB Time Constant :" << tau_TPB
        << "\n  UV Absorption Length :" << uv_absorption
        << "\n  Normalization :" << normalization
        << "\n  Time Offset :" << time_offset
        << "\n  Constant Offset :" << const_offset
        << "\n  N Sodium Events :" << n_sodium_events
        << "\n  Z Offset :" << z_offset
        << "\n  Light Profile Type :" << LArLightProfileTypeToName.at(light_profile_type)
        << " ]";
    return os;
}


