#ifndef AnalyticLightYieldGenerator_H_INCLUDED
#define AnalyticLightYieldGenerator_H_INCLUDED

#include "dataclasses/Utility.h"

#include <string>
#include <iostream>
#include <sstream>

static const unsigned analyticlightyieldgenerator_version_ = 0;

class AnalyticLightYieldGenerator : public I3FrameObject {
    public:
    
    enum class LArLightProfileType : int8_t {
        Unknown = 0,
        Full = 1,
        Simplified = 2,
    };

    float Rs;
    float Rt;
    float tau_s;
    float tau_t;
    float tau_rec;
    float tau_TPB;
    float normalization;
    float time_offset;

    float n_sodium_events;
    float uv_absorption;
    float z_offset;
    LArLightProfileType light_profile_type; 

    SET_LOGGER("AnalyticLightYieldGenerator");

    bool operator==(const AnalyticLightYieldGenerator& rhs) const {
        return Rs == rhs.Rs
            && Rt == rhs.Rt
            && tau_s == rhs.tau_s
            && tau_t == rhs.tau_t
            && tau_rec == rhs.tau_rec
            && tau_TPB == rhs.tau_TPB
            && uv_absorption == rhs.uv_absorption
            && normalization == rhs.normalization
            && n_sodium_events == rhs.n_sodium_events
            && z_offset == rhs.z_offset
            && light_profile_type == rhs.light_profile_type;
    }


  AnalyticLightYieldGenerator(float Rs_ = 0, float Rt_ = 0, float tau_s_ = 0, float tau_t_ = 0, float tau_rec_ = 0, float tau_TPB_ = 0,
          float uv_absorption_ = 0, float normalization_ = 0, float n_sodium_events_ = 0, float z_offset_ = 0,
          LArLightProfileType light_profile_type_ = AnalyticLightYieldGenerator::LArLightProfileType::Unknown):
        Rs(Rs_), Rt(Rt_), tau_s(tau_s_), tau_t(tau_t_), tau_rec(tau_rec_), tau_TPB(tau_TPB_),
        uv_absorption(uv_absorption_), normalization(normalization_), n_sodium_events(n_sodium_events_), z_offset(z_offset_), light_profile_type(light_profile_type_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
        static const std::unordered_map<AnalyticLightYieldGenerator::LArLightProfileType, std::string> LArLightProfileTypeToName;
 
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>analyticlightyieldgenerator_version_)
            log_fatal("Attempting to read version %u from file but running version %u of analyticlightyieldgenerator class.",
                    version,analyticlightyieldgenerator_version_);
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
        ar & make_nvp("Rs",Rs);
        ar & make_nvp("Rt",Rt);
        ar & make_nvp("tau_s",tau_s);
        ar & make_nvp("tau_t",tau_t);
        ar & make_nvp("tau_rec",tau_rec);
        ar & make_nvp("tau_TPB",tau_TPB);
        ar & make_nvp("normalization",normalization);
        ar & make_nvp("time_offset", time_offset);

        ar & make_nvp("n_sodium_events",n_sodium_events);
        ar & make_nvp("uv_absorption",uv_absorption);
        ar & make_nvp("z_offset",z_offset);
        ar & make_nvp("light_profile_type",light_profile_type);
    }

};

I3_CLASS_VERSION(AnalyticLightYieldGenerator,analyticlightyieldgenerator_version_);

std::ostream& operator<<(std::ostream&, const AnalyticLightYieldGenerator&);

I3_POINTER_TYPEDEFS(AnalyticLightYieldGenerator);
#endif
