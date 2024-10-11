#ifndef DetectorResponseConfig_H_INCLUDED
#define DetectorResponseConfig_H_INCLUDED

#include <icetray/I3DefaultName.h>
#include "icetray/I3FrameObject.h"
#include "dataclasses/Utility.h"

#include <string>
#include <iostream>
#include <sstream>

static const unsigned detectorresponseconfig_version_ = 1;
class DetectorResponseConfig : public I3FrameObject {
public:
    double rayleigh_scattering_length_;
    double pmt_tpb_qe_;
    double endcap_tpb_qe_;
    double side_tpb_qe_;
    double pmt_tpb_thickness_;
    double endcap_tpb_thickness_;
    double side_tpb_thickness_;
    double tpb_abs_tau_;
    double tpb_abs_norm_;

    DetectorResponseConfig() = default;
    virtual ~DetectorResponseConfig() override = default;

    std::ostream& Print(std::ostream&) const;

    bool operator==(const DetectorResponseConfig& rhs) const;
private:
    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm);
std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm);

template <class Archive>
void DetectorResponseConfig::save(Archive& ar, unsigned version) const {
    if (version>detectorresponseconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of DetectorResponseConfig class.",version,detectorresponseconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("rayleigh_scattering_length", rayleigh_scattering_length_);
    ar & make_nvp("pmt_tpb_qe", pmt_tpb_qe_);
    ar & make_nvp("endcap_tpb_qe", endcap_tpb_qe_);
    ar & make_nvp("side_tpb_qe", side_tpb_qe_);
    ar & make_nvp("pmt_tpb_thickness", pmt_tpb_thickness_);
    ar & make_nvp("endcap_tpb_thickness", endcap_tpb_thickness_);
    ar & make_nvp("side_tpb_thickness", side_tpb_thickness_);
    ar & make_nvp("tpb_abs_tau", tpb_abs_tau_);
    ar & make_nvp("tpb_abs_norm", tpb_abs_norm_);
}

template <class Archive>
void DetectorResponseConfig::load(Archive& ar, unsigned version) {
    if (version>detectorresponseconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of DetectorResponseConfig class.",version,detectorresponseconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("rayleigh_scattering_length", rayleigh_scattering_length_);
    ar & make_nvp("pmt_tpb_qe", pmt_tpb_qe_);
    ar & make_nvp("endcap_tpb_qe", endcap_tpb_qe_);
    ar & make_nvp("side_tpb_qe", side_tpb_qe_);
    ar & make_nvp("pmt_tpb_thickness", pmt_tpb_thickness_);
    ar & make_nvp("endcap_tpb_thickness", endcap_tpb_thickness_);
    ar & make_nvp("side_tpb_thickness", side_tpb_thickness_);
    if (detectorresponseconfig_version_ == 0){
        tpb_abs_tau_ = 0.0;
        tpb_abs_norm_ = 0.0;
    } else {
        ar & make_nvp("tpb_abs_tau", tpb_abs_tau_);
        ar & make_nvp("tpb_abs_norm", tpb_abs_norm_);
    }
}

I3_CLASS_VERSION(DetectorResponseConfig, detectorresponseconfig_version_);
I3_POINTER_TYPEDEFS(DetectorResponseConfig);
I3_DEFAULT_NAME(DetectorResponseConfig);

#endif // DetectorResponseConfig_H_INCLUDED

