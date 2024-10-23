#include "simclasses/DetectorResponseConfig.h"

#include "icetray/I3FrameObject.h"

std::ostream & DetectorResponseConfig::Print(std::ostream& oss) const {
    oss << "DetectorResponseConfig: "
        << "\n  Rayleigh Scattering Length :" << rayleigh_scattering_length_
        << "\n  UV Absorption Length :" << uv_absorption_length_ 
        << "\n  PMT TPB QE :" << pmt_tpb_qe_
        << "\n  Endcap TPB QE :" << endcap_tpb_qe_
        << "\n  Side TPB QE :" << side_tpb_qe_
        << "\n  PMT TPB Thickness :" << pmt_tpb_thickness_
        << "\n  Endcap TPB Thickness :" << endcap_tpb_thickness_
        << "\n  Side TPB Thickness :" << side_tpb_thickness_
        << "\n  TPB Absorption Tau :" << tpb_abs_tau_
        << "\n  TPB Absorption Norm :" << tpb_abs_norm_
        << "\n  TPB Absorption Scale :" << tpb_abs_scale_
        << "\n  Mie GG :" << mie_gg_
        << "\n  Mie Ratio :" << mie_ratio_;
    return oss;
}

bool DetectorResponseConfig::operator==(const DetectorResponseConfig& rhs) const {
    return (rayleigh_scattering_length_ == rhs.rayleigh_scattering_length_ &&
            uv_absorption_length_ == rhs.uv_absorption_length_ &&
            pmt_tpb_qe_ == rhs.pmt_tpb_qe_ &&
            endcap_tpb_qe_ == rhs.endcap_tpb_qe_ &&
            side_tpb_qe_ == rhs.side_tpb_qe_ &&
            pmt_tpb_thickness_ == rhs.pmt_tpb_thickness_ &&
            endcap_tpb_thickness_ == rhs.endcap_tpb_thickness_ &&
            side_tpb_thickness_ == rhs.side_tpb_thickness_ &&
            tpb_abs_tau_ == rhs.tpb_abs_tau_ &&
            tpb_abs_norm_ == rhs.tpb_abs_norm_ &&
            tpb_abs_scale_ == rhs.tpb_abs_scale_ &&
            mie_gg_ == rhs.mie_gg_ &&
            mie_ratio_ == rhs.mie_ratio_);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm) {
    return bcm.Print(oss);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm) {
    return bcm.Print(oss);
}

I3_SERIALIZABLE(DetectorResponseConfig);
