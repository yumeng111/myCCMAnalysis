#include "simclasses/DetectorResponseConfig.h"

#include "icetray/I3FrameObject.h"

std::ostream & DetectorResponseConfig::Print(std::ostream& oss) const {
    oss << "DetectorResponseConfig: "
        << "\n  Rayleigh Scattering Length :" << rayleigh_scattering_length_ 
        << "\n  PMT TPB QE :" << pmt_tpb_qe_ 
        << "\n  Endcap TPB QE :" << endcap_tpb_qe_ 
        << "\n  Side TPB QE :" << side_tpb_qe_ 
        << "\n  PMT TPB Thickness :" << pmt_tpb_thickness_ 
        << "\n  Endcap TPB Thickness :" << endcap_tpb_thickness_ 
        << "\n  Side TPB Thickness :" << side_tpb_thickness_;
    return oss;
}

bool DetectorResponseConfig::operator==(const DetectorResponseConfig& rhs) const {
    return (rayleigh_scattering_length_ == rhs.rayleigh_scattering_length_ &&
            pmt_tpb_qe_ == rhs.pmt_tpb_qe_ &&
            endcap_tpb_qe_ == rhs.endcap_tpb_qe_ &&
            side_tpb_qe_ == rhs.side_tpb_qe_ &&
            pmt_tpb_thickness_ == rhs.pmt_tpb_thickness_ &&
            endcap_tpb_thickness_ == rhs.endcap_tpb_thickness_ &&
            side_tpb_thickness_ == rhs.side_tpb_thickness_);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm) {
    return bcm.Print(oss);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm) {
    return bcm.Print(oss);
}

I3_SERIALIZABLE(DetectorResponseConfig);
