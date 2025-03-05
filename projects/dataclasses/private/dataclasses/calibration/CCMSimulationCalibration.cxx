#include <dataclasses/calibration/CCMSimulationCalibration.h>
#include <ostream>
#include <icetray/serialization.h>

CCMSimulationCalibration::CCMSimulationCalibration(){}

void CCMSimulationCalibration::SetPMTEfficiencies(I3MapPMTKeyDouble x){
    PMTEfficiencies = x;
}
I3MapPMTKeyDouble CCMSimulationCalibration::GetPMTEfficiencies() const {
    return PMTEfficiencies;
}

void CCMSimulationCalibration::SetLatePulseMu(I3MapPMTKeyDouble x){
    LatePulseMu = x;
}
I3MapPMTKeyDouble CCMSimulationCalibration::GetLatePulseMu() const {
    return LatePulseMu;
}

void CCMSimulationCalibration::SetLatePulseSigma(I3MapPMTKeyDouble x){
    LatePulseSigma = x;
}
I3MapPMTKeyDouble CCMSimulationCalibration::GetLatePulseSigma() const {
    return LatePulseSigma;
}

void CCMSimulationCalibration::SetLatePulseScale(I3MapPMTKeyDouble x){
    LatePulseScale = x;
}
I3MapPMTKeyDouble CCMSimulationCalibration::GetLatePulseScale() const {
    return LatePulseScale;
}

void CCMSimulationCalibration::SetPMTSPEMu(I3MapPMTKeyDouble x){
    PMTSPEMu = x;
}
I3MapPMTKeyDouble CCMSimulationCalibration::GetPMTSPEMu() const {
    return PMTSPEMu;
}

void CCMSimulationCalibration::SetPMTSPESigma(I3MapPMTKeyDouble x){
    PMTSPESigma = x;
}
I3MapPMTKeyDouble CCMSimulationCalibration::GetPMTSPESigma() const {
    return PMTSPESigma;
}

void CCMSimulationCalibration::SetRs(double x){
    Rs = x;
}
double CCMSimulationCalibration::GetRs() const {
    return Rs;
}

void CCMSimulationCalibration::SetRt(double x){
    Rt = x;
}
double CCMSimulationCalibration::GetRt() const {
    return Rt;
}

void CCMSimulationCalibration::SetTauS(double x){
    tau_s = x;
}
double CCMSimulationCalibration::GetTauS() const {
    return tau_s;
}

void CCMSimulationCalibration::SetTauT(double x){
    tau_t = x;
}
double CCMSimulationCalibration::GetTauT() const {
    return tau_t;
}

void CCMSimulationCalibration::SetTauOther(double x){
    tau_other = x;
}
double CCMSimulationCalibration::GetTauOther() const {
    return tau_other;
}

template <class Archive>
void CCMSimulationCalibration::save(Archive& ar, unsigned version) const {
    if (version > ccmsimulationcalibration_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMSimulationCalibration class.", version, ccmsimulationcalibration_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTEfficiencies",PMTEfficiencies);
    ar & make_nvp("LatePulseMu",LatePulseMu);
    ar & make_nvp("LatePulseSigma",LatePulseSigma);
    ar & make_nvp("LatePulseScale",LatePulseScale);
    ar & make_nvp("PMTSPEMu",PMTSPEMu);
    ar & make_nvp("PMTSPESigma",PMTSPESigma);
    ar & make_nvp("Rs",Rs);
    ar & make_nvp("Rt",Rt);
    ar & make_nvp("tau_s",tau_s);
    ar & make_nvp("tau_t",tau_t);
    ar & make_nvp("tau_other",tau_other);
}

template <class Archive>
void CCMSimulationCalibration::load(Archive& ar, unsigned version) {
    if (version > ccmsimulationcalibration_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMSimulationCalibration class.", version, ccmsimulationcalibration_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTEfficiencies",PMTEfficiencies);
    ar & make_nvp("LatePulseMu",LatePulseMu);
    ar & make_nvp("LatePulseSigma",LatePulseSigma);
    ar & make_nvp("LatePulseScale",LatePulseScale);
    ar & make_nvp("PMTSPEMu",PMTSPEMu);
    ar & make_nvp("PMTSPESigma",PMTSPESigma);
    ar & make_nvp("Rs",Rs);
    if (version == 0){
        Rt = 1.0 - Rs;
    } else {
        ar & make_nvp("Rt",Rt);
    }
    ar & make_nvp("tau_s",tau_s);
    ar & make_nvp("tau_t",tau_t);
    ar & make_nvp("tau_other",tau_other);
}


std::ostream& operator<<(std::ostream& os, const CCMSimulationCalibration& pe) {
    return(pe.Print(os));
}

std::ostream& CCMSimulationCalibration::Print(std::ostream& os) const{
    os << "[ CCMSimulationCalibration::"
        << "\n  PMTEfficiencies :" << PMTEfficiencies
        << "\n  LatePulseMu :" << LatePulseMu
        << "\n  LatePulseSigma :" << LatePulseSigma
        << "\n  LatePulseScale :" << LatePulseScale
        << "\n  PMTSPEMu :" << PMTSPEMu
        << "\n  PMTSPESigma :" << PMTSPESigma
        << "\n  Rs :" << Rs
        << "\n  Rt :" << Rt
        << "\n  Tau_s :" << tau_s
        << "\n  Tau_t :" << tau_t
        << "\n  Tau_other :" << tau_other
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(CCMSimulationCalibration);


