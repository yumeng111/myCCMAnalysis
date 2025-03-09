/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef CCMSimulationCalibration_H_INCLUDED
#define CCMSimulationCalibration_H_INCLUDED

#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <dataclasses/Utility.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Vector.h>
#include <icetray/I3FrameObject.h>

#include <string>
#include <iostream>
#include <sstream>

#include <icetray/CCMPMTKey.h>

static const unsigned ccmsimulationcalibration_version_ = 2;

class CCMSimulationCalibration: public I3FrameObject {
public:
    // things we apply to our simulation after the fact
    I3MapPMTKeyDouble PMTEfficiencies;
    I3MapPMTKeyDouble LatePulseMu;
    I3MapPMTKeyDouble LatePulseSigma;
    I3MapPMTKeyDouble LatePulseScale;
    I3MapPMTKeyDouble PMTSPEMu;
    I3MapPMTKeyDouble PMTSPESigma;
    double Rs;
    double Rt;
    double tau_s;
    double tau_t;
    double tau_other;

    double uv_absorption_a;
    double uv_absorption_b;
    double uv_absorption_d;
    double uv_absorption_scaling;

    SET_LOGGER("CCMSimulationCalibration");

    bool operator==(const CCMSimulationCalibration& rhs) const {
        return PMTEfficiencies == rhs.PMTEfficiencies
            && LatePulseMu == rhs.LatePulseMu
            && LatePulseSigma == rhs.LatePulseSigma
            && LatePulseScale == rhs.LatePulseScale
            && PMTSPEMu == rhs.PMTSPEMu
            && PMTSPESigma == rhs.PMTSPESigma
            && Rs == rhs.Rs
            && Rt == rhs.Rt
            && tau_s == rhs.tau_s
            && tau_t == rhs.tau_t
            && tau_other == rhs.tau_other
            && uv_absorption_a == rhs.uv_absorption_a
            && uv_absorption_b == rhs.uv_absorption_b
            && uv_absorption_d == rhs.uv_absorption_d
            && uv_absorption_scaling == rhs.uv_absorption_scaling;
    }

    CCMSimulationCalibration();

    std::ostream& Print(std::ostream&) const;

    void SetPMTEfficiencies(I3MapPMTKeyDouble x);
    I3MapPMTKeyDouble GetPMTEfficiencies() const;

    void SetLatePulseMu(I3MapPMTKeyDouble x);
    I3MapPMTKeyDouble GetLatePulseMu() const;

    void SetLatePulseSigma(I3MapPMTKeyDouble x);
    I3MapPMTKeyDouble GetLatePulseSigma() const;

    void SetLatePulseScale(I3MapPMTKeyDouble x);
    I3MapPMTKeyDouble GetLatePulseScale() const;

    void SetPMTSPEMu(I3MapPMTKeyDouble x);
    I3MapPMTKeyDouble GetPMTSPEMu() const;

    void SetPMTSPESigma(I3MapPMTKeyDouble x);
    I3MapPMTKeyDouble GetPMTSPESigma() const;

    void SetRs(double x);
    double GetRs() const;
    
    void SetRt(double x);
    double GetRt() const;

    void SetTauS(double x);
    double GetTauS() const;

    void SetTauT(double x);
    double GetTauT() const;

    void SetTauOther(double x);
    double GetTauOther() const;

    void SetUVAbsorptionA(double x);
    double GetUVAbsorptionA() const;

    void SetUVAbsorptionB(double x);
    double GetUVAbsorptionB() const;

    void SetUVAbsorptionD(double x);
    double GetUVAbsorptionD() const;

    void SetUVAbsorptionScaling(double x);
    double GetUVAbsorptionScaling() const;

friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

};

I3_CLASS_VERSION(CCMSimulationCalibration,ccmsimulationcalibration_version_);

std::ostream& operator<<(std::ostream&, const CCMSimulationCalibration&);

I3_POINTER_TYPEDEFS(CCMSimulationCalibration);
#endif
