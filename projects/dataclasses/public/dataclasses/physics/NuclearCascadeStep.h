#ifndef NuclearCascadeStep_H_INCLUDED
#define NuclearCascadeStep_H_INCLUDED

#include "dataclasses/I3Vector.h"

#include <iostream>

#include <icetray/CCMPMTKey.h>

static const unsigned nuclearcascadestep_version_ = 0;

struct NuclearCascadeStep : public I3FrameObject {
    int initial_level_index;
    double initial_level_energy;
    double initial_level_spin; //2J
    int initial_level_parity; // +1 or -1

    int final_level_index;
    double final_level_energy;
    double final_level_spin;
    int final_level_parity;

    double gamma_energy;

    double T12_ns;  // Half-life in ns
    double tau_ns;  // Mean lifetime in ns

    double sampled_delay_ns;
    double cumulative_time_ns;

    SET_LOGGER("NuclearCascadeStep");

    bool operator==(const NuclearCascadeStep& rhs) const {
        return std::tie(
                    initial_level_index, initial_level_energy, initial_level_spin, initial_level_parity,
                    final_level_index, final_level_energy, final_level_spin, final_level_parity,
                    gamma_energy, T12_ns, tau_ns, sampled_delay_ns, cumulative_time_ns)
               ==
               std::tie(
                    rhs.initial_level_index, rhs.initial_level_energy, rhs.initial_level_spin, rhs.initial_level_parity,
                    rhs.final_level_index, rhs.final_level_energy, rhs.final_level_spin, rhs.final_level_parity,
                    rhs.gamma_energy, rhs.T12_ns, rhs.tau_ns, rhs.sampled_delay_ns, rhs.cumulative_time_ns);
    }

    std::ostream& Print(std::ostream&) const;

    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>nuclearcascadestep_version_)
            log_fatal("Attempting to read version %u from file but running version %u of nuclearcascadestep class.",
                    version,nuclearcascadestep_version_);
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
        ar & make_nvp("initial_level_index", initial_level_index);
        ar & make_nvp("initial_level_energy", initial_level_energy);
        ar & make_nvp("initial_level_spin", initial_level_spin);
        ar & make_nvp("initial_level_parity", initial_level_parity);
        ar & make_nvp("final_level_index", final_level_index);
        ar & make_nvp("final_level_energy", final_level_energy);
        ar & make_nvp("final_level_spin", final_level_spin);
        ar & make_nvp("final_level_parity", final_level_parity);
        ar & make_nvp("gamma_energy", gamma_energy);
        ar & make_nvp("T12_ns", T12_ns);
        ar & make_nvp("tau_ns", tau_ns);
        ar & make_nvp("sampled_delay_ns", sampled_delay_ns);
        ar & make_nvp("cumulative_time_ns", cumulative_time_ns);
    }
};

I3_CLASS_VERSION(NuclearCascadeStep, nuclearcascadestep_version_);

typedef I3Vector<NuclearCascadeStep> NuclearCascadeStepSeries;

std::ostream& operator<<(std::ostream&, const NuclearCascadeStep&);

I3_POINTER_TYPEDEFS(NuclearCascadeStep);
I3_POINTER_TYPEDEFS(NuclearCascadeStepSeries);
#endif
