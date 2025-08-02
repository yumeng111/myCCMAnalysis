#include "dataclasses/physics/NuclearCascadeStep.h"

#include <string>
#include <ostream>
#include <sstream>

#include "icetray/I3Units.h"

I3_SERIALIZABLE(NuclearCascadeStepSeries);

std::ostream& operator<<(std::ostream& os, const NuclearCascadeStep& pe) {
    return(pe.Print(os));
}

std::ostream& NuclearCascadeStep::Print(std::ostream& os) const {
    os << "[ Nuclear Cascade Step::\n";
    os << "  From Level [" << (initial_level_index >= 0 ? std::to_string(initial_level_index) : "UNKNOWN") << "] "
        << initial_level_energy / I3Units::keV << " keV, "
        << "J=" << initial_level_spin
        << (initial_level_parity > 0 ? "+" : "-") << "\n";

    os << "  To Level [" << (final_level_index >= 0 ? std::to_string(final_level_index) : "UNKNOWN") << "] "
        << final_level_energy / I3Units::keV << " keV, "
        << "J=" << final_level_spin
        << (final_level_parity > 0 ? "+" : "-") << "\n";

    os << "  Gamma Energy = " << gamma_energy / I3Units::keV << " keV\n";
    os << "  T1/2 = " << T12_ns << " ns\n";
    os << "  Tau = " << tau_ns << " ns\n";
    os << "  Sampled Delay = " << sampled_delay_ns << " ns\n";
    os << "  Cumulative time = " << cumulative_time_ns << " ns\n";
    os << " ]";
    return os;
}


