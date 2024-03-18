#include <icetray/serialization.h>
#include "dataclasses/physics/CCMFP3Summary.h"

std::ostream & CCMFP3Gamma::Print(std::ostream& oss) const {
    oss << "[Gamma:\n"
        << "      StartTime : " << gamma_start_time << '\n'
        << "        EndTime : " << gamma_end_time << '\n'
        << "       PeakTime : " << gamma_peak_time << '\n'
        << "   PeakADCCount : " << gamma_peak_value << '\n'
        << "       Integral : " << gamma_integral << '\n'
        << "     Derivative : " << gamma_derivative << '\n'
        << "  2ndDerivative : " << gamma_second_derivative << '\n'
        << "   LocalAverage : " << gamma_local_average << "\n]";

    return oss;
}

std::ostream& operator<<(std::ostream& oss, CCMFP3Gamma const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, CCMFP3Gamma & bcm) {
    return(bcm.Print(oss));
}

std::ostream & CCMFP3Summary::Print(std::ostream& oss) const {
    oss << "[FP3Summary:\n"
        << "      WaveformLength : " << fp3_waveform_length << '\n'
        << "            Baseline : " << fp3_baseline << '\n'
        << "      BaselineStdDev : " << fp3_baseline_stddev << '\n'
        << "       NumNoisePeaks : " << fp3_num_noise_peaks << '\n'
        << "   AverageNoiseLevel : " << fp3_average_noise_level << '\n'
        << "       MaxNoiseLevel : " << fp3_max_noise_level << '\n'
        << "    NeutronStartTime : " << fp3_neutron_start_time << '\n'
        << "      NeutronEndTime : " << fp3_neutron_end_time << '\n'
        << "   NeutronDerivative : " << fp3_neutron_derivative << '\n'
        << "Neutron2ndDerivative : " << fp3_neutron_second_derivative << '\n'
        << " NeutronLocalAverage : " << fp3_neutron_local_average << '\n'
        << "     NeutronPeakTime : " << fp3_neutron_peak_time << '\n'
        << "    NeutronPeakValue : " << fp3_neutron_peak_value << '\n'
        << "     NeutronIntegral : " << fp3_neutron_integral << '\n'
        << "              Gammas :\n"
        << "              [\n";
        for(CCMFP3Gamma const & gamma : fp3_gammas) {
    oss << "                [Gamma:\n"
        << "                    StartTime : " << gamma.gamma_start_time << '\n'
        << "                      EndTime : " << gamma.gamma_end_time << '\n'
        << "                     PeakTime : " << gamma.gamma_peak_time << '\n'
        << "                 PeakADCCount : " << gamma.gamma_peak_value << '\n'
        << "                     Integral : " << gamma.gamma_integral << '\n'
        << "                   Derivative : " << gamma.gamma_derivative << '\n'
        << "                2ndDerivative : " << gamma.gamma_second_derivative << '\n'
        << "                 LocalAverage : " << gamma.gamma_local_average
        << "                ]\n";
        }
        oss << "              ]\n"
            << "]";

    return oss;
}

std::ostream& operator<<(std::ostream& oss, CCMFP3Summary const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, CCMFP3Summary & bcm) {
    return(bcm.Print(oss));
}

I3_SERIALIZABLE(CCMFP3Gamma);
I3_SERIALIZABLE(CCMFP3Summary);
