#ifndef CCMFP3Summary_H_INCLUDED
#define CCMFP3Summary_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>

static const unsigned ccmfp3gamma_version_ = 1;
class CCMFP3Gamma {
public:
    // The following times are measured in ns from the beginning of the DAQ window

    // The beginning of gamma waveform rise above baseline+noise
    double gamma_start_time;

    // The time at which the gamma waveform falls below the baseline+noise
    double gamma_end_time;

    // The time at which the gamma waveform is at its maximum
    double gamma_peak_time;


    // The maximum value of the gamma waveform
    // Measured in ADC counts above the baseline
    double gamma_peak_value;

    // The integral of the gamma waveform between the beginning and end times
    // i.e. the sum of values in the interval [gamma_start_time, gamma_end_time] multiplied by the bin width
    // Measured in ADC * ns
    double gamma_integral;

    // The derivative of the gamma waveform at the peak
    double gamma_derivative;

    // The second derivative of the gamma waveform at the peak
    double gamma_second_derivative;

    // The local average of the gamma waveform at the peak
    double gamma_local_average;

public:

    CCMFP3Gamma() :
        gamma_start_time(std::numeric_limits<double>::quiet_NaN()),
        gamma_end_time(std::numeric_limits<double>::quiet_NaN()),
        gamma_peak_time(std::numeric_limits<double>::quiet_NaN()),
        gamma_peak_value(std::numeric_limits<double>::quiet_NaN()),
        gamma_integral(std::numeric_limits<double>::quiet_NaN()),
        gamma_derivative(std::numeric_limits<double>::quiet_NaN()),
        gamma_second_derivative(std::numeric_limits<double>::quiet_NaN()),
        gamma_local_average(std::numeric_limits<double>::quiet_NaN())
    {}

    std::ostream& Print(std::ostream&) const;

    bool operator==(const CCMFP3Gamma& rhs) const {
        return std::tie(
                gamma_start_time,
                gamma_end_time,
                gamma_peak_time,
                gamma_peak_value,
                gamma_integral,
                gamma_derivative,
                gamma_second_derivative,
                gamma_local_average
                )
            == std::tie(
                rhs.gamma_start_time,
                rhs.gamma_end_time,
                rhs.gamma_peak_time,
                rhs.gamma_peak_value,
                rhs.gamma_integral,
                rhs.gamma_derivative,
                rhs.gamma_second_derivative,
                rhs.gamma_local_average
                );
    }
private:
    friend class icecube::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned version);
};

static const unsigned ccmfp3summary_version_ = 1;
class CCMFP3Summary : public I3FrameObject {
public:
    size_t fp3_waveform_length;
    double fp3_baseline;
    double fp3_baseline_stddev;

    size_t fp3_num_noise_peaks;
    double fp3_average_noise_level;
    double fp3_max_noise_level;

    double fp3_neutron_start_time;
    double fp3_neutron_end_time;
    double fp3_neutron_derivative;
    double fp3_neutron_second_derivative;
    double fp3_neutron_local_average;
    double fp3_neutron_peak_time;
    double fp3_neutron_peak_value;
    double fp3_neutron_integral;

    std::vector<CCMFP3Gamma> fp3_gammas;
public:
    CCMFP3Summary() :
        fp3_waveform_length(0),
        fp3_baseline(std::numeric_limits<double>::quiet_NaN()),
        fp3_baseline_stddev(std::numeric_limits<double>::quiet_NaN()),
        fp3_num_noise_peaks(0),
        fp3_average_noise_level(std::numeric_limits<double>::quiet_NaN()),
        fp3_max_noise_level(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_start_time(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_end_time(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_derivative(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_second_derivative(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_local_average(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_peak_time(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_peak_value(std::numeric_limits<double>::quiet_NaN()),
        fp3_neutron_integral(std::numeric_limits<double>::quiet_NaN())
    {}

    std::ostream& Print(std::ostream&) const;

    bool operator==(const CCMFP3Summary& rhs) const {
        return std::tie(
                fp3_waveform_length,
                fp3_baseline,
                fp3_baseline_stddev,
                fp3_num_noise_peaks,
                fp3_average_noise_level,
                fp3_max_noise_level,
                fp3_neutron_start_time,
                fp3_neutron_end_time,
                fp3_neutron_derivative,
                fp3_neutron_second_derivative,
                fp3_neutron_local_average,
                fp3_neutron_peak_time,
                fp3_neutron_peak_value,
                fp3_neutron_integral,
                fp3_gammas
                )
            == std::tie(
                rhs.fp3_waveform_length,
                rhs.fp3_baseline,
                rhs.fp3_baseline_stddev,
                rhs.fp3_num_noise_peaks,
                rhs.fp3_average_noise_level,
                rhs.fp3_max_noise_level,
                rhs.fp3_neutron_start_time,
                rhs.fp3_neutron_end_time,
                rhs.fp3_neutron_derivative,
                rhs.fp3_neutron_second_derivative,
                rhs.fp3_neutron_local_average,
                rhs.fp3_neutron_peak_time,
                rhs.fp3_neutron_peak_value,
                rhs.fp3_neutron_integral,
                rhs.fp3_gammas
                );
    }
private:
    friend class icecube::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned version);
};

std::ostream& operator<<(std::ostream& oss, CCMFP3Gamma const & fp3);
std::ostream& operator<<(std::ostream& oss, CCMFP3Gamma & fp3);
std::ostream& operator<<(std::ostream& oss, CCMFP3Summary const & fp3);
std::ostream& operator<<(std::ostream& oss, CCMFP3Summary & fp3);

template <class Archive>
void CCMFP3Gamma::serialize(Archive& ar, unsigned version) {
    if (version>ccmfp3gamma_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMFP3Gamma class.",version,ccmfp3gamma_version_);
    ar & make_nvp("gamma_start_time", gamma_start_time);
    ar & make_nvp("gamma_end_time", gamma_end_time);
    ar & make_nvp("gamma_peak_time", gamma_peak_time);
    ar & make_nvp("gamma_peak_value", gamma_peak_value);
    ar & make_nvp("gamma_integral", gamma_integral);
    ar & make_nvp("gamma_derivative", gamma_derivative);
    ar & make_nvp("gamma_second_derivative", gamma_second_derivative);
    ar & make_nvp("gamma_local_average", gamma_local_average);
}

template <class Archive>
void CCMFP3Summary::serialize(Archive& ar, unsigned version) {
    if (version>ccmfp3summary_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMFP3Summary class.",version,ccmfp3summary_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));

    ar & make_nvp("fp3_waveform_length", fp3_waveform_length);
    ar & make_nvp("fp3_baseline", fp3_baseline);
    ar & make_nvp("fp3_baseline_stddev", fp3_baseline_stddev);
    ar & make_nvp("fp3_num_noise_peaks", fp3_num_noise_peaks);
    ar & make_nvp("fp3_average_noise_level", fp3_average_noise_level);
    ar & make_nvp("fp3_max_noise_level", fp3_max_noise_level);
    ar & make_nvp("fp3_neutron_start_time", fp3_neutron_start_time);
    ar & make_nvp("fp3_neutron_end_time", fp3_neutron_end_time);
    ar & make_nvp("fp3_neutron_derivative", fp3_neutron_derivative);
    ar & make_nvp("fp3_neutron_second_derivative", fp3_neutron_second_derivative);
    ar & make_nvp("fp3_neutron_local_average", fp3_neutron_local_average);
    ar & make_nvp("fp3_neutron_peak_time", fp3_neutron_peak_time);
    ar & make_nvp("fp3_neutron_peak_value", fp3_neutron_peak_value);
    ar & make_nvp("fp3_neutron_integral", fp3_neutron_integral);
    ar & make_nvp("fp3_gammas", fp3_gammas);
}

I3_CLASS_VERSION(CCMFP3Gamma, ccmfp3gamma_version_);
I3_POINTER_TYPEDEFS(CCMFP3Gamma);
I3_DEFAULT_NAME(CCMFP3Gamma);

I3_CLASS_VERSION(CCMFP3Summary, ccmfp3summary_version_);
I3_POINTER_TYPEDEFS(CCMFP3Summary);
I3_DEFAULT_NAME(CCMFP3Summary);

#endif // CCMFP3Summary_H_INCLUDED
