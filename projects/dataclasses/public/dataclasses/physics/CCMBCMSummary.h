#ifndef CCMBCMSummary_H_INCLUDED
#define CCMBCMSummary_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>

static const unsigned ccmbcmsummary_version_ = 1;
class CCMBCMSummary : public I3FrameObject {
public:
    // The following times are measured in ns from the beginning of the DAQ window

    // The beginning of BCM waveform rise above baseline+noise
    double bcm_start_time;

    // The time at which the BCM waveform falls below the baseline+noise
    double bcm_end_time;

    // The time at which the BCM waveform is at its maximum
    double bcm_peak_time;

    // The maximum value of the BCM waveform
    // Measured in ADC counts above the baseline
    double bcm_peak_value;

    // The integral of the BCM waveform between the beginning and end times
    // i.e. the sum of values in the interval [bcm_start_time, bcm_end_time] multiplied by the bin width
    // Measured in ADC * ns
    double bcm_integral;

    // The baseline value of the BCM waveform before the beam pulse
    // Measured in ADC counts
    double bcm_baseline;

    // The standard deviation of the BCM baseline in the measured region (from 2000ns before the peak to the peak)
    // Actually computed as the median absolute deviation from the median
    // Measured in ADC counts
    double bcm_baseline_stddev;
public:
    CCMBCMSummary() :
        bcm_start_time(std::numeric_limits<double>::quiet_NaN()),
        bcm_end_time(std::numeric_limits<double>::quiet_NaN()),
        bcm_peak_time(std::numeric_limits<double>::quiet_NaN()),
        bcm_peak_value(std::numeric_limits<double>::quiet_NaN()),
        bcm_integral(std::numeric_limits<double>::quiet_NaN()),
        bcm_baseline(std::numeric_limits<double>::quiet_NaN()),
        bcm_baseline_stddev(std::numeric_limits<double>::quiet_NaN())
    {}

    bool operator==(const CCMBCMSummary& rhs) const {
        return std::tie(
                bcm_start_time,
                bcm_end_time,
                bcm_peak_time,
                bcm_peak_value,
                bcm_integral,
                bcm_baseline,
                bcm_baseline_stddev)
            == std::tie(
                rhs.bcm_start_time,
                rhs.bcm_end_time,
                rhs.bcm_peak_time,
                rhs.bcm_peak_value,
                rhs.bcm_integral,
                rhs.bcm_baseline,
                rhs.bcm_baseline_stddev);
    }
private:
    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

template <class Archive>
void CCMBCMSummary::save(Archive& ar, unsigned version) const {
    if (version>ccmbcmsummary_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMBCMSummary class.",version,ccmbcmsummary_version_);
    ar & make_nvp("bcm_start_time", bcm_start_time);
    ar & make_nvp("bcm_end_time", bcm_end_time);
    ar & make_nvp("bcm_peak_time", bcm_peak_time);
    ar & make_nvp("bcm_peak_value", bcm_peak_value);
    ar & make_nvp("bcm_integral", bcm_integral);
    ar & make_nvp("bcm_baseline", bcm_baseline);
    ar & make_nvp("bcm_baseline_stddev", bcm_baseline_stddev);
}

template <class Archive>
void CCMBCMSummary::load(Archive& ar, unsigned version) {
    if (version>ccmbcmsummary_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMBCMSummary class.",version,ccmbcmsummary_version_);
    ar & make_nvp("bcm_start_time", bcm_start_time);
    ar & make_nvp("bcm_end_time", bcm_end_time);
    ar & make_nvp("bcm_peak_time", bcm_peak_time);
    ar & make_nvp("bcm_peak_value", bcm_peak_value);
    ar & make_nvp("bcm_integral", bcm_integral);
    ar & make_nvp("bcm_baseline", bcm_baseline);
    ar & make_nvp("bcm_baseline_stddev", bcm_baseline_stddev);
}

I3_CLASS_VERSION(CCMBCMSummary, ccmbcmsummary_version_);
I3_POINTER_TYPEDEFS(CCMBCMSummary);
I3_DEFAULT_NAME(CCMBCMSummary);

#endif // CCMBCMSummary_H_INCLUDED
