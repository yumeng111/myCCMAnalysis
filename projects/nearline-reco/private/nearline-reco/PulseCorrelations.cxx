#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <limits>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>

#include <boost/histogram.hpp>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMPMTKeyPair.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <icetray/I3Int.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3String.h>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

namespace {
typedef std::tuple<CCMPMTKey, CCMRecoPulse> PMTKeyPulsePair;
typedef std::vector<PMTKeyPulsePair> PMTKeyPulseVector;

class KahanAccumulator {
    double sum = 0;
    double c = 0;
public:
    void Add(double x) {
        double const y = x - c;
        double const t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    operator double() const { return sum; }
    operator double&() { return sum; }
    double operator+=(double x) { Add(x); return sum; }
    double operator-=(double x) { Add(-x); return sum; }
    double Sum() const { return sum; }
    double Compensation() const { return c; }
};

/*
def online_weighted_covariance(data1, data2, data3):
    meanx = meany = 0
    wsum = wsum2 = 0
    C = 0
    for x, y, w in zip(data1, data2, data3):
        wsum += w
        wsum2 += w * w
        dx = x - meanx
        meanx += (w / wsum) * dx
        meany += (w / wsum) * (y - meany)
        C += w * dx * (y - meany)

    population_covar = C / wsum
    # Bessel's correction for sample variance
    # Frequency weights
    sample_frequency_covar = C / (wsum - 1)
    # Reliability weights
    sample_reliability_covar = C / (wsum - wsum2 / wsum)
*/

class OnlineCovariance {
    KahanAccumulator meanx;
    KahanAccumulator meany;
    KahanAccumulator wsum;
    KahanAccumulator wsum2;
    KahanAccumulator C;
public:
    void AddValue(double x, double y, double w) {
        wsum += w;
        wsum2 += w * w;
        double dx = x - meanx;
        meanx += (w / wsum) * dx;
        meany += (w / wsum) * (y - meany);
        C += w * dx * (y - meany);
    }
    double Covariance() {
        return C / wsum;
    }
    double SampleCovariance() {
        return C / (wsum - 1);
    }
    double ReliabilityCovar() {
        return C / (wsum - wsum2 / wsum);
    }
    void Add(OnlineCovariance const & other) {
        C += other.C + (other.meanx - meanx) * (other.meany - meany) * wsum * other.wsum / (wsum + other.wsum);
        wsum += other.wsum;
        wsum2 += other.wsum2;
        meanx += (other.wsum / wsum) * (other.meanx - meanx);
        meany += (other.wsum / wsum) * (other.meany - meany);
    }
    template<typename Container>
    void Fill(Container & c) const {
        c.push_back(meanx.Sum());
        c.push_back(meanx.Compensation());
        c.push_back(meany.Sum());
        c.push_back(meany.Compensation());
        c.push_back(wsum.Sum());
        c.push_back(wsum.Compensation());
        c.push_back(wsum2.Sum());
        c.push_back(wsum2.Compensation());
        c.push_back(C.Sum());
        c.push_back(C.Compensation());
     }
};
}

class PulseCorrelations: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    std::string pulses_name_;
    std::string output_prefix_;

    double time_window_ns_;
    bool pulse_correlations_;

    I3Map<std::tuple<CCMPMTKey, CCMPMTKey>, OnlineCovariance> pulse_covariance_;
    I3Map<std::tuple<CCMPMTKey, CCMPMTKey>, OnlineCovariance> charge_covariance_;

public:
    void Geometry(I3FramePtr frame);
    PulseCorrelations(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(PulseCorrelations);

std::tuple<std::tuple<CCMPMTKey, CCMPMTKey>, std::tuple<double, double>> present_key_pair(CCMPMTKey const & present, CCMPMTKey const & not_present) {
    if(present < not_present) {
        return {{present, not_present}, {1.0, 0.0}};
    } else {
        return {{not_present, present}, {0.0, 1.0}};
    }
}

std::tuple<std::tuple<CCMPMTKey, CCMPMTKey>, std::tuple<double, double>> charge_key_pair(PMTKeyPulsePair const & p0, PMTKeyPulsePair const & p1) {
    CCMPMTKey const & k0 = std::get<0>(p0);
    CCMPMTKey const & k1 = std::get<0>(p1);
    double c0 = std::get<1>(p0).GetCharge();
    double c1 = std::get<1>(p1).GetCharge();
    if(k0 < k1) {
        return {{k0, k1}, {c0, c1}};
    } else {
        return {{k1, k0}, {c1, c0}};
    }
}

PulseCorrelations::PulseCorrelations(const I3Context& context) : I3Module(context),
    geo_seen(false), geometry_name_("") {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("InputPulsesName", "Name of the input pulses", std::string("WavedeformPulses"));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
    AddParameter("TimeWindow", "Time window for considering correlations in ns", double(10.0));
    AddParameter("PulseCorrelations", "ComputePulseCorrelations", bool(false));
}

void PulseCorrelations::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("InputPulsesName", pulses_name_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("TimeWindow", time_window_ns_);
    GetParameter("PulseCorrelations", pulse_correlations_);
}

void PulseCorrelations::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    if(geo_seen) {
    }
    PushFrame(frame);
}

void PulseCorrelations::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_name_);
    if(not pulses) {
        log_fatal("Could not find %s in the DAQ frame.", pulses_name_.c_str());
    }

    PMTKeyPulseVector pulse_list;

    for(CCMRecoPulseSeriesMap::const_iterator i = pulses->begin(); i != pulses->end(); i++) {
        for(CCMRecoPulse const & pulse: i->second) {
            pulse_list.push_back(PMTKeyPulsePair(i->first, pulse));
        }
    }

    std::sort(pulse_list.begin(), pulse_list.end(), [](auto const & t0, auto const & t1){return std::get<1>(t0).GetTime() < std::get<1>(t1).GetTime();});

    std::set<std::tuple<CCMPMTKey, CCMPMTKey>> pending;
    std::map<CCMPMTKey, double> time_of_last_pulse;
    std::map<CCMPMTKey, size_t> pulses_in_window;

    for(PMTKeyPulseVector::const_iterator i = pulse_list.begin(), j = pulse_list.begin(); j != pulse_list.end(); j++) {
        // Skip over pulses not in the window anymore
        while(std::get<1>(*j).GetTime() - std::get<1>(*i).GetTime() > time_window_ns_) {
            CCMPMTKey const & i_key = std::get<0>(*i);
            pulses_in_window[i_key] -= 1;
            i++;
        }

        CCMPMTKey const & j_key = std::get<0>(*j);
        float const & j_time = std::get<1>(*j).GetTime();

        // Add new pulse to the window and record its time
        if(time_of_last_pulse.find(j_key) == time_of_last_pulse.end()) {
            pulses_in_window[j_key] = 1;
        } else {
            if(pulse_correlations_) {
                double previous_j_time = time_of_last_pulse[j_key];
                if(j_time - previous_j_time > time_window_ns_) {
                    // Add entries for pulses not present in the window
                    for(std::pair<CCMPMTKey const, double> const & p : time_of_last_pulse) {
                        CCMPMTKey const & p_key = p.first;
                        std::tuple<std::tuple<CCMPMTKey, CCMPMTKey>, std::tuple<double, double>> kp = present_key_pair(j_key, p_key);
                        std::tuple<CCMPMTKey, CCMPMTKey> const & key = std::get<0>(kp);
                        std::set<std::tuple<CCMPMTKey, CCMPMTKey>>::iterator it = pending.find(key);
                        if(it != pending.end())
                            pending.erase(it);
                        if(pulses_in_window[p_key] == 0) {
                            double p_time = time_of_last_pulse[p_key];
                            if(p_time - previous_j_time < time_window_ns_)
                                continue;
                            std::tuple<double, double> const & present_not_present = std::get<1>(kp);
                            double weight = std::min(j_time - p_time, p_time - previous_j_time) / time_window_ns_;
                            pulse_covariance_[key].AddValue(std::get<0>(present_not_present), std::get<1>(present_not_present), weight);
                        }
                    }
                } else {
                    for(std::pair<CCMPMTKey const, double> const & p : time_of_last_pulse) {
                        CCMPMTKey const & p_key = p.first;
                        if(pulses_in_window[p_key] == 0) {
                            std::tuple<std::tuple<CCMPMTKey, CCMPMTKey>, std::tuple<double, double>> kp = present_key_pair(j_key, p_key);
                            std::tuple<CCMPMTKey, CCMPMTKey> const & key = std::get<0>(kp);
                            pending.insert(key);
                        }
                    }
                }
            }

            pulses_in_window[j_key] += 1;
        }
        time_of_last_pulse[j_key] = std::get<1>(*j).GetTime();

        // Add entries for pulse pairs
        for(PMTKeyPulseVector::const_iterator k = i; k != j; ++k) {
            std::tuple<std::tuple<CCMPMTKey, CCMPMTKey>, std::tuple<double, double>> kp = charge_key_pair(*j, *k);
            std::tuple<CCMPMTKey, CCMPMTKey> const & key = std::get<0>(kp);
            std::tuple<double, double> const & charge = std::get<1>(kp);
            charge_covariance_[key].AddValue(std::get<0>(charge), std::get<1>(charge), 1.0);
            if(pulse_correlations_)
                pulse_covariance_[key].AddValue(1.0, 1.0, 1.0);
        }
    }

    // Iterate over pending pairs and add them as empty entries to the covariance matrix
    if(pulse_correlations_)
        for(std::tuple<CCMPMTKey, CCMPMTKey> const & key : pending) {
            CCMPMTKey const & j_key = std::get<0>(key);
            CCMPMTKey const & p_key = std::get<1>(key);
            double const & j_time = time_of_last_pulse.at(j_key);
            double const & p_time = time_of_last_pulse.at(p_key);
            std::tuple<std::tuple<CCMPMTKey, CCMPMTKey>, std::tuple<double, double>> kp;
            if(p_time < j_time) {
                kp = present_key_pair(j_key, p_key);
            } else {
                kp = present_key_pair(p_key, j_key);
            }
            present_key_pair(j_key, p_key);
            //std::tuple<CCMPMTKey, CCMPMTKey> const & key = std::get<0>(kp);
            std::tuple<double, double> const & present_not_present = std::get<1>(kp);
            double weight = std::abs(j_time - p_time) / time_window_ns_;
            pulse_covariance_[key].AddValue(std::get<0>(present_not_present), std::get<1>(present_not_present), weight);
        }

    PushFrame(frame);
}

void PulseCorrelations::Finish() {
    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::Physics);

    if(pulse_correlations_) {
        boost::shared_ptr<I3Map<CCMPMTKeyPair, std::vector<double>>> pulse_covariance_out = boost::make_shared<I3Map<CCMPMTKeyPair, std::vector<double>>>();
        for(I3Map<std::tuple<CCMPMTKey, CCMPMTKey>, OnlineCovariance>::const_iterator i = pulse_covariance_.begin(); i != pulse_covariance_.end(); i++) {
            std::tuple<CCMPMTKey, CCMPMTKey> const & key = i->first;
            CCMPMTKeyPair key_(std::get<0>(key), std::get<1>(key));
            OnlineCovariance const & cov = i->second;
            pulse_covariance_out->operator[](key_) = std::vector<double>();
            std::vector<double> & c = pulse_covariance_out->operator[](key_);
            cov.Fill(c);
        }
        frame->Put(output_prefix_ + "PulseCovariance", pulse_covariance_out);
    }

    boost::shared_ptr<I3Map<CCMPMTKeyPair, std::vector<double>>> charge_covariance_out = boost::make_shared<I3Map<CCMPMTKeyPair, std::vector<double>>>();
    for(I3Map<std::tuple<CCMPMTKey, CCMPMTKey>, OnlineCovariance>::const_iterator i = charge_covariance_.begin(); i != charge_covariance_.end(); i++) {
        std::tuple<CCMPMTKey, CCMPMTKey> const & key = i->first;
        CCMPMTKeyPair key_(std::get<0>(key), std::get<1>(key));
        OnlineCovariance const & cov = i->second;
        charge_covariance_out->operator[](key_) = std::vector<double>();
        std::vector<double> & c = charge_covariance_out->operator[](key_);
        cov.Fill(c);
    }
    frame->Put(output_prefix_ + "ChargeCovariance", charge_covariance_out);

    PushFrame(frame);
}
