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

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <icetray/I3Int.h>
#include <dataclasses/physics/CCMEventHeader.h>
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

typedef std::tuple<CCMPMTKey, CCMRecoPulse> PMTKeyPulsePair;
typedef std::vector<PMTKeyPulsePair> PMTKeyPulseVector;

namespace {
// Template function for computing the sum of a container
template<typename T>
double sum(T const & data) {
    return std::accumulate(data.begin(), data.end(), 0.0);
}
}

template<typename T>
void compute_charge_unfolding(std::vector<double> const & data, double beta_s, double beta_t, double alpha, double charge_target, T & total, T & singlet, T & triplet) {
    total[0] = data[0];

    double running_singlet = alpha * data[0];
    double running_triplet = (1 - alpha) * data[0];

    singlet[0] = running_singlet;
    triplet[0] = running_triplet;

    std::deque<double> singlet_window;
    std::deque<double> triplet_window;
    std::deque<double> data_window;

    singlet_window.push_back(running_singlet);
    triplet_window.push_back(running_triplet);
    data_window.push_back(data[0]);

    for(size_t i=1; i<data.size(); ++i) {
        running_singlet *= beta_s;
        running_triplet *= beta_t;

        if(data[i] > 0) {
            double S_pred = sum(singlet_window);
            double T_pred = sum(triplet_window);
            double Data = sum(data_window);

            double dLdgamma = 1.0 -
                (1.0 / beta_t - 1.0)
                /
                (std::pow(1.0 / beta_t, data_window.size() + 1.0) - 1.0);

            double gamma = Data - (S_pred + T_pred) / dLdgamma;
            gamma = std::max(0.0, gamma);

            double k = data[i];
            double b = running_singlet + running_triplet;
            double c = std::pow(beta_t, data_window.size() + 1.0);

            double new_charge = k - b - gamma * c;
            new_charge = std::max(0.0, new_charge);

            total[i] = new_charge;
            running_singlet += alpha * total[i];
            running_triplet += (1 - alpha) * total[i];
            singlet[i] = running_singlet;
            triplet[i] = running_triplet;
        }

        singlet_window.push_back(running_singlet);
        triplet_window.push_back(running_triplet);
        data_window.push_back(data[i]);

        if(i+1 == data.size())
            continue;
        double next_data = data[i+1];
        double target = charge_target - next_data;
        double lost = 0.0;
        double tot_charge = sum(data_window);
        while(true) {
            if(data_window.size() <= 1)
                break;
            if(tot_charge - data_window.front() < target)
                break;

            lost = triplet_window.front();
            tot_charge -= data_window.front();

            singlet_window.pop_front();
            triplet_window.pop_front();
            data_window.pop_front();
        }
        if(lost > 0.0) {
            for(size_t j=0; j<data_window.size(); ++j) {
                lost *= beta_t;
                singlet_window[j] -= lost;
                triplet_window[j] -= lost;
            }
        }
    }
}

class ChargeUnfolding: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;

    bool split_by_pmt_;

    double tau_s_;
    double tau_t_;
    double singlet_triplet_ratio_;
    double delta_t_;

    double charge_target_;

    double alpha_;
    double beta_s_;
    double beta_t_;

    std::string input_prefix_;
    std::string output_prefix_;

    bool check_masked_pulses_;
    bool check_raw_pulses_;
    std::string raw_pulses_name_;
    std::string pulses_mask_name_;

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

public:
    void Geometry(I3FramePtr frame);
    ChargeUnfolding(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
};

I3_MODULE(ChargeUnfolding);

ChargeUnfolding::ChargeUnfolding(const I3Context& context) : I3Module(context),
    geo_seen(false), geometry_name_("") {
    I3Vector<double> default_time_windows;
    default_time_windows.push_back(90.0);
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("SplitByPMT", "Split the output by PMT", false);
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("TauSinglet", "Time constant for singlet light", 8.13);
    AddParameter("TauTriplet", "Time constant for triplet light", 743);
    AddParameter("SingletTripletRatio", "Ratio of singlet to triplet light", 0.33);
    AddParameter("ChargeTarget", "Target charge for  triplet light inference", 1.0);
    AddParameter("BinWidth", "Width of the time bins", 2.0);
    AddParameter("InputPulsesMaskName", "Name of the input pulses mask", std::string(""));
    AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
    AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
}

void ChargeUnfolding::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("SplitByPMT", split_by_pmt_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("TauSinglet", tau_s_);
    GetParameter("TauTriplet", tau_t_);
    GetParameter("SingletTripletRatio", singlet_triplet_ratio_);
    GetParameter("ChargeTarget", charge_target_);
    GetParameter("BinWidth", delta_t_);
    GetParameter("InputPulsesMaskName", pulses_mask_name_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);

    double norm_s = singlet_triplet_ratio_ / tau_s_;
    double norm_t = (1.0 - singlet_triplet_ratio_) / tau_t_;
    alpha_ = norm_s / (norm_s + norm_t);
    beta_s_ = exp(-delta_t_ / tau_s_);
    beta_t_ = exp(-delta_t_ / tau_t_);

    check_masked_pulses_ = false;
    check_raw_pulses_ = (raw_pulses_name_ != "");
    check_masked_pulses_ = (pulses_mask_name_ != "");
    if(not check_masked_pulses_) {
        pulses_mask_name_ = input_prefix_ + "EventPulses";
        if(not check_raw_pulses_) {
            check_masked_pulses_ = true;
        }
    }
}

void ChargeUnfolding::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    pmt_keys.clear();
    if(geo_seen) {
        std::set<CCMOMGeo::OMType> allowed_pmt_types(pmt_types.begin(), pmt_types.end());
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            if(allowed_pmt_types.count(p.second.omtype) == 0)
                continue;
            pmt_keys.insert(p.first);
        }
    }
    PushFrame(frame);
}

void ChargeUnfolding::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    bool raw_pulses = false;
    CCMRecoPulseSeriesMapConstPtr pulses;
    if(check_raw_pulses_) {
        pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(raw_pulses_name_);
        raw_pulses = bool(pulses);
    }
    if(check_masked_pulses_) {
        CCMRecoPulseSeriesMapMaskConstPtr mask = frame->Get<CCMRecoPulseSeriesMapMaskConstPtr>(pulses_mask_name_);
        if(mask) {
            raw_pulses = false;
            pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_mask_name_);
        }
    }

    if(not pulses) {
        std::stringstream ss;
        ss << "Could not find ";
        if(check_raw_pulses_) {
            ss << raw_pulses_name_;
            if(check_masked_pulses_) {
                ss << " or ";
            }
        }
        if(check_masked_pulses_) {
            ss << pulses_mask_name_;
        }
        ss << " in the DAQ frame.";
        log_fatal("%s", ss.str().c_str());
        PushFrame(frame);
        return;
    }

    PMTKeyPulseVector pulse_list;
    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin();
            i != pulses->end(); i++) {
        if(pmt_keys.count(i->first) == 0)
            continue;
        for(CCMRecoPulse const & pulse: i->second) {
            pulse_list.push_back(PMTKeyPulsePair(i->first, pulse));
        }
    }

    std::sort(pulse_list.begin(), pulse_list.end(), [](auto const & t0, auto const & t1){return std::get<1>(t0).GetTime() < std::get<1>(t1).GetTime();});

    double min_time = std::get<1>(pulse_list.front()).GetTime();
    double max_time = std::get<1>(pulse_list.back()).GetTime();
    min_time = std::floor(min_time / delta_t_) * delta_t_;
    max_time = std::ceil(max_time / delta_t_) * delta_t_;
    double end_time = max_time + delta_t_;

    std::vector<double> data((max_time - min_time) / delta_t_, 0.0);
    for (PMTKeyPulseVector::const_iterator i = pulse_list.begin(); i != pulse_list.end(); ++i) {
        double time = std::get<1>(*i).GetTime();
        size_t bin = (time - min_time) / delta_t_;
        data[bin] += std::get<1>(*i).GetCharge();
    }

    I3VectorDoublePtr raw = boost::make_shared<I3VectorDouble>(data.begin(), data.end());
    I3VectorDoublePtr total = boost::make_shared<I3VectorDouble>(data.size(), 0.0);

    I3VectorDoublePtr singlet = boost::make_shared<I3VectorDouble>(data.size(), 0.0);
    I3VectorDoublePtr triplet = boost::make_shared<I3VectorDouble>(data.size(), 0.0);

    compute_charge_unfolding(data, beta_s_, beta_t_, alpha_, charge_target_, *total, *singlet, *triplet);

    frame->Put(output_prefix_ + "ChargeUnfoldingMinTime", boost::make_shared<I3Double>(min_time));

    frame->Put(output_prefix_ + "ChargeUnfoldingTotalRaw", raw);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotal", total);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalSinglet", singlet);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalTriplet", triplet);

    if(not split_by_pmt_) {
        PushFrame(frame);
        return;
    }

    I3MapPMTKeyVectorDoublePtr raw_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr total_map = boost::make_shared<I3MapPMTKeyVectorDouble>();

    I3MapPMTKeyVectorDoublePtr singlet_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr triplet_map = boost::make_shared<I3MapPMTKeyVectorDouble>();

    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin();
            i != pulses->end(); i++) {
        if(pmt_keys.count(i->first) == 0)
            continue;
        data = std::vector<double>((max_time - min_time) / delta_t_, 0.0);
        for(CCMRecoPulse const & pulse : i->second) {
            double time = pulse.GetTime();
            size_t bin = (time - min_time) / delta_t_;
            data[bin] += pulse.GetCharge();
        }

        raw_map->insert(std::make_pair(i->first, std::vector<double>(data.begin(), data.end())));
        total_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        singlet_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        triplet_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        compute_charge_unfolding(data, beta_s_, beta_t_, alpha_, charge_target_, total_map->at(i->first), singlet_map->at(i->first), triplet_map->at(i->first));
    }

    frame->Put(output_prefix_ + "ChargeUnfoldingRaw", raw_map);
    frame->Put(output_prefix_ + "ChargeUnfolding", total_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingSinglet", singlet_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingTriplet", triplet_map);

    PushFrame(frame);
}

