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

template<typename T>
void compute_charge_unfolding(std::vector<double> const & data, double beta_s, double beta_t, double alpha, T & total, T & running) {
    total[0] = data[0];
    running[0] = data[0];

    double running_singlet = alpha * data[0];
    double running_triplet = (1 - alpha) * data[0];

    for(size_t i=1; i<data.size(); ++i) {
        running_singlet *= beta_s;
        running_triplet *= beta_t;

        total[i] = data[i] - running_singlet - running_triplet;

        running_singlet += alpha * total[i];
        running_triplet += (1 - alpha) * total[i];
        running[i] = running_singlet + running_triplet;
    }
}

template<typename T>
void compute_charge_unfolding_positive(std::vector<double> const & data, double beta_s, double beta_t, double alpha, T & total, T & running) {
    total[0] = data[0];

    double running_singlet = alpha * data[0];
    double running_triplet = (1 - alpha) * data[0];

    for(size_t i=1; i<data.size(); ++i) {
        running_singlet *= beta_s;
        running_triplet *= beta_t;

        double running_total = running_singlet + running_triplet;

        total[i] = data[i] - running_total;
        total[i] = std::max(0.0, total[i]);
        running_singlet += alpha * total[i];
        running_triplet += (1 - alpha) * total[i];
        running[i] = running_singlet + running_triplet;
    }
}

template<typename T>
void compute_charge_unfolding_bayesian(std::vector<double> const & data, double beta_s, double beta_t, double alpha, T & total, T & running) {
    total[0] = data[0];

    double running_singlet = alpha * data[0];
    double running_triplet = (1 - alpha) * data[0];

    for(size_t i=1; i<data.size(); ++i) {
        running_singlet *= beta_s;
        running_triplet *= beta_t;

        double a = data[i] + running_singlet + running_triplet;
        double b = 2;
        double r = a;
        double p = b / (1.0 + b);
        double mode;

        if(r > 1) {
            mode = (r - 1) * (1 - p) / p;
        } else {
            mode = 0;
        }

        total[i] = data[i] - mode;
        running_singlet += alpha * total[i];
        running_triplet += (1 - alpha) * total[i];
        running[i] = running_singlet + running_triplet;
    }
}

template<typename T>
void compute_charge_unfolding_bayesian_positive(std::vector<double> const & data, double beta_s, double beta_t, double alpha, T & total, T & running) {
    total[0] = data[0];

    double running_singlet = alpha * data[0];
    double running_triplet = (1 - alpha) * data[0];
    running[0] = running_singlet;

    for(size_t i=1; i<data.size(); ++i) {
        running_singlet *= beta_s;
        running_triplet *= beta_t;

        double a = data[i] + running_singlet + running_triplet;
        double b = 2;
        double r = a;
        double p = b / (1.0 + b);
        double mode;

        if(r > 1) {
            mode = (r - 1) * (1 - p) / p;
        } else {
            mode = 0;
        }

        total[i] = data[i] - mode;
        total[i] = std::max(0.0, total[i]);
        running_singlet += alpha * total[i];
        running_triplet += (1 - alpha) * total[i];
        running[i] = running_singlet + running_triplet;
    }
}

class ChargeUnfolding: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;


    double tau_s_;
    double tau_t_;
    double alpha_;
    double delta_t_;

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
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("TauSinglet", "Time constant for singlet light", 8.13);
    AddParameter("TauTriplet", "Time constant for triplet light", 743);
    AddParameter("SingletTripletRatio", "Ratio of singlet to triplet light", 0.97);
    AddParameter("BinWidth", "Width of the time bins", 2.0);
    AddParameter("InputPulsesMaskName", "Name of the input pulses mask", std::string(""));
    AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
    AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
}

void ChargeUnfolding::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("TauSinglet", tau_s_);
    GetParameter("TauTriplet", tau_t_);
    GetParameter("SingletTripletRatio", alpha_);
    GetParameter("BinWidth", delta_t_);
    GetParameter("InputPulsesMaskName", pulses_mask_name_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);

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
    I3VectorDoublePtr total_positive = boost::make_shared<I3VectorDouble>(data.size(), 0.0);
    I3VectorDoublePtr total_bayesian = boost::make_shared<I3VectorDouble>(data.size(), 0.0);
    I3VectorDoublePtr total_bayesian_positive = boost::make_shared<I3VectorDouble>(data.size(), 0.0);

    I3VectorDoublePtr running = boost::make_shared<I3VectorDouble>(data.size(), 0.0);
    I3VectorDoublePtr running_positive = boost::make_shared<I3VectorDouble>(data.size(), 0.0);
    I3VectorDoublePtr running_bayesian = boost::make_shared<I3VectorDouble>(data.size(), 0.0);
    I3VectorDoublePtr running_bayesian_positive = boost::make_shared<I3VectorDouble>(data.size(), 0.0);

    compute_charge_unfolding(data, beta_s_, beta_t_, alpha_, *total, *running);
    compute_charge_unfolding_positive(data, beta_s_, beta_t_, alpha_, *total_positive, *running_positive);
    compute_charge_unfolding_bayesian(data, beta_s_, beta_t_, alpha_, *total_bayesian, *running_bayesian);
    compute_charge_unfolding_bayesian_positive(data, beta_s_, beta_t_, alpha_, *total_bayesian_positive, *running_bayesian_positive);

    I3MapPMTKeyVectorDoublePtr raw_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr total_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr total_positive_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr total_bayesian_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr total_bayesian_positive_map = boost::make_shared<I3MapPMTKeyVectorDouble>();

    I3MapPMTKeyVectorDoublePtr running_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr running_positive_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr running_bayesian_map = boost::make_shared<I3MapPMTKeyVectorDouble>();
    I3MapPMTKeyVectorDoublePtr running_bayesian_positive_map = boost::make_shared<I3MapPMTKeyVectorDouble>();

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
        total_positive_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        total_bayesian_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        total_bayesian_positive_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));

        running_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        running_positive_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        running_bayesian_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));
        running_bayesian_positive_map->insert(std::make_pair(i->first, std::vector<double>(data.size(), 0.0)));

        compute_charge_unfolding(data, beta_s_, beta_t_, alpha_, total_map->at(i->first), running_map->at(i->first));
        compute_charge_unfolding_positive(data, beta_s_, beta_t_, alpha_, total_positive_map->at(i->first), running_positive_map->at(i->first));
        compute_charge_unfolding_bayesian(data, beta_s_, beta_t_, alpha_, total_bayesian_map->at(i->first), running_bayesian_map->at(i->first));
        compute_charge_unfolding_bayesian_positive(data, beta_s_, beta_t_, alpha_, total_bayesian_positive_map->at(i->first), running_bayesian_positive_map->at(i->first));
    }

    frame->Put(output_prefix_ + "ChargeUnfoldingTotalRaw", raw);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotal", total);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalPositive", total_positive);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalBayesian", total_bayesian);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalBayesianPositive", total_bayesian_positive);

    frame->Put(output_prefix_ + "ChargeUnfoldingTotalRunning", running);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalRunningPositive", running_positive);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalRunningBayesian", running_bayesian);
    frame->Put(output_prefix_ + "ChargeUnfoldingTotalRunningBayesianPositive", running_bayesian_positive);

    frame->Put(output_prefix_ + "ChargeUnfoldingRaw", raw_map);
    frame->Put(output_prefix_ + "ChargeUnfolding", total_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingPositive", total_positive_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingBayesian", total_bayesian_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingBayesianPositive", total_bayesian_positive_map);

    frame->Put(output_prefix_ + "ChargeUnfoldingRunning", running_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingRunningPositive", running_positive_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingRunningBayesian", running_bayesian_map);
    frame->Put(output_prefix_ + "ChargeUnfoldingRunningBayesianPositive", running_bayesian_positive_map);

    PushFrame(frame);
}

