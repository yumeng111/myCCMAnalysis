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
#include <string>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3ConditionalModule.h>
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

class IntervalChargeSumQ: public I3ConditionalModule {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    I3Vector<double> time_windows_;
    std::string input_prefix_;
    std::string output_prefix_;

    std::string raw_pulses_name_;

    I3Vector<CCMOMGeo::OMType> pmt_types;
    std::set<CCMPMTKey> pmt_keys;

    I3Vector<std::string> bad_keys_names_ = {"BadKeys"};
    std::set<CCMPMTKey> bad_keys;

    public:
    void Geometry(I3FramePtr frame);
    IntervalChargeSumQ(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
};

I3_MODULE(IntervalChargeSumQ);

IntervalChargeSumQ::IntervalChargeSumQ(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_(""), pmt_types(I3Vector<CCMOMGeo::OMType>{CCMOMGeo::OMType::CCM8inCoated, CCMOMGeo::OMType::CCM8inUncoated}) {
    I3Vector<double> default_time_windows;
    default_time_windows.push_back(90.0);
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("TimeWindows", "Time window for charge estimate", default_time_windows);
    AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
    AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
    AddParameter("BadKeysNames", "List of bad keys to exclude from event finding", bad_keys_names_);
}

void IntervalChargeSumQ::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("TimeWindows", time_windows_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("BadKeysNames", bad_keys_names_);

    std::sort(time_windows_.begin(), time_windows_.end());

    if(raw_pulses_name_.empty()) {
        raw_pulses_name_ = input_prefix_ + "EventPulses";
    }
}

void IntervalChargeSumQ::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    // check if frame has list of bad keys to exclude
    bad_keys.clear();
    for(std::string const & bad_keys_name : bad_keys_names_) {
        if(frame->Has(bad_keys_name)) {
            boost::shared_ptr<const I3Vector<CCMPMTKey>> keys = frame->Get<boost::shared_ptr<const I3Vector<CCMPMTKey>>>(bad_keys_name);
            for(CCMPMTKey const & k : *keys) {
                bad_keys.insert(k);
            }
        }
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    pmt_keys.clear();
    if(geo_seen) {
        std::set<CCMOMGeo::OMType> allowed_pmt_types(pmt_types.begin(), pmt_types.end());
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            if(allowed_pmt_types.count(p.second.omtype) == 0)
                continue;
            if(bad_keys.count(p.first) != 0)
                continue;
            pmt_keys.insert(p.first);
        }
    }
    PushFrame(frame);
}

void IntervalChargeSumQ::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    CCMRecoPulseSeriesMapConstPtr input_pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(raw_pulses_name_);
    boost::shared_ptr<I3Map<I3ParticleID, CCMRecoPulseSeriesMap>const> input_by_particle = frame->Get<boost::shared_ptr<I3Map<I3ParticleID, CCMRecoPulseSeriesMap>const>> (raw_pulses_name_);
    std::map<I3ParticleID, std::reference_wrapper<CCMRecoPulseSeriesMap const>> inputs;

    std::vector<I3VectorDoublePtr> output_charges;
    std::vector<boost::shared_ptr<I3Map<I3ParticleID, std::vector<double>>>> output_by_particle;
    std::map<I3ParticleID, std::vector<std::reference_wrapper<std::vector<double>>>> outputs;

    I3VectorDoubleConstPtr event_start_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "EventStartTimes");
    I3VectorDoubleConstPtr event_end_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "EventEndTimes");

    if(input_pulses) {
        // Map the pulses into the "inputs" object
        inputs.insert({I3ParticleID(), *input_pulses});

        // Allocating memory for the final output
        for(auto const & time_window : time_windows_) {
            I3VectorDoublePtr event_charge = boost::make_shared<I3VectorDouble>(event_start_times->size(), 0.0);
            output_charges.push_back(event_charge);
        }

        // Map the final output locations into the "outputs" object
        std::vector<std::reference_wrapper<std::vector<double>>> & output = outputs[I3ParticleID()];
        for(size_t i=0; i<time_windows_.size(); ++i) {
            output.push_back(*output_charges[i]);
        }
    } else if(input_by_particle) {
        // Map the pulses for each particle into the "inputs" object
        for(auto const & [particleID, pulses] : *input_by_particle) {
            inputs.insert({particleID, pulses});
        }

        // Allocating memory for the final output
        for(auto const & time_window : time_windows_) {
            boost::shared_ptr<I3Map<I3ParticleID, std::vector<double>>> map_of_event_charge = boost::make_shared<I3Map<I3ParticleID, std::vector<double>>>();
            for(auto const & [particleID, pulses] : *input_by_particle) {
                (*map_of_event_charge)[particleID] = std::vector<double>(event_start_times->size(), 0.0);
            }
            output_by_particle.push_back(map_of_event_charge);
        }

        // Map the final output locations into the "outputs" object
        for(auto const & [particleID, pulses] : *input_by_particle) {
            std::vector<std::reference_wrapper<std::vector<double>>> & output = outputs[particleID];
            for(size_t i=0; i<time_windows_.size(); ++i) {
                output.push_back((*output_by_particle[i])[particleID]);
            }
        }
    } else {
        std::stringstream ss;
        ss << "Could not find ";
        ss << raw_pulses_name_;
        ss << " in the DAQ frame.";
        log_fatal("%s", ss.str().c_str());
        PushFrame(frame);
        return;
    }

    for(auto const & [particleID, pulses_ref] : inputs) {
        CCMRecoPulseSeriesMap const & pulses = pulses_ref.get();
        std::vector<std::reference_wrapper<std::vector<double>>> & output = outputs[particleID];
        PMTKeyPulseVector pulse_list;
        for(auto const & [pmt_key, pulse_series] : pulses) {
            if(pmt_keys.count(pmt_key) == 0)
                continue;
            for(CCMRecoPulse const & pulse : pulse_series) {
                pulse_list.push_back(PMTKeyPulsePair(pmt_key, pulse));
            }
        }

        std::sort(pulse_list.begin(), pulse_list.end(), [](auto const & t0, auto const & t1){return std::get<1>(t0).GetTime() < std::get<1>(t1).GetTime();});

        size_t event_idx = 0;
        for (PMTKeyPulseVector::const_iterator i = pulse_list.begin(); i != pulse_list.end() and event_idx < event_start_times->size(); ++i) {
            if(event_start_times->at(event_idx) > std::get<1>(*i).GetTime()) {
                continue;
            }
            size_t time_bin = 0;
            double total_charge = 0.0;
            double start_time = event_start_times->at(event_idx);
            double end_time = event_end_times->at(event_idx);
            for (PMTKeyPulseVector::const_iterator j = i; j != pulse_list.end(); ++j) {
                double time = std::get<1>(*j).GetTime();
                if(time - start_time > time_windows_[time_bin] or time > end_time) {
                    output[time_bin].get().at(event_idx) = total_charge;
                    ++time_bin;
                    if(time_bin == time_windows_.size()) {
                        break;
                    }
                }
                total_charge += std::get<1>(*j).GetCharge();
            }
            if(time_bin < time_windows_.size())
                output[time_bin].get().at(event_idx) = total_charge;
            time_bin = 0;
            ++event_idx;
        }
    }

    if(input_pulses) {
        for(size_t i = 0; i < time_windows_.size(); ++i) {
            std::string output_key = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_[i])) + "NS";
            frame->Put(output_key, output_charges[i]);
        }
    } else if(input_by_particle) {
        for(size_t i = 0; i < time_windows_.size(); ++i) {
            std::string output_key = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_[i])) + "NS";
            frame->Put(output_key, output_by_particle[i]);
        }
    }

    PushFrame(frame);
}

class IntervalChargeSumP: public I3ConditionalModule {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    I3Vector<double> time_windows_;
    std::string input_prefix_;
    std::string output_prefix_;

    std::string raw_pulses_name_;

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

    I3Vector<std::string> bad_keys_names_ = {"BadKeys"};
    std::set<CCMPMTKey> bad_keys;

  public:
    void Geometry(I3FramePtr frame);
    IntervalChargeSumP(const I3Context&);
    void Configure();
    void Physics(I3FramePtr frame);
};

I3_MODULE(IntervalChargeSumP);

IntervalChargeSumP::IntervalChargeSumP(const I3Context& context)
  : I3ConditionalModule(context),
    geo_seen(false),
    geometry_name_("") {

    I3Vector<double> default_time_windows;
    default_time_windows.push_back(90.0);

    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("TimeWindows", "Time windows for charge estimate", default_time_windows);
    AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
    AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
    AddParameter("BadKeysNames", "List of bad keys to exclude from event finding", bad_keys_names_);
}

void IntervalChargeSumP::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("TimeWindows", time_windows_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("BadKeysNames", bad_keys_names_);

    std::sort(time_windows_.begin(), time_windows_.end());

    if(raw_pulses_name_.empty()) {
        raw_pulses_name_ = input_prefix_ + "EventPulses";
    }
}

void IntervalChargeSumP::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    // check if frame has list of bad keys to exclude
    bad_keys.clear();
    for(std::string const & bad_keys_name : bad_keys_names_) {
        if(frame->Has(bad_keys_name)) {
            boost::shared_ptr<const I3Vector<CCMPMTKey>> keys = frame->Get<boost::shared_ptr<const I3Vector<CCMPMTKey>>>(bad_keys_name);
            for(CCMPMTKey const & k : *keys) {
                bad_keys.insert(k);
            }
        }
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    pmt_keys.clear();
    if(geo_seen) {
        std::set<CCMOMGeo::OMType> allowed_pmt_types(pmt_types.begin(), pmt_types.end());
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            if(allowed_pmt_types.count(p.second.omtype) == 0)
                continue;
            if(bad_keys.count(p.first) != 0)
                continue;
            pmt_keys.insert(p.first);
        }
    }
    PushFrame(frame);
}

void IntervalChargeSumP::Physics(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    // If we already have per-event Q-frame outputs, downselect them to this P frame:
    //  - legacy vector<double> -> I3Double
    //  - new map<I3ParticleID, vector<double>> -> map<I3ParticleID, double>
    size_t const event_index =
        frame->Get<CCMEventHeaderConstPtr>("CCMEventHeader")->GetSubEventID();

    std::string const first_key =
        output_prefix_ + "EventCharges" + std::to_string(int(time_windows_.at(0))) + "NS";

    if(frame->Has(first_key)) {
        // Branch 1: legacy vector<double>
        if(auto vec = frame->Get<I3VectorDoubleConstPtr>(first_key)) {
            for(size_t i = 0; i < time_windows_.size(); ++i) {
                const std::string in_key  = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_[i])) + "NS";
                const std::string out_key = output_prefix_ + "EventCharge"  + std::to_string(int(time_windows_[i])) + "NS";
                auto v = frame->Get<I3VectorDoubleConstPtr>(in_key);
                if(not v) log_fatal("Expected \"%s\" but it's missing.", in_key.c_str());
                frame->Put(out_key, boost::make_shared<I3Double>(v->at(event_index)));
            }
            PushFrame(frame);
            return;
        }
        // Branch 2: new per-particle map<I3ParticleID, vector<double>>
        if(auto mp = frame->Get<boost::shared_ptr<const I3Map<I3ParticleID, std::vector<double>>>>(first_key)) {
            for(size_t i = 0; i < time_windows_.size(); ++i) {
                const std::string in_key  = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_[i])) + "NS";
                const std::string out_key = output_prefix_ + "EventCharge"  + std::to_string(int(time_windows_[i])) + "NS";

                auto m = frame->Get<boost::shared_ptr<const I3Map<I3ParticleID, std::vector<double>>>>(in_key);
                if(not m) log_fatal("Expected \"%s\" but it's missing.", in_key.c_str());

                auto out_map = boost::make_shared<I3Map<I3ParticleID, double>>();
                for(auto const & [pid, vec] : *m) {
                    if(event_index >= vec.size())
                        log_fatal("Event index %zu out of range (size=%zu) for particle map \"%s\".", event_index, vec.size(), in_key.c_str());
                    (*out_map)[pid] = vec[event_index];
                }
                frame->Put(out_key, out_map);
            }
            PushFrame(frame);
            return;
        }
        // If we get here, the type under first_key wasn't recognized; fall through to compute.
    }

    // Otherwise: compute charges now from pulses.
    // Accept either legacy pulses or new per-particle pulses.
    CCMRecoPulseSeriesMapConstPtr pulses_single =
        frame->Get<CCMRecoPulseSeriesMapConstPtr>(raw_pulses_name_);
    auto pulses_by_particle =
        frame->Get<boost::shared_ptr<I3Map<I3ParticleID, CCMRecoPulseSeriesMap> const>>(raw_pulses_name_);

    if(not pulses_single and not pulses_by_particle) {
        std::stringstream ss;
        ss << "Could not find \"" << raw_pulses_name_ << "\" (neither single-map nor per-particle) in Physics frame.";
        log_fatal("%s", ss.str().c_str());
        PushFrame(frame);
        return;
    }

    // Event window
    I3DoubleConstPtr start_time_ptr = frame->Get<I3DoubleConstPtr>(input_prefix_ + "EventStartTime");
    I3DoubleConstPtr end_time_ptr   = frame->Get<I3DoubleConstPtr>(input_prefix_ + "EventEndTime");
    if (!start_time_ptr)
        log_fatal("Could not find event start time \"%s\".", (input_prefix_ + "EventStartTime").c_str());
    if (!end_time_ptr)
        log_fatal("Could not find event end time \"%s\".", (input_prefix_ + "EventEndTime").c_str());

    const double start_time = start_time_ptr->value;
    const double end_time = end_time_ptr->value;

    auto compute_totals = [&](const CCMRecoPulseSeriesMap &pulses) -> std::vector<double> {
        std::vector<CCMRecoPulse> pulse_list;
        pulse_list.reserve(1024);

        for(auto const & [key, series] : pulses) {
            if(not pmt_keys.count(key)) continue;
            if(bad_keys.count(key)) continue;
            for(CCMRecoPulse const & p : series) {
                pulse_list.emplace_back(p);
            }
        }
        std::sort(pulse_list.begin(), pulse_list.end(),
                  [](CCMRecoPulse const & a, CCMRecoPulse const & b){ return a.GetTime() < b.GetTime(); });

        std::vector<double> totals(time_windows_.size(), 0.0);
        if(pulse_list.empty()) return totals;

        // Find first pulse >= start_time
        auto it = std::lower_bound(
            pulse_list.begin(), pulse_list.end(), start_time,
            [](CCMRecoPulse const & pulse, double t){ return pulse.GetTime() < t; });

        size_t bin = 0;
        double accum = 0.0;

        size_t size = pulse_list.size();
        for(size_t i = std::distance(pulse_list.begin(), it); i < size; ++i) {
            CCMRecoPulse const & pulse = pulse_list[i];
            double const t = pulse.GetTime();
            if(t > end_time)
                break;

            // Write out any windows we’ve just passed
            while(bin < time_windows_.size() and (t - start_time) > time_windows_[bin]) {
                totals[bin] = accum;
                ++bin;
            }
            if(bin == time_windows_.size())
                break;

            accum += pulse.GetCharge();
        }

        // Fill remaining windows with the final accumulated total
        while(bin < time_windows_.size()) {
            totals[bin] = accum;
            ++bin;
        }
        return totals;
    };

    if(pulses_single) {
        std::vector<double> const totals = compute_totals(*pulses_single);
        for(size_t i = 0; i < time_windows_.size(); ++i) {
            const std::string key = output_prefix_ + "EventCharge" + std::to_string(int(time_windows_[i])) + "NS";
            frame->Put(key, boost::make_shared<I3Double>(totals[i]));
        }
    } else {
        // Per-particle inputs -> per-particle outputs (one map per time window)
        auto const & by_pid = *pulses_by_particle;
        // Pre-compute per-particle totals
        std::map<I3ParticleID, std::vector<double>> totals_by_pid;
        for(auto const & [pid, pulses] : by_pid) {
            totals_by_pid.emplace(pid, compute_totals(pulses));
        }
        // Emit one map per window
        for(size_t w = 0; w < time_windows_.size(); ++w) {
            std::string const key = output_prefix_ + "EventCharge" + std::to_string(int(time_windows_[w])) + "NS";
            auto out_map = boost::make_shared<I3Map<I3ParticleID, double>>();
            for(auto const & [pid, totals] : totals_by_pid) {
                (*out_map)[pid] = totals[w];
            }
            frame->Put(key, out_map);
        }
    }

    PushFrame(frame);
}

