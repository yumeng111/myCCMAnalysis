#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <string>
#include <random>
#include <iostream>
#include <algorithm>

#include <icetray/I3Int.h>
#include <icetray/I3Frame.h>
#include <icetray/I3ConditionalModule.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/calibration/CCMCalibration.h>
#include <dataclasses/calibration/BaselineEstimate.h>

#include "daqtools/WaveformAccumulator.h"

class SumPulses : public I3ConditionalModule {
    std::string geometry_key_;
    std::string pulses_key_;
    std::string output_prefix_;
    bool skip_missing_;
    double bin_width_;
    double bin_offset_;

    bool geo_seen_ = false;
    CCMGeometry geo_;
    I3Vector<CCMPMTKey> allowed_pmt_keys_;
    I3Vector<CCMOMGeo::OMType> allowed_pmt_types_;
    std::vector<CCMPMTKey> pmt_keys_;
    std::string reference_time_key_;

public:
    SumPulses(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(SumPulses);

void SumPulses::ProcessFrame(I3FramePtr frame) {
    CCMRecoPulseSeriesMapConstPtr pulses =
        frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_key_);
    if(!pulses)
        log_fatal("Couldn't find '%s' in the frame!",
                pulses_key_.c_str());

    WaveformAccumulator summed_waveform;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        if(pulses->find(pmt_key) == pulses->end()) {
            if(skip_missing_) {
                log_debug("Skipping missing PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        pulses_key_.c_str());
                continue;
            } else {
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        pulses_key_.c_str());
            }
        }
        CCMRecoPulseSeries const & pulse_series = pulses->at(pmt_key);
        if(pulse_series.size() == 0) {
            log_debug("Skipping empty pulse series for PMT (%i/%u) in '%s'",
                    pmt_key.GetRegion(), pmt_key.GetSensor(),
                    pulses_key_.c_str());
            continue;
        }
        // Assume the pulses are sorted in time
        CCMRecoPulse const & pulse0 = pulse_series.front();
        CCMRecoPulse const & pulse1 = pulse_series.back();
        double t0 = pulse0.GetTime();
        double t1 = pulse1.GetTime();
        int32_t min_binning_index = (t0 - bin_offset_) / bin_width_;
        int32_t reference_index = -min_binning_index;
        int32_t max_binning_index = (t1 - bin_offset_) / bin_width_;
        size_t length = max_binning_index - min_binning_index + 1;
        std::vector<double> waveform(length, 0.0);
        for(CCMRecoPulse const & pulse : pulse_series) {
            double t = pulse.GetTime();
            int32_t binning_index = (t - bin_offset_) / bin_width_;
            waveform[binning_index + reference_index] += pulse.GetCharge();
        }

        summed_waveform.AddWaveform(waveform, reference_index);
    }

    std::deque<double> deque_samples = summed_waveform.GetSummedWaveform();
    std::deque<unsigned int> deque_counts = summed_waveform.GetCounts();
    int32_t fixed_position = summed_waveform.GetFixedPosition();

    boost::shared_ptr<I3Vector<double>> p_samples = boost::make_shared<I3Vector<double>>(deque_samples.begin(), deque_samples.end());
    boost::shared_ptr<I3Vector<uint32_t>> p_counts = boost::make_shared<I3Vector<uint32_t>>(deque_counts.begin(), deque_counts.end());
    boost::shared_ptr<I3Int> p_fixed_position = boost::make_shared<I3Int>(fixed_position);

    frame->Put((output_prefix_).c_str(), p_samples);
    frame->Put((output_prefix_ + "Counts").c_str(), p_counts);
    frame->Put((output_prefix_ + "FixedPosition").c_str(), p_fixed_position);
    frame->Put((output_prefix_ + "BinWidth").c_str(), boost::make_shared<I3Double>(bin_width_));
    frame->Put((output_prefix_ + "BinOffset").c_str(), boost::make_shared<I3Double>(bin_offset_));
}

SumPulses::SumPulses(const I3Context& context) : I3ConditionalModule(context) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("Pulses", "Name of pulse series to use", "WavedeformPulses");
    AddParameter("OutputPrefix", "Prefix for the module output", std::string("SummedPulses"));
    AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>());
    AddParameter("PMTTypes", "PMTKeys to run over", I3Vector<CCMOMGeo::OMType>{CCMOMGeo::OMType::CCM8inCoated, CCMOMGeo::OMType::CCM8inUncoated});
    AddParameter("SkipMissingInformation", "Skip frames that are missing information?", bool(false));
    AddParameter("BinWidth", "Width of the bins in ns", double(2.0));
    AddParameter("BinOffset", "Offset of the bins in ns", double(0.0));
}

void SumPulses::Configure() {
    GetParameter("CCMGeometryName", geometry_key_);
    GetParameter("Pulses", pulses_key_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("PMTKeys", allowed_pmt_keys_);
    GetParameter("PMTTypes", allowed_pmt_types_);
    GetParameter("SkipMissingInformation", skip_missing_);
    GetParameter("BinWidth", bin_width_);
    GetParameter("BinOffset", bin_offset_);

    // Preemtively sort the allowed_pmt_keys_ so they're ready for the Geometry function
    if(allowed_pmt_keys_.size() > 0)
        std::sort(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end());
}

void SumPulses::Geometry(I3FramePtr frame) {
    // Assumes allowed_pmt_keys_ is already sorted
    CCMGeometryConstPtr geo = frame->Get<CCMGeometryConstPtr>(geometry_key_);
    if (!geo)
        log_fatal("Couldn't find '%s' in the frame!",
                geometry_key_.c_str());
    geo_ = *geo;
    geo_seen_ = true;

    // Copy and sort the pmt keys from the geometry
    std::vector<CCMPMTKey> geo_keys; geo_keys.reserve(geo_.pmt_channel_map.size());
    for(std::pair<CCMPMTKey const, uint32_t> const & p : geo_.pmt_channel_map)
        geo_keys.push_back(p.first);
    std::sort(geo_keys.begin(), geo_keys.end());

    // Assume we're doing all PMTs if none are specified
    if(allowed_pmt_keys_.size() == 0) {
        if(allowed_pmt_types_.size() == 0) {
            pmt_keys_ = geo_keys;
        } else {
            // Clear the final pmt key list
            pmt_keys_.clear();
            // Fill the final pmt key list with the intersection of specified pmt types and those available in the geometry
            for(CCMPMTKey const & pmt_key : geo_keys) {
                if(std::find(allowed_pmt_types_.begin(), allowed_pmt_types_.end(), geo_.pmt_geo.at(pmt_key).omtype) != allowed_pmt_types_.end())
                    pmt_keys_.push_back(pmt_key);
            }
        }
    } else {
        // Clear the final pmt key list
        pmt_keys_.clear();
        // Fill the final pmt key list with the intersection of specified pmt keys and those available in the geometry
        // If allowed_pmt_keys_ is not sorted then this will fail horribly
        std::set_intersection(geo_keys.begin(), geo_keys.end(), allowed_pmt_keys_.begin(), allowed_pmt_keys_.end(), std::back_inserter(pmt_keys_));

        if(pmt_keys_.size() < allowed_pmt_keys_.size()) {
            std::vector<CCMPMTKey> missing_keys;
            std::set_difference(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end(), pmt_keys_.begin(), pmt_keys_.end(), std::back_inserter(missing_keys));
            std::stringstream ss;
            ss << "Some specified CCMPMTKeys are not present in the geometry:";
            for(CCMPMTKey const & key : missing_keys)
                ss << " " << key;
            log_warn(ss.str().c_str());
        }
    }

    PushFrame(frame);
}

void SumPulses::DAQ(I3FramePtr frame) {
    ProcessFrame(frame);
    PushFrame(frame);
}
