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
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <thread>
#include <cmath>
#include <functional>

#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/BaselineEstimate.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include "daqtools/WaveformDerivative.h"
#include "daqtools/WaveformSmoother.h"

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

class CCMDerivativePulseFinder: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string ccm_waveform_name_;
    std::string ccm_trigger_name_;
    std::string nim_pulses_name_;
    bool throw_without_input_;
    std::string output_name_;

    size_t num_threads;
    std::vector<std::tuple<size_t, size_t>> thread_ranges_;
    I3Vector<CCMPMTKey> allowed_pmt_keys_;
    std::vector<CCMPMTKey> pmt_keys_;

    double pulse_initial_deriv_threshold;
    size_t pulse_max_width;
    size_t pulse_min_width;
    double pulse_min_height;
    double pulse_min_deriv_magnitude;
    double pulse_min_integral;

    bool use_raw_waveforms_;

    CCMGeometry geo;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
public:
    CCMDerivativePulseFinder(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
    template<typename Waveform>
    static std::vector<CCMRecoPulse> CheckForPulse(Waveform & deriv, size_t & start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold, uint16_t & failure_reason);
    template<typename Waveform>
    static void FindPulses(
        CCMRecoPulseSeries & output,
        Waveform & deriv,
        double deriv_threshold,
        size_t max_pulse_width,
        size_t min_pulse_width,
        double min_pulse_height,
        double min_deriv_magnitude,
        double min_integral,
        double nim_pulse_time);
    static void ApplyTimeOffset(CCMRecoPulseSeries & output, double offset);
    static void FindPulsesThread(
        std::tuple<size_t, size_t> thread_range,
        std::vector<CCMPMTKey> const & pmt_keys,
        I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map,
        CCMGeometry const & geo,
        I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>> const & nim_pulses,
        boost::shared_ptr<CCMWaveformDoubleSeries const> const & cal_waveforms,
        boost::shared_ptr<CCMWaveformUInt16Series const> const & raw_waveforms,
        bool use_raw_waveforms,
        std::vector<std::reference_wrapper<CCMRecoPulseSeries>> & output_pulses_references,
        double pulse_initial_deriv_threshold,
        size_t pulse_max_width,
        size_t pulse_min_width,
        double pulse_min_height,
        double pulse_min_deriv_magnitude,
        double pulse_min_integral
	); 
};

I3_MODULE(CCMDerivativePulseFinder);

CCMDerivativePulseFinder::CCMDerivativePulseFinder(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NumThreads", "Number of worker threads to use for pulse fitting", (size_t)(0));
    AddParameter("NIMPulsesName", "Key for NIMPulses", std::string("NIMPulses"));
    AddParameter("CCMWaveformsName", "Key for the input CCMWaveformDoubleSeries or CCMWaveformUInt16Series", std::string("CCMCalibratedWaveforms"));
    AddParameter("ThrowWithoutInput", "Whether to throw an error when there is no input", false);

    AddParameter("InitialDerivativeThreshold", "Initial positive derivative threshold for a pulse", double(0.3));
    AddParameter("MinPulseWidth", "Minimum width for defining a pulse", size_t(5));
    AddParameter("MaxPulseWidth", "Maxiumum width for defining a pulse", size_t(100));
    AddParameter("MinPulseHeight", "Minimum height for defining a pulse", double(5.0));
    AddParameter("MinPulseDerivativeMagnitude", "Minimum derivative magnitude for defining a pulse", double(0.325));
    AddParameter("MinPulseIntegral", "Minimum integral for defining a pulse", double(10.0));
    AddParameter("ExpectRawWaveforms", "Whether the input is the raw waveforms (CCMWaveformUInt16Series), or the calibrated waveforms (CCMWaveformDoubleSeries", bool(false));
	AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>());
    AddParameter("OutputName", "Key for output CCMRecoPulseSeriesMap", std::string("DerivativePulses"));
}

void CCMDerivativePulseFinder::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NumThreads", num_threads);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("CCMWaveformsName", ccm_waveform_name_);
    GetParameter("ThrowWithoutInput", throw_without_input_);
    GetParameter("InitialDerivativeThreshold", pulse_initial_deriv_threshold);
    GetParameter("MinPulseWidth", pulse_min_width);
    GetParameter("MaxPulseWidth", pulse_max_width);
    GetParameter("MinPulseHeight", pulse_min_height);
    GetParameter("MinPulseDerivativeMagnitude", pulse_min_deriv_magnitude);
    GetParameter("MinPulseIntegral", pulse_min_integral);
    GetParameter("ExpectRawWaveforms", use_raw_waveforms_);
	GetParameter("PMTKeys", allowed_pmt_keys_);
    GetParameter("OutputName", output_name_);
}


void CCMDerivativePulseFinder::Geometry(I3FramePtr frame) {
    // std::cout << "CCMDerivativePulseFinder::Geometry" << std::endl;
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    geo = frame->Get<CCMGeometry const>(geometry_name_);
    // std::cout << "geo.pmt_channel_map.size() == " << geo.pmt_channel_map.size() << std::endl;
    pmt_channel_map_ = geo.pmt_channel_map;
    // std::cout << "pmt_channel_map_.size() == " << pmt_channel_map_.size() << std::endl;
    I3Map<CCMPMTKey, CCMOMGeo> const & pmt_geo_ = geo.pmt_geo;
    std::set<CCMPMTKey> allowed_pmt_keys(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end());
    bool filter_pmts = allowed_pmt_keys.size() > 0;
    for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : pmt_geo_) {
        if(filter_pmts and allowed_pmt_keys.find(p.first) == allowed_pmt_keys.end()) {
            continue;
        }
        if(p.second.omtype == CCMOMGeo::OMType::CCM8inCoated or p.second.omtype == CCMOMGeo::OMType::CCM8inUncoated or p.second.omtype == CCMOMGeo::OMType::CCM1in) {
            pmt_keys_.push_back(p.first);
        }
    }
    size_t min_channels_per_thread = pmt_keys_.size() / num_threads;
    size_t max_channels_per_thread = min_channels_per_thread + 1;
    size_t num_threads_with_max = pmt_keys_.size() - num_threads * min_channels_per_thread;
    size_t channel_start = 0;
    for(size_t i=0; i<num_threads; ++i) {
        size_t n_channels = min_channels_per_thread;
        if(i < num_threads_with_max) {
            n_channels = max_channels_per_thread;
        }
        thread_ranges_.emplace_back(channel_start, channel_start + n_channels);
        channel_start += n_channels;
    }
    geo_seen = true;
    PushFrame(frame);
}

void CCMDerivativePulseFinder::DAQ(I3FramePtr frame) {
    // std::cout << "CCMDerivativePulseFinder::DAQ" << std::endl;
    if(not geo_seen) {
        log_fatal("No Geometry frame seen befor DAQ frame!");
    }

    if(not frame->Has(ccm_waveform_name_)) {
        if(throw_without_input_) {
            log_fatal("Input waveform name %s not present", ccm_waveform_name_.c_str());
        } else {
            log_warn("Input waveform name %s not present", ccm_waveform_name_.c_str());
        }
        PushFrame(frame);
        return;
    }
    if(not frame->Has(nim_pulses_name_)) {
        if(throw_without_input_) {
            log_fatal("Input NIMPulses name %s not present", nim_pulses_name_.c_str());
        } else {
            log_warn("Input NIMPulses named %s not present", nim_pulses_name_.c_str());
        }
        PushFrame(frame);
        return;
    }

    // let's read in our waveform
    boost::shared_ptr<const CCMWaveformDoubleSeries> cal_waveforms;
    boost::shared_ptr<const CCMWaveformUInt16Series> raw_waveforms;
    if(use_raw_waveforms_) {
        raw_waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>(ccm_waveform_name_);
        if(raw_waveforms == nullptr) {
            if(throw_without_input_) {
                log_fatal("I3FrameObject of type \"CCMWaveformUInt16Series\" does not exist under key \"%s\"", ccm_waveform_name_.c_str());
            } else {
                log_warn("I3FrameObject of type \"CCMWaveformUInt16Series\" does not exist under key \"%s\"", ccm_waveform_name_.c_str());
            }
            PushFrame(frame);
            return;
        }
    } else {
        cal_waveforms = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>(ccm_waveform_name_);
        if(cal_waveforms == nullptr) {
            if(throw_without_input_) {
                log_fatal("I3FrameObject of type \"CCMWaveformDoubleSeries\" does not exist under key \"%s\"", ccm_waveform_name_.c_str());
            } else {
                log_warn("I3FrameObject of type \"CCMWaveformDoubleSeries\" does not exist under key \"%s\"", ccm_waveform_name_.c_str());
            }
            PushFrame(frame);
            return;
        }
    }
    boost::shared_ptr<I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>> const> nim_pulses = frame->Get<boost::shared_ptr<I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>> const>>(nim_pulses_name_);
    if(nim_pulses == nullptr) {
        if(throw_without_input_) {
            log_fatal("I3FrameObject of type \"I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>>\" does not exist under key \"%s\"", nim_pulses_name_.c_str());
        } else {
            log_warn("I3FrameObject of type \"I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>>\" does not exist under key \"%s\"", nim_pulses_name_.c_str());
        }
        PushFrame(frame);
        return;
    }

    if(num_threads == 0) {
        num_threads = pmt_keys_.size();
    }

    std::vector<std::thread> threads;
    threads.reserve(thread_ranges_.size());

    boost::shared_ptr<CCMRecoPulseSeriesMap> output(new CCMRecoPulseSeriesMap);
    for(size_t i = 0; i < pmt_keys_.size(); ++i) {
        output->operator[](pmt_keys_[i]) = CCMRecoPulseSeries();
    }

    std::vector<std::reference_wrapper<CCMRecoPulseSeries>> output_pulses_references;
    output_pulses_references.reserve(pmt_keys_.size());
    for(size_t i = 0; i < pmt_keys_.size(); ++i) {
        output_pulses_references.emplace_back(output->at(pmt_keys_[i]));
    }

    // loop over each channel in waveforms
    for(size_t i=0; i<num_threads; ++i) {
        threads.emplace_back(
            FindPulsesThread,
            thread_ranges_[i],
            std::cref(pmt_keys_),
            std::cref(pmt_channel_map_),
            std::cref(geo),
            std::cref(*nim_pulses),
            std::cref(cal_waveforms),
            std::cref(raw_waveforms),
            use_raw_waveforms_,
            std::ref(output_pulses_references),
            pulse_initial_deriv_threshold,
            pulse_max_width,
            pulse_min_width,
            pulse_min_height,
            pulse_min_deriv_magnitude,
            pulse_min_integral
        );
    }

    for(size_t i=0; i<num_threads; ++i) {
        threads[i].join();
    }

    frame->Put(output_name_, output);
    PushFrame(frame);
}

void CCMDerivativePulseFinder::FindPulsesThread(
	std::tuple<size_t, size_t> thread_range,
	std::vector<CCMPMTKey> const & pmt_keys,
	I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map,
    CCMGeometry const & geo,
    I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>> const & nim_pulses,
	boost::shared_ptr<CCMWaveformDoubleSeries const> const & cal_waveforms,
	boost::shared_ptr<CCMWaveformUInt16Series const> const & raw_waveforms,
    bool use_raw_waveforms,
    std::vector<std::reference_wrapper<CCMRecoPulseSeries>> & output_pulses_references,
    double pulse_initial_deriv_threshold,
    size_t pulse_max_width,
    size_t pulse_min_width,
    double pulse_min_height,
    double pulse_min_deriv_magnitude,
    double pulse_min_integral
	) { 
    for(size_t i=std::get<0>(thread_range); i<std::get<1>(thread_range); ++i) {
        CCMPMTKey key = pmt_keys[i];
        uint32_t channel = pmt_channel_map.at(key);
        CCMTriggerKey trigger_key =  geo.trigger_copy_map.at(key);
        I3Vector<NIMLogicPulse> const & trigger_nim_pulses = nim_pulses.at(trigger_key);
        if(trigger_nim_pulses.size() == 0) {
            log_warn("No board timing information available!");
            continue;
        }
        double nim_pulse_time = trigger_nim_pulses[0].GetNIMPulseTime();
        if(use_raw_waveforms) {
            WaveformSmoother smoother = WaveformSmoother(raw_waveforms->at(channel).GetWaveform().begin(), raw_waveforms->at(channel).GetWaveform().end(), 2.0, 12.0);
            CCMDerivativePulseFinder::FindPulses<WaveformSmoother>(output_pulses_references.at(i),
                    smoother,
                    pulse_initial_deriv_threshold,
                    pulse_max_width,
                    pulse_min_width,
                    pulse_min_height,
                    pulse_min_deriv_magnitude,
                    pulse_min_integral,
                    nim_pulse_time);
        } else {
            WaveformDerivative<double> deriv = WaveformDerivative<double>(cal_waveforms->at(channel).GetWaveform().begin(), cal_waveforms->at(channel).GetWaveform().end(), 2.0);
            CCMDerivativePulseFinder::FindPulses<WaveformDerivative<double>>(output_pulses_references.at(i),
                    deriv,
                    pulse_initial_deriv_threshold,
                    pulse_max_width,
                    pulse_min_width,
                    pulse_min_height,
                    pulse_min_deriv_magnitude,
                    pulse_min_integral,
                    nim_pulse_time);
        }
        ApplyTimeOffset(output_pulses_references.at(i), -nim_pulse_time);
    }
}

template<typename Waveform>
void CCMDerivativePulseFinder::FindPulses(
        CCMRecoPulseSeries & output,
        Waveform & deriv,
        double deriv_threshold,
        size_t max_pulse_width,
        size_t min_pulse_width,
        double min_pulse_height,
        double min_deriv_magnitude,
        double min_integral,
        double nim_pulse_time) {
    // Determine the maximum starting sample for a pulse
    // We need at least 5 samples after the starting pulse
    size_t N = deriv.Size();
    // Define the maximum sample to consider for the start of a pulse
    // A pulse needs at least 5 samples, so subtract off 5 from the max sample we can look at
    if(N > 5)
        N -= 5;
    size_t pulse_first_index = 0;
    size_t pulse_last_index = 0;

    std::map<uint16_t, size_t> failure_counts;
    // Reset the smoother position to the start of the waveform
    deriv.Reset();
    // nim_pulse_time is in ns, so divide by 2 to get bins and subtract 250 to get -500. Add 250 to get +500.
    for(size_t i=0; i<N; ++i) {//Loop from start to end of the time region desired.
        uint16_t failure_reason = 0;
        // Check the condition for the beginning of a pulse
        if(deriv.Derivative() > deriv_threshold) {
            // std::cout << "Derivative above threshold" << std::endl;
            // Get the last index of the found pulse, 0 if no pulse is found
            std::vector<CCMRecoPulse> pulses = CheckForPulse(deriv, i, max_pulse_width, min_pulse_width, min_pulse_height, min_deriv_magnitude, min_integral, failure_reason);
            if(pulses.size() == 0) {
                if(failure_counts.find(failure_reason) == failure_counts.end()) {
                    failure_counts[failure_reason] = 0;
                }
                failure_counts[failure_reason] += 1;
            }
            for(size_t j=0; j<pulses.size(); ++j) {
                output.push_back(pulses[j]);
                pulse_first_index = pulses[j].GetTime() / 2;
                pulse_last_index = pulse_first_index + pulses[j].GetWidth() / 2;
                i = pulse_last_index;
                deriv.Reset(pulse_last_index);
            }
        }
        deriv.Next();
    }
    /*
    for(std::pair<uint16_t const, size_t> failure : failure_counts) {
        std::cout << "Failure mode: b";
        for(size_t i=0; i<5; ++i) {
            if(failure.first & (0x1 <<i)) {
                std::cout << "1";
            } else {
                std::cout << "0";
            }
        }
        std::cout << " appeared " << failure.second << " times" << std::endl;
    }
    */
}

template<typename Waveform>
std::vector<CCMRecoPulse> CCMDerivativePulseFinder::CheckForPulse(Waveform & deriv, size_t & start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold, uint16_t & failure_reason) {
    deriv.Reset(start_idx);
    // 0 before start of pulse checking [checking for positive derivative]
    // 1 derivative is positive (in rising edge of pulse) [checking for negative derivative]
    // 2 derivative is negative (in falling edge of pulse) [checking for positive derivative]
    // 3 derivative is positive (recovering from droop) [done]
    double baseline = deriv.Baseline();
    int state = 1;
    double max_value = deriv.Value() - baseline;
    double max_abs_derivative = deriv.Derivative();
    double integral = 0;
    bool found_pulse = false;
    size_t final_idx = start_idx;
    bool state_2_reached_zero = false;
    size_t N = std::min(deriv.Size(), start_idx + max_samples);
    // std::cout << "state = 0" << std::endl;
    for(size_t i=start_idx; i<N; ++i) {
        max_value = std::max(max_value, deriv.Value() - baseline);
        max_abs_derivative = std::max(max_abs_derivative, std::abs(deriv.Derivative()));
        integral += deriv.Value() - baseline;
        if(state % 2) {
            if(deriv.Derivative() < 0) {
                state += 1;
                // std::cout << "state = " << state << std::endl;
            }
        } else {
            if(state == 2) {
                // std::cout << "Deriv == " << deriv.Value() << std::endl;
                if(deriv.Value() - baseline <= 0.1) {
                    state_2_reached_zero = true;
                    // std::cout << "Reached zero!" << std::endl;
                }
                if(not state_2_reached_zero) {
                    deriv.Next();
                    continue;
                }
            }
            if(deriv.Derivative() > 0) {
                state += 1;
                // std::cout << "state = " << state << std::endl;
            }
        }
        if(state >= 3) {
            // std::cout << "Found pulse" << std::endl;
            found_pulse = true;
            final_idx = i;
            break;
        }
        deriv.Next();
    }
    if(found_pulse) {
        // std::cout << "Found a pulse" << std::endl;
        if(max_value >= value_threshold) {
            // std::cout  << "Max value  : passed" << std::endl;
        } else {
            // std::cout  << "Max value  : failed" << std::endl;
            failure_reason |= (0x1 << 1);
        }
        if(max_abs_derivative >= derivative_threshold) {
            // std::cout  << "Max abs der: passed" << std::endl;
        } else {
            // std::cout  << "Max abs der: failed" << std::endl;
            failure_reason |= (0x1 << 2);
        }
        if(integral >= integral_threshold) {
            // std::cout  << "Integral   : passed" << std::endl;
        } else {
            // std::cout  << "Integral   : failed" << std::endl;
            failure_reason |= (0x1 << 3);
        }
        if(final_idx - start_idx >= min_length) {
            // std::cout  << "Min Length : passed" << std::endl;
        } else {
            // std::cout  << "Min Length : failed" << std::endl;
            failure_reason |= (0x1 << 4);
        }
    } else {
        // std::cout << "No pulse found" << std::endl;
        failure_reason |= 0x1;
    }
    if(found_pulse
            and (max_value >= value_threshold)
            and (max_abs_derivative >= derivative_threshold)
            and (integral >= integral_threshold)
            and (final_idx - start_idx >= min_length)
            ) {
        CCMRecoPulse pulse;
        pulse.SetTime(double(start_idx) * 2.0);
        pulse.SetCharge(integral);
        pulse.SetWidth((double(final_idx) - double(start_idx)) * 2.0);
        return {pulse};
    } else {
        if(not found_pulse) {
            deriv.Reset(start_idx);
        } else {
            start_idx = final_idx;
        }
        return {};
    }
}

void CCMDerivativePulseFinder::ApplyTimeOffset(CCMRecoPulseSeries & output, double offset) {
    for(CCMRecoPulse & pulse : output) {
        pulse.SetTime(pulse.GetTime() + offset);
    }
}
