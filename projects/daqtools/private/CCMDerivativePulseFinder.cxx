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
#include <daqtools/WaveformDerivative.h>

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

class CCMDerivativePulseFinder: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string ccm_waveform_name_;
    std::string ccm_trigger_name_;
    std::string nim_pulses_name_;
    bool throw_without_input_;
    std::string output_name_;

    double pulse_initial_deriv_threshold;
    size_t pulse_max_width;
    size_t pulse_min_width;
    double pulse_min_height;
    double pulse_min_deriv_magnitude;
    double pulse_min_integral;

    CCMGeometry geo;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
public:
    CCMDerivativePulseFinder(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
    static std::vector<CCMRecoPulse> CheckForPulse(WaveformDerivative & deriv, size_t start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold, uint16_t & failure_reason);
    static void FindPulses(
        CCMRecoPulseSeries & output,
        WaveformDerivative & deriv,
        double deriv_threshold,
        size_t max_pulse_width,
        size_t min_pulse_width,
        double min_pulse_height,
        double min_deriv_magnitude,
        double min_integral);
    static void ApplyTimeOffset(CCMRecoPulseSeries & output, double offset);
};

I3_MODULE(CCMDerivativePulseFinder);

CCMDerivativePulseFinder::CCMDerivativePulseFinder(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMPulses", std::string("NIMPulses"));
    AddParameter("CCMCalibratedWaveformsName", "Key for the input CCMWaveformDoubleSeries", std::string("CCMCalibratedWaveforms"));
    AddParameter("ThrowWithoutInput", "Whether to throw an error when there is no input", true);

    AddParameter("InitialDerivativeThreshold", "Initial positive derivative threshold for a pulse", double(0.3));
    AddParameter("MinPulseWidth", "Minimum width for defining a pulse", size_t(5));
    AddParameter("MaxPulseWidth", "Maxiumum width for defining a pulse", size_t(100));
    AddParameter("MinPulseHeight", "Minimum height for defining a pulse", double(5.0));
    AddParameter("MinPulseDerivativeMagnitude", "Minimum derivative magnitude for defining a pulse", double(0.325));
    AddParameter("MinPulseIntegral", "Minimum integral for defining a pulse", double(10.0));
    AddParameter("OutputName", "Key for output CCMRecoPulseSeriesMap", std::string("DerivativePulses"));
}

void CCMDerivativePulseFinder::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("CCMCalibratedWaveformsName", ccm_waveform_name_);
    GetParameter("ThrowWithoutInput", throw_without_input_);
    GetParameter("InitialDerivativeThreshold", pulse_initial_deriv_threshold);
    GetParameter("MinPulseWidth", pulse_min_width);
    GetParameter("MaxPulseWidth", pulse_max_width);
    GetParameter("MinPulseHeight", pulse_min_height);
    GetParameter("MinPulseDerivativeMagnitude", pulse_min_deriv_magnitude);
    GetParameter("MinPulseIntegral", pulse_min_integral);
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
        return;
    }
    if(not frame->Has(nim_pulses_name_)) {
        if(throw_without_input_) {
            log_fatal("Input NIMPulses name %s not present", nim_pulses_name_.c_str());
        } else {
            log_warn("Input NIMPulses named %s not present", nim_pulses_name_.c_str());
        }
        return;
    }

    // let's read in our waveform
    boost::shared_ptr<const CCMWaveformDoubleSeries> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>(ccm_waveform_name_);
    if(waveforms == nullptr) {
        log_fatal("I3FrameObject of type \"CCMWaveformDoubleSeries\" does not exist under key \"%s\"", ccm_waveform_name_.c_str());
    }
    boost::shared_ptr<I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>> const> nim_pulses = frame->Get<boost::shared_ptr<I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>> const>>(nim_pulses_name_);
    if(waveforms == nullptr) {
        log_fatal("I3FrameObject of type \"I3Map<CCMTriggerKey, I3Vector<NIMLogicPulse>>\" does not exist under key \"%s\"", nim_pulses_name_.c_str());
    }

    std::vector<CCMPMTKey> pmt_keys;
    // std::cout << "pmt_channel_map_.size() == " << pmt_channel_map_.size() << std::endl;
    pmt_keys.reserve(pmt_channel_map_.size());
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        pmt_keys.push_back(p.first);
    }
    // std::cout << "Have " << pmt_keys.size() << " pmt keys" << std::endl;

    boost::shared_ptr<CCMRecoPulseSeriesMap> output(new CCMRecoPulseSeriesMap);

    // loop over each channel in waveforms
    for(size_t i=0; i<pmt_keys.size(); ++i) {
        CCMPMTKey key = pmt_keys[i];
        // std::cout << "Looking at " << key << std::endl;
        uint32_t channel = pmt_channel_map_.at(key);
        CCMTriggerKey trigger_key =  geo.trigger_copy_map.at(key);
        I3Vector<NIMLogicPulse> const & trigger_nim_pulses = nim_pulses->at(trigger_key);
        if(trigger_nim_pulses.size() == 0) {
            log_fatal("No board timing information available!");
        }
        double nim_pulse_time = trigger_nim_pulses[0].GetNIMPulseTime();
        WaveformDerivative deriv(waveforms->at(channel).GetWaveform().begin(), waveforms->at(channel).GetWaveform().end(), 2.0);
        FindPulses(output->operator[](key),
                deriv,
                pulse_initial_deriv_threshold,
                pulse_max_width,
                pulse_min_width,
                pulse_min_height,
                pulse_min_deriv_magnitude,
                pulse_min_integral);
        std::cout << key << " has " << output->operator[](key).size() << " pulses" << std::endl;
        ApplyTimeOffset(output->operator[](key), -nim_pulse_time);
    }

    frame->Put(output_name_, output);
    PushFrame(frame);
}

void CCMDerivativePulseFinder::FindPulses(
        CCMRecoPulseSeries & output,
        WaveformDerivative & deriv,
        double deriv_threshold,
        size_t max_pulse_width,
        size_t min_pulse_width,
        double min_pulse_height,
        double min_deriv_magnitude,
        double min_integral) {
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
    for(size_t i=0; i<N; ++i) {
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
}

std::vector<CCMRecoPulse> CCMDerivativePulseFinder::CheckForPulse(WaveformDerivative & deriv, size_t start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold, uint16_t & failure_reason) {
    deriv.Reset(start_idx);
    // 0 before start of pulse checking [checking for positive derivative]
    // 1 derivative is positive (in rising edge of pulse) [checking for negative derivative]
    // 2 derivative is negative (in falling edge of pulse) [checking for positive derivative]
    // 3 derivative is positive (recovering from droop) [done]
    int state = 1;
    double max_value = deriv.Value();
    double max_abs_derivative = deriv.Derivative();
    double integral = 0;
    bool found_pulse = false;
    size_t final_idx = start_idx;
    bool state_2_reached_zero = false;
    size_t N = std::min(deriv.Size(), start_idx + max_samples);
    // std::cout << "state = 0" << std::endl;
    for(size_t i=start_idx; i<N; ++i) {
        max_value = std::max(max_value, deriv.Value());
        max_abs_derivative = std::max(max_abs_derivative, std::abs(deriv.Derivative()));
        integral += deriv.Value();
        if(state % 2) {
            if(deriv.Derivative() < 0) {
                state += 1;
                // std::cout << "state = " << state << std::endl;
            }
        } else {
            if(state == 2) {
                // std::cout << "Deriv == " << deriv.Value() << std::endl;
                if(deriv.Value() <= 0.1) {
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
        deriv.Reset(start_idx);
        return {};
    }
}

void CCMDerivativePulseFinder::ApplyTimeOffset(CCMRecoPulseSeries & output, double offset) {
    for(CCMRecoPulse & pulse : output) {
        pulse.SetTime(pulse.GetTime() + offset);
    }
}
