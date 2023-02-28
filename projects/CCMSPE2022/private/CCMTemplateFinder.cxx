#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <nlopt.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibraiton.h>
#include <dataclasses/calibration/I3DomCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"



// Fits the SPE template shape to find a waveform template 
// Collects statistics about these regions

class CCMTemplateFinder : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string daq_config_name_;
    std::string waveforms_name_;
    std::string peak_amplitude_name_;
    std::string peak_time_name_;
    std::string peak_length_name_;

    std::string template_output_name_;

    std::vector<I3FramePtr> cached_frames;

    int fit_length;
    int min_amplitude;;
    int max_amplitude;;
    size_t n_daq_frames = 0;
    

    std::vector<WindowStats> GetPeaks(CCMWaveformUInt16 const & wf, double derivative_threshold, size_t min_window_size);
    void AddPeaks(I3FramePtr frame);

public:
    CCMTemplateFinder(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMTemplateFinder);

CCMTemplateFinder::CCMTemplateFinder(const I3Context& context) : I3Module(context),
    geometry_name_(""), daq_config_name_(""), waveforms_name_(""), peak_amplitude_name_(""), peak_time_name_(""), peak_length_name_(""), template_output_name_(""),
    fit_length(10), min_amplitude(20), max_amplitude(40), cached_frames(),
    thresholds_tuned_(false) {

    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMAnalysis::Binary::CCMDAQConfig>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("CCMPeakAmplitudesName", "Key to output vector of peak amplitudes", std::string("PeakAmplitudes"));
    AddParameter("CCMPeakLengthsName", "Key to output vector of peak lengths", std::string("PeakLengths"));
    AddParameter("CCMPeakTimesName", "Key to output vector of peak times", std::string("PeakTimes"));
    AddParameter("TemplateOutputName", "The output key of the baseline fit.", std::string("SPETemplate"));
    AddParameter("AbsoluteTimeName", "Key for absoluting timing information of frame", std::string("FirstTriggerTime"));
}

void CCMTemplateFinder::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMDAQConfigName", daq_config_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("CCMPeakAmplitudesName", peak_amplitude_name_);
    GetParameter("CCMPeakLengthsName", peak_lengths_name_);
    GetParameter("CCMPeakTimesName", peak_times_name_);
    GetParameter("TemplateOutputName", template_output_name_);
    GetParameter("AbsoluteTimeName", abs_time_name_);
}


std::vector<WindowStats> CCMTemplateFinder::GetPeaks(CCMWaveformUInt16 const & wf, double derivative_threshold, size_t min_window_size) {
    WaveformSource<uint16_t> wf_source(wf);
    MultiplicativeFilter<uint16_t, double> m_filter(1.0, &wf_source); // don't make negative yet
    BufferView<double, double> buffer_m(&m_filter);
    ExponentialFilter<double, double> exp_filter(10.0, 2.0, &buffer_m);
    std::deque<double> smoothing_kernel = {1,1,1};
    SmoothingFilter<double, double> smoothing_filter(1, smoothing_kernel, &exp_filter);
    BackDerivativeFilter<double, double> derivative_filter(2.0, &smoothing_filter);
    DerivativeWindowCapture<double, double, double> window(derivative_filter, buffer_m, derivative_threshold);

    std::vector<WindowStats> peaks;

    try {
        while(true) {
            bool have_window = window.Next();
            if(have_window) {
                std::tuple<std::deque<double>::const_iterator, std::deque<double>::const_iterator, size_t> window_buffer = window.GetBuffer();
                std::deque<double>::const_iterator begin = std::get<0>(window_buffer);
                std::deque<double>::const_iterator end = std::get<1>(window_buffer);
                size_t size = std::get<2>(window_buffer);
                if(size >= min_window_size) {
                    peaks.emplace_back(begin, end);
                    peaks.back().time = window.pos - size;// + wf.GetStartTime();
                }
            }
        }
    } catch(SourceExhausted const & e) {
    }

    return peaks;
}

double CCMTemplateFinder::PeakLossFunction(const std::vector<double> &x,
                                           const std::vector<double> &grad,
                                           void * f_data) {
}

SPETemplate CCMTemplateFinder::FitPeak(int64_t const & length,
                                       int64_t const & time,
                                       std::vector<double> const & waveform,
                                       int64_t pre_window = 3) {
  
  std::vector<double> 
}

void CCMTemplateFidner::AddTemplates(I3framePtr frame) {
    
    
    // let's grab some stuff from the frame
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series>(waveforms_name_);
    size_t size = waveforms.size();
    I3Vector<I3Vector<int>> const & window_amplitudes = frame->Get<I3Vector<I3Vector<int>>>(peak_amplitude_name_);
    I3Vector<I3Vector<int64_t>> const & window_lengths = frame->Get<I3Vector<I3Vector<int>>>(peak_length_name_);
    I3Vector<I3Vector<int64_t>> const & window_times = frame->Get<I3Vector<I3Vector<int>>>(peak_time_name_);
    
    // a vector storing the template for each channel
    boost::shared_ptr<I3Vector<I3Vector<SPETemplate>>> channel_templates(new I3Vector<I3Vector<SPETemplate>>(size));
     
    for(size_t ich = 0; ich < waveforms.size(); ++ich) {
      for(size_t ipk = 0; ipk < window_amplitudes[ich].size(); ++ipk) {
        FitPeak(window_amplitudes[ich][ipk],
                window_lengths[ich][ipk],
                window_times[ich][ipk],
      }
      
    }
}

void CCMTemplateFinder::AddPeaks(I3FramePtr frame) {
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series>(waveforms_name_);
    int64_t frame_time = frame->Get<I3PODHolder<int64_t>>(abs_time_name_).value;
    size_t size = waveforms.size();
    boost::shared_ptr<I3Vector<I3Vector<int>>> window_peaks(new I3Vector<I3Vector<int>>(size));
    boost::shared_ptr<I3Vector<I3Vector<int64_t>>> window_times(new I3Vector<I3Vector<int64_t>>(size));
    boost::shared_ptr<I3Vector<I3Vector<int64_t>>> window_lengths(new I3Vector<I3Vector<int64_t>>(size));

    for(size_t i=0; i<waveforms.size(); ++i) {
        std::vector<WindowStats> peaks =
            GetPeaks(waveforms[i], initial_derivative_threshold_, minimum_sample_length_);
        I3Vector<int> window_peak;
        I3Vector<int64_t> window_time;
        I3Vector<int64_t> window_length;
        for(size_t j=0; j<peaks.size(); ++j) {
            window_peak.push_back(peaks[j].peak);
            window_time.push_back(peaks[j].time);
            window_length.push_back(peaks[j].len);
        }
        window_peaks->operator[](i) = window_peak;
        window_times->operator[](i) = window_time;
        window_lengths->operator[](i) = window_length;
    }

    frame->Put("PeakAmplitudes", window_peaks);
    frame->Put("PeakTimes", window_times);
    frame->Put("PeakLengths", window_lengths);
}

void CCMTemplateFinder::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() == I3Frame::Geometry) {
        PushFrame(frame);
        return;
    }

    if(frame->GetStop() != I3Frame::DAQ) {
        PushFrame(frame);
        return;
    }

    n_daq_frames += 1;
    AddPeaks(frame);
    PushFrame(frame);

}

void CCMTemplateFinder::Finish() {
}
