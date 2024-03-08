
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

#include <icetray/ctpl.h>
#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

struct ElectronicsCorrectionJob {
    std::atomic<bool> running = false;
    std::thread thread;
    size_t thread_index = 0;
    I3FramePtr frame = nullptr;
    size_t frame_index = 0;
};

struct ElectronicsCorrectionResult {
    I3FramePtr frame = nullptr;
    bool done = false;
};

class  ElectronicsCorrection: public I3Module {
    std::exception_ptr teptr = nullptr;
    std::string geometry_name_;
    std::string pmt_channel_map_name_;
    bool geo_seen;
    bool calib_seen;
    std::string ccm_calibration_name_;
    std::string ccm_waveforms_name_;
    std::string baseline_estimates_name_;
    std::string output_name_;
    std::vector<double> default_droop_tau_;
    CCMPMTCalibration default_calib_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    CCMCalibration ccm_calibration_;

    bool remove_waveforms;
    size_t num_threads;
    size_t max_cached_frames;

    // initialize droop correction time constant map
    double delta_t = 2.0;
    //void BoxFilter(std::vector<double> const & samples, std::vector<double> & box_filter_results);
    //void OutlierFilter(std::vector<double> const & samples, std::vector<double> & outlier_filter_results);
    // void ECorrection(std::vector<uint16_t> const & samples, double const & mode, const CCMPMTCalibration& calibration,  std::vector<double> & electronics_correction_samples);
    size_t frame_index = 0;
    size_t min_frame_idx = 0;
    std::deque<ElectronicsCorrectionJob *> free_jobs;
    std::deque<ElectronicsCorrectionJob *> running_jobs;
    std::deque<ElectronicsCorrectionResult> results;
public:
    ElectronicsCorrection(const I3Context&);
    void Configure();
    void Process();
    void Finish();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
    void Calibration(I3FramePtr frame);
};

I3_MODULE( ElectronicsCorrection);

ElectronicsCorrection:: ElectronicsCorrection(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false), calib_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("PMTChannelMapName", "Key for PMTChannelMap", std::string(""));
        AddParameter("CCMCalibrationName", "Key for CCMCalibration", std::string("CCMCalibration"));
        AddParameter("CCMWaveformsName", "Key for input CCMWaveforms", std::string("CCMWaveforms"));
        AddParameter("BaselineEstimatesName", "Key for input BaselineEstimates", std::string("BaselineEstimates"));
        AddParameter("DefaultDroopTau", "The default droop time constant (in ns) to use if there is no entry in the calibratioin", std::vector((double)1200, (double)5000));
        AddParameter("OutputName", "Key to save output CCMWaveformDoubleSeries to", std::string("CCMCalibratedWaveforms"));
        AddParameter("NumThreads", "Number of worker threads to use for baseline estimation", (size_t)(0));
        AddParameter("MaxCachedFrames", "The maximum number of frames this module is allowed to have cached", (size_t)(1000));
        AddParameter("RemoveWaveforms", "Remove the input waveforms?", bool(false));
}


void  ElectronicsCorrection::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTChannelMapName", pmt_channel_map_name_);
    GetParameter("CCMCalibrationName", ccm_calibration_name_);
    GetParameter("CCMWaveformsName", ccm_waveforms_name_);
    GetParameter("BaselineEstimatesName", baseline_estimates_name_);
    GetParameter("DefaultDroopTau", default_droop_tau_);
    GetParameter("OutputName", output_name_);
    default_calib_.SetDroopTimeConstant(default_droop_tau_);
    GetParameter("NumThreads", num_threads);
    GetParameter("MaxCachedFrames", max_cached_frames);
    GetParameter("RemoveWaveforms", remove_waveforms);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
}


void  ElectronicsCorrection::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    if(pmt_channel_map_name_ != "") {
        if(not frame->Has(pmt_channel_map_name_)) {
            log_fatal("Could not find I3Map<CCMPMTKey, uint32_t> object with the key named \"%s\" in the Geometry frame.", pmt_channel_map_name_.c_str());
        }
        boost::shared_ptr<I3Map<CCMPMTKey, uint32_t> const> pmt_channel_map = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, uint32_t> const>>(pmt_channel_map_name_);
        boost::shared_ptr<I3Map<CCMPMTKey, int> const> pmt_channel_map_alt = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, int> const>>(pmt_channel_map_name_);
        if(pmt_channel_map != nullptr) {
            pmt_channel_map_ = *pmt_channel_map;
        } else if(pmt_channel_map_alt != nullptr) {
            for(std::pair<CCMPMTKey const, int> const & it : *pmt_channel_map_alt) {
                pmt_channel_map_[it.first] = uint32_t(it.second);
            }
        } else {
            log_fatal("Could not find I3Map<CCMPMTKey, uint32_t> object with the key named \"%s\" in the Geometry frame.", pmt_channel_map_name_.c_str());
        }
    } else {
        CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
        pmt_channel_map_ = geo.pmt_channel_map;
    }
    geo_seen = true;
    PushFrame(frame);
}

void ElectronicsCorrection::Calibration(I3FramePtr frame) {
    if(not frame->Has(ccm_calibration_name_)) {
        log_fatal("Could not find CCMCalibration object with the key named \"%s\" in the Calibration frame.", ccm_calibration_name_.c_str());
    }
    ccm_calibration_ = frame->Get<CCMCalibration>(ccm_calibration_name_);
    calib_seen = true;
    PushFrame(frame);
}

void OutlierFilter(std::vector<double> const & samples, std::vector<double> & outlier_filter_results) {
    double delta_tau = 20;
    double prev_tau = 2.0;
    double next_tau = 2.0;

    // first let's find the mode of the first 800 bins of our wf as the starting value
    std::vector<double> starting_samples(samples.begin(), samples.begin() + std::min(size_t(800), samples.size()));
    std::sort(starting_samples.begin(), starting_samples.end());
    double value = robust_stats::Mode(starting_samples.begin(), starting_samples.end());

    // now let's loop over the waveform
    for(size_t wf_it = 0; wf_it < samples.size(); ++wf_it) {
        double delta = samples[wf_it] - value;
        double e = std::fabs(delta) / delta_tau;

        if(wf_it > 0) {
            e += std::fabs(samples[wf_it] - samples[wf_it - 1]) / prev_tau;
        }

        if(wf_it + 1 < samples.size()) {
            e += std::fabs(samples[wf_it] - samples[wf_it + 1]) / next_tau;
        }
        delta *= std::exp(-e);
        value += delta;
        outlier_filter_results[wf_it] = value;
    }
}

void BoxFilter(std::vector<double> const & samples, std::vector<double> & box_filter_results) {
    for(size_t i = 0; i < samples.size(); ++i) {
        // special logic for first element
        if(i == 0)
            box_filter_results[i] = (samples[i] + samples[i+1] + samples[i+2])/3;
        else if (i == samples.size() - 1)
            // special logic for last element
            box_filter_results[i] = (samples[i] + samples[i-1] + samples[i-2])/3;
        else
            box_filter_results[i] = (samples[i-1] + samples[i] + samples[i+1])/3;
    }
}

void ECorrection(std::vector<uint16_t> const & samples, double const & mode, const CCMPMTCalibration& calibration,  std::vector<double> & electronics_correction_samples, std::vector<double> default_droop_tau_, double delta_t) {

    // let's start off by inverting samples and subtracting off the baseline
    std::vector<double> inv_wf_min_baseline(samples.size());

    for (size_t wf_it = 0; wf_it < samples.size(); ++wf_it) {
        inv_wf_min_baseline[wf_it] = -1 * (double(samples[wf_it]) + mode); // mode is negative and samples is positive!
    }

    size_t N = samples.size();

    // now let's run inv_wf_min_baseline through our box filter
    std::vector<double> box_filter_results(N);
    BoxFilter(inv_wf_min_baseline, box_filter_results);

    // now droop correction time
    std::vector<double> tau = calibration.GetDroopTimeConstant();
    if(std::isnan(tau.at(0))) {
        tau = default_droop_tau_;
    }

    // now loop over tau

    double C;
    double B;
    double A;
    double S;
    double X;

    // place to store our results from droop correcting once
    std::vector<double> first_droop_correction_results(N);
    electronics_correction_samples.resize(N);

    std::vector<double *> inputs = {box_filter_results.data(), first_droop_correction_results.data()};
    std::vector<double *> outputs = {first_droop_correction_results.data(), electronics_correction_samples.data()};

    for (size_t tau_it = 0; tau_it < 2; ++tau_it) {
        double * input = inputs[tau_it];
        double * output = outputs[tau_it];

        // define our numbers on the first iteration
        C = std::exp(-delta_t/tau.at(tau_it));
        B = (1.0 - C);
        A = tau.at(tau_it) / delta_t * B;
        S = 0.0;

        X = double(input[0]) / A + B * S;
        output[0] = X;

        for(size_t i=1; i<N; ++i) {
            S = X + C * S;
            X = input[i] / A + B * S;
            output[i] = X;
        }
    }

    // now let's subtract off our outlier filter
    // place to store results of the outlier filter
    std::vector<double> outlier_filter_results(samples.size());
    OutlierFilter(electronics_correction_samples, outlier_filter_results);

    // now let's subtract off
    for(size_t it = 0; it < samples.size(); ++it) {
        electronics_correction_samples[it] -= outlier_filter_results[it];
    }
}

int ProcessWaveform(
        std::pair<CCMPMTKey const, BaselineEstimate> const & it,
        I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map_,
        CCMCalibration const & calibration,
        CCMWaveformUInt16Series const & waveforms,
        CCMWaveformDoubleSeries & electronics_corrected_wf,
        CCMPMTCalibration const & default_calib_,
        std::vector<double> default_droop_tau_,
        double delta_t) {
    CCMPMTKey key = it.first;
    BaselineEstimate value = it.second;
    double mode = value.baseline; // baseline mode is negative fyi
    uint32_t channel = pmt_channel_map_.at(key);

    std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);

    CCMWaveformUInt16 const & waveform = waveforms.at(channel);

    // get the vector of samples from the CCMWaveform object;
    std::vector<uint16_t> const & samples = waveform.GetWaveform();

    // no work to be done if there is not waveform available
    if (samples.size() == 0) {
        return 0;
    }

    // place to store our waveform after electronics correction
    std::vector<double> & electronics_correction_samples = electronics_corrected_wf.at(channel).GetWaveform();
    electronics_correction_samples.reserve(samples.size());

    // let's pass it to electronics correction module
    // this module first droop corrects
    // then substracts off the outlier filter
    if(calib == calibration.pmtCal.cend()) {
        ECorrection(samples, mode, default_calib_, electronics_correction_samples, default_droop_tau_, delta_t);
    } else {
        ECorrection(samples, mode, calib->second, electronics_correction_samples, default_droop_tau_, delta_t);
    }
    return 0;
}

void ElectronicsCorrection::Process() {
  if (inbox_)
    log_trace("%zu frames in inbox", inbox_->size());

  I3FramePtr frame = PopFrame();

  if(!frame or frame->GetStop() == I3Frame::DAQ) {
    DAQ(frame);
    return;
  }

  if(frame->GetStop() == I3Frame::Physics && ShouldDoPhysics(frame))
    {
      Physics(frame);
    }
  else if(frame->GetStop() == I3Frame::Geometry && ShouldDoGeometry(frame))
    Geometry(frame);
  else if(frame->GetStop() == I3Frame::Calibration && ShouldDoCalibration(frame))
    Calibration(frame);
  else if(frame->GetStop() == I3Frame::DetectorStatus && ShouldDoDetectorStatus(frame))
    DetectorStatus(frame);
  else if(frame->GetStop() == I3Frame::Simulation && ShouldDoSimulation(frame))
    Simulation(frame);
  else if(frame->GetStop() == I3Frame::DAQ && ShouldDoDAQ(frame)) {
    DAQ(frame);
  } else if(ShouldDoOtherStops(frame))
    OtherStops(frame);
}

void FrameThread(std::atomic<bool> & running, I3Frame * frame, CCMCalibration const & ccm_calibration_, std::string const & ccm_waveforms_name_, std::string const & baseline_estimates_name_, std::string const & output_name_, bool geo_seen, bool calib_seen, I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map_, CCMPMTCalibration const & default_calib_, std::vector<double> const & default_droop_tau_, double delta_t, bool remove_waveforms) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen before DAQ frame! Did you forget to include a geometry file?");
    }
    if(not calib_seen) {
        log_fatal("No Calibration frame seen before DAQ frame! Did you forget to inlcude a calibration file?");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    CCMCalibration const & calibration = ccm_calibration_;

    boost::shared_ptr<CCMWaveformUInt16Series const> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>(ccm_waveforms_name_);
    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const> baseline_mode = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const>>(baseline_estimates_name_);
    boost::shared_ptr<I3Map<uint32_t, BaselineEstimate> const> baseline_mode_by_channel = frame->Get<boost::shared_ptr<I3Map<uint32_t, BaselineEstimate> const>>(baseline_estimates_name_);
    if(waveforms == nullptr) {
        log_fatal("No CCMWaveformUInt16Series under key name \"%s\"", ccm_waveforms_name_.c_str());
    }

    size_t size = waveforms->size();
    // a vector storing the electronics-corrected wfs for each channel
    boost::shared_ptr<CCMWaveformDoubleSeries> electronics_corrected_wf = boost::make_shared<CCMWaveformDoubleSeries>(size);
    if(baseline_mode != nullptr) {
        // loop over each pmt
        // We loop over pmt keys in baseline_estimates because a baseline estimate is required for the electronics correction
        // and the baseline_estimates object should have its pmt keys derived from the same geometry
        for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : *baseline_mode) {
            ProcessWaveform(
                it,
                pmt_channel_map_,
                calibration,
                std::cref(*waveforms),
                std::ref(*electronics_corrected_wf),
                default_calib_,
                default_droop_tau_,
                delta_t
            );
        }
    } else if(baseline_mode_by_channel != nullptr) {
        std::map<uint32_t, CCMPMTKey> reverse_channel_map;
        for(std::pair<CCMPMTKey const, uint32_t> const & it : pmt_channel_map_) {
            reverse_channel_map[it.second] = it.first;
        }
        // loop over each pmt
        // We loop over pmt keys in baseline_estimates because a baseline estimate is required for the electronics correction
        // and the baseline_estimates object should have its pmt keys derived from the same geometry
        for(std::pair<uint32_t const, BaselineEstimate> const & it : *baseline_mode_by_channel) {
            std::pair<CCMPMTKey const, BaselineEstimate> fake_it = std::make_pair(reverse_channel_map.at(it.first), it.second);
            ProcessWaveform(
                fake_it,
                pmt_channel_map_,
                calibration,
                std::cref(*waveforms),
                std::ref(*electronics_corrected_wf),
                default_calib_,
                default_droop_tau_,
                delta_t
            );
        }
    } else {
        log_fatal("No I3Map<CCMPMTKey, BaselineEstimate> or I3Vector<BaselineEstimate> under key name \"%s\"", baseline_estimates_name_.c_str());
    }

    double total_adc_count = 0.0;
    for(size_t i=0; i<electronics_corrected_wf->size(); ++i) {
        std::vector<double> const & wf = electronics_corrected_wf->at(i).GetWaveform();
        for(size_t j=0; j<wf.size(); ++j) {
            total_adc_count += wf.at(j);
        }
    }

    if(remove_waveforms) {
        frame->Delete(ccm_waveforms_name_);
    }
    frame->Put(output_name_, electronics_corrected_wf);
    frame->Put(output_name_ + "TotalADC", boost::make_shared<I3Double>(total_adc_count));
    running.store(false);
}

void RunFrameThread(ElectronicsCorrectionJob * job,
        CCMCalibration const & ccm_calibration_,
        std::string const & ccm_waveforms_name_,
        std::string const & baseline_estimates_name_,
        std::string const & output_name_,
        bool geo_seen,
        bool calib_seen,
        I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map_,
        CCMPMTCalibration const & default_calib_,
        std::vector<double> const & default_droop_tau_,
        double delta_t,
        bool remove_waveforms) {

    job->running.store(true);

    job->thread = std::thread(FrameThread,
        std::ref(job->running),
        job->frame.get(),
        std::cref(ccm_calibration_),
        std::cref(ccm_waveforms_name_),
        std::cref(baseline_estimates_name_),
        std::cref(output_name_),
        geo_seen,
        calib_seen,
        std::cref(pmt_channel_map_),
        std::cref(default_calib_),
        std::cref(default_droop_tau_),
        delta_t,
        remove_waveforms
    );
}

void ElectronicsCorrection::DAQ(I3FramePtr frame) {
	while(true) {
		// Check if any jobs have finished
		for(int i=int(running_jobs.size())-1; i>=0; --i) {
			if (teptr) {
				try{
					std::rethrow_exception(teptr);
				}
				catch(const std::exception &ex)
				{
					std::cerr << "Thread exited with exception: " << ex.what() << "\n";
				}
			}
			if(not running_jobs[i]->running.load()) {
				ElectronicsCorrectionJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                results[job->frame_index - min_frame_idx].done = true;
            } else {
                ElectronicsCorrectionJob * job = running_jobs[i];
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                PushFrame(results[i].frame);
                results[i].frame = nullptr;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_frame_idx += results_done;
        }

        if(not frame)
            break;

        // Attempt to queue up a new job for the frame
        ElectronicsCorrectionJob * job = nullptr;

        if(free_jobs.size() > 0) {
            job = free_jobs.front();
            job->running.store(false);
            free_jobs.pop_front();
        } else if(running_jobs.size() < num_threads) {
            job = new ElectronicsCorrectionJob();
            job->running.store(false);
            job->thread_index = running_jobs.size();
        }

        if(job != nullptr and results.size() < max_cached_frames) {
            job->running.store(true);
            running_jobs.push_back(job);
            job->frame = frame;
            job->frame_index = frame_index;
            results.emplace_back();
            results.back().frame = frame;
            results.back().done = false;
            frame_index += 1;
            RunFrameThread(job,
                ccm_calibration_,
                ccm_waveforms_name_,
                baseline_estimates_name_,
                output_name_,
                geo_seen,
                calib_seen,
                pmt_channel_map_,
                default_calib_,
                default_droop_tau_,
                delta_t,
                remove_waveforms);
            break;
        } else if(job != nullptr) {
            free_jobs.push_back(job);
        }
    }
}

void ElectronicsCorrection::Finish() {
    while(running_jobs.size() > 0) {
        // Check if any jobs have finished
        for(int i=int(running_jobs.size())-1; i>=0; --i) {
            if(not running_jobs[i]->running.load()) {
                ElectronicsCorrectionJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                results[job->frame_index - min_frame_idx].done = true;
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                PushFrame(results[i].frame);
                results[i].frame = nullptr;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_frame_idx += results_done;
        }
    }
}


