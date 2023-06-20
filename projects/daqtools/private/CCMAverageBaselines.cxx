// this modules reads in the interpolated baselines saved to the frame and finds the average across +- 25 frames
// it also does a wee bit of peak finding using cuts on the derivative
// so it returns average baselines in regions where we want to fit for SPEs
// and returns 0 for all over samples

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

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Orientation.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"


class OnlineRobustStatsBatched {
    std::deque<std::vector<double>> buffer;
    std::multiset<double> sorted_samples;
public:
    OnlineRobustStatsBatched() {}

    void AddValue(double x) {
        buffer.push_back({x});
        sorted_samples.insert(x);
    }

    void AddValues(std::vector<double> const & x) {
        buffer.push_back(x);
        sorted_samples.insert(x.begin(), x.end());
    }

    void AddValues(I3Vector<double> const & x) {
        buffer.push_back(x);
        sorted_samples.insert(x.begin(), x.end());
    }

    void RemoveValue() {
        if(buffer.size() == 0)
            return;
        {
            std::vector<double> const & to_remove = buffer.front();
            for(double const & x : to_remove) {
                sorted_samples.erase(x);
            }
        }
        buffer.pop_front();
    }   

    double Median() {
        const size_t N = sorted_samples.size();
        const size_t half = N / 2;
        std::multiset<double>::iterator it = sorted_samples.begin();
        for(size_t i=0; i<(half-1); ++i)
            ++it;

        // Odd count: return middle
        if(N % 2) {
            ++it;
            return *it;
        }

        // Even count: return average of middle two.
        double ret = *it;
        ++it;
        ret += *it;
        ret /= 2;
        return ret;
    }

    double Mode() {
        return robust_stats::Mode(sorted_samples.begin(), sorted_samples.end());
    }

    double Stddev(double median) {
        return robust_stats::MedianAbsoluteDeviation(
                sorted_samples.begin(),
                sorted_samples.end(),
                median);
    }
};


class CCMAverageBaselines : public I3Module {
    std::deque<I3FramePtr> frame_cache_;
    size_t n_daq_seen = 0;
    size_t n_to_cache = 25;
    bool geo_seen;
    std::string geometry_name_;
    std::map<CCMPMTKey, OnlineRobustStatsBatched> baseline_cache_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void EstimateCurrentBaselines(I3FramePtr frame);
    void RemoveBaselinesFromCache(I3FramePtr frame);
    void AddBaselinesToCache(I3FramePtr frame);
    void AddFrameToCache(I3FramePtr frame);
    void RemoveFrameFromCache(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
    void Process();
    void ProcessWaveform(CCMWaveformUInt16 const & waveform, I3Vector<double> const & baselines, I3Vector<double> const & derivs, I3Vector<double>& smoothed_wf);
public:
    CCMAverageBaselines(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(CCMAverageBaselines);

CCMAverageBaselines::CCMAverageBaselines(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
}


void CCMAverageBaselines::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
}


void CCMAverageBaselines::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);

    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;
    for(std::pair<CCMPMTKey const, uint32_t> p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        baseline_cache_.insert({key, OnlineRobustStatsBatched()});
    }

    geo_seen = true;
    PushFrame(frame);
}

void CCMAverageBaselines::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Vector<I3Vector<double>> baselines = frame->Get<I3Vector<I3Vector<double>>>("Baselines");
    I3Vector<I3Vector<double>> derivs = frame->Get<I3Vector<I3Vector<double>>>("Derivatives");

    // Place to store average baselines
    boost::shared_ptr<I3Vector<I3Vector<double>>> averageBaselines (new I3Vector<I3Vector<double>>(waveforms->size()));
    
    // loop over each channel / each waveform
    for(size_t i=0; i<waveforms->size(); ++i) {

        // get the CCMWaveform object from the vector
        CCMWaveformUInt16 const & waveform = waveforms->at(i);
        I3Vector<double> const & baseline = baselines[i];
        I3Vector<double> const & deriv = derivs[i];
        
	// pass the waveform to the function that computes the baselines
	ProcessWaveform(waveform, baseline, deriv, averageBaselines->at(i));
    }
    
 //   frame->Put("AverageBaselines", averageBaselines);
 //   PushFrame(frame);
}

void CCMAverageBaselines::EstimateCurrentBaselines(I3FramePtr frame){
    // estimate baseline mode and add to frame
    
    boost::shared_ptr<I3Map<CCMPMTKey, double>> baseline_mode = boost::make_shared<I3Map<CCMPMTKey, double>>();
    
    for(std::pair<CCMPMTKey const, uint32_t>  const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
	double mode = baseline_cache_[key].Mode();
        baseline_mode->emplace(std::make_pair(key, mode));
    }

    frame->Put("BaselineMode", baseline_mode);
    std::cout << "baseline mode  = " << baseline_mode << std::endl;
}


void CCMAverageBaselines::RemoveBaselinesFromCache(I3FramePtr frame) {
    // time to remove baselines from cache

    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
        baseline_cache_[key].RemoveValue();
    }
}


void CCMAverageBaselines::AddBaselinesToCache(I3FramePtr frame) {
    // need to get baselines from frame
    // then save to a map associated with key
    I3Vector<I3Vector<double>> const & baselines = frame->Get<I3Vector<I3Vector<double>> const>("Baselines");
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
        I3Vector<double> const & baseline_per_channel = baselines[channel];
	baseline_cache_[key].AddValues(baseline_per_channel);
    }
}

void CCMAverageBaselines::AddFrameToCache(I3FramePtr frame){
    frame_cache_.push_back(frame);
}

void CCMAverageBaselines::RemoveFrameFromCache(I3FramePtr frame){
    frame_cache_.pop_front();
}


void CCMAverageBaselines::Process(){
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() == I3Frame::Geometry) {
       Geometry(frame);
    }

    if(frame->GetStop() != I3Frame::DAQ) {
        PushFrame(frame);
        return;
    }
    
    ++n_daq_seen;

    // cache frame and baselines
    AddFrameToCache(frame);
    AddBaselinesToCache(frame);

    if(n_daq_seen < n_to_cache) {
        return;
    } else if(n_daq_seen == n_to_cache){
        I3FramePtr cached_frame;
        for(size_t i=0; i<frame_cache_.size(); ++i) {
            // Get the frame from the cache
            cached_frame = frame_cache_[i];
            // add baseline mode to frame
            EstimateCurrentBaselines(cached_frame);
            // push frame to next module
            PushFrame(cached_frame);
        }
        // pop one old frame from frame cache
        RemoveFrameFromCache(frame);
        // pop baselines in old frame from baseline cache
        RemoveBaselinesFromCache(frame);
        return;
    }

    // add baseline mode to frame
    EstimateCurrentBaselines(frame);
    // push frame to next module
    PushFrame(frame);

    // pop one old frame from frame cache
    RemoveFrameFromCache(frame);
    // pop baselines in old frame from baseline cache
    RemoveBaselinesFromCache(frame);
}


//void CCMAverageBaselines::FindRegionofSPEs(I3Vector<double> const & samples,I3Vector<double> const & derivs,I3Vector<double>& averageBaseline) {
//  // ok so we have average baselines and derivs
//  // let's find regions where we want to look for SPEs
//  // and fill all other regions with averageBaseline = 0
//  
//  double deriv_threshold = 0.3;
//  double wf_threshold_lower = 20;
//  double wf_threshold_upper = 40;
//  
//  // first let's loop over deriv
//  for(size_t deriv_idx = 0; deriv_idx<derivs.size(); ++deriv_idx){
//     
//     if(derivs[deriv_idx] < deriv_threshold && derivs[deriv_idx] > (-1*deriv_threshold) ){
//       // within threhold in the derivs
//
//       if(samples[deriv_idx] < (-1*averageBaseline[deriv_idx]+wf_threshold_upper) && samples[deriv_idx] > (-1*averageBaseline[deriv_idx]+wf_threshold_lower)){
//         // within threhold in the samples
//	 // want to also do a check on timing 
//	 // plan to add that in later
//
//	 //no change to average baseline
//	 continue; 
//         
//       }
//
//       else{
//         averageBaseline[deriv_idx] == 0;
//       }
//     
//     }
//
//    else{
//         averageBaseline[deriv_idx] == 0;
//       }
//  
//  }
//
//}


void CCMAverageBaselines::ProcessWaveform(CCMWaveformUInt16 const & waveform, I3Vector<double> const & baseline, I3Vector<double> const & deriv, I3Vector<double>& averageBaselines) {
    // get the vector of samples from the CCMWaveform object;
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();
    I3Vector<double> samples_double(samples.size());
   
    if (samples.size() == 0) {
    	I3Vector<double> average_baselines_empty (samples.size(), 0.0);
	averageBaselines = average_baselines_empty;
    	return;
    }

    for (size_t i = 0; i < samples.size(); ++i) {
    samples_double[i] = static_cast<double>(samples[i]);
    } 

//    FindRegionofSPEs(samples_double, deriv, averageBaselines);
}

void CCMAverageBaselines::Finish() {
}
