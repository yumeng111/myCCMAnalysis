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
#include <dataclasses/calibration/BaselineEstimate.h>

typedef std::tuple<CCMPMTKey, CCMRecoPulse> PMTKeyPulsePair;
typedef std::vector<PMTKeyPulsePair> PMTKeyPulseVector;

class BaselineCheck: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    std::string pulses_name_;
    std::string baseline_name_;

    std::string output_prefix_;

    std::map<CCMPMTKey, std::ofstream> output_files_;

public:
    void Geometry(I3FramePtr frame);
    BaselineCheck(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
};

I3_MODULE(BaselineCheck);

BaselineCheck::BaselineCheck(const I3Context& context) : I3Module(context),
    geo_seen(false), geometry_name_("") {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("InputPulsesName", "Name of the input pulses", std::string("WavedeformPulses"));
        AddParameter("BaselineEstimatesName", "Name of the output baseline estimates", std::string("BaselineEstimates"));
        AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
}

void BaselineCheck::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("InputPulsesName", pulses_name_);
    GetParameter("BaselineEstimatesName", baseline_name_);
    GetParameter("OutputPrefix", output_prefix_);
}

void BaselineCheck::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
        std::stringstream ss;
        ss << output_prefix_ << "_baselines_pmt_" << p.first.GetRegion() << "_" << p.first.GetSensor();
        ss << "_" << static_cast<uint64_t>(p.first.GetSubsensor()) << ".txt";
        std::string output_name = ss.str();
        output_files_.emplace(p.first, std::ofstream(output_name.c_str()));
        output_files_.at(p.first) << "# num_pulses total_charge baseline baseline_stddev" << std::endl;
    }
    PushFrame(frame);
}

void BaselineCheck::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }
    
    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_name_);
    if(not pulses) {
        log_fatal("Could not find %s in the DAQ frame.", pulses_name_.c_str());
    }
    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const> baselines = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const>>(baseline_name_);
    if(not baselines) {
        log_fatal("Could not find %s in the DAQ frame.", baseline_name_.c_str());
    }

    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin(); i != pulses->end(); i++) {
        CCMPMTKey const & pmt_key = i->first;
        double total_charge = 0.0;
        for(CCMRecoPulse const & pulse: i->second) {
            total_charge += pulse.GetCharge();
        }
        double baseline = baselines->at(pmt_key).baseline;
        double num_pulses = i->second.size();
        double baseline_stddev = baselines->at(pmt_key).stddev;

        output_files_.at(pmt_key) << num_pulses << " " << total_charge << " " << baseline << " " << baseline_stddev << std::endl;
    }
    PushFrame(frame);
}
