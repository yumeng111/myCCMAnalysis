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

typedef std::tuple<CCMPMTKey, CCMRecoPulse> PMTKeyPulsePair;
typedef std::vector<PMTKeyPulsePair> PMTKeyPulseVector;

std::vector<double> np_logspace(double log_low, double log_high, size_t n_edges, double base=10.0) {
    assert(n_edges > 1);
    if(n_edges == 1)
        return {std::pow(base, log_low)};
    if(n_edges == 2)
        return {std::pow(base, log_low), std::pow(base, log_high)};

    double width = (log_high - log_low) / (n_edges - 1);
    std::vector<double> bin_edges;
    bin_edges.reserve(n_edges);

    for(size_t i=0; i<n_edges; ++i) {
        bin_edges.push_back(std::pow(base, log_low + i * width));
    }

    return bin_edges;
}

class CountPulses: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    std::string pulses_name_;

    std::string output_prefix_;

    double log_low_;
    double log_high_;
    size_t n_edges_;
    double base_;

    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> charge_hists_;
    boost::shared_ptr<I3Map<CCMPMTKey, unsigned int>> num_triggers_;

public:
    void Geometry(I3FramePtr frame);
    CountPulses(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
    static inline void AddEntryLog(std::vector<double> & hist, double charge, double log_low, double log_high, size_t n_edges, double base=10.0);
};

I3_MODULE(CountPulses);

void CountPulses::AddEntryLog(std::vector<double> & hist, double charge, double log_low, double log_high, size_t n_edges, double base) {
    if(charge <= 0.0) {
        return;
    }
    if(charge < std::pow(base, log_low)) {
        return;
    }
    if(charge >= std::pow(base, log_high)) {
        return;
    }
    double width = (log_high - log_low) / (n_edges - 1);
    size_t bin = std::floor((std::log(charge) / std::log(base) - log_low) / width);
    hist.at(bin) += 1.0;
}

CountPulses::CountPulses(const I3Context& context) : I3Module(context),
    geo_seen(false), geometry_name_("") {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("InputPulsesName", "Name of the input pulses", std::string("WavedeformPulses"));
        AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
        AddParameter("LogLow", "Low edge of the log binning", double(-2.0));
        AddParameter("LogHigh", "High edge of the log binning", double(9.0));
        AddParameter("NEdges", "Number of edges in the log binning", size_t(10*11+1));
        AddParameter("Base", "Base of the log binning", double(10.0));
}

void CountPulses::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("InputPulsesName", pulses_name_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("LogLow", log_low_);
    GetParameter("LogHigh", log_high_);
    GetParameter("NEdges", n_edges_);
    GetParameter("Base", base_);

    charge_hists_ = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    num_triggers_ = boost::make_shared<I3Map<CCMPMTKey, unsigned int>>();
}

void CountPulses::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    if(geo_seen) {
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            charge_hists_->insert(std::make_pair(p.first, std::vector<double>(n_edges_, 0.0)));
            num_triggers_->insert(std::make_pair(p.first, 0));
        }
    }
    PushFrame(frame);
}

void CountPulses::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }
    
    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_name_);
    if(not pulses) {
        log_fatal("Could not find %s in the DAQ frame.", pulses_name_.c_str());
    }

    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin(); i != pulses->end(); i++) {
        std::vector<double> & charge_hist = charge_hists_->at(i->first);
        for(CCMRecoPulse const & pulse: i->second) {
            AddEntryLog(charge_hist, pulse.GetCharge(), log_low_, log_high_, n_edges_, base_);
        }
    }

    PushFrame(frame);
}

void CountPulses::Finish() {
    std::string output_prefix = output_prefix_;
    if(output_prefix.empty()) {
        output_prefix = pulses_name_;
    }

    std::vector<double> bin_edges = np_logspace(log_low_, log_high_, n_edges_, base_);
    for(std::pair<CCMPMTKey const, std::vector<double>> const & p : *charge_hists_) {
        std::stringstream ss;
        ss << output_prefix << "_pmt_" << p.first.GetRegion() << " " << p.first.GetSensor() << " " << p.first.GetSubsensor() << ".txt";
        std::string output_name = ss.str();
        std::ofstream output_file(output_name.c_str());
        if(not output_file) {
            log_fatal("Could not open file \"%s\" for writing.", output_name.c_str());
        }
        std::vector<double> const & charge_hist = p.second;
        output_file << bin_edges.at(0);
        for(size_t bin_idx=1; bin_idx<n_edges_; ++bin_idx) {
            output_file << " " << bin_edges.at(bin_idx);
        }

        output_file << charge_hist.at(0);
        for(size_t bin_idx=1; (bin_idx + 1) <n_edges_; ++bin_idx) {
            output_file << " " << charge_hist.at(bin_idx);
        }
        output_file << std::endl;
    }
}
