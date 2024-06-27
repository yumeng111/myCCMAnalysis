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

#include <boost/histogram.hpp>

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

namespace {
std::vector<double> np_linspace(double low, double high, size_t n_edges) {
    assert(n_edges > 1);
    if(n_edges == 1)
        return {low};
    if(n_edges == 2)
        return {low, high};

    double width = (high - low) / (n_edges - 1);
    std::vector<double> bin_edges;
    bin_edges.reserve(n_edges);

    for(size_t i=0; i<n_edges; ++i) {
        bin_edges.push_back(low + i * width);
    }

    return bin_edges;
}

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
}

class PulseChargeTimeHistogram: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    std::string pulses_name_;

    std::string output_prefix_;

    double time_low_;
    double time_high_;
    size_t time_n_edges_;

    double charge_log_low_;
    double charge_log_high_;
    size_t charge_n_edges_;
    double charge_base_;

    using storage_t = boost::histogram::dense_storage<boost::histogram::accumulators::weighted_sum<double>>;
    using time_axis_t =
        boost::histogram::axis::regular<double>;
    using charge_axis_t =
        boost::histogram::axis::regular<double, boost::histogram::axis::transform::log>;

    storage_t storage_;
    time_axis_t time_axis_;
    charge_axis_t charge_axis_;

    using axes_t = std::tuple<
        time_axis_t,
        charge_axis_t
    >;
    using hist_t = boost::histogram::histogram<axes_t, storage_t>;

    using time_axes_t = std::tuple<
        time_axis_t
    >;
    using charge_axes_t = std::tuple<
        charge_axis_t
    >;
    using time_hist_t = boost::histogram::histogram<time_axes_t, storage_t>;
    using charge_hist_t = boost::histogram::histogram<charge_axes_t, storage_t>;


    I3Map<CCMPMTKey, hist_t> hists_;
    I3Map<CCMPMTKey, time_hist_t> time_hists_;
    I3Map<CCMPMTKey, charge_hist_t> charge_hists_;
    I3Map<CCMPMTKey, time_hist_t> time_chargeW_hists_;

    hist_t new_hist() {
        return boost::histogram::make_histogram_with(
            storage_,
            time_axis_,
            charge_axis_
        );
    }
    time_hist_t new_time_hist() {
        return boost::histogram::make_histogram_with(
            storage_,
            time_axis_
        );
    }
    charge_hist_t new_charge_hist() {
        return boost::histogram::make_histogram_with(
            storage_,
            charge_axis_
        );
    }

public:
    void Geometry(I3FramePtr frame);
    PulseChargeTimeHistogram(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();

    /*
    static void SaveHists(
        std::string fname,
        size_t n_time_bins,
        size_t n_charge_bins,
        hist_t const & hist,
        time_hist_t const & time_hist,
        charge_hist_t const & charge_hist,
        time_hist_t const & time_charge_hist);
        */
    static void SaveBins(
        std::string fname,
        std::vector<double> const & time_bins,
        std::vector<double> const & charge_bins);
};

I3_MODULE(PulseChargeTimeHistogram);

template<typename Hist>
void Add1DHists(size_t nbins, Hist & hist, Hist const & other) {
    for(size_t i=0; i<nbins; ++i) {
        hist.at(i) += other.at(i);
    }
}

template<typename Hist>
void Add2DHists(size_t nbins_1, size_t nbins_2, Hist & hist, Hist const & other) {
    for(size_t i=0; i<nbins_1; ++i) {
        for(size_t j=0; j<nbins_2; ++j) {
            hist.at(i, j) += other.at(i, j);
        }
    }
}

template<typename Hist>
void Save1DHist(std::ostream & os, size_t nbins, Hist const & hist) {
    size_t first_nonzero = std::max(nbins, size_t(1))-1;
    for(size_t i=0; i<nbins; ++i) {
        if(hist.at(i).value() > 0.0) {
            first_nonzero = i;
            break;
        }
    }
    size_t last_nonzero = 0;
    for(size_t i=nbins; i>0; --i) {
        if(hist.at(i-1).value() > 0.0) {
            last_nonzero = i-1;
            break;
        }
    }
    os << first_nonzero << ":";
    for(size_t i=first_nonzero; i<=last_nonzero; ++i) {
        os << " " << hist.at(i).value();
    }
}

template<typename Hist>
void Save1DHistVar(std::ostream & os, size_t nbins, Hist const & hist) {
    size_t first_nonzero = std::max(nbins, size_t(1))-1;
    for(size_t i=0; i<nbins; ++i) {
        if(hist.at(i).variance() > 0.0) {
            first_nonzero = i;
            break;
        }
    }
    size_t last_nonzero = 0;
    for(size_t i=nbins; i>0; --i) {
        if(hist.at(i-1).variance() > 0.0) {
            last_nonzero = i-1;
            break;
        }
    }
    os << first_nonzero << ":";
    for(size_t i=first_nonzero; i<=last_nonzero; ++i) {
        os << " " << hist.at(i).variance();
    }
}

template<typename Hist>
void Save2DHist(std::ostream & os, size_t nbins_1, size_t nbins_2, Hist const & hist) {
    size_t first_nonzero = std::max(nbins_1, size_t(1))-1;
    bool found_first = false;
    for(size_t i=0; i<nbins_1; ++i) {
        for(size_t j=0; j<nbins_2; ++j) {
            if(hist.at(i, j).value() > 0.0) {
                first_nonzero = i;
                found_first = true;
                break;
            }
        }
        if(found_first)
            break;
    }
    size_t last_nonzero = 0;
    bool found_last = false;
    for(size_t i=nbins_1; i>0; --i) {
        for(size_t j=0; j<nbins_2; ++j) {
            if(hist.at(i-1, j).value() > 0.0) {
                last_nonzero = i-1;
                found_last = true;
                break;
            }
        }
        if(found_last)
            break;
    }
    os << first_nonzero << ":";
    for(size_t i=first_nonzero; i<=last_nonzero; ++i) {
        os << std::endl;
        for(size_t j=0; j<nbins_2; ++j) {
            if(j > 0)
                os << " ";
            os << hist.at(i, j).value();
        }
    }
}

template<typename Hist>
void Save1DHist(std::vector<double> & dest, size_t nbins, Hist const & hist) {
    dest.clear();
    dest.reserve(nbins);
    for(size_t i=0; i<nbins; ++i) {
        dest.push_back(hist.at(i).value());
    }
}

template<typename Hist>
void Save1DHistVar(std::vector<double> & dest, size_t nbins, Hist const & hist) {
    dest.clear();
    dest.reserve(nbins);
    for(size_t i=0; i<nbins; ++i) {
        dest.push_back(hist.at(i).variance());
    }
}

template<typename Hist>
void Save1DHist(I3Vector<double> & dest, size_t nbins, Hist const & hist) {
    dest.clear();
    dest.reserve(nbins);
    for(size_t i=0; i<nbins; ++i) {
        dest.push_back(hist.at(i).value());
    }
}

template<typename Hist>
void Save1DHistVar(I3Vector<double> & dest, size_t nbins, Hist const & hist) {
    dest.clear();
    dest.reserve(nbins);
    for(size_t i=0; i<nbins; ++i) {
        dest.push_back(hist.at(i).variance());
    }
}

template<typename Hist>
void Save2DHist(I3Vector<I3Vector<double>> & dest, size_t nbins_1, size_t nbins_2, Hist const & hist) {
    dest.clear();
    dest.reserve(nbins_1);
    for(size_t i=0; i<nbins_1; ++i) {
        dest.emplace_back();
        dest.back().reserve(nbins_2);
        for(size_t j=0; j<nbins_2; ++j) {
            dest.back().push_back(hist.at(i, j).value());
        }
    }
}

template<typename Hist>
void Save2DHist(std::vector<std::vector<double>> & dest, size_t nbins_1, size_t nbins_2, Hist const & hist) {
    dest.clear();
    dest.reserve(nbins_1);
    for(size_t i=0; i<nbins_1; ++i) {
        dest.emplace_back();
        dest.back().reserve(nbins_2);
        for(size_t j=0; j<nbins_2; ++j) {
            dest.back().push_back(hist.at(i, j).value());
        }
    }
}

void SaveValues(std::ostream & os, std::vector<double> const & bin_edges) {
    for(size_t i=0; i<bin_edges.size(); ++i) {
        if(i > 0)
            os << " ";
        os << bin_edges.at(i);
    }
}

void PulseChargeTimeHistogram::SaveBins(
        std::string fname,
        std::vector<double> const & time_bins,
        std::vector<double> const & charge_bins) {
    std::ofstream output_file(fname);
    if(not output_file) {
        log_fatal("Could not open file \"%s\" for writing.", fname.c_str());
    }

    output_file << "# time_bin_edges:" << std::endl;
    SaveValues(output_file, time_bins);
    output_file << std::endl;

    output_file << "# charge_bin_edges:" << std::endl;
    SaveValues(output_file, charge_bins);
    output_file << std::endl;
}

template<typename hist_t, typename time_hist_t, typename charge_hist_t>
void SaveHists(
        std::vector<std::vector<double>> & hist_dest,
        std::vector<double> & time_hist_dest,
        std::vector<double> & charge_hist_dest,
        std::vector<double> & time_charge_hist_dest,
        std::vector<double> & time_charge_hist_var_dest,
        size_t n_time_bins,
        size_t n_charge_bins,
        hist_t const & hist,
        time_hist_t const & time_hist,
        charge_hist_t const & charge_hist,
        time_hist_t const & time_charge_hist) {
    Save1DHist(time_hist_dest, n_time_bins, time_hist);
    Save1DHist(charge_hist_dest, n_charge_bins, charge_hist);
    Save1DHist(time_charge_hist_dest, n_time_bins, time_charge_hist);
    Save1DHistVar(time_charge_hist_var_dest, n_time_bins, time_charge_hist);
    Save2DHist(hist_dest, n_time_bins, n_charge_bins, hist);
}

template<typename hist_t, typename time_hist_t, typename charge_hist_t>
void SaveHists(
        I3Vector<I3Vector<double>> & hist_dest,
        I3Vector<double> & time_hist_dest,
        I3Vector<double> & charge_hist_dest,
        I3Vector<double> & time_charge_hist_dest,
        I3Vector<double> & time_charge_hist_var_dest,
        size_t n_time_bins,
        size_t n_charge_bins,
        hist_t const & hist,
        time_hist_t const & time_hist,
        charge_hist_t const & charge_hist,
        time_hist_t const & time_charge_hist) {
    Save1DHist(time_hist_dest, n_time_bins, time_hist);
    Save1DHist(charge_hist_dest, n_charge_bins, charge_hist);
    Save1DHist(time_charge_hist_dest, n_time_bins, time_charge_hist);
    Save1DHistVar(time_charge_hist_var_dest, n_time_bins, time_charge_hist);
    Save2DHist(hist_dest, n_time_bins, n_charge_bins, hist);
}

/*
void PulseChargeTimeHistogram::SaveHists(
        std::string fname,
        size_t n_time_bins,
        size_t n_charge_bins,
        hist_t const & hist,
        time_hist_t const & time_hist,
        charge_hist_t const & charge_hist,
        time_hist_t const & time_charge_hist) {
    std::ofstream output_file(fname);
    if(not output_file) {
        log_fatal("Could not open file \"%s\" for writing.", fname.c_str());
    }

    output_file << "# time_hist:" << std::endl;
    Save1DHist(output_file, n_time_bins, time_hist);
    output_file << std::endl;

    output_file << "# charge_hist:" << std::endl;
    Save1DHist(output_file, n_charge_bins, charge_hist);
    output_file << std::endl;

    output_file << "# charge_weighted_time_hist:" << std::endl;
    Save1DHist(output_file, n_time_bins, time_charge_hist);
    output_file << std::endl;

    output_file << "# charge_weighted_time_hist_variance:" << std::endl;
    Save1DHistVar(output_file, n_time_bins, time_charge_hist);
    output_file << std::endl;

    output_file << "# time_charge_hist" << std::endl;
    Save2DHist(output_file, n_time_bins, n_charge_bins, hist);
    output_file << std::endl;
}
*/

PulseChargeTimeHistogram::PulseChargeTimeHistogram(const I3Context& context) : I3Module(context),
    geo_seen(false), geometry_name_("") {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("InputPulsesName", "Name of the input pulses", std::string("WavedeformPulses"));
        AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
        AddParameter("TimeLow", "Low edge of the log binning", double(-16000));
        AddParameter("TimeHigh", "High edge of the log binning", double(16000));
        AddParameter("TimeNEdges", "Number of edges in the log binning", size_t(16000 + 1));
        AddParameter("ChargeLogLow", "Low edge of the log binning", double(-2.0));
        AddParameter("ChargeLogHigh", "High edge of the log binning", double(9.0));
        AddParameter("ChargeNEdges", "Number of edges in the log binning", size_t(10*11+1));
        AddParameter("ChargeBase", "Base of the log binning", double(10.0));
}

void PulseChargeTimeHistogram::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("InputPulsesName", pulses_name_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("TimeLow", time_low_);
    GetParameter("TimeHigh", time_high_);
    GetParameter("TimeNEdges", time_n_edges_);
    GetParameter("ChargeLogLow", charge_log_low_);
    GetParameter("ChargeLogHigh", charge_log_high_);
    GetParameter("ChargeNEdges", charge_n_edges_);
    GetParameter("ChargeBase", charge_base_);

    time_axis_ = time_axis_t(time_n_edges_, time_low_, time_high_);
    // charge_axis_ = charge_axis_t(charge_n_edges_, log(charge_base_) * charge_log_low_, log(charge_base_) * charge_log_high_);
    charge_axis_ = charge_axis_t(charge_n_edges_, std::pow(charge_base_, charge_log_low_), std::pow(charge_base_, charge_log_high_));
}

void PulseChargeTimeHistogram::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    if(geo_seen) {
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            hists_.insert(
                    std::make_pair(
                        p.first,
                        new_hist()
                        )
                    );
            time_hists_.insert(
                    std::make_pair(
                        p.first,
                        new_time_hist()
                        )
                    );
            charge_hists_.insert(
                    std::make_pair(
                        p.first,
                        new_charge_hist()
                        )
                    );
            time_chargeW_hists_.insert(
                    std::make_pair(
                        p.first,
                        new_time_hist()
                        )
                    );
        }
    }
    PushFrame(frame);
}

void PulseChargeTimeHistogram::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_name_);
    if(not pulses) {
        log_fatal("Could not find %s in the DAQ frame.", pulses_name_.c_str());
    }

    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin(); i != pulses->end(); i++) {
        hist_t & hist = hists_.at(i->first);
        time_hist_t & time_hist = time_hists_.at(i->first);
        charge_hist_t & charge_hist = charge_hists_.at(i->first);
        time_hist_t & time_chargeW_hist = time_chargeW_hists_.at(i->first);
        for(CCMRecoPulse const & pulse: i->second) {
            hist(boost::histogram::weight(1.0), pulse.GetTime(), pulse.GetCharge());
            time_hist(boost::histogram::weight(1.0), pulse.GetTime());
            charge_hist(boost::histogram::weight(1.0), pulse.GetCharge());
            time_chargeW_hist(boost::histogram::weight(pulse.GetCharge()), pulse.GetTime());
        }
    }

    PushFrame(frame);
}

void PulseChargeTimeHistogram::Finish() {
    std::string output_prefix = output_prefix_;
    if(output_prefix.empty()) {
        output_prefix = pulses_name_;
    }

    std::vector<double> charge_bin_edges = np_logspace(charge_log_low_, charge_log_high_, charge_n_edges_, charge_base_);
    std::vector<double> time_bin_edges = np_linspace(time_low_, time_high_, time_n_edges_);
    size_t ntime = time_bin_edges.size();
    size_t ncharge = charge_bin_edges.size();

    hist_t tot_hist = new_hist();
    time_hist_t tot_time_hist = new_time_hist();
    charge_hist_t tot_charge_hist = new_charge_hist();
    time_hist_t tot_time_chargeW_hist = new_time_hist();

    //std::stringstream ss;
    //ss << output_prefix << "_binning.txt";

    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::Physics);

    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> time_hists = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> charge_hists = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> time_chargeW_hists = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> time_chargeW_hists_var = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<std::vector<double>>>> time_charge_hists = boost::make_shared<I3Map<CCMPMTKey, std::vector<std::vector<double>>>>();

    //std::string binning_output_name = ss.str();

    //SaveBins(binning_output_name, time_bin_edges, charge_bin_edges);

    for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
        Add2DHists(ntime, ncharge, tot_hist, hists_.at(p.first));
        Add1DHists(ntime, tot_time_hist, time_hists_.at(p.first));
        Add1DHists(ncharge, tot_charge_hist, charge_hists_.at(p.first));
        Add1DHists(ntime, tot_time_chargeW_hist, time_chargeW_hists_.at(p.first));

        time_charge_hists->insert(std::make_pair(p.first, std::vector<std::vector<double>>()));
        time_hists->insert(std::make_pair(p.first, std::vector<double>()));
        charge_hists->insert(std::make_pair(p.first, std::vector<double>()));
        time_chargeW_hists->insert(std::make_pair(p.first, std::vector<double>()));
        time_chargeW_hists_var->insert(std::make_pair(p.first, std::vector<double>()));
        SaveHists(
            time_charge_hists->at(p.first),
            time_hists->at(p.first),
            charge_hists->at(p.first),
            time_chargeW_hists->at(p.first),
            time_chargeW_hists_var->at(p.first),
            ntime,
            ncharge,
            hists_.at(p.first),
            time_hists_.at(p.first),
            charge_hists_.at(p.first),
            time_chargeW_hists_.at(p.first)
        );
    }

    boost::shared_ptr<I3Vector<I3Vector<double>>> tot_time_charge_output = boost::make_shared<I3Vector<I3Vector<double>>>();
    boost::shared_ptr<I3Vector<double>> tot_time_output = boost::make_shared<I3Vector<double>>();
    boost::shared_ptr<I3Vector<double>> tot_charge_output = boost::make_shared<I3Vector<double>>();
    boost::shared_ptr<I3Vector<double>> tot_time_chargeW_output = boost::make_shared<I3Vector<double>>();
    boost::shared_ptr<I3Vector<double>> tot_time_chargeW_output_var = boost::make_shared<I3Vector<double>>();

    SaveHists(
        *tot_time_charge_output,
        *tot_time_output,
        *tot_charge_output,
        *tot_time_chargeW_output,
        *tot_time_chargeW_output_var,
        ntime,
        ncharge,
        tot_hist,
        tot_time_hist,
        tot_charge_hist,
        tot_time_chargeW_hist);

    frame->Put(output_prefix + "TimeBinEdges", boost::make_shared<I3VectorDouble>(time_bin_edges.begin(), time_bin_edges.end()));
    frame->Put(output_prefix + "ChargeBinEdges", boost::make_shared<I3VectorDouble>(charge_bin_edges.begin(), charge_bin_edges.end()));
    frame->Put(output_prefix + "TimeChargeHists", tot_time_charge_output);
    frame->Put(output_prefix + "TimeHists", tot_time_output);
    frame->Put(output_prefix + "ChargeHists", tot_charge_output);
    frame->Put(output_prefix + "TimeChargeWHists", tot_time_chargeW_output);
    frame->Put(output_prefix + "TimeChargeWHistsVar", tot_time_chargeW_output_var);

    PushFrame(frame);
}
