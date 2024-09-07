#ifndef G4CCMReadout_h
#define G4CCMReadout_h

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3Map.h"
#include "simclasses/PhotonSummary.h"

#include <vector>
#include <string>
#include <sstream>

#include <G4SystemOfUnits.hh>

class G4CCMReadout {
public:
    G4CCMReadout() = default;
    G4CCMReadout(size_t n_threads) : readout(n_threads), primaries(n_threads), edep_trees(n_threads, nullptr)
    {}
    ~G4CCMReadout() = default;

    struct SingleEventReadout {
        I3MCTreePtr edep_tree;
        CCMMCPESeriesMapPtr mcpeseries;
        PhotonSummarySeriesPtr photon_summary_series;
        boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map;
    };

    typedef std::map<int, SingleEventReadout> SingleThreadReadout;

    typedef std::deque<SingleThreadReadout> Readout;

    void SetNumberOfThreads(size_t n_threads) { readout.resize(n_threads); }

    SingleThreadReadout & GetReadout(size_t thread_id) { return readout.at(thread_id); }

    void AddEntry(size_t thread_id, int event_id, I3MCTreePtr edep_tree, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map) {
        SingleThreadReadout & thread_readout = GetReadout(thread_id);
        SingleThreadReadout::iterator it = thread_readout.find(event_id);
        // Add a blank entry if it doesn't exist
        if (it != thread_readout.end()) {
            thread_readout.insert(it, {event_id, {edep_tree, nullptr, photon_summary_series, photon_summary_series_map}});
        } else {
            CCMMCPESeriesMapPtr mcpeseries = it->second.mcpeseries;
            UpdateMCPESeries(mcpeseries, edep_tree, photon_summary_series, photon_summary_series_map);
            it->second = {edep_tree, mcpeseries, nullptr, nullptr};
        }
    }

    void AddEntry(size_t thread_id, int event_id, CCMMCPESeriesMapPtr mcpeseries) {
        SingleThreadReadout & thread_readout = GetReadout(thread_id);
        SingleThreadReadout::iterator it = thread_readout.find(event_id);
        // Add a blank entry if it doesn't exist
        if (it != thread_readout.end()) {
            thread_readout.insert(it, {event_id, {nullptr, mcpeseries, nullptr, nullptr}});
        } else {
            I3MCTreePtr edep_tree = it->second.edep_tree;
            PhotonSummarySeriesPtr photon_summary_series = it->second.photon_summary_series;
            boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map = it->second.photon_summary_series_map;
            UpdateMCPESeries(mcpeseries, edep_tree, photon_summary_series, photon_summary_series_map);
            it->second = {edep_tree, mcpeseries, nullptr, nullptr};
        }
    }

    void AddEntry(size_t thread_id, int event_id, I3MCTreePtr edep_tree, CCMMCPESeriesMapPtr mcpeseries, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map) {
        SingleThreadReadout & thread_readout = GetReadout(thread_id);
        UpdateMCPESeries(mcpeseries, edep_tree, photon_summary_series, photon_summary_series_map);
        thread_readout[event_id] = {edep_tree, mcpeseries, nullptr, nullptr};
    }

    I3Particle GetPrimary(size_t i) const { return primaries.at(i); }
    void SetPrimary(size_t i, I3Particle primary) { primaries.at(i) = primary; }

    I3MCTreePtr GetEdepTree(size_t i) const { return edep_trees.at(i); }
    void SetEdepTree(size_t i, I3MCTreePtr edep_tree) { edep_trees.at(i) = edep_tree; }

private:
    Readout readout;
    std::vector<I3Particle> primaries;
    std::vector<I3MCTreePtr> edep_trees;
};

#endif // G4CCMReadout_h
