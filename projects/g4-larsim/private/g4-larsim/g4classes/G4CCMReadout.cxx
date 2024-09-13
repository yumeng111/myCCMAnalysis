#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3Map.h"
#include "simclasses/PhotonSummary.h"
#include "simclasses/CCMMCPE.h"

#include "g4-larsim/g4classes/G4CCMReadout.h"

#include <vector>
#include <string>
#include <sstream>

#include <G4SystemOfUnits.hh>

const std::unordered_map<PhotonSummary::PhotonSource, CCMMCPE::PhotonSource> G4CCMReadout::PhotonSummarytoCCMMCPEPhotonSource = {
    {PhotonSummary::PhotonSource::Unknown, CCMMCPE::PhotonSource::Unknown},
    {PhotonSummary::PhotonSource::Scintillation, CCMMCPE::PhotonSource::Scintillation},
    {PhotonSummary::PhotonSource::Cerenkov, CCMMCPE::PhotonSource::Cerenkov}
};


G4CCMReadout::G4CCMReadout(size_t n_threads) : n_threads_(n_threads), readout(n_threads), primaries(), edep_trees() {}

void G4CCMReadout::SetInput(std::vector<I3Particle> primaries, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> edep_trees) {
    this->primaries = primaries;
    this->mcpeseries = mcpeseries;
    this->edep_trees = edep_trees;
}

void G4CCMReadout::SetNumberOfThreads(size_t n_threads) {
    n_threads_ = n_threads;
    readout.resize(n_threads_);
}

G4CCMReadout::SingleThreadReadout & G4CCMReadout::GetReadout(size_t thread_id) { return readout.at(thread_id); }

void G4CCMReadout::AddEntry(size_t thread_id, int event_id, I3MCTreePtr edep_tree, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map) {
    SingleThreadReadout & thread_readout = GetReadout(thread_id);
    SingleThreadReadout::iterator it = thread_readout.find(event_id);
    // Add a blank entry if it doesn't exist
    if (it == thread_readout.end()) {
        thread_readout.insert(it, {event_id, {edep_tree, nullptr, photon_summary_series, photon_summary_series_map}});
    } else {
        CCMMCPESeriesMapPtr mcpeseries = it->second.mcpeseries;
        UpdateMCPESeries(mcpeseries, photon_summary_series, photon_summary_series_map);
        it->second = {edep_tree, mcpeseries, nullptr, nullptr};
    }
}

void G4CCMReadout::AddEntry(size_t thread_id, int event_id, CCMMCPESeriesMapPtr mcpeseries) {
    SingleThreadReadout & thread_readout = GetReadout(thread_id);
    SingleThreadReadout::iterator it = thread_readout.find(event_id);
    // Add a blank entry if it doesn't exist
    if (it == thread_readout.end()) {
        thread_readout.insert(it, {event_id, {nullptr, mcpeseries, nullptr, nullptr}});
    } else {
        I3MCTreePtr edep_tree = it->second.edep_tree;
        PhotonSummarySeriesPtr photon_summary_series = it->second.photon_summary_series;
        boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map = it->second.photon_summary_series_map;
        UpdateMCPESeries(mcpeseries, photon_summary_series, photon_summary_series_map);
        it->second = {edep_tree, mcpeseries, nullptr, nullptr};
    }
}

void G4CCMReadout::AddEntry(size_t thread_id, int event_id, I3MCTreePtr edep_tree, CCMMCPESeriesMapPtr mcpeseries, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map) {
    SingleThreadReadout & thread_readout = GetReadout(thread_id);
    UpdateMCPESeries(mcpeseries, photon_summary_series, photon_summary_series_map);
    thread_readout[event_id] = {edep_tree, mcpeseries, photon_summary_series, photon_summary_series_map};
}

I3Particle G4CCMReadout::GetPrimary(size_t i) const { return primaries.at(i); }
void G4CCMReadout::SetPrimary(size_t i, I3Particle primary) { primaries.at(i) = primary; }

I3MCTreePtr G4CCMReadout::GetMCTree(size_t i) const { return edep_trees.at(i); }
void G4CCMReadout::SetMCTree(size_t i, I3MCTreePtr edep_tree) { edep_trees.at(i) = edep_tree; }

CCMMCPESeriesMapPtr G4CCMReadout::GetMCPESeries(size_t i) const { return mcpeseries.at(i); }
void G4CCMReadout::SetMCPESeries(size_t i, CCMMCPESeriesMapPtr mcpeseries) { this->mcpeseries.at(i) = mcpeseries; }

void G4CCMReadout::Reset() {
    readout.clear();
    primaries.clear();
    edep_trees.clear();

    readout = Readout(n_threads_);
    primaries = std::vector<I3Particle>();
    edep_trees = std::vector<I3MCTreePtr>();
}

void G4CCMReadout::UpdateMCPESeries(CCMMCPESeriesMapPtr mcpeseries, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map) {
    // Iterate over PMTs in source map
    for(CCMMCPESeriesMap::iterator it = mcpeseries->begin(); it != mcpeseries->end(); ++it) {
        // Reference to the PE series
        CCMMCPESeries & pe_series = it->second;

        // Iterate backwards over the vector of CCMMCPE in the source map for this PMT
        for(int i=pe_series.size()-1; i>=0; --i) {
            CCMMCPE & pe = pe_series[i];

            // Check if the photon has everything properly recorded
            I3Map<int, size_t>::iterator it_map = photon_summary_series_map->find(pe.track_id);
            // If not trash it as an edge case
            if(it_map == photon_summary_series_map->end()) {
                pe_series.erase(pe_series.begin() + i);
                continue;
            }

            // Grab the summary information for this photon track
            PhotonSummary const & this_photon_summary = photon_summary_series->at(it_map->second);

            // Update the destination CCMMCPE with the summary information
            pe.g4_time = this_photon_summary.g4_time;
            pe.calculated_time = this_photon_summary.calculated_time;
            pe.g4_distance_uv = this_photon_summary.g4_distance_uv;
            pe.g4_distance_visible = this_photon_summary.g4_distance_visible;
            pe.calculated_distance_uv = this_photon_summary.calculated_distance_uv;
            pe.calculated_distance_visible = this_photon_summary.calculated_distance_visible;
            pe.photon_source = PhotonSummarytoCCMMCPEPhotonSource.at(this_photon_summary.photon_source);
        }

        std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.g4_time < b.g4_time; });
    }
}
