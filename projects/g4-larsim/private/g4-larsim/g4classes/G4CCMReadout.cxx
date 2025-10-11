#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/I3Map.h"
#include "simclasses/PhotonSummary.h"
#include "simclasses/CCMMCPE.h"

#include "g4-larsim/g4classes/G4CCMReadout.h"

#include <vector>
#include <sstream>

#include <G4SystemOfUnits.hh>

const std::unordered_map<PhotonSummary::PhotonSource, CCMMCPE::PhotonSource> G4CCMReadout::PhotonSummarytoCCMMCPEPhotonSource = {
    {PhotonSummary::PhotonSource::Unknown, CCMMCPE::PhotonSource::Unknown},
    {PhotonSummary::PhotonSource::Scintillation, CCMMCPE::PhotonSource::Scintillation},
    {PhotonSummary::PhotonSource::Cherenkov, CCMMCPE::PhotonSource::Cherenkov}
};


G4CCMReadout::G4CCMReadout(size_t n_threads) : n_threads_(n_threads) {}

void G4CCMReadout::SetInput(std::vector<I3Particle> primaries) {
    this->primaries = primaries;
}

void G4CCMReadout::SetOutputs(std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> edep_trees, std::vector<I3MCTreePtr> veto_edep_trees, std::vector<I3MCTreePtr> inner_edep_trees, std::vector<I3VectorI3ParticlePtr> veto_edep_vector, std::vector<I3VectorI3ParticlePtr> inner_edep_vector) {
    std::vector<size_t> sizes = {mcpeseries.size(), edep_trees.size(), veto_edep_trees.size(), inner_edep_trees.size(), veto_edep_vector.size(), inner_edep_vector.size()};

    size_t n_events = 0;
    for(size_t i = 0; i < sizes.size(); ++i) {
        if(n_events == 0 and sizes[i] > 0) {
            n_events = sizes[i];
            continue;
        }
        if(sizes[i] != 0 and sizes[i] != n_events) {
            std::stringstream ss;
            ss << "Mismatch in number of events between outputs: " << n_events << " vs " << sizes[i];
            throw std::runtime_error(ss.str());
        }
    }

    if(n_events == 0) {
        throw std::runtime_error("No events found in outputs");
    }

    for(size_t i = 0; i < n_events; ++i) {
        if(edep_trees[i].get() == nullptr) {
            std::stringstream ss;
            ss << "edep_trees[" << i << "] is null";
            throw std::runtime_error(ss.str());
        }
    }

    this->mcpeseries = mcpeseries;
    this->edep_trees = edep_trees;
    this->veto_edep_trees = veto_edep_trees;
    this->inner_edep_trees = inner_edep_trees;
    this->veto_edep_vector = veto_edep_vector;
    this->inner_edep_vector = inner_edep_vector;

    this->photon_summary_map = std::vector<boost::shared_ptr<I3Map<int, std::tuple<ParentInfo, PhotonSummary>>>>(n_events, nullptr);
    this->photon_source_map = std::vector<boost::shared_ptr<I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>>>(n_events, nullptr);

    this->tracking_logged = std::vector<bool>(n_events, false);
    this->detailed_photon_tracking = std::vector<bool>(n_events, false);
    this->pmts_logged = std::vector<bool>(n_events, false);
}

void G4CCMReadout::SetNumberOfThreads(size_t n_threads) {
    n_threads_ = n_threads;
}

void G4CCMReadout::LogPMTResult(int event_id, CCMMCPESeriesMapPtr mcpeseries) {
    this->mcpeseries[event_id] = mcpeseries;
    this->pmts_logged[event_id] = true;

    if(this->tracking_logged[event_id]) {
        if(this->detailed_photon_tracking[event_id]) {
            boost::shared_ptr<I3Map<int, std::tuple<ParentInfo, PhotonSummary>>> photon_summary_map = this->photon_summary_map[event_id];
            UpdateMCPESeries(mcpeseries, photon_summary_map);
        } else {
            boost::shared_ptr<I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>> photon_source_map = this->photon_source_map[event_id];
            UpdateMCPESeries(mcpeseries, photon_source_map);
        }
    }
}

I3Particle G4CCMReadout::GetPrimary(size_t i) const {
    return primaries.at(i);
}

CCMMCPESeriesMapPtr G4CCMReadout::GetMCPESeries(size_t i) const { return mcpeseries.at(i); }
I3MCTreePtr G4CCMReadout::GetEDepMCTree(size_t i) const { return edep_trees.at(i); }
I3MCTreePtr G4CCMReadout::GetVetoEDepMCTree(size_t i) const { return veto_edep_trees.at(i); }
I3MCTreePtr G4CCMReadout::GetInnerEDepMCTree(size_t i) const { return inner_edep_trees.at(i); }
I3VectorI3ParticlePtr G4CCMReadout::GetVetoEDepVector(size_t i) const { return veto_edep_vector.at(i); }
I3VectorI3ParticlePtr G4CCMReadout::GetInnerEDepVector(size_t i) const { return inner_edep_vector.at(i); }
I3MCTreePtr G4CCMReadout::GetVolumeEDepMCTree(size_t i, VolumeType volume) const {
    switch(volume) {
        case VolumeType::Detector:
            return edep_trees.at(i);
        case VolumeType::Veto:
            return veto_edep_trees.at(i);
        case VolumeType::Inner:
            return inner_edep_trees.at(i);
    }
    throw std::runtime_error("Invalid volume type");
}
I3VectorI3ParticlePtr G4CCMReadout::GetVolumeEDepVector(size_t i, VolumeType volume) const {
    switch(volume) {
        case VolumeType::Veto:
            return veto_edep_vector.at(i);
        case VolumeType::Inner:
            return inner_edep_vector.at(i);
        case VolumeType::Detector:
            throw std::runtime_error("No edep vector for full detector volume");
    }
    throw std::runtime_error("Invalid volume type");
}

void G4CCMReadout::Reset() {
    primaries.clear();
    edep_trees.clear();
    veto_edep_trees.clear();
    inner_edep_trees.clear();
    veto_edep_vector.clear();
    inner_edep_vector.clear();
    mcpeseries.clear();
    photon_summary_map.clear();
    tracking_logged.clear();
    detailed_photon_tracking.clear();
    pmts_logged.clear();
}

void G4CCMReadout::UpdateMCPESeries(CCMMCPESeriesMapPtr mcpeseries, boost::shared_ptr<I3Map<int, std::tuple<G4CCMReadout::ParentInfo, PhotonSummary::PhotonSource>>> photon_source_map) {
    // Iterate over PMTs in source map
    for(CCMMCPESeriesMap::iterator it = mcpeseries->begin(); it != mcpeseries->end(); ++it) {
        // Reference to the PE series
        CCMMCPESeries & pe_series = it->second;

        // Iterate backwards over the vector of CCMMCPE in the source map for this PMT
        for(int i=pe_series.size()-1; i>=0; --i) {
            CCMMCPE & pe = pe_series[i];

            // Check if the photon has everything properly recorded
            typename I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>::const_iterator it_map = photon_source_map->find(pe.track_id);
            // If not trash it as an edge case
            if(it_map == photon_source_map->end()) {
                pe_series.erase(pe_series.begin() + i);
                continue;
            }

            // Grab the summary information for this photon track
            PhotonSummary::PhotonSource const & this_photon_source = std::get<1>(it_map->second);

            // Update the destination CCMMCPE with the summary information
            pe.photon_source = PhotonSummarytoCCMMCPEPhotonSource.at(this_photon_source);
        }
        std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.time < b.time; });
    }
}

void G4CCMReadout::UpdateMCPESeries(CCMMCPESeriesMapPtr mcpeseries, boost::shared_ptr<I3Map<int, std::tuple<G4CCMReadout::ParentInfo, PhotonSummary>>> photon_summary_map) {
    // Iterate over PMTs in source map
    for(CCMMCPESeriesMap::iterator it = mcpeseries->begin(); it != mcpeseries->end(); ++it) {
        // Reference to the PE series
        CCMMCPESeries & pe_series = it->second;

        // Iterate backwards over the vector of CCMMCPE in the source map for this PMT
        for(int i=pe_series.size()-1; i>=0; --i) {
            CCMMCPE & pe = pe_series[i];

            // Check if the photon has everything properly recorded
            typename I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::const_iterator it_map = photon_summary_map->find(pe.track_id);
            // If not trash it as an edge case
            if(it_map == photon_summary_map->end()) {
                pe_series.erase(pe_series.begin() + i);
                continue;
            }

            // Grab the summary information for this photon track
            PhotonSummary const & this_photon_summary = std::get<1>(it_map->second);

            // Update the destination CCMMCPE with the summary information
            pe.distance_uv = this_photon_summary.distance_uv;
            pe.original_wavelength = this_photon_summary.original_wavelength;
            pe.distance_visible = this_photon_summary.distance_visible;
            pe.photon_source = PhotonSummarytoCCMMCPEPhotonSource.at(this_photon_summary.photon_source);
            pe.n_photons_per_wls = this_photon_summary.n_photons_per_wls;
            pe.wls_loc = this_photon_summary.wls_loc;
        }
        std::sort(it->second.begin(), it->second.end(), [](const CCMMCPE& a, const CCMMCPE& b) { return a.time < b.time; });
    }
}

void G4CCMReadout::LogTrackingResult(int event_id, boost::shared_ptr<I3Map<int, std::tuple<G4CCMReadout::ParentInfo, PhotonSummary::PhotonSource>>> photon_source_map) {
    this->photon_source_map[event_id] = photon_source_map;
    this->tracking_logged[event_id] = true;
    this->detailed_photon_tracking[event_id] = false;

    if(this->pmts_logged[event_id]) {
        CCMMCPESeriesMapPtr mcpeseries = this->mcpeseries[event_id];
        UpdateMCPESeries(mcpeseries, photon_source_map);
    }
}

void G4CCMReadout::LogTrackingResult(int event_id, boost::shared_ptr<I3Map<int, std::tuple<G4CCMReadout::ParentInfo, PhotonSummary>>> photon_summary_map) {
    this->photon_summary_map[event_id] = photon_summary_map;
    this->tracking_logged[event_id] = true;
    this->detailed_photon_tracking[event_id] = true;

    if(this->pmts_logged[event_id]) {
        CCMMCPESeriesMapPtr mcpeseries = this->mcpeseries[event_id];
        UpdateMCPESeries(mcpeseries, photon_summary_map);
    }
}

