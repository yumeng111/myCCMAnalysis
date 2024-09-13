#ifndef G4CCMReadout_h
#define G4CCMReadout_h

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3Map.h"
#include "simclasses/PhotonSummary.h"
#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>
#include <sstream>

#include <G4SystemOfUnits.hh>

class G4CCMReadout {
public:
    struct SingleEventReadout {
        I3MCTreePtr edep_tree = nullptr;
        CCMMCPESeriesMapPtr mcpeseries = nullptr;
        PhotonSummarySeriesPtr photon_summary_series = nullptr;
        boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map = nullptr;
    };

    using SingleThreadReadout = std::map<int, SingleEventReadout>;
    using Readout = std::vector<SingleThreadReadout>;

    G4CCMReadout() = default;
    G4CCMReadout(size_t n_threads);
    ~G4CCMReadout() = default;

    static void UpdateMCPESeries(CCMMCPESeriesMapPtr mcpeseries, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map);

    void SetInput(std::vector<I3Particle> primaries, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> edep_trees);
    void SetNumberOfThreads(size_t n_threads);

    SingleThreadReadout & GetReadout(size_t thread_id);

    void AddEntry(size_t thread_id, int event_id, I3MCTreePtr edep_tree, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map);
    void AddEntry(size_t thread_id, int event_id, CCMMCPESeriesMapPtr mcpeseries);
    void AddEntry(size_t thread_id, int event_id, I3MCTreePtr edep_tree, CCMMCPESeriesMapPtr mcpeseries, PhotonSummarySeriesPtr photon_summary_series, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map);
    I3Particle GetPrimary(size_t i) const;
    void SetPrimary(size_t i, I3Particle primary);
    I3MCTreePtr GetMCTree(size_t i) const;
    CCMMCPESeriesMapPtr GetMCPESeries(size_t i) const;
    void SetMCTree(size_t i, I3MCTreePtr edep_tree);
    void SetMCPESeries(size_t i, CCMMCPESeriesMapPtr mcpeseries);
    void Reset();

private:
    size_t n_threads_ = 0;

    Readout readout;
    std::vector<I3Particle> primaries;
    std::vector<CCMMCPESeriesMapPtr> mcpeseries;
    std::vector<I3MCTreePtr> edep_trees;

    static const std::unordered_map<PhotonSummary::PhotonSource, CCMMCPE::PhotonSource> PhotonSummarytoCCMMCPEPhotonSource;
};

#endif // G4CCMReadout_h
