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
    enum class VolumeType {
        Detector,
        Veto,
        Inner
    };

    G4CCMReadout() = default;
    G4CCMReadout(size_t n_threads);
    ~G4CCMReadout() = default;

    static void UpdateMCPESeries(CCMMCPESeriesMapPtr mcpeseries, boost::shared_ptr<I3Map<int, PhotonSummary>> photon_summary_map, bool full_photon_tracking);

    void SetInput(std::vector<I3Particle> primaries);
    void SetOutputs(std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> edep_trees, std::vector<I3MCTreePtr> veto_edep_trees, std::vector<I3MCTreePtr> inner_edep_trees, std::vector<I3VectorI3ParticlePtr> veto_edep_vector, std::vector<I3VectorI3ParticlePtr> inner_edep_vector);
    void SetNumberOfThreads(size_t n_threads);

    I3Particle GetPrimary(size_t i) const;
    CCMMCPESeriesMapPtr GetMCPESeries(size_t i) const;
    I3MCTreePtr GetEDepMCTree(size_t i) const;
    I3MCTreePtr GetVetoEDepMCTree(size_t i) const;
    I3MCTreePtr GetInnerEDepMCTree(size_t i) const;
    I3VectorI3ParticlePtr GetVetoEDepVector(size_t i) const;
    I3VectorI3ParticlePtr GetInnerEDepVector(size_t i) const;

    I3MCTreePtr GetVolumeEDepMCTree(size_t i, VolumeType volume) const;
    I3VectorI3ParticlePtr GetVolumeEDepVector(size_t i, VolumeType volume) const;

    void Reset();

    void LogTrackingResult(int event_id, boost::shared_ptr<I3Map<int, PhotonSummary>> photon_summary_map, bool full_photon_tracking);
    void LogPMTResult(int event_id, CCMMCPESeriesMapPtr mcpeseries, bool full_photon_tracking);

private:
    size_t n_threads_ = 0;

    std::vector<I3Particle> primaries;

    std::vector<CCMMCPESeriesMapPtr> mcpeseries;

    std::vector<I3MCTreePtr> edep_trees;
    std::vector<I3MCTreePtr> veto_edep_trees;
    std::vector<I3MCTreePtr> inner_edep_trees;
    std::vector<I3VectorI3ParticlePtr> veto_edep_vector;
    std::vector<I3VectorI3ParticlePtr> inner_edep_vector;

    std::vector<boost::shared_ptr<I3Map<int, PhotonSummary>>> photon_summary_map;

    std::vector<bool> tracking_logged;
    std::vector<bool> pmts_logged;

    static const std::unordered_map<PhotonSummary::PhotonSource, CCMMCPE::PhotonSource> PhotonSummarytoCCMMCPEPhotonSource;
};

#endif // G4CCMReadout_h
