#include <dataclasses/physics/PhotonYieldSummary.h>
#include <ostream>

I3_SERIALIZABLE(PhotonYieldSummarySeriesMap);
I3_SERIALIZABLE(PhotonYieldSummarySeries);

const std::unordered_map<PhotonYieldSummary::PhotonSource, std::string> PhotonYieldSummary::photonSourceToProcessName = {{PhotonYieldSummary::PhotonSource::Unknown, "Unknown"},
                                                                                                      {PhotonYieldSummary::PhotonSource::Vertex, "Vertex"},
                                                                                                      {PhotonYieldSummary::PhotonSource::TPBFoil, "TPBFoil"}};
std::ostream& operator<<(std::ostream& os, const PhotonYieldSummary& pe) {
    return(pe.Print(os));
}

std::ostream& PhotonYieldSummary::Print(std::ostream& os) const{
    os << "[ PhotonYieldSummary::"
        << "\n  Time :" << time
        << "\n  Yield :" << yield
        << "\n PhotonSource :" << photonSourceToProcessName.at(photon_source) 
        << " ]";
    return os;
}


