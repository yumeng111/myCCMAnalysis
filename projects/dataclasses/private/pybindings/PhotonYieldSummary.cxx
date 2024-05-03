#include <vector>
#include <iostream>

#include <dataclasses/physics/PhotonYieldSummary.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>


namespace bp = boost::python;

std::string to_str(PhotonYieldSummary const & theSummary) {
    std::ostringstream oss;
    oss << theSummary << std::flush;
    return oss.str();
}


void register_PhotonYieldSummary() {
    std::string(*convert_to_str)(PhotonYieldSummary const &) = to_str;
    bp::class_<PhotonYieldSummary, bp::bases<I3FrameObject>, boost::shared_ptr<PhotonYieldSummary> >("PhotonYieldSummary")

    .def_readwrite("time",&PhotonYieldSummary::time)
    .def_readwrite("yield",&PhotonYieldSummary::yield)
    .def_readwrite("photon_source",&PhotonYieldSummary::photon_source)
    .def(bp::dataclass_suite<PhotonYieldSummary>())
    ;
    boost::python::enum_<PhotonYieldSummary::PhotonSource>("PhotonSource")
      .value("Unknown", PhotonYieldSummary::PhotonSource::Unknown)
      .value("Vertex", PhotonYieldSummary::PhotonSource::Vertex)
      .value("TPBFoil", PhotonYieldSummary::PhotonSource::TPBFoil)
      .export_values()
    ;
    register_pointer_conversions<PhotonYieldSummary>();

}
