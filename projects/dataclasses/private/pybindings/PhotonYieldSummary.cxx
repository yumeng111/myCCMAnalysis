#include <dataclasses/physics/PhotonYieldSummary.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_PhotonYieldSummary() {
    {
    scope photonyield_scope =
        class_<PhotonYieldSummary, boost::shared_ptr<PhotonYieldSummary> >("PhotonYieldSummary")
        .def(dataclass_suite<PhotonYieldSummary>())
        .def(init<float, float, I3Position, I3Direction, PhotonYieldSummary::PhotonSource >())
        .def_readwrite("time",&PhotonYieldSummary::time)
        .def_readwrite("yield",&PhotonYieldSummary::yield)
        .def_readwrite("photon_source",&PhotonYieldSummary::photon_source)
        ;
    enum_<PhotonYieldSummary::PhotonSource>("PhotonSource")
      .value("Unknown", PhotonYieldSummary::PhotonSource::Unknown)
      .value("Vertex", PhotonYieldSummary::PhotonSource::Vertex)
      .value("TPBFoil", PhotonYieldSummary::PhotonSource::TPBFoil)
      .export_values()
      ;
        
    }


    class_<PhotonYieldSummarySeries, PhotonYieldSummarySeriesPtr, bases<I3FrameObject>>("PhotonYieldSummarySeries")
        .def(dataclass_suite<PhotonYieldSummarySeries>())
        ;
    register_pointer_conversions<PhotonYieldSummarySeries>();

    class_<PhotonYieldSummarySeriesMap, PhotonYieldSummarySeriesMapPtr, bases<I3FrameObject> >("PhotonYieldSummarySeriesMap")
        .def(dataclass_suite<PhotonYieldSummarySeriesMap>())
        ;
    register_pointer_conversions<PhotonYieldSummarySeriesMap>();
}
