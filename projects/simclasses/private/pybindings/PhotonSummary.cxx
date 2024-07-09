#include <vector>
#include <iostream>

#include <simclasses/PhotonSummary.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

void register_PhotonSummary() {
    {
    scope mcpe_scope =
        class_<PhotonSummary, boost::shared_ptr<PhotonSummary> >("PhotonSummary")
        .def(dataclass_suite<PhotonSummary>())
        .def(init<size_t, size_t, float, I3Position, PhotonSummary::CreationProcess >())
        .def_readwrite("parent_id",&PhotonSummary::parent_id)
        .def_readwrite("track_id",&PhotonSummary::track_id)
        .def_readwrite("time",&PhotonSummary::time)
        .def_readwrite("position",&PhotonSummary::position)
        .def_readwrite("creation_process",&PhotonSummary::creation_process)
        ;
    enum_<PhotonSummary::CreationProcess>("CreationProcess")
      .value("Unknown", PhotonSummary::CreationProcess::Unknown)
      .value("Scintillation", PhotonSummary::CreationProcess::Scintillation)
      .value("Cerenkov", PhotonSummary::CreationProcess::Cerenkov)
      .value("OpWLS", PhotonSummary::CreationProcess::OpWLS)
      .export_values()
      ;
        
    }


    class_<PhotonSummarySeries, PhotonSummarySeriesPtr, bases<I3FrameObject>>("PhotonSummarySeries")
        .def(dataclass_suite<PhotonSummarySeries>())
        ;
    register_pointer_conversions<PhotonSummarySeries>();

}

//namespace bp = boost::python;
//
//std::string to_str(PhotonSummary const & theSummary) {
//    std::ostringstream oss;
//    oss << theSummary << std::flush;
//    return oss.str();
//}
//
//void register_PhotonSummary() {
//    std::string(*convert_to_str)(PhotonSummary const &) = to_str;
//    bp::class_<PhotonSummary, bp::bases<I3FrameObject>, boost::shared_ptr<PhotonSummary> >("PhotonSummary")
//
//    .def_readwrite("parent_id",&PhotonSummary::parent_id)
//    .def_readwrite("track_id",&PhotonSummary::track_id)
//    .def_readwrite("time",&PhotonSummary::time)
//    .def_readwrite("position",&PhotonSummary::position)
//    .def_readwrite("creation_process",&PhotonSummary::creation_process)
//    .def(bp::dataclass_suite<PhotonSummary>())
//    ;
//    boost::python::enum_<PhotonSummary::CreationProcess>("CreationProcess")
//      .value("Unknown", PhotonSummary::CreationProcess::Unknown)
//      .value("Scintillation", PhotonSummary::CreationProcess::Scintillation)
//      .value("Cerenkov", PhotonSummary::CreationProcess::Cerenkov)
//      .value("OpWLS", PhotonSummary::CreationProcess::OpWLS)
//      .export_values()
//    ;
//    register_pointer_conversions<PhotonSummary>();
//    
//    class_<PhotonSummarySeries, PhotonSummarySeriesPtr, bases<I3FrameObject>>("PhotonSummarySeries")
//        .def(dataclass_suite<PhotonSummarySeries>())
//        ;
//    register_pointer_conversions<PhotonSummarySeries>();
//
//}
