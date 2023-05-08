#include <vector>
#include <iostream>

#include <dataclasses/physics/CCMBCMSummary.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

namespace bp = boost::python;

std::string to_str(CCMBCMSummary const & theBCMSummary) {
    std::ostringstream oss;
    oss << theBCMSummary << std::flush;
    return oss.str();
}

void register_CCMBCMSummary() {

    std::string(*convert_to_str)(CCMBCMSummary const &) = to_str;
    bp::class_<CCMBCMSummary, bp::bases<I3FrameObject>, boost::shared_ptr<CCMBCMSummary> >("CCMBCMSummary")

    .def_readwrite("start_time", &CCMBCMSummary::bcm_start_time)
    .def_readwrite("end_time", &CCMBCMSummary::bcm_end_time)
    .def_readwrite("peak_time", &CCMBCMSummary::bcm_peak_time)
    .def_readwrite("peak_value", &CCMBCMSummary::bcm_peak_value)
    .def_readwrite("integral", &CCMBCMSummary::bcm_integral)
    .def_readwrite("baseline", &CCMBCMSummary::bcm_baseline)
    .def_readwrite("baseline_stddev", &CCMBCMSummary::bcm_baseline_stddev)
    .def(bp::dataclass_suite<CCMBCMSummary>())
    ;

    register_pointer_conversions<CCMBCMSummary>();
}
