#include <vector>
#include <iostream>

#include <dataclasses/physics/CCMBCMSummary.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

std::string to_str(CCMBCMSummary const & theBCMSummary) {
    std::ostringstream oss;
    oss << theBCMSummary << std::flush;
    return oss.str();
}

void register_CCMSummary() {
    {
        scope waveform_scope =
            class_<CCMBCMSummary, boost::shared_ptr<CCMBCMSummary> >("CCMBCMSummary")
            .def(copy_suite<CCMBCMSummary>())
            .add_property("start_time", &CCMBCMSummary::bcm_start_time)
            .add_property("end_time", &CCMBCMSummary::bcm_end_time)
            .add_property("peak_time", &CCMBCMSummary::bcm_peak_time)
            .add_property("peak_value", &CCMBCMSummary::bcm_peak_value)
            .add_property("integral", &CCMBCMSummary::bcm_integral)
            .add_property("baseline", &CCMBCMSummary::bcm_baseline)
            .add_property("baseline_stddev", &CCMBCMSummary::bcm_baseline_stddev)

            .def(self == self)
            .def(dataclass_suite<CCMBCMSummary>())
            .def("__str__", to_str)
            ;
    }
}
