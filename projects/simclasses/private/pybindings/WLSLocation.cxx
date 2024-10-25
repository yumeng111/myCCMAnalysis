#include <vector>
#include <iostream>

#include <simclasses/WLSLocation.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

void register_WLSLocation() {
    {
    scope mcpe_scope =
        class_<WLSLocation, boost::shared_ptr<WLSLocation> >("WLSLocation")
        .def(dataclass_suite<WLSLocation>())
        .def(init<WLSLocation::WLSLoc>())
        .def_readwrite("wls_loc",&WLSLocation::wls_loc)
        ;
    enum_<WLSLocation::WLSLoc>("WLSLocation")
      .value("Unknown", WLSLocation::WLSLoc::Unknown)
      .value("PMT", WLSLocation::WLSLoc::PMT)
      .value("FoilTop", WLSLocation::WLSLoc::FoilTop)
      .value("FoilBottom", WLSLocation::WLSLoc::FoilBottom)
      .value("FoilSides", WLSLocation::WLSLoc::FoilSides)
      .export_values()
      ;
    }


    class_<WLSLocationSeries, WLSLocationSeriesPtr, bases<I3FrameObject>>("WLSLocationSeries")
        .def(dataclass_suite<WLSLocationSeries>())
        ;
    register_pointer_conversions<WLSLocationSeries>();

}

