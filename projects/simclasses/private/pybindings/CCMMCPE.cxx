//
//   Copyright (c) 2013   Alex Olivas
//

#include <simclasses/CCMMCPE.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_CCMMCPE() {
{
    scope mcpe_scope =
        class_<CCMMCPE, boost::shared_ptr<CCMMCPE> >("CCMMCPE")
        .def(dataclass_suite<CCMMCPE>())
        .def(init<>())
        .def(init<uint32_t>())
        .def(init<uint32_t,double>())
        .def_readwrite("time",&CCMMCPE::time)
        .def_readwrite("npe",&CCMMCPE::npe)
        .def_readonly("ID",&CCMMCPE::ID)
        ;
}

    class_<CCMMCPESeries, CCMMCPESeriesPtr>("CCMMCPESeries")
        .def(dataclass_suite<CCMMCPESeries>())
        ;

    class_<CCMMCPESeriesMap, CCMMCPESeriesMapPtr, bases<I3FrameObject> >("CCMMCPESeriesMap")
        .def(dataclass_suite<CCMMCPESeriesMap>())
        ;
    register_pointer_conversions<CCMMCPESeriesMap>();
}
