
#include <dataclasses/I3Int16.h>
#include <icetray/I3PODHolder.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/operator_suite.hpp>

namespace bp = boost::python;

static std::string 
I3Int16_prettyprint(const I3Int16& d)
{
  std::ostringstream oss;
  oss << "I3Int16(" << d.value << ")";
  return oss.str();
}

bool I3Int16_bool(const I3Int16& v)
{
    return v.value != int16_t(0);
}

void register_I3Int16() {
    bp::class_<I3Int16, bp::bases<I3FrameObject>, boost::shared_ptr<I3Int16> >("I3Int16",
      "A serializable int16_t. Can compare directly with numeric types.\n\
  Note that python assignment is by reference, creating two links to one object.")
      .def(bp::init<>())
      .def(bp::init<int16_t>())
      .def(bp::init<const I3Int16&>())
      .def_readwrite("value", &I3Int16::value)
      .def("__repr__",I3Int16_prettyprint)
      .def_pickle(bp::boost_serializable_pickle_suite<I3Int16>())
      .def(bp::operator_suite<I3Int16>())
      .def(bp::operator_int_suite<I3Int16>())
      .def(bp::operator_float_suite<I3Int16>())
      .def("__nonzero__", I3Int16_bool)
      .def("__bool__", I3Int16_bool)      
      .def( freeze() )
      ;

    register_pointer_conversions<I3Int16>();
}
