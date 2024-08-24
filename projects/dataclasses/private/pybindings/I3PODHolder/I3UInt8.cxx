
#include <dataclasses/I3UInt8.h>
#include <icetray/I3PODHolder.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/operator_suite.hpp>

namespace bp = boost::python;

static std::string 
I3UInt8_prettyprint(const I3UInt8& d)
{
  std::ostringstream oss;
  oss << "I3UInt8(" << d.value << ")";
  return oss.str();
}

bool I3UInt8_bool(const I3UInt8& v)
{
    return v.value != uint8_t(0);
}

void register_I3UInt8() {
    bp::class_<I3UInt8, bp::bases<I3FrameObject>, boost::shared_ptr<I3UInt8> >("I3UInt8",
      "A serializable uint8_t. Can compare directly with numeric types.\n\
  Note that python assignment is by reference, creating two links to one object.")
      .def(bp::init<>())
      .def(bp::init<uint8_t>())
      .def(bp::init<const I3UInt8&>())
      .def_readwrite("value", &I3UInt8::value)
      .def("__repr__",I3UInt8_prettyprint)
      .def_pickle(bp::boost_serializable_pickle_suite<I3UInt8>())
      .def(bp::operator_suite<I3UInt8>())
      .def(bp::operator_int_suite<I3UInt8>())
      .def(bp::operator_float_suite<I3UInt8>())
      .def("__nonzero__", I3UInt8_bool)
      .def("__bool__", I3UInt8_bool)      
      .def( freeze() )
      ;

    register_pointer_conversions<I3UInt8>();
}
