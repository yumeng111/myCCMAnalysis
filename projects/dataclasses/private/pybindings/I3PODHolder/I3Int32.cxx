
#include <dataclasses/I3Int32.h>
#include <icetray/I3PODHolder.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/operator_suite.hpp>

namespace bp = boost::python;

static std::string 
I3Int32_prettyprint(const I3Int32& d)
{
  std::ostringstream oss;
  oss << "I3Int32(" << d.value << ")";
  return oss.str();
}

bool I3Int32_bool(const I3Int32& v)
{
    return v.value != int32_t(0);
}

void register_I3Int32() {
    bp::class_<I3Int32, bp::bases<I3FrameObject>, boost::shared_ptr<I3Int32> >("I3Int32",
      "A serializable int32_t. Can compare directly with numeric types.\n\
  Note that python assignment is by reference, creating two links to one object.")
      .def(bp::init<>())
      .def(bp::init<int32_t>())
      .def(bp::init<const I3Int32&>())
      .def_readwrite("value", &I3Int32::value)
      .def("__repr__",I3Int32_prettyprint)
      .def_pickle(bp::boost_serializable_pickle_suite<I3Int32>())
      .def(bp::operator_suite<I3Int32>())
      .def(bp::operator_int_suite<I3Int32>())
      .def(bp::operator_float_suite<I3Int32>())
      .def("__nonzero__", I3Int32_bool)
      .def("__bool__", I3Int32_bool)      
      .def( freeze() )
      ;

    register_pointer_conversions<I3Int32>();
}
