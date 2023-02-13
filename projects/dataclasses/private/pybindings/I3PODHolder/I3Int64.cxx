
#include <icetray/I3PODHolder.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/operator_suite.hpp>

namespace bp = boost::python;
typedef I3PODHolder<int64_t> I3Int64;

static std::string 
I3Int64_prettyprint(const I3Int64& d)
{
  std::ostringstream oss;
  oss << "I3Int64(" << d.value << ")";
  return oss.str();
}

bool I3Int64_bool(const I3Int64& v)
{
    return v.value != int64_t(0);
}

void register_I3Int64() {
    bp::class_<I3Int64, bp::bases<I3FrameObject>, boost::shared_ptr<I3Int64> >("I3Int64",
      "A serializable int64_t. Can compare directly with numeric types.\n\
  Note that python assignment is by reference, creating two links to one object.")
      .def(bp::init<>())
      .def(bp::init<int64_t>())
      .def(bp::init<const I3Int64&>())
      .def_readwrite("value", &I3Int64::value)
      .def("__repr__",I3Int64_prettyprint)
      .def_pickle(bp::boost_serializable_pickle_suite<I3Int64>())
      .def(bp::operator_suite<I3Int64>())
      .def(bp::operator_int_suite<I3Int64>())
      .def(bp::operator_float_suite<I3Int64>())
      .def("__nonzero__", I3Int64_bool)
      .def("__bool__", I3Int64_bool)      
      .def( freeze() )
      ;

    register_pointer_conversions<I3Int64>();
}
