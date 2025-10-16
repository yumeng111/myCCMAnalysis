
#include <dataclasses/I3Map.h>
#include <dataclasses/physics/I3ParticleID.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_I3MapI3ParticleIDVectorDouble()
{
    class_<I3MapI3ParticleIDVectorDouble, bases<I3FrameObject>,
        boost::shared_ptr<I3MapI3ParticleIDVectorDouble> >("I3MapI3ParticleIDVectorDouble")
        .def(dataclass_suite<I3MapI3ParticleIDVectorDouble >())
    ;

    register_pointer_conversions<I3MapI3ParticleIDVectorDouble>();
}

