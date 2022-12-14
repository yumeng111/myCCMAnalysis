
#include <icetray/serialization.h>
#include <dataclasses/I3Vector.h>
#include <string>

I3_SERIALIZABLE(I3VectorBool);
I3_SERIALIZABLE(I3VectorChar);

I3_SERIALIZABLE(I3VectorShort);
I3_SERIALIZABLE(I3VectorUShort);
I3_SERIALIZABLE(I3VectorInt);
I3_SERIALIZABLE(I3VectorUnsignedInt); // this is anomalous... but *dont* change it
                                      // until we work in multiple-class-registration
                                      // names.
I3_SERIALIZABLE(I3VectorInt64);
I3_SERIALIZABLE(I3VectorUInt64);

I3_SERIALIZABLE(I3VectorFloat);
I3_SERIALIZABLE(I3VectorDouble);

I3_SERIALIZABLE(I3VectorString);
I3_SERIALIZABLE(I3VectorOMKey);
I3_SERIALIZABLE(I3VectorTankKey);
I3_SERIALIZABLE(I3VectorModuleKey);
// This is just a typedef for an int:
//I3_SERIALIZABLE(I3VectorStationKey);

I3_SERIALIZABLE(I3VectorDoubleDouble);
I3_SERIALIZABLE(I3VectorI3Position);

I3_SERIALIZABLE(I3VectorI3VectorUShort);
I3_EXPORT_EXTRA_KEY(I3VectorI3VectorUShort, CCMAnalysis::Binary::I3VectorI3VectorUInt16)
I3_SERIALIZABLE(I3VectorI3VectorUInt);
I3_SERIALIZABLE(I3VectorI3VectorUInt64);

I3_SERIALIZABLE(I3VectorI3VectorShort);
I3_SERIALIZABLE(I3VectorI3VectorInt);
I3_SERIALIZABLE(I3VectorI3VectorInt64);
I3_EXPORT_EXTRA_KEY(I3VectorI3VectorInt64, I3Vector<I3Vector<int64_t>>);

