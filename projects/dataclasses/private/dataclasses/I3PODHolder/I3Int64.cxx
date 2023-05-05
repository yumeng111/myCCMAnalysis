#include <dataclasses/I3Int64.h>
#include <icetray/I3PODHolder.h>
#include <icetray/serialization.h>

I3_SERIALIZABLE(I3Int64);
I3_EXPORT_EXTRA_KEY(I3Int64, I3PODHolder<long>);
