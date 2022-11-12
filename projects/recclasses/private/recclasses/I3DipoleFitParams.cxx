#include "recclasses/I3DipoleFitParams.h"
#include "recclasses/Utility.h"
#include <icetray/serialization.h>

I3DipoleFitParams::~I3DipoleFitParams() {}

template <class Archive>
void I3DipoleFitParams::serialize(Archive& archive, unsigned version)
{
  archive & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  archive & make_nvp("magnet", magnet);
  archive & make_nvp("magnetX", magnetX);
  archive & make_nvp("magnetY", magnetY);
  archive & make_nvp("magnetZ", magnetZ);
  archive & make_nvp("ampSum", ampSum);
  archive & make_nvp("nHits", nHits);
  archive & make_nvp("nPairs", nPairs);
  archive & make_nvp("maxAmp", maxAmp);
}

I3_SERIALIZABLE(I3DipoleFitParams);

bool I3DipoleFitParams::operator==(const I3DipoleFitParams& other) const
{
  return nan_aware_equality(magnet, other.magnet) &&
         nan_aware_equality(magnetX, other.magnetX) &&
         nan_aware_equality(magnetY, other.magnetY) &&
         nan_aware_equality(magnetZ, other.magnetZ) &&
         nan_aware_equality(ampSum, other.ampSum) &&
         nHits == other.nHits &&
         nPairs == other.nPairs &&
         nan_aware_equality(maxAmp, other.maxAmp);
}

std::ostream& operator<<(std::ostream& os,
                         const I3DipoleFitParams& x){
  return(x.Print(os));
}

std::ostream& I3DipoleFitParams::Print(std::ostream& os) const
{
  os << "[ I3DipoleFitParams::"
        "\n  magnet : " << magnet <<
        "\n  magnetX: " << magnetX <<
        "\n  magnetY: " << magnetY <<
        "\n  magnetZ: " << magnetZ <<
        "\n  ampSum : " << ampSum <<
        "\n  nHits  : " << nHits <<
        "\n  nPairs : " << nPairs <<
        "\n  maxAmp : " << maxAmp <<
        " ]";
  return os;
}
