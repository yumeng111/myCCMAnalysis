/**
 * @brief implementation of the I3FiniteCuts class
 *
 * @file I3FiniteCuts.cxx
 * @version $Revision$
 * @date $Date$
 * @author Sebastian Euler <sebastian.euler@icecube.wisc.edu>
 *
 * The serialization method for the data class
 */

#include "recclasses/I3FiniteCuts.h"
#include "recclasses/Utility.h"
#include <icetray/serialization.h>

I3FiniteCuts::~I3FiniteCuts() {}

template <class Archive>
void I3FiniteCuts::serialize(Archive& ar, unsigned version){
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("Length_",Length);
  ar & make_nvp("Lend_",Lend);
  ar & make_nvp("Lstart_",Lstart);
  ar & make_nvp("Sdet_",Sdet);
  ar & make_nvp("finiteCut_",finiteCut);
  ar & make_nvp("DetectorLength_",DetectorLength);
  
  if (version > 1){
    log_warn( "too high version number" );
  }
}

void I3FiniteCuts::Reset(){
  Sdet = NAN;
  finiteCut = NAN;
  Length = NAN;
  Lend = NAN;
  Lstart = NAN;
  DetectorLength = NAN;
}

std::ostream& operator<<(std::ostream& os, const I3FiniteCuts& fc)
{
  return(fc.Print(os));
}

std::ostream& I3FiniteCuts::Print(std::ostream& os) const
{
  os << "[ I3FiniteCuts Length : " << Length << std::endl
     << "          endFraction : " << Lend << std::endl
     << "        startFraction : " << Lstart  << std::endl
     << "                 Sdet : " << Sdet << std::endl
     << "            finiteCut : " << finiteCut  << std::endl
     << "       DetectorLength : " << DetectorLength << std::endl
     << "]" ;
  return os;
}

I3_CLASS_VERSION(I3FiniteCuts, 1);
I3_SERIALIZABLE(I3FiniteCuts);

bool I3FiniteCuts::operator==(const I3FiniteCuts& other) const
{
  return
    nan_aware_equality(Length, other.Length) &&
    nan_aware_equality(Lend, other.Lend) &&
    nan_aware_equality(Lstart, other.Lstart) &&
    nan_aware_equality(Sdet, other.Sdet) &&
    nan_aware_equality(finiteCut, other.finiteCut) &&
    nan_aware_equality(DetectorLength, other.DetectorLength);
}
