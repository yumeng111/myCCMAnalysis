#include <icetray/serialization.h>
#include "recclasses/I3PortiaEvent.h"

I3PortiaEvent::~I3PortiaEvent() {}

template <class Archive>
void I3PortiaEvent::serialize(Archive& ar, unsigned version)
{

  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));

   ar & make_nvp("BestNPE",           BestNPE_);
   ar & make_nvp("AtwdNPE",           AtwdNPE_);
   ar & make_nvp("FadcNPE",           FadcNPE_);

   ar & make_nvp("TotalNch",          TotalNch_);
   ar & make_nvp("AtwdNch",           AtwdNch_);
   ar & make_nvp("FadcNch",           FadcNch_);

  ar & make_nvp("FirstOMKey",         FirstOMKey_);
  ar & make_nvp("LastOMKey",          LastOMKey_);

  if(version >=2){ // variables used in the IC40 analysis
    ar & make_nvp("BestNPEbtw",           BestNPEbtw_);
    ar & make_nvp("AtwdNPEbtw",           AtwdNPEbtw_);
    ar & make_nvp("FadcNPEbtw",           FadcNPEbtw_);

    ar & make_nvp("TotalNchbtw",          TotalNchbtw_);
    ar & make_nvp("AtwdNchbtw",           AtwdNchbtw_);
    ar & make_nvp("FadcNchbtw",           FadcNchbtw_);

    ar & make_nvp("LargestOMKey",          largestOMKey_);
    ar & make_nvp("FirstOMKeybtw",         FirstOMKeybtw_);
    ar & make_nvp("LastOMKeybtw",          LastOMKeybtw_);
  }

  //for future use
  /*
  if(version > 0)
  ar & make_nvp("LCBit", LCBit_);
  else LCBit_ = false;

  if(version > 1)
  ar & make_nvp("pmtGain", PMTGain_);
  else PMTGain_ = 1;
  */

}

std::ostream& I3PortiaEvent::Print(std::ostream& os) const{
  os << "I3PortiaEvent:\n"
     << "      BestNPE: " << BestNPE_ << '\n'
     << "      ATWDNPE: " << AtwdNPE_ << '\n'
     << "      FADCNPE: " << FadcNPE_ << '\n'
     << "     TotalNch: " << TotalNch_ << '\n'
     << "      ATWDNch: " << AtwdNch_ << '\n'
     << "      FADCNch: " << FadcNch_ << '\n'
     << "      FirstOM: " << FirstOMKey_ << '\n'
     << "       LastOM: " << LastOMKey_ << '\n'
     << "   BestNPEbtw: " << BestNPEbtw_ << '\n'
     << "   ATWDNPEbtw: " << AtwdNPEbtw_ << '\n'
     << "   FADCNPEbtw: " << FadcNPEbtw_ << '\n'
     << "  TotalNchbtw: " << TotalNchbtw_ << '\n'
     << "   ATWDNchbtw: " << AtwdNchbtw_ << '\n'
     << "   FADCNchbtw: " << FadcNchbtw_ << '\n'
     << "    LargestOM: " << largestOMKey_ << '\n'
     << "   FirstOMbtw: " << FirstOMKeybtw_ << '\n'
     << "    LastOMbtw: " << LastOMKeybtw_ << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const I3PortiaEvent& e){
  return(e.Print(os));
}

I3_SERIALIZABLE(I3PortiaEvent);
