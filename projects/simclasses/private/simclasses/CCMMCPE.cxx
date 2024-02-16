#include <simclasses/CCMMCPE.h>
#include <ostream>

I3_SERIALIZABLE(CCMMCPESeriesMap);

std::ostream& operator<<(std::ostream& os, const CCMMCPE& pe) {
  return(pe.Print(os));
}

std::ostream& CCMMCPE::Print(std::ostream& os) const{
  os << "[ CCMMCPE::"
     << "\n  Time :" << time
     << "\n  NPE  :" << npe
     << "\n  " << ID
     << " ]";
  return os;
}

