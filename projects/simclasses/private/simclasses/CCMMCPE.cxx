#include <simclasses/CCMMCPE.h>
#include <ostream>

I3_SERIALIZABLE(CCMMCPESeriesMap);

std::ostream& operator<<(std::ostream& os, const CCMMCPE& pe) {
    return(pe.Print(os));
}

std::ostream& CCMMCPE::Print(std::ostream& os) const{
    os << "[ CCMMCPE::"
        << "\n  Time :" << time
        << "\n  Wavelngth :" << wavelength
        << "\n  Position :" << position 
        << "\n  Direction :" << direction 
        //std::ostream & operator<<(std::ostream & os, CCMMCPE::PhotonSource const & p) {
        //    if (CCMMCPE::PhotonSourceNames.find(p) != CCMMCPE::PhotonSourceNames.end())
        //        << "\n  PhotonSource :" << CCMMCPE::PhotonSourceNames.at(p);
        << " ]";
    return os;
}


