#include "dataclasses/physics/HESodiumEvent.h"
#include <icetray/serialization.h>
#include <ostream>

std::ostream& HESodiumEvent::Print(std::ostream& os) const{
    os << "[ HESodiumEvent::"
        << "\n 1275kev Photon Vertex :" << photon_vertex
        << "\n 511kev Vertex :" << electron_vertex
        << "\n 511kev Vertex :" << positron_vertex
        << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& oss, HESodiumEvent const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, HESodiumEvent & bcm) {
    return(bcm.Print(oss));
}

I3_SERIALIZABLE(HESodiumEventSeries);
I3_SERIALIZABLE(HESodiumEvent);
