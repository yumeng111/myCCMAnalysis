#include <icetray/serialization.h>
#include "dataclasses/physics/CCMBCMSummary.h"

std::ostream & CCMBCMSummary::Print(std::ostream& oss) const {
    oss << "[BCMSummary:\n"
        << "      StartTime : " << bcm_start_time << '\n'
        << "        EndTime : " << bcm_end_time << '\n'
        << "       PeakTime : " << bcm_peak_time << '\n'
        << "   PeakADCCount : " << bcm_peak_value << '\n'
        << "       Integral : " << bcm_integral << '\n'
        << "       Baseline : " << bcm_baseline << '\n'
        << " BaselineStdDev : " << bcm_baseline_stddev << "\n]";
    return oss;
}

std::ostream& operator<<(std::ostream& oss, CCMBCMSummary const & bcm) {
    return(bcm.Print(oss));
}

I3_SERIALIZABLE(CCMBCMSummary);
