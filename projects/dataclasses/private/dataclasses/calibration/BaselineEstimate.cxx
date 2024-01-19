#include "dataclasses/calibration/BaselineEstimate.h"

std::ostream& operator<<(std::ostream& oss, BaselineEstimate const & c) {
    oss << "[ BaselineEstimate  :: " << std::endl
        << "        Baseline : " << c.baseline << std::endl
        << "          StdDev : " << c.stddev << std::endl
        << " TargetNumFrames : " << c.target_num_frames << std::endl
        << "       NumFrames : " << c.num_frames << std::endl
        << "      NumSamples : " << c.num_samples << std::endl
        << "]" ;
    return oss;
}

std::ostream& operator<<(std::ostream& oss, BaselineEstimateMap const & m) {
    oss << "[ BaselineEstimateMap :: " << std::endl;
    BaselineEstimateMap::const_iterator iter = m.begin();
    for (;iter != m.end();iter++)
    {
        oss << "  " << iter->first << " : " << iter->second << std::endl;
    }
    oss << "]" ;
    return oss;
}

std::ostream& operator<<(std::ostream& oss, BaselineEstimateChannelMap const & m) {
    oss << "[ BaselineEstimateChannelMap :: " << std::endl;
    BaselineEstimateChannelMap::const_iterator iter = m.begin();
    for (;iter != m.end();iter++)
    {
        oss << "  " << iter->first << " : " << iter->second << std::endl;
    }
    oss << "]" ;
    return oss;
}

I3_SERIALIZABLE(BaselineEstimate);
I3_SERIALIZABLE(BaselineEstimateMap);
I3_SERIALIZABLE(BaselineEstimateChannelMap);
