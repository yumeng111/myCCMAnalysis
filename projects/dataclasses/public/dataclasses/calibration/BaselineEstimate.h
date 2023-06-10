#ifndef BaselineEstimate_H
#define BaselineEstimate_H

#include <cstddef>
#include <sstream>

#include <boost/shared_ptr.hpp>
#include "icetray/serialization.h"
#include "icetray/CCMPMTKey.h"
#include "icetray/I3FrameObject.h"
#include "dataclasses/I3Map.h"


static const unsigned baselineestimate_version_ = 0;

struct BaselineEstimate : public I3FrameObject {
    double baseline;
    double stddev;
    size_t target_num_frames;
    size_t num_frames;
    size_t num_samples;
    BaselineEstimate() {}
    BaselineEstimate(double baseline, double stddev, size_t target_num_frames, size_t num_frames, size_t num_samples) :
        baseline(baseline), stddev(stddev), target_num_frames(target_num_frames), num_frames(num_frames), num_samples(num_samples) {}

    template <class Archive>
    void serialize(Archive& ar, unsigned version) {
        ar & make_nvp("baseline", baseline);
        ar & make_nvp("stddev", stddev);
        ar & make_nvp("target_num_frames", target_num_frames);
        ar & make_nvp("num_frames", num_frames);
        ar & make_nvp("num_samples", num_samples);
    }

};

typedef I3Map<CCMPMTKey, BaselineEstimate> BaselineEstimateMap;
I3_POINTER_TYPEDEFS(BaselineEstimateMap);

I3_CLASS_VERSION(BaselineEstimate, baselineestimate_version_);
I3_POINTER_TYPEDEFS(BaselineEstimate);

std::ostream& operator<<(std::ostream& oss, const BaselineEstimate& c);
std::ostream& operator<<(std::ostream& oss, const BaselineEstimateMap& m);

#endif // BaselineEstimate_H
