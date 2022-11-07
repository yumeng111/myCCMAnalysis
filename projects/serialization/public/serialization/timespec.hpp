#ifndef I3_SERIALIZATION_TIMESPEC_HPP
#define I3_SERIALIZATION_TIMESPEC_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// timespec serialization

#include <time.h>
#include <boost/config.hpp>

#include <serialization/nvp.hpp>
#include <serialization/split_free.hpp>

namespace icecube {
namespace serialization {

template<class Archive>
inline void serialize(
    Archive & ar,
    timespec & t,
    const unsigned int file_version
){
    icecube::serialization::split_free(ar, t, file_version);
}

template<class Archive>
inline void save(
    Archive & ar,
    timespec const & t,
    const unsigned int /* file_version */
){
    const int64_t tv_sec = int64_t(t.tv_sec);
    const int64_t tv_nsec = int64_t(t.tv_nsec);
    ar << icecube::serialization::make_nvp("tv_sec" , tv_sec);
    ar << icecube::serialization::make_nvp("tv_nsec", tv_nsec);
}

template<class Archive>
inline void load(
    Archive & ar,
    timespec & t,
    const unsigned int /* file_version */
){
    int64_t tv_sec;
    int64_t tv_nsec;
    ar >> icecube::serialization::make_nvp("tv_sec" , tv_sec);
    ar >> icecube::serialization::make_nvp("tv_nsec", tv_nsec);
    t.tv_sec = time_t(tv_sec);
    t.tv_nsec = long(tv_nsec);
}

} // serialization
} // namespace icecube


#endif // I3_SERIALIZATION_TIMESPEC_HPP
