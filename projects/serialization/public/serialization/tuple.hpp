#ifndef I3_SERIALIZATION_TUPLE_HPP
#define I3_SERIALIZATION_TUPLE_HPP

#if defined(_MSC_VER)
# pragma once
#endif

#include <tuple>

namespace icecube {
namespace serialization {

template<uint N>
struct Serialize
{
    template<class Archive, typename... Args>
    static void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
    {
        ar & std::get<N-1>(t);
        Serialize<N-1>::serialize(ar, t, version);
    }
};

template<>
struct Serialize<0>
{
    template<class Archive, typename... Args>
    static void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
    {
        (void) ar;
        (void) t;
        (void) version;
    }
};

template<class Archive, typename... Args>
void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
{
    Serialize<sizeof...(Args)>::serialize(ar, t, version);
}

}
}

#endif // I3_SERIALIZATION_TUPLE_HPP
