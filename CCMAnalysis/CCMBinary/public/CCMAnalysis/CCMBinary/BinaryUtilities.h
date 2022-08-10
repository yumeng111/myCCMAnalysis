#ifndef CCMAnalysis_BinaryUtilities_H
#define CCMAnalysis_BinaryUtilities_H

#include <string>
#include <vector>
#include <iostream>

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

inline std::ostream & write_binary(std::ostream & os, uint32_t const & x);
inline std::istream & read_binary(std::istream & is, uint32_t & s);
inline std::ostream & write_binary(std::ostream & os, float const & x);
inline std::istream & read_binary(std::istream & is, float & s);
inline std::ostream & write_binary(std::ostream & os, timespec const & x);
inline std::istream & read_binary(std::istream & is, timespec & s);
inline std::ostream & write_binary(std::ostream & os, std::string const & s);
inline std::istream & read_binary(std::istream & is, std::string & s);
inline std::ostream & write_binary(std::ostream & os, ChannelHeader const & header);
inline std::istream & read_binary(std::istream & is, ChannelHeader & header);
inline std::ostream & write_binary(std::ostream & os, DigitizerBoard const & board);
inline std::istream & read_binary(std::istream & is, DigitizerBoard & board);
inline std::ostream & write_binary(std::ostream & os, CCMDAQConfig const & config);
inline std::istream & read_binary(std::istream & is, CCMDAQConfig & config);
inline std::ostream & write_binary(std::ostream & os, CCMDigitizerReadout const & digi_readout);
inline std::istream & read_binary(std::istream & is, CCMDigitizerReadout const & digi_readout);
inline std::ostream & write_binary(std::ostream & os, CCMTriggerReadout const & trigger_readout);
inline std::istream & read_binary(std::istream & is, CCMTriggerReadout & trigger_readout);

template<typename T>
inline std::ostream & write_binary(std::ostream & os, std::vector<T> const & v);

template<typename T>
inline std::istream & read_binary(std::istream & is, std::vector<T> & v);

template<typename T>
inline std::ostream & write_binary_contiguous_vector(std::ostream & os, std::vector<T> const & v);

template<typename T>
inline std::istream & read_binary_contiguous_vector(std::istream & is, std::vector<T> & v);

template<>
inline std::ostream & write_binary(std::ostream & os, std::vector<uint32_t> const & v);
template<>
inline std::istream & read_binary(std::istream & is, std::vector<uint32_t> & v);

template<>
inline std::ostream & write_binary(std::ostream & os, std::vector<uint16_t> const & v);
template<>
inline std::istream & read_binary(std::istream & is, std::vector<uint16_t> & v);

template<>
inline std::ostream & write_binary(std::ostream & os, std::vector<float> const & v);
template<>
inline std::istream & read_binary(std::istream & is, std::vector<float> & v);

#include "CCMAnalysis/CCMBinary/BinaryUtilities.cxx"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.tcc"

#endif // CCMAnalysis_BinaryUtilities_H
