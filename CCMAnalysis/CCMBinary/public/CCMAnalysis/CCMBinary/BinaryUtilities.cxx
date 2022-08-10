#ifndef CCMAnalysis_BinaryUtilities_CXX
#define CCMAnalysis_BinaryUtilities_CXX

#include <string>
#include <vector>
#include <iostream>

#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

std::ostream & write_binary(std::ostream & os, uint32_t const & x) {
    constexpr size_t size = sizeof(uint32_t);
    os.write((char*)&x, size);
    return os;
}

std::istream & read_binary(std::istream & is, uint32_t & s) {
    constexpr size_t size = sizeof(uint32_t);
    is.read((char*)&s, size);
    return is;
}

std::ostream & write_binary(std::ostream & os, float const & x) {
    constexpr size_t size = sizeof(float);
    os.write((char*)&x, size);
    return os;
}

std::istream & read_binary(std::istream & is, float & s) {
    constexpr size_t size = sizeof(float);
    is.read((char*)&s, size);
    return is;
}

std::ostream & write_binary(std::ostream & os, std::string const & s) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t string_size = s.size();
    os.write((char*)&string_size, size_size);
    os.write(s.c_str(), string_size);
    return os;
}

std::istream & read_binary(std::istream & is, std::string & s) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t string_size;
    is.read((char*)&string_size, size_size);
    char * str_contents = new char[string_size + 1];
    is.read(str_contents, string_size);
    str_contents[string_size] = '\0';
    s = str_contents;
    return is;
}

std::ostream & write_binary(std::ostream & os, ChannelHeader const & header) {
    if(header.version == 0) {
        write_binary(os, header.version);
        write_binary(os, header.physical_board_id);
        write_binary(os, header.board_serial_number);
        write_binary(os, header.physical_channel_id);
        write_binary(os, header.caen_optical_link_number);
        write_binary(os, header.caen_optical_link_board_number);
        write_binary(os, header.caen_channel_number);
    } else {
        throw std::runtime_error("Can only write ChannelHeader version <= 0");
    }
    return os;
}

std::istream & read_binary(std::istream & is, ChannelHeader & header) {
    read_binary(is, header.version);
    if(header.version == 0) {
        read_binary(is, header.physical_board_id);
        read_binary(is, header.board_serial_number);
        read_binary(is, header.physical_channel_id);
        read_binary(is, header.caen_optical_link_number);
        read_binary(is, header.caen_optical_link_board_number);
        read_binary(is, header.caen_channel_number);
    } else {
        throw std::runtime_error("Can only read ChannelHeader version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, DigitizerBoard const & board) {
    if(board.version == 0) {
        write_binary(os, board.version);
        write_binary(os, board.physical_board_id);
        write_binary(os, board.board_serial_number);
        write_binary(os, board.caen_optical_link_number);
        write_binary(os, board.caen_optical_link_board_number);
        write_binary(os, board.channels);
    } else {
        throw std::runtime_error("Can only write DigitizerBoard version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, DigitizerBoard & board) {
    read_binary(is, board.version);
    if(board.version == 0) {
        read_binary(is, board.physical_board_id);
        read_binary(is, board.board_serial_number);
        read_binary(is, board.caen_optical_link_number);
        read_binary(is, board.caen_optical_link_board_number);
        read_binary(is, board.channels);
    } else {
        throw std::runtime_error("Can only read DigitizerBoard version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMDAQConfig const & config) {
    if(config.version == 0) {
        write_binary(os, config.version);
        write_binary(os, config.machine_identifier);
        write_binary(os, config.num_digitizer_boards);
        write_binary(os, config.num_channels);
        write_binary(os, config.num_samples);
        write_binary(os, config.ring_buffer_size);
        write_binary(os, config.trigger_percent_after);
        write_binary(os, config.digitizer_boards);
    } else {
        throw std::runtime_error("Can only write CCMDAQConfig version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMDAQConfig & config) {
    read_binary(is, config.version);
    if(config.version == 0) {
        read_binary(is, config.machine_identifier);
        read_binary(is, config.num_digitizer_boards);
        read_binary(is, config.num_channels);
        read_binary(is, config.num_samples);
        read_binary(is, config.ring_buffer_size);
        read_binary(is, config.trigger_percent_after);
        read_binary(is, config.digitizer_boards);
    } else {
        throw std::runtime_error("Can only read CCMDAQConfig version <= 0");
    }
    return is;
}

#endif // CCMAnalysis_BinaryUtilities_CXX
