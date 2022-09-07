#ifndef CCMAnalysis_BinaryUtilities_CXX
#define CCMAnalysis_BinaryUtilities_CXX

#include <string>
#include <vector>
#include <cstring>
#include <sstream>
#include <iostream>

#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

namespace CCMAnalysis {
namespace Binary {

inline std::ostream & write_binary(std::ostream & os, uint32_t const & x) {
    constexpr size_t size = sizeof(uint32_t);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, uint32_t & x) {
    constexpr size_t size = sizeof(uint32_t);
    is.read((char*)&x, size);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, uint16_t const & x) {
    constexpr size_t size = sizeof(uint16_t);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, uint16_t & x) {
    constexpr size_t size = sizeof(uint16_t);
    is.read((char*)&x, size);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, float const & x) {
    constexpr size_t size = sizeof(float);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, float & x) {
    constexpr size_t size = sizeof(float);
    is.read((char*)&x, size);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, char const * x, size_t n) {
    os.write(x, n);
    return os;
}

inline std::istream & read_binary(std::istream & is, char * x, size_t n) {
    is.read(x, n);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, timespec const & x) {
    constexpr size_t size = sizeof(timespec);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, timespec & s) {
    constexpr size_t size = sizeof(timespec);
    is.read((char*)&s, size);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, std::string const & s) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t string_size = s.size();
    os.write((char*)&string_size, size_size);
    os.write(s.c_str(), string_size);
    return os;
}

inline std::istream & read_binary(std::istream & is, std::string & s) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t string_size;
    is.read((char*)&string_size, size_size);
    char * str_contents = new char[string_size + 1];
    is.read(str_contents, string_size);
    str_contents[string_size] = '\0';
    s = str_contents;
    return is;
}

inline std::ostream & write_binary(std::ostream & os, ChannelHeader const & header) {
    if(header.version == 0) {
        write_binary(os, header.version);
        write_binary(os, header.physical_board_id);
        write_binary(os, header.board_serial_number);
        write_binary(os, header.physical_channel_type);
        write_binary(os, header.physical_channel_id);
        write_binary(os, header.caen_optical_link_number);
        write_binary(os, header.caen_optical_link_board_number);
        write_binary(os, header.caen_channel_number);
    } else {
        throw std::runtime_error("Can only write ChannelHeader version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, ChannelHeader & header) {
    read_binary(is, header.version);
    if(header.version == 0) {
        read_binary(is, header.physical_board_id);
        read_binary(is, header.board_serial_number);
        read_binary(is, header.physical_channel_type);
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

inline std::ostream & write_binary(std::ostream & os, CCMDigitizerReadout const & digi_readout) {
    if(digi_readout.version == 0) {
        write_binary(os, digi_readout.version);
        write_binary(os, digi_readout.channel_sizes);
        write_binary(os, digi_readout.channel_masks);
        write_binary(os, digi_readout.channel_temperatures);
        write_binary(os, digi_readout.board_event_numbers);
        write_binary(os, digi_readout.board_times);
        write_binary(os, digi_readout.samples);
    } else {
        throw std::runtime_error("Can only write CCMDigitizerReadout version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMDigitizerReadout & digi_readout) {
    read_binary(is, digi_readout.version);
    if(digi_readout.version == 0) {
        read_binary(is, digi_readout.channel_sizes);
        read_binary(is, digi_readout.channel_masks);
        read_binary(is, digi_readout.channel_temperatures);
        read_binary(is, digi_readout.board_event_numbers);
        read_binary(is, digi_readout.board_times);
        read_binary(is, digi_readout.samples);
    } else {
        throw std::runtime_error("Can only read CCMDigitizerReadout version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMTriggerReadout const & trigger_readout) {
    if(trigger_readout.version == 0) {
        write_binary(os, trigger_readout.version);
        write_binary(os, trigger_readout.event_number);
        write_binary(os, trigger_readout.computer_time);
        write_binary(os, trigger_readout.digitizer_readout);
    } else {
        throw std::runtime_error("Can only write CCMTriggerReadout version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMTriggerReadout & trigger_readout) {
    read_binary(is, trigger_readout.version);
    if(trigger_readout.version == 0) {
        read_binary(is, trigger_readout.event_number);
        read_binary(is, trigger_readout.computer_time);
        read_binary(is, trigger_readout.digitizer_readout);
    } else {
        throw std::runtime_error("Can only read CCMTriggerReadout version <= 0");
    }
    return is;
}

inline std::ostream & write_magic_number(std::ostream & os) {
    write_binary(os, lexical_magic_number, magic_size);
    return os;
}

inline std::istream & read_magic_number(std::istream & is, bool & result) {
    char file_identifier[magic_size];
    read_binary(is, file_identifier, magic_size);
    result = std::strncmp(file_identifier, lexical_magic_number, magic_size);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMData const & data) {
    if(data.version == 0) {
        write_magic_number(os);
        write_binary(os, data.version);
        write_binary(os, data.daq_config);
        write_binary(os, data.trigger_readout);
    } else {
        throw std::runtime_error("Can only write CCMData version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMData & data) {
    bool has_magic_number = false;
    read_magic_number(is, has_magic_number);
    if(not has_magic_number) {
        std::stringstream ss;
        ss << "Binary file must begin with the magic number '";
        ss << lexical_magic_number;
        ss << "'";
        throw std::runtime_error(ss.str());
    }
    read_binary(is, data.version);
    if(data.version == 0) {
        read_binary(is, data.daq_config);
        read_binary(is, data.trigger_readout);
    } else {
        throw std::runtime_error("Can only read CCMData version <= 0");
    }
    return is;
}

} // namespace Binary
} // namespace CCMAnalsysis

#endif // CCMAnalysis_BinaryUtilities_CXX
