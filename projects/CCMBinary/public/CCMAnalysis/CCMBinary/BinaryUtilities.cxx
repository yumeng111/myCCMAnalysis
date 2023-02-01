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

inline std::ostream & write_binary(std::ostream & os, double const & x) {
    constexpr size_t size = sizeof(double);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, double & x) {
    constexpr size_t size = sizeof(double);
    is.read((char*)&x, size);
    return is;
}

inline std::ostream & write_binary(std::ostream & os, long double const & x) {
    constexpr size_t size = sizeof(long double);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, long double & x) {
    constexpr size_t size = sizeof(long double);
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

inline std::ostream & write_binary(std::ostream & os, struct timespec const & x) {
    constexpr size_t size = sizeof(struct timespec);
    os.write((char*)&x, size);
    return os;
}

inline std::istream & read_binary(std::istream & is, struct timespec & s) {
    constexpr size_t size = sizeof(struct timespec);
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
        write_binary(os, board.trigger_out_type);
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
        read_binary(is, board.trigger_out_type);
    } else {
        throw std::runtime_error("Can only read DigitizerBoard version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMDAQMachineConfig const & config) {
    if(config.version == 0) {
        write_binary(os, config.version);
        write_binary(os, config.machine_identifier);
        write_binary(os, config.num_digitizer_boards);
        write_binary(os, config.num_channels);
        write_binary(os, config.num_samples);
        write_binary(os, config.trigger_percent_after);
        write_binary(os, config.trigger_time_tolerance);
        write_binary(os, config.missed_trigger_tolerance);
        write_binary(os, config.offset_estimate_min_triggers);
        write_binary(os, config.offset_estimate_abs_error_threshold);
        write_binary(os, config.offset_estimate_rel_error_threshold);
        write_binary(os, config.offset_estimate_tau);
    } else {
        throw std::runtime_error("Can only write CCMDAQMachineConfig version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMDAQMachineConfig & config) {
    read_binary(is, config.version);
    if(config.version == 0) {
        read_binary(is, config.machine_identifier);
        read_binary(is, config.num_digitizer_boards);
        read_binary(is, config.num_channels);
        read_binary(is, config.num_samples);
        read_binary(is, config.trigger_percent_after);
        read_binary(is, config.trigger_time_tolerance);
        read_binary(is, config.missed_trigger_tolerance);
        read_binary(is, config.offset_estimate_min_triggers);
        read_binary(is, config.offset_estimate_abs_error_threshold);
        read_binary(is, config.offset_estimate_rel_error_threshold);
        read_binary(is, config.offset_estimate_tau);
    } else {
        throw std::runtime_error("Can only read CCMDAQMachineConfig version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMDAQConfig const & config) {
    if(config.version == 0) {
        write_binary(os, config.version);
        write_binary(os, config.machine_configurations);
        write_binary(os, config.digitizer_boards);
    } else {
        throw std::runtime_error("Can only write CCMDAQConfig version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMDAQConfig & config) {
    read_binary(is, config.version);
    if(config.version == 0) {
        read_binary(is, config.machine_configurations);
        read_binary(is, config.digitizer_boards);
    } else {
        throw std::runtime_error("Can only read CCMDAQConfig version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMTrigger const & digi_readout) {
    if(digi_readout.version == 0) {
        write_binary(os, digi_readout.version);
        write_binary(os, digi_readout.channel_sizes);
        write_binary(os, digi_readout.channel_masks);
        write_binary(os, digi_readout.channel_temperatures);
        write_binary(os, digi_readout.board_event_numbers);
        write_binary(os, digi_readout.board_times);
        write_binary(os, digi_readout.board_computer_times);
    } else {
        throw std::runtime_error("Can only write CCMTrigger version <= 0");
    }
    return os;
}

inline std::istream & read_binary(std::istream & is, CCMTrigger & digi_readout) {
    read_binary(is, digi_readout.version);
    if(digi_readout.version == 0) {
        read_binary(is, digi_readout.channel_sizes);
        read_binary(is, digi_readout.channel_masks);
        read_binary(is, digi_readout.channel_temperatures);
        read_binary(is, digi_readout.board_event_numbers);
        read_binary(is, digi_readout.board_times);
        read_binary(is, digi_readout.board_computer_times);
    } else {
        throw std::runtime_error("Can only read CCMTrigger version <= 0");
    }
    return is;
}

inline std::ostream & write_binary(std::ostream & os, CCMTriggerReadout const & trigger_readout) {
    if(trigger_readout.version == 0) {
        write_binary(os, trigger_readout.version);
        write_binary(os, trigger_readout.event_number);
        write_binary(os, trigger_readout.computer_time);
        write_binary(os, trigger_readout.triggers);
        write_binary(os, trigger_readout.samples);
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
        read_binary(is, trigger_readout.triggers);
        read_binary(is, trigger_readout.samples);
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
    result = std::strncmp(file_identifier, lexical_magic_number, magic_size) == 0;
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

inline std::istream & read_config(std::istream & is, CCMDAQConfig & config) {
    bool has_magic_number = false;
    read_magic_number(is, has_magic_number);
    if(not has_magic_number) {
        std::stringstream ss;
        ss << "Binary file must begin with the magic number '";
        ss << lexical_magic_number;
        ss << "'";
        throw std::runtime_error(ss.str());
    }
    uint32_t version = 0;
    read_binary(is, version);
    if(version == 0) {
        read_binary(is, config);
    } else {
        throw std::runtime_error("Can only read CCMData version <= 0");
    }
    return is;
}

inline std::istream & read_size(std::istream & is, size_t & size) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t read_size;
    is.read((char*)&read_size, size_size);
    size = size_t(read_size);
    return is;
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

inline int32_t subtract_times(uint32_t t1, uint32_t t0) {
    // Compute the difference of two unsigned time counters with one overflow bit

    // Cover the trivial case
    if(t1 == t0) {
        return 0;
    }

    // We know the 32nd bit is an overflow bit so we must mask it out
    bool overflow_0 = t0 & (0x1 << 31);
    bool overflow_1 = t1 & (0x1 << 31);
    bool different_overflow = overflow_0 ^ overflow_1;
    int32_t time_0 = uint32_t(t0 & 0x7FFFFFFF);
    int32_t time_1 = uint32_t(t1 & 0x7FFFFFFF);

    int32_t diff;
    if(different_overflow) {
        // We know which time is greater
        // But one time has crossed the overflow boundary
        if(overflow_0) {
            // Convert t1 to a negative number in the frame of t0's zero then subtract
            diff = (time_1 - 0x7FFFFFFE) - time_0;
        } else {
            // Convert t0 to a negative number in the frame of t1's zero then subtract
            diff = time_1 - int32_t(int32_t(time_0) - int32_t(0x7FFFFFFE));
        }
    } else {
        // We don't know which time is greater
        // Use the difference with the smallest magnitude

        // Cover the trivial case
        int32_t diff_0 = time_1 - time_0;

        // Cover the case where the difference crosses the overflow boundary
        int32_t diff_1;
        if(time_0 < time_1) {
            // Convert t1 to a negative number in the frame of t0's zero then subtract
            diff_1 = (time_1 - 0x7FFFFFFE) - time_0;
        } else {
            // Convert t0 to a negative number in the frame of t1's zero then subtract
            diff_1 = time_1 - (time_0 - 0x7FFFFFFE);
        }
        if(std::abs(diff_0) < std::abs(diff_1)) {
            diff = diff_0;
        } else {
            diff = diff_1;
        }
    }
    return diff;
}


} // namespace Binary
} // namespace CCMAnalsysis

#endif // CCMAnalysis_BinaryUtilities_CXX
