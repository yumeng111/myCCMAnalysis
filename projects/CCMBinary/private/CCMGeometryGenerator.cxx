#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Orientation.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

class CCMGeometryGenerator : public I3Module {
    I3Map<std::string, std::string> pmt_positions_by_id;
public:
    CCMGeometryGenerator(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMGeometryGenerator);

namespace detail {
    std::string tolower(std::string const & s) {
        std::string r(s);
        std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c){return std::tolower(c);});
        return r;
    }
}

CCMGeometryGenerator::CCMGeometryGenerator(const I3Context& context) : I3Module(context) {
    AddParameter("PMTPositionsByID",
            "Map from PMTID strings to PMTPosition strings",
            I3Map<std::string, std::string>());
}

void CCMGeometryGenerator::Configure() {
    I3Map<std::string, std::string> positions;
    GetParameter("PMTPositionsByID", positions);
    for(auto const & p : positions) {
        pmt_positions_by_id.insert({detail::tolower(p.first), detail::tolower(p.second)});
    }
}

namespace detail {
    double default_start_pmt_number = 1;
    int default_num_pmts = 24;
    int default_middle_pmt_row = 3;

    double angle_from_cylinder_pmt_number(
            int pmt_number,
            double start_pmt_number = default_start_pmt_number,
            int num_pmts = default_num_pmts) {
        return double(pmt_number - start_pmt_number) / double(num_pmts) * 2.0 * M_PI;
    }

    std::pair<double, double> xy_position_from_cylinder_pmt_number(
            int pmt_number,
            double radius,
            double start_pmt_number = default_start_pmt_number,
            int num_pmts = default_num_pmts,
            double angular_offset = 0.0) {
        double angle = angle_from_cylinder_pmt_number(pmt_number, start_pmt_number, num_pmts);
        angle += angular_offset;
        return {radius * cos(angle), radius * sin(angle)};
    }

    /*
    double default_pmt_row_separation = 23.11;
    double default_row_direction = -1.0;
    double z_from_row_number(
            int pmt_row,
            int middle_pmt_row = default_middle_pmt_row,
            double pmt_row_separation = default_pmt_row_separation,
            int row_direction = default_row_direction) {
        return (pmt_row - pmt_middle_row) * row_direction * pmt_row_separation;
    }
    */

    // Z positions of the rows of pmts
    std::map<int, double> pmt_region_z_positions = {
        {-2, 75.00},
        {-1, 65.00},
        {0,  58.00},
        {1,  46.22},
        {2,  23.11},
        {3,   0.00},
        {4, -23.11},
        {5, -46.22},
        {6, -58.00},
        {7, -65.00},
        {8, -75.00},
    };

    // The radii of rings of pmts
    std::map<int, double> ring_radii = {
        {0,  96.0}, // Outer wall
        {1,  20.0}, // Inner ring on top/bottom
        {2,  42.0}, //
        {3,  64.2}, //
        {4,  85.5}, // Outer ring on top/bottom
        {5, 115.0}  // Veto PMT ring
    };

    // The number of possible pmt positions in each ring (not all positions will be filled)
    std::map<int, int> ring_pmt_pos_count = {
        {0, 24}, // Outer wall
        {1,  5}, // Inner ring on top/bottom
        {2, 10}, //
        {3, 15}, //
        {4, 20}, // Outer ring on top/bottom
    };

    std::set<std::tuple<int, int, int>> cap_uncoated_pmts = {
        // Ring, col, row
        {1,  1, 0},
        {1,  3, 0},
        {2,  3, 0},
        {2,  8, 0},
        {3,  1, 0},
        {3,  8, 0},
        {4,  5, 0},
        {4, 15, 0},
        {1,  3, 6},
        {1,  5, 6},
        {2,  2, 6},
        {2,  7, 6},
        {3,  7, 6},
        {3, 14, 6},
        {4,  4, 6},
        {4, 14, 6},
    };

    I3Position get_pmt_cap_position(int pmt_row, int ring_number, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
        double z = pmt_region_z_positions[pmt_row];
        std::pair<double, double> xy = xy_position_from_cylinder_pmt_number(
                pmt_number, 
                ring_radii[ring_number], 
                starting_pmt_number, 
                ring_pmt_pos_count[ring_number],
                angular_offset);
        I3Position pos(xy.first, xy.second, z);
        return pos;
    }

    I3Position get_pmt_wall_position(int pmt_row, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
        return get_pmt_cap_position(pmt_row, 0, pmt_number, starting_pmt_number, angular_offset);
    }

    I3Position get_pmt_veto_position(int pmt_row, int pmt_number, double starting_pmt_number = 0.5, double angular_offset = 0.0) {
        int pmt_veto_ring = 5;
        return get_pmt_cap_position(pmt_row, pmt_veto_ring, pmt_number, starting_pmt_number, angular_offset);
    }

    I3Orientation get_pmt_wall_orientation(int pmt_row, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
        I3Position pos = get_pmt_wall_position(pmt_row, pmt_number, starting_pmt_number, angular_offset);
        I3Position dir = I3Position(0, 0, 0) - pos;
        dir.Normalize();
        I3Position up(0, 0, 1);
        return I3Orientation(dir.GetX(), dir.GetY(), dir.GetZ(), up.GetX(), up.GetY(), up.GetZ());
    }

    I3Orientation get_pmt_cap_orientation(int pmt_row, int ring_number, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
        I3Position pos = get_pmt_cap_position(pmt_row, ring_number, pmt_number, starting_pmt_number, angular_offset);
        I3Position up = I3Position(0, 0, 0) - pos;
        up.Normalize();
        double z_dir = pos.GetZ() > 0 ? -1 : 1;
        I3Position dir(0, 0, z_dir);
        return I3Orientation(dir.GetX(), dir.GetY(), dir.GetZ(), up.GetX(), up.GetY(), up.GetZ());
    }

    I3Orientation get_pmt_veto_orientation(int pmt_row, int pmt_number, double up_down, double starting_pmt_number = 0.5, double angular_offset = 0.0) {
        I3Position pos = get_pmt_veto_position(pmt_row, pmt_number, starting_pmt_number, angular_offset);
        I3Position up = I3Position(0, 0, 0) - pos;
        up.Normalize();
        I3Position dir(0, 0, up_down);
        dir.Normalize();
        return I3Orientation(dir.GetX(), dir.GetY(), dir.GetZ(), up.GetX(), up.GetY(), up.GetZ());
    }

    std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> ParsePMT8inPosition(std::string position_string) {
        I3Position position;
        I3Orientation orientation;
        CCMPMTKey key;
        CCMOMGeo::OMType omtype = CCMOMGeo::OMType::CCM8inCoated;
        if(position_string.size() < 4)
            throw std::runtime_error("PMT position string must be at least 4 characters long: " + position_string);
        if(position_string[0] == 'c') {
            size_t r_pos = position_string.find("r");
            if(r_pos == std::string::npos)
                throw std::runtime_error("PMT position starts with \"C\" but does not contain \"R\": " + position_string);
            size_t n_row_chars = r_pos - 1;
            bool on_caps = n_row_chars > 2;
            int row;
            int pmt_number = 0;
            if(on_caps) {
                int col = std::atoi(position_string.substr(2, n_row_chars - 1).c_str());
                int ring = std::atoi(position_string.substr(1, 1).c_str());
                row = std::atoi(position_string.substr(r_pos + 1, std::string::npos).c_str());
                position = get_pmt_cap_position(row, ring, col);
                orientation = get_pmt_cap_orientation(row, ring, col);
                for(std::pair<int, int> const & p : ring_pmt_pos_count)
                    if(p.first < ring and p.first > 0)
                        pmt_number += p.second;
                pmt_number += col;
                std::tuple<int, int, int> cap_coated_key = {ring, col, row};
                if(cap_uncoated_pmts.count(cap_coated_key) > 0)
                    omtype = CCMOMGeo::OMType::CCM8inUncoated;
            } else {
                int col = std::atoi(position_string.substr(1, r_pos - 1).c_str());
                row = std::atoi(position_string.substr(r_pos + 1, std::string::npos).c_str());
                position = get_pmt_wall_position(row, col);
                orientation = get_pmt_wall_orientation(row, col);
                pmt_number = col;
                if( (row ==  2 and col % 4 == 1) or
                    (row ==  1 and col % 4 == 3) or
                    (row == -1 and col % 4 == 0) or
                    (row == -2 and col % 4 == 2)) {
                    omtype = CCMOMGeo::OMType::CCM8inUncoated;
                }
            }
            key = CCMPMTKey(row, pmt_number);
        } else {
            throw std::runtime_error("PMT position must start with \"C\": " + position_string);
        }
        return {position, orientation, key, omtype};
    }

/* Regions:
 * Row: 8 VT and VB on bottom
 * Row: 7 VCB on bottom
 * Row: 6 cylinder bottom
 * Row: 5 cylinder wall bottom row
 * Row: 4 cylinder wall
 * Row: 3 cylinder wall middle row
 * Row: 2 cylinder wall
 * Row: 1 cylinder wall top row
 * Row: 0 cylinder top
 * Row: -1 VCT on top
 * */
    std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> ParsePMT1inPosition(std::string position_string) {
        I3Position position;
        I3Orientation orientation;
        CCMPMTKey key;
        CCMOMGeo::OMType omtype = CCMOMGeo::OMType::CCM1in;
        if(position_string[0] == 'v') {
            int row;
            double up_down;
            size_t char_pos;
            if((char_pos = position_string.find("vct")) == 0) {
                char_pos += 3; // Skip the VCT chars
                row = -1;
                up_down = -1; // Pointing down
            } else if((char_pos = position_string.find("vcb")) == 0) {
                char_pos += 3; // Skip the VCB chars
                row = 7;
                up_down = 1; // Pointing up
            } else if((char_pos = position_string.find("vt")) == 0) {
                char_pos += 2; // Skip the VT chars
                row = 8;
                up_down = 1; // Pointing up
            } else if((char_pos = position_string.find("vb")) == 0) {
                char_pos += 2; // Skip the VB chars
                row = 8;
                up_down = -1; // Pointing down
            } else {
                throw std::runtime_error("PMT position string must start with \"VCT\", \"VCB\", \"VT\", or \"VB\": " + position_string);
            }
            int pmt_number = std::atoi(position_string.substr(char_pos, std::string::npos).c_str());
            position = get_pmt_veto_position(row, pmt_number);
            orientation = get_pmt_veto_orientation(row, pmt_number, up_down);
            key = CCMPMTKey(row, pmt_number);
        } else {
            throw std::runtime_error("PMT position string must start with \"V\": " + position_string);
        }
        return {position, orientation, key, omtype};
    }
}

void CCMGeometryGenerator::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() != I3Frame::Geometry) {
        PushFrame(frame);
        return;
    }

    CCMAnalysis::Binary::CCMDAQConfig config = frame->Get<CCMAnalysis::Binary::CCMDAQConfig>("CCMDAQConfig");

    size_t absolute_idx = 0;
    for(size_t board_idx = 0; board_idx < config.digitizer_boards.size(); ++board_idx) {
        CCMAnalysis::Binary::DigitizerBoard const & board = config.digitizer_boards[board_idx];
        size_t n_channels = board.channels.size();
        for(size_t channel_idx = 0; channel_idx < n_channels; ++channel_idx, ++absolute_idx) {
            CCMAnalysis::Binary::ChannelHeader const & channel = board.channels[channel_idx];
            std::string type = detail::tolower(channel.physical_channel_type);
            std::string id = detail::tolower(channel.physical_channel_id);
            I3Position position;
            I3Orientation orientation;
            CCMPMTKey pmt_key;
            CCMTriggerKey trigger_key;
            CCMOMGeo::OMType omtype;
            bool is_trigger = false;
            bool is_sensor = false;
/* Regions:
 * Row: 8 VT and VB on bottom
 * Row: 7 VCB on bottom
 * Row: 6 cylinder bottom
 * Row: 5 cylinder wall bottom row
 * Row: 4 cylinder wall
 * Row: 3 cylinder wall middle row
 * Row: 2 cylinder wall
 * Row: 1 cylinder wall top row
 * Row: 0 cylinder top
 * Row: -1 VCT on top
 * */
            if(pmt_positions_by_id.count(id) == 0) {
                std::stringstream ss;
                ss << "PMT ID \"" << id << "\" not found in supplied PMTPositionByID map";
                throw std::runtime_error(ss.str());
            }

            if(type == detail::tolower("PMT 1in")) {
                is_sensor = true;
                std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> p = detail::ParsePMT1inPosition(id);
                position = std::get<0>(p);
                orientation = std::get<1>(p);
                pmt_key = std::get<2>(p);
                omtype = std::get<3>(p);
            } else if(type == detail::tolower("PMT 8in")) {
                is_sensor = true;
                std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> p = detail::ParsePMT8inPosition(id);
                position = std::get<0>(p);
                orientation = std::get<1>(p);
                pmt_key = std::get<2>(p);
                omtype = std::get<3>(p);
            } else if(type == detail::tolower("EJ")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::EJ301;
            } else if(type == detail::tolower("EJ301")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::EJ301;
            } else if(type == detail::tolower("Flight Path 3 Monitor")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::EJ301;
            } else if(type == detail::tolower("Trigger Copy")) {
                is_trigger = true;
            } else if(type == detail::tolower("BCM Monitor")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::BeamCurrentMonitor;
            } else if(type == detail::tolower("BEAM Trigger")) {
                is_trigger = true;
            } else if(type == detail::tolower("STROBE Trigger")) {
                is_trigger = true;
            } else if(type == detail::tolower("LEDTOP Trigger")) {
                is_trigger = true;
            } else if(type == detail::tolower("LEDBOTTOM Trigger")) {
                is_trigger = true;
            } else if(type == detail::tolower("Cosmic Watch Trigger")) {
                is_trigger = true;
            }
        }
    }

    // Put the output objects in the frame and push the frame
    // frame->Put("CCMDigitalReadout", std::get<0>(readout), I3Frame::DAQ);
    // frame->Put("CCMTriggers", std::get<1>(readout), I3Frame::DAQ);
    // ++counter;
}

void CCMGeometryGenerator::Finish() {
    // log_notice_stream(
    //     "Merged " << offsets.size() << " DAQ streams into " << counter << " separate triggers. Encountered "
    //     << counter-incomplete_counter << " complete triggers and " << incomplete_counter << " incomplete triggers");
}
