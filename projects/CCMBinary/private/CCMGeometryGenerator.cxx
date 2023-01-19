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
    CCMAnalysis::Binary::CCMDAQConfig cached_config;
    size_t n_configs_seen = 0;;
    bool squash_duplicate_configs;
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

CCMGeometryGenerator::CCMGeometryGenerator(const I3Context& context) : I3Module(context), n_configs_seen(0) {
    AddParameter("PMTPositionsByID",
            "Map from PMTID strings to PMTPosition strings",
            I3Map<std::string, std::string>());
    AddParameter("SquashDuplicateConfigs",
            "Prevent duplicate CCMDAQConfig objects from producing additional geometries.",
            bool(true));
}

void CCMGeometryGenerator::Configure() {
    I3Map<std::string, std::string> positions;
    GetParameter("PMTPositionsByID", positions);
    for(auto const & p : positions) {
        pmt_positions_by_id.insert({detail::tolower(p.first), detail::tolower(p.second)});
    }
    GetParameter("SquashDuplicateConfigs", squash_duplicate_configs);
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
        {5, 101.6},  // Veto VT and VB rings
        {6, 111.6}  // Veto VCT and VCB rings
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
        {3,  6, 6},
        {3, 13, 6},
        {4,  2, 6},
        {4, 12, 6},
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
        int pmt_veto_ring;
        if(pmt_row >= -1 and pmt_row <= 7)
            pmt_veto_ring = 5;
        else
            pmt_veto_ring = 6;
        return get_pmt_cap_position(pmt_row, pmt_veto_ring, pmt_number, starting_pmt_number, angular_offset);
    }

    I3Orientation get_pmt_wall_orientation(int pmt_row, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
        I3Position pos = get_pmt_wall_position(pmt_row, pmt_number, starting_pmt_number, angular_offset);
        I3Position dir = I3Position(0, 0, 0) - pos;
        dir.SetZ(0);
        dir.Normalize();
        I3Position up(0, 0, 1);
        return I3Orientation(dir.GetX(), dir.GetY(), dir.GetZ(), up.GetX(), up.GetY(), up.GetZ());
    }

    I3Orientation get_pmt_cap_orientation(int pmt_row, int ring_number, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
        I3Position pos = get_pmt_cap_position(pmt_row, ring_number, pmt_number, starting_pmt_number, angular_offset);
        I3Position up = I3Position(0, 0, 0) - pos;
        up.SetZ(0);
        up.Normalize();
        double z_dir = pos.GetZ() > 0 ? -1 : 1;
        I3Position dir(0, 0, z_dir);
        return I3Orientation(dir.GetX(), dir.GetY(), dir.GetZ(), up.GetX(), up.GetY(), up.GetZ());
    }

    I3Orientation get_pmt_veto_orientation(int pmt_row, int pmt_number, double up_down, int inward, double starting_pmt_number = 0.5, double angular_offset = 0.0) {
        I3Position pos = get_pmt_veto_position(pmt_row, pmt_number, starting_pmt_number, angular_offset);
        I3Position dir;
        I3Position up;
        if(inward){
          dir = I3Position(0, 0, 0) - pos;
          dir.SetZ(0);
          dir.Normalize();
          up = I3Position(0,0,1);
        }
        else {
          up = I3Position(0, 0, 0) - pos;
          up.SetZ(0);
          up.Normalize();
          dir = I3Position(0, 0, up_down);
          dir.Normalize();
        }
        return I3Orientation(dir.GetX(), dir.GetY(), dir.GetZ(), up.GetX(), up.GetY(), up.GetZ());
    }

    std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> ParsePMT8inPosition(std::string position_string) {
        I3Position position;
        I3Orientation orientation;
        CCMPMTKey key;
        CCMOMGeo::OMType omtype = CCMOMGeo::OMType::CCM8inCoated;
        if(position_string.size() < 4)
            throw std::runtime_error("PMT position string must be at least 4 characters long: \"" + position_string + "\"");
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
                for(std::pair<const int, int> const & p : ring_pmt_pos_count)
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
 * 10: BCM
 * 9: External sensors
 * Row: 8 VB on bottom
 * Row: 7 VCB on bottom
 * Row: 6 cylinder bottom
 * Row: 5 cylinder wall bottom row
 * Row: 4 cylinder wall
 * Row: 3 cylinder wall middle row
 * Row: 2 cylinder wall
 * Row: 1 cylinder wall top row
 * Row: 0 cylinder top
 * Row: -1 VCT on top
 * Row: -2 VT on top
 * */
    std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> ParsePMT1inPosition(std::string position_string) {
        I3Position position;
        I3Orientation orientation;
        CCMPMTKey key;
        CCMOMGeo::OMType omtype = CCMOMGeo::OMType::CCM1in;
        if(position_string[0] == 'v') {
            int row;
            double up_down;
            int inward;
            size_t char_pos;
            if((char_pos = position_string.find("vct")) == 0) {
                char_pos += 3; // Skip the VCT chars
                row = -1;
                up_down = -1; // Pointing down
                inward = 0; // Not pointing inward
            } else if((char_pos = position_string.find("vcb")) == 0) {
                char_pos += 3; // Skip the VCB chars
                row = 7;
                up_down = 1; // Pointing up
                inward = 0; // Not pointing inward
            } else if((char_pos = position_string.find("vt")) == 0) {
                char_pos += 2; // Skip the VT chars
                row = -2;
                up_down = 0; // Pointing inward
                inward = 1; // Pointing inward
            } else if((char_pos = position_string.find("vb")) == 0) {
                char_pos += 2; // Skip the VB chars
                row = 8;
                up_down = 0; // Pointing inward
                inward = 1; // Pointing inward
            } else {
                throw std::runtime_error("PMT position string must start with \"VCT\", \"VCB\", \"VT\", or \"VB\": " + position_string);
            }
            int pmt_number = std::atoi(position_string.substr(char_pos, std::string::npos).c_str());
            position = get_pmt_veto_position(row, pmt_number);
            orientation = get_pmt_veto_orientation(row, pmt_number, up_down, inward);
            key = CCMPMTKey(row, pmt_number);
        } else {
            throw std::runtime_error("PMT position string must start with \"V\": " + position_string);
        }
        return {position, orientation, key, omtype};
    }

    CCMTriggerKey ParseTriggerCopy(std::string copy_string) {
        std::string prefix = detail::tolower("Trigger Copy ");
        size_t char_pos = copy_string.find(prefix);
        if(char_pos == std::string::npos) {
            throw std::runtime_error("Trigger Copy ID must begin with \"Trigger Copy \". Saw: " + copy_string);
        }
        char_pos += prefix.size();

        int trigger_copy_number = std::atoi(copy_string.substr(char_pos, std::string::npos).c_str());
        if(trigger_copy_number < 0)
            throw std::runtime_error("Trigger copy number must be positive: " + copy_string);

        return CCMTriggerKey(CCMTriggerKey::TriggerType::BoardTriggerCopy, size_t(trigger_copy_number));
    }
}

void CCMGeometryGenerator::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() != I3Frame::Geometry) {
        PushFrame(frame);
        return;
    }


    CCMAnalysis::Binary::CCMDAQConfig config = frame->Get<CCMAnalysis::Binary::CCMDAQConfig>("CCMDAQConfig");
    n_configs_seen += 1;
    bool is_duplicate_config = config == cached_config;
    cached_config = config;
    if(n_configs_seen > 1 and is_duplicate_config and squash_duplicate_configs)
        return;

    boost::shared_ptr<CCMGeometry> geometry = boost::make_shared<CCMGeometry>();

    unsigned int fp3_monitor_count = 0;
    unsigned int ej_idx = 2;
    unsigned int bcm_idx = 1;
// enum TriggerType {UnknownType = 0, BoardTriggerCopy = 10, BeamTrigger = 20, StrobeTrigger = 30, CosmicTrigger = 40, LaserTrigger = 50, LEDTopTrigger = 60, LEDBottomTrigger = 70};
    unsigned int beam_trigger_count = 0;
    unsigned int strobe_trigger_count = 0;
    unsigned int cosmic_trigger_count = 0;
    unsigned int laser_trigger_count = 0;
    unsigned int led_top_trigger_count = 0;
    unsigned int led_bottom_trigger_count = 0;

    size_t absolute_idx = 0;
    for(size_t board_idx = 0; board_idx < config.digitizer_boards.size(); ++board_idx) {
        CCMAnalysis::Binary::DigitizerBoard const & board = config.digitizer_boards[board_idx];
        size_t n_channels = board.channels.size();
        std::vector<CCMPMTKey> board_pmt_keys;
        board_pmt_keys.reserve(n_channels);
        CCMTriggerKey board_trigger_copy;
        bool board_has_trigger_copy = false;
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

            if(type == detail::tolower("PMT 1in")) {
                is_sensor = true;
                if(pmt_positions_by_id.count(id) == 0) {
                    if(id == "") {
                        log_warn("PMT ID is empty but channel type indicates this is a 1in PMT!");
                        continue;
                    } else {
                        std::stringstream ss;
                        ss << "PMT ID \"" << id << "\" not found in supplied PMTPositionByID map";
                        throw std::runtime_error(ss.str());
                    }
                }
                std::string position_string = pmt_positions_by_id[id];
                std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> p = detail::ParsePMT1inPosition(position_string);
                position = std::get<0>(p);
                orientation = std::get<1>(p);
                pmt_key = std::get<2>(p);
                omtype = std::get<3>(p);
            } else if(type == detail::tolower("PMT 8in")) {
                is_sensor = true;
                if(pmt_positions_by_id.count(id) == 0) {
                    if(id == "") {
                        log_warn("PMT ID is empty but channel type indicates this is an 8in PMT!");
                        continue;
                    } else {
                        std::stringstream ss;
                        ss << "PMT ID \"" << id << "\" not found in supplied PMTPositionByID map";
                        throw std::runtime_error(ss.str());
                    }
                }
                std::string position_string = pmt_positions_by_id[id];
                std::tuple<I3Position, I3Orientation, CCMPMTKey, CCMOMGeo::OMType> p = detail::ParsePMT8inPosition(position_string);
                position = std::get<0>(p);
                orientation = std::get<1>(p);
                pmt_key = std::get<2>(p);
                omtype = std::get<3>(p);
            } else if(type == detail::tolower("EJ")
                   or type == detail::tolower("EJ301")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::EJ301;
                position = I3Position(0,0, -1000);
                orientation = I3Orientation(0,0,1,1,0,0);
                unsigned int sensor_number = ej_idx;
                int region = 9;
                ej_idx += 1;
                pmt_key = CCMPMTKey(region, sensor_number);
            } else if(type == detail::tolower("Flight Path 3 Monitor")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::EJ301;
                position = I3Position(0, 0, -1000);
                orientation = I3Orientation(0,0,1,1,0,0);
                unsigned int sensor_number = 0;
                int region = 9;
                if(fp3_monitor_count > 0) {
                    sensor_number = ej_idx;
                    ej_idx += 1;
                } else {
                    sensor_number = 1;
                    fp3_monitor_count += 1;
                }
                pmt_key = CCMPMTKey(region, sensor_number);
            } else if (type == detail::tolower("CsI Crystal Single PMT")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::CsISinglePMT;
                position = I3Position(0, 0, -1000);
                orientation = I3Orientation(0,0,1,1,0,0);
                unsigned int sensor_number = 10;
                int region = 9;
                pmt_key = CCMPMTKey(region, sensor_number);
            } else if(type == detail::tolower("Trigger Copy")) {
                is_trigger = true;
                trigger_key = detail::ParseTriggerCopy(id);
                board_trigger_copy = trigger_key;
                board_has_trigger_copy = true;
            } else if(type == detail::tolower("BCM Monitor")) {
                is_sensor = true;
                omtype = CCMOMGeo::OMType::BeamCurrentMonitor;
                unsigned int sensor_number = bcm_idx;
                int region = 10;
                bcm_idx += 1;
                pmt_key = CCMPMTKey(region, sensor_number);
            } else if(type == detail::tolower("BEAM Trigger")) {
                is_trigger = true;
                beam_trigger_count += 1;
                trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::BeamTrigger, beam_trigger_count);
            } else if(type == detail::tolower("STROBE Trigger")) {
                is_trigger = true;
                strobe_trigger_count += 1;
                trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::StrobeTrigger, strobe_trigger_count);
            } else if(type == detail::tolower("LEDTOP Trigger")) {
                is_trigger = true;
                led_top_trigger_count += 1;
                trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::LEDTopTrigger, led_top_trigger_count);
            } else if(type == detail::tolower("LEDBOTTOM Trigger")) {
                is_trigger = true;
                led_bottom_trigger_count += 1;
                trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::LEDBottomTrigger, led_bottom_trigger_count);
            } else if(type == detail::tolower("Cosmic Watch Trigger")) {
                is_trigger = true;
                cosmic_trigger_count += 1;
                trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::CosmicTrigger, cosmic_trigger_count);
            } else if(type == detail::tolower("LASER Trigger")) {
                is_trigger = true;
                laser_trigger_count += 1;
                trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::LaserTrigger, laser_trigger_count);
            } else if(type == detail::tolower("")) {
                continue;
            } else {
                throw std::runtime_error("Unknown channel type: \"" + channel.physical_channel_type + "\"");
            }

            if(is_trigger) {
                geometry->trigger_channel_map.insert({trigger_key, absolute_idx});
            } else if(is_sensor) {
                geometry->pmt_channel_map.insert({pmt_key, absolute_idx});
                CCMOMGeo om;
                om.position = position;
                om.orientation = orientation;
                om.omtype = omtype;
                geometry->pmt_geo.insert({pmt_key, om});
                board_pmt_keys.push_back(pmt_key);
            } else {
                throw std::runtime_error("Channel must be either trigger or sensor: " + channel.physical_channel_type);
            }
        }

        if(board_has_trigger_copy) {
            for(CCMPMTKey const & k : board_pmt_keys) {
                geometry->trigger_copy_map.insert({k, board_trigger_copy});
            }
        }
    }

    frame->Put("CCMGeometry", geometry);
    PushFrame(frame);
}

void CCMGeometryGenerator::Finish() {
    // log_notice_stream(
    //     "Merged " << offsets.size() << " DAQ streams into " << counter << " separate triggers. Encountered "
    //     << counter-incomplete_counter << " complete triggers and " << incomplete_counter << " incomplete triggers");
}
