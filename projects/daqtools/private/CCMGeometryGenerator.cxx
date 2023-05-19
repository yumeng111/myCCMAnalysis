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

namespace position_mappings {
    I3Map<std::string, std::string> mapping_2022 = std::map<std::string, std::string>({
            {"FB0784", "C101R0"},
            {"FB0792", "C401R0"},
            {"FB0831", "C105R0"},
            {"FB0791", "C210R0"},
            {"FB0790", "C201R0"},
            {"FB0783", "C103R0"},
            {"FB0793", "C315R0"},
            {"FB0787", "C301R0"},
            {"FB0835", "C102R0"},
            {"FB0789", "C314R0"},
            {"FB0796", "C419R0"},
            {"FB0797", "C313R0"},
            {"FB0788", "C209R0"},
            {"FB0728", "C417R0"},
            {"FB0827", "C104R0"},
            {"FB0692", "C312R0"},
            {"FB0707", "C207R0"},
            {"FB0778", "C208R0"},
            {"FB0785", "C415R0"},
            {"FB0698", "C311R0"},
            {"FB0703", "C302R0"},
            {"FB0718", "C305R0"},
            {"FB0706", "C303R0"},
            {"FB0723", "C204R0"},
            {"FB0716", "C306R0"},
            {"FB0780", "C405R0"},
            {"FB0708", "C304R0"},
            {"FB0711", "C403R0"},
            {"FB0833", "C202R0"},
            {"FB0779", "C203R0"},
            {"FB0667", "C407R0"},
            {"FB0665", "C307R0"},
            {"FB0715", "C205R0"},
            {"FB0721", "C409R0"},
            {"FB0782", "C308R0"},
            {"FB0668", "C309R0"},
            {"FB0663", "C411R0"},
            {"FB0671", "C206R0"},
            {"FB0794", "C310R0"},
            {"FB0839", "C413R0"},
            {"FB0795", "C2R3"},
            {"FB0830", "C2R4"},
            {"FB0836", "C2R2"},
            {"FB0852", "C2R5"},
            {"FB0786", "C2R1"},
            {"FB0832", "C1R3"},
            {"FB0838", "C1R4"},
            {"FB0834", "C1R5"},
            {"FB0837", "C1R1"},
            {"FB0687", "C1R2"},
            {"FB0777", "C3R3"},
            {"FB0695", "C3R5"},
            {"FB0775", "C3R1"},
            {"FB0846", "C3R4"},
            {"FB0803", "C3R2"},
            {"FB0680", "C5R2"},
            {"FB0709", "C5R5"},
            {"FB0851", "C5R1"},
            {"FB0727", "C5R4"},
            {"FB0712", "C5R3"},
            {"FB0696", "C24R5"},
            {"FB0674", "C24R2"},
            {"FB0853", "C24R4"},
            {"FB0731", "C24R1"},
            {"FB0705", "C24R3"},
            {"FB0772", "C23R4"},
            {"FB0773", "C23R3"},
            {"FB0768", "C23R1"},
            {"FB0855", "C23R2"},
            {"FB0845", "C23R5"},
            {"FB0713", "C4R5"},
            {"FB0681", "C4R3"},
            {"FB0798", "C4R4"},
            {"FB0774", "C4R2"},
            {"FB0771", "C4R1"},
            {"FB0766", "C8R1"},
            {"FB0729", "C8R5"},
            {"FB0807", "C8R4"},
            {"FB0676", "C8R3"},
            {"FB0679", "C8R2"},
            {"FB0724", "C10R1"},
            {"FB0673", "C9R4"},
            {"FB0677", "C9R5"},
            {"FB0678", "C9R2"},
            {"FB0808", "C10R3"},
            {"FB0813", "C10R2"},
            {"FB0802", "C10R5"},
            {"FB0675", "C9R3"},
            {"FB0804", "C9R1"},
            {"FB0815", "C10R4"},
            {"FB0840", "C11R1"},
            {"FB0809", "C11R3"},
            {"FB0799", "C11R2"},
            {"FB0776", "C11R4"},
            {"FB0844", "C11R5"},
            {"FB0812", "C13R5"},
            {"FB0841", "C12R1"},
            {"FB0842", "C12R3"},
            {"FB0849", "C12R2"},
            {"FB0800", "C13R1"},
            {"FB0848", "C13R2"},
            {"FB0847", "C13R4"},
            {"FB0769", "C12R5"},
            {"FB0805", "C12R4"},
            {"FB0843", "C13R3"},
            {"FB0810", "C14R1"},
            {"FB0816", "C14R2"},
            {"FB0806", "C14R5"},
            {"FB0814", "C14R4"},
            {"FB0739", "C14R3"},
            {"FB0736", "C15R3"},
            {"FB0801", "C15R2"},
            {"FB0738", "C15R1"},
            {"FB0743", "C15R4"},
            {"FB0662", "C15R5"},
            {"FB0664", "C16R3"},
            {"FB0714", "C16R2"},
            {"FB07122", "C16R5"},
            {"FB0761", "C16R4"},
            {"FB0666", "C16R1"},
            {"FB0620", "C19R1"},
            {"FB0719", "C19R5"},
            {"FB0717", "C19R4"},
            {"FB0749", "C19R2"},
            {"FB0720", "C19R3"},
            {"FB0735", "C20R1"},
            {"FB0669", "C20R2"},
            {"FB0734", "C20R5"},
            {"FB0762", "C20R4"},
            {"FB0737", "C20R3"},
            {"FB0741", "C21R3"},
            {"FB0811", "C21R4"},
            {"FB0817", "C21R5"},
            {"FB0701", "C21R2"},
            {"FB0757", "C21R1"},
            {"FB0616", "C22R2"},
            {"FB0617", "C22R4"},
            {"FB0619", "C22R1"},
            {"FB0740", "C22R3"},
            {"FB0759", "C22R5"},
            {"FB0825", "C18R2"},
            {"FB0672", "C18R1"},
            {"FB0730", "C18R4"},
            {"FB0732", "C17R3"},
            {"FB0820", "C18R3"},
            {"FB0821", "C17R2"},
            {"FB0818", "C17R5"},
            {"FB0748", "C18R5"},
            {"FB0726", "C17R4"},
            {"FB0763", "C17R1"},
            {"FB0618", "C6R4"},
            {"FB0744", "C6R5"},
            {"FB0612", "C6R2"},
            {"FB0755", "C7R2"},
            {"FB0472", "C7R3"},
            {"FB0469", "C6R3"},
            {"FB0613", "C7R4"},
            {"FB0620", "C6R1"},
            {"FB0614", "C7R5"},
            {"FB0482", "C7R1"},
            {"FB0694", "C410R6"},
            {"FB0628", "C308R6"},
            {"FB0690", "C206R6"},
            {"FB0691", "C309R6"},
            {"FB0747", "C412R6"},
            {"FB0765", "C310R6"},
            {"FB0682", "C207R6"},
            {"FB0688", "C416R6"},
            {"FB0824", "C313R6"},
            {"FB0828", "C209R6"},
            {"FB0686", "C208R6"},
            {"FB0684", "C312R6"},
            {"FB0688", "C311R6"},
            {"FB0829", "C105R6"},
            {"FB0750", "C104R6"},
            {"FB0699", "C101R6"},
            {"FB0685", "C314R6"},
            {"FB0756", "C418R6"},
            {"FB0693", "C414R6"},
            {"FB0700", "C315R6"},
            {"FB0781", "C210R6"},
            {"FB0733", "C420R6"},
            {"FB0697", "C201R6"},
            {"FB0683", "C301R6"},
            {"FB0651", "C202R6"},
            {"FB0588", "C404R6"},
            {"FB0645", "C303R6"},
            {"FB0597", "C302R6"},
            {"FB0526", "C402R6"},
            {"FB0636", "C203R6"},
            {"FB0584", "C304R6"},
            {"FB0626", "C102R6"},
            {"FB0691", "C103R6"},
            {"FB0468", "C307R6"},
            {"FB0819", "C204R6"},
            {"FB0615", "C406R6"},
            {"FB0755", "C205R6"},
            {"FB0850", "C408R6"},
            {"FB0440", "C306R6"},
            {"FB0471", "C305R6"},
            {"V15", "VCB10"},
            {"V29", "VT10"},
            {"V34", "VCB6"},
            {"V11", "VCT9"},
            {"V33", "VCB12"},
            {"V27", "VCT11"},
            {"V28", "VCT13"},
            {"V8", "VT12"},
            {"V28", "VCB2"},
            {"V6", "VCT7"},
            {"V16", "VCT15"},
            {"V31", "VT18"},
            {"V25", "VT14"},
            {"V30", "VT16"},
            {"V36", "VT4"},
            {"V35", "VT6"},
            {"V16", "VCT17"},
            {"V30", "VT20"},
            {"V5", "VCT19"},
            {"V14", "VCT21"},
            {"V10", "VCT5"},
            {"V12", "VCT3"},
            {"V13", "VCT1"},
            {"V2", "VCT23"},
            {"V37", "VT24"},
            {"V38", "VT22"},
            {"V41", "VCB20"},
            {"V39", "VCB18"},
            {"V128", "VB9"},
            {"V44", "VB23"},
            {"V127", "VCB24"},
            {"V88", "VT2"},
            {"V129", "VCB22"},
            {"VU", "VB17"},
            {"VS", "VB15"},
            {"V55", "VB13"},
            {"V-Omega", "VB1"},
            {"V125", "VCB4"},
            {"V126", "VB5"},
    });
}

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
            position_mappings::mapping_2022);
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
        {5, 24},  // Veto VT and VB rings
        {6, 24}  // Veto VCT and VCB rings
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
                if( (row ==  1 and col % 4 == 1) or
                    (row ==  2 and col % 4 == 3) or
                    (row ==  4 and col % 4 == 0) or
                    (row ==  5 and col % 4 == 2)) {
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
            // Special case to handle missing trigger copy in 2022 data
            if(board.physical_board_id == "physical_board_18" and channel_idx == 15 and id == "" and type == "") {
                type = detail::tolower("Trigger Copy");
                id = detail::tolower("Trigger Copy 18");
            }
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
                // Special check for board 17 since its trigger copy was mislabeled in 2022 as "Trigger Copy 15"
                if(board.physical_board_id == "physical_board_17")
                    trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::BoardTriggerCopy, 17);
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
                log_fatal(("Unknown channel type: \"" + channel.physical_channel_type + "\"").c_str());
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
        } else {
            log_fatal(("Board " + board.physical_board_id + " does not have a trigger copy!").c_str());
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
