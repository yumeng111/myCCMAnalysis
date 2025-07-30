#ifndef MARLEY_SIMULATOR_H
#define MARLEY_SIMULATOR_H

#include <string>
#include <vector>

#include <marley/Generator.hh>
#include <marley/FileManager.hh>

#include "icetray/I3Frame.h"
#include "icetray/I3ConditionalModule.h"
#include "dataclasses/physics/I3MCTree.h"
#include "phys-services/I3RandomService.h"


struct NuclearCascadeStep {
    int initial_level_index;
    double initial_level_energy;
    double initial_level_spin; //2J
    int initial_level_parity; // +1 or -1

    int final_level_index;
    double final_level_energy;
    double final_level_spin;
    int final_level_parity;

    double gamma_energy;

    double T12_ns;  // Half-life in ns
    double tau_ns;  // Mean lifetime in ns

    double sampled_delay_ns;
    double cumulative_time_ns;
};

struct LevelInfo {
    int level_index;
    double energy;
    double spin; //2J
    int parity;

    double T12_ns;  // Half-life in ns
    double tau_ns;  // Mean lifetime in ns

    struct Transition {
        double gamma_energy;
        double branching_ratio;  //called RI in MARLEY
        int final_level_index;
    };

    std::vector<Transition> transitions;

    bool operator<(double const energy) const {
        return this->energy < energy;
    }
};


class FileManager_ : public ::marley::FileManager {
friend class MarleySimulator;
    static std::string get_search_path() {
        return ::marley::FileManager::default_search_path_;
    }
    static void set_search_path(std::string const & search_path) {
        ::marley::FileManager::default_search_path_ = search_path;
    }
    static void set_marley_path(std::string const & marley_path) {
        std::string default_search_path;
        default_search_path = marley_path + "/data";
        default_search_path += ':' + marley_path + "/data/react";
        default_search_path += ':' + marley_path + "/data/structure";
        set_search_path(default_search_path);
    }
};


class MarleySimulator : public I3ConditionalModule {
public:
    MarleySimulator(const I3Context& context);
    ~MarleySimulator() = default;

    void Configure();
    void Simulation(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void FillSimulationFrame(I3FramePtr frame);
    void Finish();
private:
    bool seen_s_frame_;

    std::string output_mc_tree_name_;
    std::string input_mc_tree_name_;
    std::string marley_search_path_;
    // Add marley_generator_as member of the class
    marley::Generator marley_generator_;

    bool enable_gamma_time_offset_;
    bool save_levels_file_;
    std::string levels_filename_;
    std::vector<LevelInfo>::iterator ClosestLevel(std::vector<LevelInfo> & levels, double energy);
    void LoadK40Transitions(const std::string& filename);
    void AdjustGammaTimes(I3MCTreePtr mcTree, I3FramePtr frame);

    std::vector<LevelInfo> levels_map_;

    double SampleDelay(double mean_lifetime_ns);
    I3RandomServicePtr rng_;

    size_t frames_with_cascade = 0;
    size_t frames_with_metastable = 0;
    size_t frames_in_total = 0;
    std::vector<double> metastable_energy_diffs;

    SET_LOGGER("MarleySimulator");
};

#endif
