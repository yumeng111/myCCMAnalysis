#ifndef MARLEY_SIMULATOR_H
#define MARLEY_SIMULATOR_H

#include <string>
#include <vector>

#include <marley/Generator.hh>

#include "icetray/I3Frame.h"
#include "icetray/I3ConditionalModule.h"
#include "dataclasses/physics/I3MCTree.h"
#include "phys-services/I3RandomService.h"


struct NuclearCascadeStep {
    int initial_level_index;
    double initial_level_energy_keV;
    double initial_level_spin; //2J
    int initial_level_parity; // +1 or -1

    int final_level_index;
    double final_level_energy_keV;
    double final_level_spin;
    int final_level_parity;

    double gamma_energy_keV;

    double T12_ns;  // Half-life in ns
    double tau_ns;  // Mean lifetime in ns

    double sampled_delay_ns;
    double cumulative_time_ns;
};

struct LevelInfo {
    int level_index;
    double energy_keV;
    double spin; //2J
    int parity;

    double T12_ns;  // Half-life in ns
    double tau_ns;  // Mean lifetime in ns

    struct Transition {
        double gamma_energy_keV;
        double branching_ratio;  //called RI in MARLEY
        int final_level_index;
    };

    std::vector<Transition> transitions;

    bool operator<(double const energy) const {
        double keV = energy / I3Units::keV;
        return energy_keV < energy;
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
    LevelInfo & ClosestLevel(std::vector<LevelInfo> & levels, double energy);
    void LoadK40Transitions(const std::string& filename);
    void AdjustGammaTimes(I3MCTreePtr mcTree, I3FramePtr frame);

    std::vector<LevelInfo> levels_map_;

    double SampleDelay(double mean_lifetime_ns);
    I3RandomServicePtr rng_;

    size_t frames_with_cascade = 0;
    size_t frames_with_metastable = 0;
    size_t frames_in_total = 0;
    std::vector<double> metastable_energy_diffs;

    // Tolerance for matching energy sums from marley to real levels real levels from .dat
    static constexpr double match_tolerance_keV = 5.0;


    SET_LOGGER("MarleySimulator");
};

#endif
