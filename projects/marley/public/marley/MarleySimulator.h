#ifndef MARLEY_SIMULATOR_H
#define MARLEY_SIMULATOR_H

#include <icetray/I3ConditionalModule.h>
#include <icetray/I3Frame.h>
#include "phys-services/I3RandomService.h"

#include <string>
#include <vector>

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
};


// Marley forward-declaration
namespace marley {
    class Generator;
}

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

    void LoadK40Transitions(const std::string& filename);
    void AdjustGammaTimes(I3MCTreePtr mcTree, I3FramePtr frame);

    std::map<int, LevelInfo> levels_map_;

    double SampleDelay(double mean_lifetime_ns);
    I3RandomServicePtr rng_;

    size_t frames_with_cascade = 0;
    size_t frames_with_unknowns = 0;
    size_t frames_with_metastable = 0;
    size_t frames_with_corrections_ = 0;
    size_t corrections_total_ = 0;

    SET_LOGGER("MarleySimulator");
};

#endif
