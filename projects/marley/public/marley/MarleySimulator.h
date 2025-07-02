#ifndef MARLEY_SIMULATOR_H
#define MARLEY_SIMULATOR_H

#include <icetray/I3ConditionalModule.h>
#include <icetray/I3Frame.h>

#include <string>
#include <vector>

struct K40Transition {
    double initial_energy;
    double final_energy;
    double gamma_energy;
    double relative_intensity;
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

private:
    bool seen_s_frame_;

    std::string output_mc_tree_name_;
    std::string input_mc_tree_name_;
    std::string marley_search_path_;
    // Add marley_generator_as member of the class
    marley::Generator marley_generator_;

    std::vector<K40Transition> k40_transitions_;
    std::vector<double> energy_levels_;

    void LoadK40Transitions(const std::string& filename);
    void PrintLevelsAndTransitions(const std::string& outfilename);

    SET_LOGGER("MarleySimulator");
};

#endif
