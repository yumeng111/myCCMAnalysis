// standard library stuff

#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>

#include "dataclasses/I3Double.h"
#include "dataclasses/I3Int32.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/physics/I3Particle.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3Bool.h"

#include "phys-services/I3RandomService.h"
#include "phys-services/I3GSLRandomService.h"

#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/MarleySimulator.h"

#include <fstream>
#include <sstream>
#include <regex> //searches for patterns
#include <cmath> //math functions
#include <random> //for the seed
#include <exception>

I3_MODULE(MarleySimulator);

MarleySimulator::MarleySimulator(const I3Context& context) : I3ConditionalModule(context),
    output_mc_tree_name_("I3MCTree"), input_mc_tree_name_("SIRENMarleyInjectionTree") {
        std::random_device rd;
        unsigned int default_seed = rd() & 0x7FFFFFFF;

        AddParameter("MarleySearchPath", "Search path for Marley environment file", marley_search_path_);
        AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
        AddParameter("InputMCTreeName", "Name of the MCTree in the frame.", input_mc_tree_name_);
        AddParameter("RandomSeed", "Seed for the random number generator in SampleDelay", static_cast<int>(default_seed));
        AddParameter("EnableGammaTimeOffset",
                 "True: applies time offsets to gamma rays from nuclear cascade. "
                 "False: leaves gamma times untouched.",
                 true);
        AddParameter("SaveLevelsFile", "Name of file to which level information will be saved. An empty string (default)"
                                       "will prevent the file from being produced.", std::string(""));
    }

void MarleySimulator::Configure() {
    GetParameter("MarleySearchPath", marley_search_path_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
    GetParameter("InputMCTreeName", input_mc_tree_name_);

    int random_seed;
    GetParameter("RandomSeed",random_seed);
    log_debug("Using random seed for Sample Delay in gamma times = %d", random_seed);
    rng_ = boost::make_shared<I3GSLRandomService>(random_seed);

    GetParameter("EnableGammaTimeOffset", enable_gamma_time_offset_);
    log_debug("EnableGammaTimeOffset = %s", enable_gamma_time_offset_ ? "True" : "False");

    GetParameter("SaveLevelsFile", levels_filename_);
    save_levels_file_ = levels_filename_ != std::string("");

    setenv("MARLEY", "", 0);
    setenv("MARLEY_SEARCH_PATH", marley_search_path_.c_str(), 0);

    // String that contains the contents of the config file (from examples/config/annotated.js, we only will need the reaction, in this case ve40arCC, but there are more.)
    std::string marley_config = R"({reactions: ["ve40ArCC_Bhattacharya2009.react"]})";
    // Now convert the string into a JSON object of MARLEY (Using include/marley/JSON.hh)
    marley::JSON marley_json = marley::JSON::load(marley_config);
    // Pass the json object to the marley::JSONConfig constructor to create a config object (as in examples/marg4/src/MarleyPrimaryGeneratorAction.cc)
    marley::JSONConfig config(marley_json);
    // Call the config.create_generator function to get a marley::Generator object
    marley_generator_ = config.create_generator();

    //Loads all the info of levels and transitions for K40
    //This works for the version with MARLEY-SIREN-CCMAnalysis shared framework
    //The data files are in local/share/marley/structure
    const char* marley_path = std::getenv("MARLEY_SEARCH_PATH");
    if (!marley_path) {
        log_fatal("MARLEY_SEARCH_PATH is not set! Cannot load K.dat.");
    }

    std::string k40_file = std::string(marley_path) + "/structure/K.dat";
    log_debug("Loading K40 transitions from: %s", k40_file.c_str());
    this->LoadK40Transitions(k40_file);
}

void MarleySimulator::Simulation(I3FramePtr frame) {
    seen_s_frame_ = true;
    FillSimulationFrame(frame);
    PushFrame(frame);
}

void MarleySimulator::DAQ(I3FramePtr frame) {

    if(not seen_s_frame_) {
        I3FramePtr sim_frame = boost::make_shared<I3Frame>(I3Frame::Simulation);
        FillSimulationFrame(sim_frame);
        PushFrame(sim_frame);
        seen_s_frame_ = true;
    }

    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(input_mc_tree_name_);

    //The script crashes if the tree does not exist
    if (!inputMCTree) {
        log_warn("Input I3MCTree not found in frame!");
        PushFrame(frame);
        return;
    }

    // look for the primary particle in the input I3MCTree
    std::vector<I3Particle const *> primaries = I3MCTreeUtils::GetPrimariesPtr(inputMCTree);

    if(primaries.size() == 0 or (primaries.size() == 1 and not primaries[0]->IsNeutrino())) {
        log_warn("Neutrino not found in input tree");
        I3MCTreePtr outputMCTree = boost::make_shared<I3MCTree>(*inputMCTree);
        frame->Put(output_mc_tree_name_, outputMCTree);
        PushFrame(frame);
        return;
    } else if(primaries.size() > 1) {
        log_fatal("Found more than one primary in the tree: %s", input_mc_tree_name_.c_str());
    }

    I3Particle const * neutrino = primaries.at(0);

    ++frames_in_total; //Count how many frames have Marley events

    //Pass the parameters to marley_generator_.create_event()
    int pdg_a = neutrino->GetPdgEncoding();
    double KEa = neutrino->GetEnergy() / I3Units::MeV;  // Energy in MeV (it was in GeV in the input Tree)
    int pdg_atom = I3Particle::ParticleType::K40Nucleus;
    std::array<double, 3> dir_vec = { neutrino->GetDir().GetX(),
                                      neutrino->GetDir().GetY(),
                                      neutrino->GetDir().GetZ() };

    // Call the function create_event from marley generator
    marley::Event ev = marley_generator_.create_event(pdg_a, KEa, pdg_atom, dir_vec);

    //Create a new I3MCTree to add the particles from Marley event
    I3MCTreePtr outputMCTree = boost::make_shared<I3MCTree>(*inputMCTree);
    outputMCTree->erase_children(neutrino->GetID());

    //Put particles from Marley into output tree
    //
    // Loop over each of the final particles (fp) in the MARLEY event
    // and save their properties into the output I3MCTree
    for ( const auto& fp : ev.get_final_particles() ) {
        I3Particle p;
        p.SetType(static_cast<I3Particle::ParticleType>(fp->pdg_code()));;
        double total_energy = fp->total_energy();
        double mass = fp->mass();
        double kinetic_energy = total_energy - mass;
        p.SetEnergy(kinetic_energy * I3Units::MeV); //kinetic_energy is in MeV. The energies in the tree are in GeV

        //Set the position of the interaction using the initial neutrino position, direction and length
        p.SetPos(neutrino->GetStopPos());

        //Calculate the Direction using the momentum
        double px = fp->px();
        double py = fp->py();
        double pz = fp->pz();
        double magnitude = std::sqrt(px * px + py * py + pz * pz);

        double dirx = px / magnitude;
        double diry = py / magnitude;
        double dirz = pz / magnitude;

        p.SetDir(dirx, diry, dirz);

        p.SetTime(neutrino->GetStopTime());

        // Add the particle as daughter of the primary neutrino in the output tree
        outputMCTree->append_child(*neutrino, p);
    }

    //Pass the outputMCTree to the function for gammas
    AdjustGammaTimes(outputMCTree, frame);
    // Now save to a frame
    frame->Put(output_mc_tree_name_, outputMCTree);
    PushFrame(frame);
}

LevelInfo & MarleySimulator::ClosestLevel(std::vector<LevelInfo> & levels, double energy) {
    // high points to the first element that is *not* ordered before the energy
    // i.e. high >= energy
    std::vector<LevelInfo>::iterator high = std::lower_bound(levels.begin(), levels.end(), energy);
    if(high == levels.begin())
        return levels.front();
    if(high == levels.end())
        return levels.back();
    std::vector<LevelInfo>::iterator low = high-1;
    double diff_high = std::abs(high->energy_keV / I3Units::keV - energy);
    double diff_low = std::abs(energy - low->energy_keV / I3Units::keV);
    if(diff_low < diff_high)
        return *low;
    else
        return *high;
}

//This method gets the transitions from the K40.dat file
//That file is an extract of the marley/data/structure/K.dat only for the K40
//and is located in the same path as the main python script
//In this method we save all the info of the levels in levels_map_
void MarleySimulator::LoadK40Transitions(const std::string& filename) {
    std::ifstream infile(filename);
    std::string line;

    //Excitation energy of level, 2*Spin, Parity, Total gammas (transitions)
    //  Raw string()
    //  ^ Start of the line
    //  \s+ One or more blank spaces
    // ([\d\.Ee+-]+) Digit, point, scientific notation E, sign
    //  \d+ One or more digits
    // [-+]? Optional sign
    // [+-] Mandatory sign (Parity)
    std::regex level_regex(R"(^\s*([\d\.Ee+-]+)\s+([-+]?\d+)\s+([+-])\s+(\d+))");

    //Gamma energy, Relative intensity (Branching ratio), Final level (to which it descends).
    std::regex gamma_regex(R"(^\s*([\d\.Ee+-]+)\s+([\d\.Ee+-]+)\s+(\d+))");

    std::vector<LevelInfo> levels_temp;

    int current_level_index = -1;

    //We need to look for the block of K40 in the K.dat file
    bool inside_K40_block = false;

    if (!infile.is_open()) {
        log_error("Could not open file %s", filename.c_str());
        return;
    }

    //Start reading the file line by line
    while (std::getline(infile, line)) {
        std::smatch match;

        //Match header lines (Z A num_levels)
        std::istringstream line_stream(line);
        int Z, A, num_levels;

        if (line_stream >> Z >> A >> num_levels){
            if(Z == 19 && A == 40) {
                inside_K40_block = true;
                log_debug("Found start of the K40 block: Z=%d A=%d N=%d", Z, A, num_levels);
            }else if (inside_K40_block){
                inside_K40_block = false;
                log_debug("Exiting block");
            }
            continue;
        }
        //Ignore everything outside the K40 block
        if(!inside_K40_block) continue;

        //Now match level line
        if(std::regex_match(line, match, level_regex)) {
            try{
                //Get Energy of the level, 2J, parity, and number of gammas (or transitions)
                double energy_MeV = std::stod(match[1]);
                double spin_times_two = std::stod(match[2]);
                std::string parity_sign = match[3];
                int num_transitions = std::stoi(match[4]);

                log_debug("Found level: %.5f MeV with %d gammas", energy_MeV, num_transitions);

                //Then save the info:
                LevelInfo lvl;
                lvl.level_index = static_cast<int>(levels_temp.size());
                lvl.energy_keV = energy_MeV * 1e3;
                lvl.spin = spin_times_two / 2.0;
                lvl.parity = (parity_sign == "+") ? +1 : -1;

                lvl.T12_ns = 0.0; //half life time
                lvl.tau_ns = 0.0; //mean life time

                current_level_index = lvl.level_index;

                levels_temp.push_back(lvl);
            } catch (const std::exception& e){
                log_warn("Failed to parse level line: '%s'. Exception: %s", line.c_str(), e.what());
                continue;
            }
        }

        //Now we look at the lines for gammas
        else if (inside_K40_block && std::regex_match(line, match, gamma_regex) && current_level_index >= 0) {
            try{
                double E_gamma_MeV = std::stod(match[1]);  //Energy of the gamma
                double branching_ratio = std::stod(match[2]); //Relative intensity
                int final_level_idx = std::stoi(match[3]);  //index of the level of de-excitation (to which it descends)

                LevelInfo::Transition t;
                t.gamma_energy_keV = E_gamma_MeV * 1e3;
                t.branching_ratio = branching_ratio;
                t.final_level_index = final_level_idx;

                levels_temp[current_level_index].transitions.push_back(t);
            } catch (const std::exception& e){
                log_warn("Failed to parse gamma line: '%s'. Error: %s", line.c_str(), e.what());
                continue;
            }
        }
    }

    // Copy levels_temp into levels_map_
    levels_map_ = levels_temp;
    log_debug("Loaded %zu K-40 levels.", levels_map_.size());
    // There are only 3 levels that will add a significant T1/2
    // All the other levels are of the order of pico or femto seconds
    // Manually set T1/2 where known from ENSDF

    // Level n = 15 (2542.79 keV)
    LevelInfo & level_15 = ClosestLevel(levels_map_, 2542.79 * I3Units::keV);
    level_15.T12_ns = 1.09; //from www.nndc.bnl.gov/
    level_15.tau_ns = level_15.T12_ns / std::log(2.0);
    log_debug("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns) ",
             level_15.level_index, level_15.energy_keV, level_15.T12_ns, level_15.tau_ns);

    // Level n = 4 (1643.64 keV)
    LevelInfo & level_4 = ClosestLevel(levels_map_, 1643.64 * I3Units::keV);
    level_4.T12_ns = 336.0; //from www.nndc.bnl.gov/
    level_4.tau_ns = level_4.T12_ns / std::log(2.0);
    log_debug("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns)",
             level_4.level_index, level_4.energy_keV, level_4.T12_ns, level_4.tau_ns);

    // Level n = 1  (29.83 keV)
    LevelInfo & level_1 = ClosestLevel(levels_map_, 29.8299 * I3Units::keV);
    level_1.T12_ns = 4.25; //from www.nndc.bnl.gov/
    level_1.tau_ns = level_1.T12_ns / std::log(2.0);
    log_debug("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns)",
             level_1.level_index, level_1.energy_keV, level_1.T12_ns, level_1.tau_ns);

    // Check that everything is ok for the first levels
    size_t N = std::min(levels_map_.size(), size_t(5));
    for (size_t i = 0; i < N; ++i) {
        const LevelInfo& lvl = levels_map_[i];
        std::string parity_init = (lvl.parity > 0) ? "+" : "-";

        log_debug("Level [%d]: %.3f keV, spin %.1f%s, T1/2 = %.3f ns",
                 lvl.level_index,
                 lvl.energy_keV,
                 lvl.spin,
                 parity_init.c_str(),
                 lvl.T12_ns);

        for (const auto& t : lvl.transitions) {
            const LevelInfo& lvl_final = levels_map_[t.final_level_index];
            std::string parity_final = (lvl_final.parity > 0) ? "+" : "-";

            log_debug("   Transition: [%d] %.3f keV J=%.1f%s → [%d] %.3f keV J=%.1f%s | Gamma = %.3f keV | BR = %.5f",
                     lvl.level_index,
                     lvl.energy_keV,
                     lvl.spin,
                     parity_init.c_str(),
                     t.final_level_index,
                     lvl_final.energy_keV,
                     lvl_final.spin,
                     parity_final.c_str(),
                     t.gamma_energy_keV,
                     t.branching_ratio);
        }
    }


    // Print totals
    size_t total_transitions = 0;
    for (LevelInfo const & level : levels_map_) {
        total_transitions += level.transitions.size();
    }

    log_debug("Finished reading K40.dat");
    log_debug("Total number of levels parsed: %zu", levels_map_.size());
    log_debug("Total number of transitions parsed: %zu", total_transitions);


    //Finally we can write all the levels and transitions into a file for reference
    // Write all levels and transitions to a file for reference
    if(save_levels_file_) {
        std::ofstream outfile(levels_filename_);
        if (!outfile.is_open()) {
            log_error("Could not open output file %s", levels_filename_.c_str());
            return;
        }

        outfile << "# K40 Nuclear Levels and Transitions\n";
        outfile << "# =================================\n\n";

        outfile << "# ENERGY LEVELS (in keV)\n";
        for (LevelInfo const & lvl : levels_map_) {
            outfile << "Level[" << lvl.level_index << "] "
                    << lvl.energy_keV << " keV, "
                    << "J=" << lvl.spin
                    << ((lvl.parity > 0) ? "+" : "-");

            if (lvl.T12_ns > 0.0) {
                outfile << ", T1/2=" << lvl.T12_ns << " ns";
            }

            outfile << "\n";
        }

        outfile << "\n# TRANSITIONS\n";
        for (LevelInfo const & lvl : levels_map_) {
            for (const auto& t : lvl.transitions) {
                const LevelInfo& lvl_final = levels_map_.at(t.final_level_index);
                outfile << "Transition: "
                        << "[" << lvl.level_index << "] "
                        << lvl.energy_keV << " keV → "
                        << "[" << lvl_final.level_index << "] "
                        << lvl_final.energy_keV << " keV, "
                        << "Gamma = " << t.gamma_energy_keV << " keV, "
                        << "BR = " << t.branching_ratio
                        << "\n";
            }
        }

        outfile.close();
        log_debug("K40 levels and transitions written to %s", levels_filename_.c_str());
        log_debug("Loaded %zu levels and %zu transitions.", levels_map_.size(), total_transitions);
    }

}


// Look for events with K40 in the I3MCTree and save save pointers to the gamma particles
// Reconstruct the cascade for each event
// Save the cascade info in a key
// Sample a time delay for those de excitations coming from levels with a tau > ns
// Add new time delay in the I3MCTree
void MarleySimulator::AdjustGammaTimes(I3MCTreePtr mcTree, I3FramePtr frame) {
    log_debug("AdjustGammaTimes called.");

    if (!enable_gamma_time_offset_) {
        log_debug("Skipping time offsets for gammas because EnableGammaTimeOffset=False.");
        return;
    }

    std::function<bool(I3Particle const &)> is_K40 = [](I3Particle const & p) -> bool {
        return p.GetType() == I3Particle::ParticleType::K40Nucleus;
    };

    std::function<bool(I3Particle const &)> is_gamma = [](I3Particle const & p) -> bool {
        return p.GetType() == I3Particle::ParticleType::Gamma;
    };

    std::vector<typename I3MCTree::fast_const_iterator> const K40_vec = I3MCTreeUtils::GetFilterPtr(mcTree, is_K40);
    bool has_K40 = K40_vec.size() > 0;
    if(K40_vec.size() > 1)
        log_fatal("Found more than one K40 in the MARLEY tree");

    if(not has_K40) {
        frame->Put("HasK40", I3BoolPtr(new I3Bool(has_K40)));
        PushFrame(frame);
        return;
    }

    I3Particle const & K40 = *K40_vec.front();
    I3Particle const & parent = *I3MCTreeUtils::GetParentPtr(mcTree, K40);
    std::vector<I3Particle*> const siblings = I3MCTreeUtils::GetDaughtersPtr(mcTree, parent);

    std::function<bool(I3Particle const *)> is_gamma_ptr = [=](I3Particle const * p) -> bool {
        return p->GetType() == I3Particle::ParticleType::Gamma;
    };

    std::vector<I3Particle*> gamma_rays; gamma_rays.reserve(siblings.size() - 1);
    std::copy_if(siblings.begin(), siblings.end(), std::back_inserter(gamma_rays), is_gamma_ptr);

    //Now we add the logic for the time delays
    //First we need to reconstruct the cascade

    log_debug("=== Begin reconstructing gamma cascade ===");

    // We start from the ground state (E=0)
    double running_energy = 0.0;

    // Vector for saving info of the cascade and then pass it to a key in the frame
    std::vector<NuclearCascadeStep> cascade_steps(gamma_rays.size());

    // Iterate in reverse order since we're working our way up from the ground state
    for(int i = gamma_rays.size()-1; i >= 0; --i) {
        I3Particle* g = gamma_rays.at(i);
        double gamma_energy = g->GetEnergy();

        //Find the matching level in the nuclear levels map
        //Remember the map wad made from the .dat file
        LevelInfo & final_level = ClosestLevel(levels_map_, running_energy);
        running_energy += gamma_energy;
        LevelInfo & initial_level = ClosestLevel(levels_map_, running_energy);

        // Update the runnning energy for the next step
        running_energy = initial_level.energy_keV * I3Units::keV;

        // Sample delay with SampleDelay function if level has lifetime > ns
        double delay_ns = 0.0;
        if(initial_level.tau_ns > 0.0) {
            delay_ns = SampleDelay(initial_level.tau_ns);
        }

        //And save the step in the cascade
        NuclearCascadeStep step;
        step.initial_level_index = initial_level.level_index;
        step.initial_level_energy_keV = initial_level.energy_keV;
        step.initial_level_spin = initial_level.spin;
        step.initial_level_parity = initial_level.parity;

        step.final_level_index = final_level.level_index;
        step.final_level_energy_keV = final_level.energy_keV;
        step.final_level_spin = final_level.spin;
        step.final_level_parity = final_level.parity;

        step.gamma_energy_keV = gamma_energy;
        step.T12_ns = initial_level.T12_ns;
        step.tau_ns = initial_level.tau_ns;
        step.sampled_delay_ns = delay_ns;

        cascade_steps.at(i) = step;

        log_debug("Cascade step: [%d] %.3f keV J=%.1f%s → [%d] %.3f keV J=%.1f%s | gamma = %.3f keV | Sampled delay = %.3f ns",
                step.initial_level_index,
                step.initial_level_energy_keV,
                step.initial_level_spin,
                (step.initial_level_parity > 0) ? "+" : "-",
                step.final_level_index,
                step.final_level_energy_keV,
                step.final_level_spin,
                (step.final_level_parity > 0) ? "+" : "-",
                step.gamma_energy_keV,
                step.sampled_delay_ns
                );

    } //end of 'for' for gammas
    //END OF LEVELS RECONSTRUCTION

    // Assign the new time to the gamma
    double cumulative_time_ns = 0.0;
    for(size_t i=0; i<gamma_rays.size(); ++i) {
        I3Particle* g = gamma_rays.at(i);
        NuclearCascadeStep & step = cascade_steps.at(i);
        cumulative_time_ns += step.sampled_delay_ns;
        step.cumulative_time_ns = cumulative_time_ns;
        double new_time = g->GetTime() + cumulative_time_ns;
        g->SetTime(new_time);
        log_debug("Gamma energy %.3f keV assigned new time %.3f ns "
                "(cumulative delay %.3f ns)",
                g->GetEnergy() / I3Units::keV,
                new_time,
                step.cumulative_time_ns);
    }

    //Save the cascade into the frame
    I3VectorStringPtr cascade_info = boost::make_shared<I3VectorString>();

    for (const auto& step : cascade_steps) {
        std::ostringstream ss;

        ss << "Step:\n";
        ss << "  From Level [" << (step.initial_level_index >= 0 ? std::to_string(step.initial_level_index) : "UNKNOWN") << "] "
            << step.initial_level_energy_keV << " keV, "
            << "J=" << step.initial_level_spin
            << (step.initial_level_parity > 0 ? "+" : "-") << "\n";

        ss << "  To Level [" << (step.final_level_index >= 0 ? std::to_string(step.final_level_index) : "UNKNOWN") << "] "
            << step.final_level_energy_keV << " keV, "
            << "J=" << step.final_level_spin
            << (step.final_level_parity > 0 ? "+" : "-") << "\n";

        ss << "  Gamma Energy = " << step.gamma_energy_keV << " keV\n";
        ss << "  T1/2 = " << step.T12_ns << " ns\n";
        ss << "  Tau = " << step.tau_ns << " ns\n";
        ss << "  Sampled Delay = " << step.sampled_delay_ns << " ns\n";
        ss << "  Cumulative time = " << step.cumulative_time_ns << " ns\n";

        cascade_info->push_back(ss.str());
    }

    // ===Some counters ===
    ++frames_with_cascade;

    // Save the energy difference if it falls into the meta-stable
    double highest_level_energy = cascade_steps.front().initial_level_energy_keV;
    double metastable_energy_diff = -1.0;

    bool has_metastable = false;
    for(const auto& step : cascade_steps) {
        // Check metastable (use index 4 directly)
        if (step.initial_level_index == 4) {
            has_metastable = true;
            break;
        }
    }

    if (has_metastable) {
        ++frames_with_metastable;

        // Calculate difference between max E and metastable
        const double metastable_energy = levels_map_.at(4).energy_keV;
        metastable_energy_diff = highest_level_energy - metastable_energy;
        metastable_energy_diffs.push_back(metastable_energy_diff);

        log_debug("Frame with metastable level. Energy diff from highest level = %.3f keV",
                metastable_energy_diff);
    }
    // Saves in the frame if has metastable level n the cascade
    frame->Put("HasMetaStable", I3BoolPtr(new I3Bool(has_metastable)));
    // And the difference in Energy from max to meta stable
    frame->Put("MarleyMetastableDeltaE",
            boost::make_shared<I3Double>(metastable_energy_diff));
    // ==End of counters=====


    frame->Put("MarleyGammaCascadeInfo", cascade_info);
    log_debug("Stored MarleyGammaCascadeInfo in frame.");
    log_debug("=== Finished reconstructing gamma cascade ===");

    //Save other important variables for analysis

    auto cumulative_times_vec = boost::make_shared<I3VectorDouble>();
    auto gamma_energies_vec = boost::make_shared<I3VectorDouble>();

    for(const auto& step : cascade_steps) {
        cumulative_times_vec->push_back(step.cumulative_time_ns);
        gamma_energies_vec->push_back(step.gamma_energy_keV);
    }

    frame->Put("MarleyGammaCumulativeTimes", cumulative_times_vec);
    frame->Put("MarleyGammaEnergies", gamma_energies_vec);

    frame->Put("HasK40", I3BoolPtr(new I3Bool(has_K40)));

} //End of AdjustGammaTimes

//This method samples the delay for a given lifetime
//using an exponential distribution.
//t= -tau * ln (1-u)
//tau in ns
//u= random number between [0,1]
double MarleySimulator::SampleDelay(double mean_lifetime_ns) {
    if (mean_lifetime_ns <= 0.0)
        return 0.0;

    double u = rng_->Uniform(0.0, 1.0);
    return -mean_lifetime_ns * std::log1p(-u);
}


void MarleySimulator::FillSimulationFrame(I3FramePtr frame) {
    I3Int32Ptr obj = boost::make_shared<I3Int32>(1);
    frame->Put("MarleyConfiguration", obj);
}

void MarleySimulator::Finish(){
    log_debug("====== MARLEY SIMULATOR SUMMARY ======");
    log_debug("Frames with Marley Events: %zu", frames_in_total);
    log_debug("Frames with MarleyGammaCascadeInfo: %zu", frames_with_cascade);
    log_debug("Frames passing through metastable level [n=4] (1.64364 MeV): %zu", frames_with_metastable);

    //I want to check if the E of these prompt gamma is around 2.74 MeV as with the solar neutrinos
    if (!metastable_energy_diffs.empty()) {
        double sum = 0.0;
        for (auto x : metastable_energy_diffs) sum += x;
        double mean_diff = sum / metastable_energy_diffs.size();
        log_debug("Average energy difference to metastable level: %.3f keV", mean_diff);
    }

    log_debug("=======================================");
}
