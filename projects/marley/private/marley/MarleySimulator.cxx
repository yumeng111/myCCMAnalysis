// standard library stuff

#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>

#include "dataclasses/I3Double.h"
#include "dataclasses/I3Int32.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"
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
    }

void MarleySimulator::Configure() {
    GetParameter("MarleySearchPath", marley_search_path_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
    GetParameter("InputMCTreeName", input_mc_tree_name_);

    int random_seed;
    GetParameter("RandomSeed",random_seed);
    log_info("Using random seed for Sample Delay in gamma times = %d", random_seed);
    rng_ = boost::make_shared<I3GSLRandomService>(random_seed);

    GetParameter("EnableGammaTimeOffset", enable_gamma_time_offset_);
    log_info("EnableGammaTimeOffset = %s", enable_gamma_time_offset_ ? "True" : "False");

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

    //Loads and saves all the info of levels and transitions for K40
    this->LoadK40Transitions("/users/marisolc/workspaces/CCM2/sources/marley/data/structure/K.dat"); //need TO DO a generic path
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
    I3Particle neutrino;
    bool found_neutrino = false;

    // Iterate over all the particles in the tree
    for (auto iter = inputMCTree->begin(); iter != inputMCTree->end(); ++iter) {
        // If the particle does not have parent, is primary
        if(not inputMCTree->parent(*iter)) {
            neutrino = *iter;  // Save the primary particle
            found_neutrino = true;
            break;
        }
    }
    if(not found_neutrino) {
        log_warn("Neutrino not found in input tree");
        I3MCTreePtr outputMCTree = boost::make_shared<I3MCTree>(*inputMCTree);
        frame->Put(output_mc_tree_name_, outputMCTree);
        PushFrame(frame);
        return;
    }

    frames_in_total++; //Count how many frames have Marley events

    //Pass the parameters to marley_generator_.create_event()
    int pdg_a = neutrino.GetPdgEncoding();
    double KEa = neutrino.GetEnergy() * 1000.0;  // Energy in MeV (it was in GeV in the input Tree)
    int pdg_atom = 1000180400;  //  40Ar
    std::array<double, 3> dir_vec = { neutrino.GetDir().GetX(),
                                      neutrino.GetDir().GetY(),
                                      neutrino.GetDir().GetZ() };

    // Call the function create_event from marley generator
    marley::Event ev = marley_generator_.create_event(pdg_a, KEa, pdg_atom, dir_vec);

    //Create a new I3MCTree to add the particles from Marley event
    I3MCTreePtr outputMCTree = boost::make_shared<I3MCTree>(*inputMCTree);
    outputMCTree->erase_children(neutrino.GetID());

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
        p.SetEnergy(kinetic_energy*1.0e-3); //kinetic_energy is in MeV. The energies in the tree are in GeV

        //Set the position of the interaction using the initial neutrino position, direction and length
        p.SetPos(neutrino.GetPos() + (neutrino.GetLength() * neutrino.GetDir()));
        //if (p.pdg_code() == 11) { //electron
        //    p.SetPos(neutrino.GetPos() + (neutrino.GetLength() * neutrino.GetDir()));

        //Calculate the Direction using the momentum
        double px = fp->px();
        double py = fp->py();
        double pz = fp->pz();
        double magnitude = std::sqrt(px * px + py * py + pz * pz);

        double dirx = px / magnitude;
        double diry = py / magnitude;
        double dirz = pz /magnitude;

        p.SetDir(dirx, diry, dirz);

        //Calculate the time when the neutrino arrives. This will be the start time for p
        const double c = 3.0e8;
        p.SetTime((neutrino.GetTime() + neutrino.GetLength() /c)*1.0e9 ); //time in s. Time in tree in ns

        // Add the particle as daughter of the primary neutrino in the output tree
        outputMCTree->append_child(neutrino, p);
    }


    //I3MCTreePtr outputMCTree = boost::make_shared<I3MCTree>(); //output is an empty tree

    //Pass the outputMCTree to the function for gammas
    AdjustGammaTimes(outputMCTree, frame);
    // Now save to a frame
    frame->Put(output_mc_tree_name_, outputMCTree);
    PushFrame(frame);
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

    std::map<int, LevelInfo> levels_temp;

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
                log_info("Found start of the K40 block: Z=%d A=%d N=%d", Z, A, num_levels);
            }else if (inside_K40_block){
                inside_K40_block = false;
                log_info("Exiting block");
            }
            continue;
        }
        //Ignore everything outside the K40 block
        if(!inside_K40_block) continue;

        //Now match level line
        if (inside_K40_block && std::regex_match(line, match, level_regex)) {
            try{
                //Get Energy of the level, 2J, parity, and number of gammas (or transitions)
                double energy_MeV = std::stod(match[1]);
                double spin_times_two = std::stod(match[2]);
                std::string parity_sign = match[3];
                int num_transitions = std::stoi(match[4]);

                log_info("Found level: %.5f MeV with %d gammas", energy_MeV, num_transitions);

                //Then save the info:
                LevelInfo lvl;
                lvl.level_index = static_cast<int>(levels_temp.size());
                lvl.energy_keV = energy_MeV * 1e3;
                lvl.spin = spin_times_two / 2.0;
                lvl.parity = (parity_sign == "+") ? +1 : -1; //I don't like this logic (??)

                lvl.T12_ns = 0.0; //half life time
                lvl.tau_ns = 0.0; //mean life time

                current_level_index = lvl.level_index;

                levels_temp[current_level_index] = lvl;
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
    log_info("Loaded %zu K-40 levels.", levels_map_.size());
    // There are only 3 levels that will add a significant T1/2
    // All the other levels are of the order of pico or femto seconds
    // Manually set T1/2 where known from ENSDF

    for (auto& kv : levels_map_) {
        LevelInfo& lvl = kv.second;

        // Level n = 15 (2542.79 keV)
        if (std::abs(lvl.energy_keV - 2542.79) < 0.02) {
            lvl.T12_ns = 1.09; //from www.nndc.bnl.gov/
            lvl.tau_ns = lvl.T12_ns / std::log(2.0);
            log_info("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns) ",
                     lvl.level_index, lvl.energy_keV, lvl.T12_ns, lvl.tau_ns);
        }

        // Level n = 4 (1643.64 keV)
        if (std::abs(lvl.energy_keV - 1643.64) < 0.02) {
            lvl.T12_ns = 336.0; //from www.nndc.bnl.gov/
            lvl.tau_ns = lvl.T12_ns / std::log(2.0);
            log_info("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns)",
                     lvl.level_index, lvl.energy_keV, lvl.T12_ns, lvl.tau_ns);
        }

        // Level n = 1  (29.83 keV)
        if (std::abs(lvl.energy_keV - 29.8299) < 0.02) {
            lvl.T12_ns = 4.25; //from www.nndc.bnl.gov/
            lvl.tau_ns = lvl.T12_ns / std::log(2.0);
            log_info("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns)",
                     lvl.level_index, lvl.energy_keV, lvl.T12_ns, lvl.tau_ns);
        }
    }

    // Check that everything is ok for the first levels
    size_t N = std::min(levels_map_.size(), size_t(5));
    for (size_t i = 0; i < N; ++i) {
        const LevelInfo& lvl = levels_map_[i];
        std::string parity_init = (lvl.parity > 0) ? "+" : "-";

        log_info("Level [%d]: %.3f keV, spin %.1f%s, T1/2 = %.3f ns",
                 lvl.level_index,
                 lvl.energy_keV,
                 lvl.spin,
                 parity_init.c_str(),
                 lvl.T12_ns);

        for (const auto& t : lvl.transitions) {
            const LevelInfo& lvl_final = levels_map_[t.final_level_index];
            std::string parity_final = (lvl_final.parity > 0) ? "+" : "-";

            log_info("   Transition: [%d] %.3f keV J=%.1f%s → [%d] %.3f keV J=%.1f%s | Gamma = %.3f keV | BR = %.5f",
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
    for (const auto& kv : levels_map_) {
        total_transitions += kv.second.transitions.size();
    }

    log_info("Finished reading K40.dat");
    log_info("Total number of levels parsed: %zu", levels_map_.size());
    log_info("Total number of transitions parsed: %zu", total_transitions);


    //Finally we can write all the levels and transitions into a file for reference
    // Write all levels and transitions to a file for reference
    std::ofstream outfile("K40_levels_and_transitions.txt");
    if (!outfile.is_open()) {
        log_error("Could not open output file K40_levels_and_transitions.txt");
        return;
    }

    outfile << "# K40 Nuclear Levels and Transitions\n";
    outfile << "# =================================\n\n";

    outfile << "# ENERGY LEVELS (in keV)\n";
    for (const auto& kv : levels_map_) {
        const LevelInfo& lvl = kv.second;

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
    for (const auto& kv : levels_map_) {
        const LevelInfo& lvl = kv.second;

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
    log_info("K40 levels and transitions written to K40_levels_and_transitions.txt");
    log_info("Loaded %zu levels and %zu transitions.", levels_map_.size(), total_transitions);

}


// Look for events with K40 in the I3MCTree and save save pointers to the gamma particles
// Reconstruct the cascade for each event
// Save the cascade info in a key
// Sample a time delay for those de excitations coming from levels with a tau > ns
// Add new time delay in the I3MCTree
void MarleySimulator::AdjustGammaTimes(I3MCTreePtr mcTree, I3FramePtr frame) {
    log_info("AdjustGammaTimes called.");

    if (!enable_gamma_time_offset_) {
        log_info("Skipping time offsets for gammas because EnableGammaTimeOffset=False.");
        return;
    }

    bool has_k40 = false;
    bool has_metastable = false;

    //Let's look into the I3MCtree
    for (auto iter = mcTree->begin(); iter != mcTree->end(); ++iter) {

        const I3Particle& particle = *iter;

        //Checking everything is ok...
        /*log_info("Particle PDG: %d, E = %.5f GeV, time = %.5f ns",
                 particle.GetPdgEncoding(),
                 particle.GetEnergy(),
                 particle.GetTime());
        */

        //Then look for interactions with K40
        if (particle.GetPdgEncoding() == 1000190400) { // K40

            //Add a flag for K40 events
            has_k40 = true;
            //Look for gammas in the tree
            std::vector<I3Particle*> gamma_candidates;
            std::vector<double> gamma_energies_keV;

            for (auto gamma_iter = mcTree->begin(); gamma_iter != mcTree->end(); ++gamma_iter) {
                // Filter only gammas
                if (gamma_iter->GetPdgEncoding() == 22) {
                    gamma_candidates.push_back(&(*gamma_iter));
                    gamma_energies_keV.push_back(gamma_iter->GetEnergy() * 1e6); // GeV -> keV
                }
            }

            if (gamma_candidates.empty()) {
                log_info("No gamma particles found in this event.");
                continue;
            }

            //Calculate the sum of all the gammas
            //I don't use gamma_total_energy but it could be helpful
            double gamma_total_energy = std::accumulate(
                gamma_energies_keV.begin(),
                gamma_energies_keV.end(),
                0.0
            );

            log_info("Collected %zu gammas. Total energy sum = %.3f keV",
                    gamma_energies_keV.size(), gamma_total_energy);


            //Now we add the logic for the time delays
            //First we need to reconstruct the cascade

            log_info("=== Begin reconstructing gamma cascade ===");

            // Tolerance for matching energy sums from marley to real levels real levels from .dat
            const double match_tolerance_keV = 5.0;

            // Gammas from Marley in the I3MCTree are in order of production
            // We need to reverse order, to start from the last gamma produced
            std::reverse(gamma_candidates.begin(), gamma_candidates.end());
            std::reverse(gamma_energies_keV.begin(), gamma_energies_keV.end());

            // We start from the ground state (E=0)
            double running_energy = 0.0;
            double cumulative_time_ns = 0.0;

            // Vector for saving info of the cascade and then pass it to a key in the frame
            std::vector<NuclearCascadeStep> cascade_steps;

            for (size_t i = 0; i < gamma_candidates.size(); ++i) {
                double gamma_energy = gamma_energies_keV[i];
                I3Particle* g = gamma_candidates[i];
                // Calculate the energy of the level we’d reach by emitting this gamma
                // Note: initial = local higher level, running = local lower level
                double initial_energy_candidate = running_energy + gamma_energy;

                //Find the matching level in the nuclear levels map
                //Remember the map wad made from the .dat file

                double min_diff_initial = 1e9;
                int initial_lvl_idx = -1;

                for (const auto& kv : levels_map_) {
                    double diff = std::abs(kv.second.energy_keV - initial_energy_candidate);
                    if (diff < min_diff_initial) {
                        min_diff_initial = diff;
                        initial_lvl_idx = kv.first;
                    }
                }

                // Check if we matched the level directly
                if (min_diff_initial > match_tolerance_keV) {
                    // If there is no match within tolerance, force it to fall into the to nearest level
                    double nearest_diff = 1e9;
                    int nearest_idx = -1;
                    for (const auto& kv : levels_map_) {
                        double diff = std::abs(kv.second.energy_keV - initial_energy_candidate);
                        if (diff < nearest_diff) {
                            nearest_diff = diff;
                            nearest_idx = kv.first;
                        }
                    }

                    if (nearest_idx >= 0) {
                        const auto& lvl_nearest = levels_map_.at(nearest_idx);
                        log_warn("No direct match found for level %.3f keV. Assigning nearest level [index %d, %.3f keV]",
                                 initial_energy_candidate,
                                 nearest_idx,
                                 lvl_nearest.energy_keV);

                        // According to the energy of the gamma
                        // There will be only 3 levels that will add a significant delay
                        // n=15 E=2.54279 MeV  T1/2=1.09 ns
                        // n=4  E=1.64364 MeV  T1/2=336  ns <--------THIS IS THE META-STABLE LEVEL
                        // n=1  E=0.0298  MeV  T1/2=4.25 ns
                        // All the other levels have half life times of pico or femto seconds
                        //-----------------------------------------------------------------------
                        //Possible transitions from these levels:
                        //15-3-1-0
                        //15-3-0
                        //15-0
                        //4-2-1-0
                        //4-2-0
                        //4-1-0

                        // Check if it's one of the important delayed levels
                        if (std::abs(lvl_nearest.energy_keV - 29.8299) < 1.0 ||
                            std::abs(lvl_nearest.energy_keV - 1643.64) < 1.0 ||
                            std::abs(lvl_nearest.energy_keV - 2542.79) < 1.0) {
                            log_warn("!!! IMPORTANT: This level has significant lifetime and may impact delays !!!");
                        } else {
                            log_warn("This level has negligible time delay, nothing to worry about.");
                        }

                        initial_lvl_idx = nearest_idx;
                    } else {
                        log_warn("Could not even assign nearest level for %.3f keV.", initial_energy_candidate);
                        initial_lvl_idx = -1;
                    }
                } //end of 'if' for matching level

                // Now that the initial level index was set:
                // Assign and save the level info
                double initial_energy = initial_energy_candidate;
                double initial_T12_ns = 0.0;
                double initial_tau_ns = 0.0;
                double initial_spin = 0.0;
                int initial_parity = 0;

                if (initial_lvl_idx >= 0) {
                    const auto& lvl = levels_map_.at(initial_lvl_idx);
                    initial_energy = lvl.energy_keV; // We'll replace energy by nuclear data
                    initial_T12_ns = lvl.T12_ns;
                    initial_tau_ns = lvl.tau_ns;
                    initial_spin = lvl.spin;
                    initial_parity = lvl.parity;
                }

                // The “final” level is the running one we are leaving behind
                double final_energy = running_energy;
                int final_lvl_idx = -1;
                double final_spin = 0.0;
                int final_parity = 0;

                // Attempt to match running_energy to a level
                double min_diff_final = 1e9;
                for (const auto& kv : levels_map_) {
                    double diff = std::abs(kv.second.energy_keV - final_energy);
                    if (diff < min_diff_final) {
                        min_diff_final = diff;
                        final_lvl_idx = kv.first;
                    }
                }

                if (min_diff_final > match_tolerance_keV) {
                    final_lvl_idx = -1;
                } else {
                    const auto& lvl_final = levels_map_.at(final_lvl_idx);
                    final_energy = lvl_final.energy_keV; // replace by nuclear data
                    final_spin = lvl_final.spin;
                    final_parity = lvl_final.parity;
                }

                // Sample delay with SampleDelay function if level has lifetime > ns
                double delay_ns = 0.0;
                if (initial_tau_ns > 0.0) {
                    delay_ns = SampleDelay(initial_tau_ns);
                    cumulative_time_ns += delay_ns;
                }


                //And save the step in the cascade
                NuclearCascadeStep step;
                step.initial_level_index = initial_lvl_idx;
                step.initial_level_energy_keV = initial_energy;
                step.initial_level_spin = initial_spin;
                step.initial_level_parity = initial_parity;

                step.final_level_index = final_lvl_idx;
                step.final_level_energy_keV = final_energy;
                step.final_level_spin = final_spin;
                step.final_level_parity = final_parity;

                step.gamma_energy_keV = gamma_energy;
                step.T12_ns = initial_T12_ns;
                step.tau_ns = initial_tau_ns;
                step.sampled_delay_ns = delay_ns;
                step.cumulative_time_ns = cumulative_time_ns;

                cascade_steps.push_back(step);

                log_info("Cascade step: [%d] %.3f keV J=%.1f%s → [%d] %.3f keV J=%.1f%s | gamma = %.3f keV | Sampled delay = %.3f ns",
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

                // Update the runnning energy for the next step
                running_energy = initial_energy;
            } //end of 'for' for gammas
            //END OF LEVELS RECONSTRUCTION


            // Assign the new time to the gamma

            // If the initial nuclear level has a non-zero lifetime,
            // we sample a random delay to simulate how long it stays excited

            double initial_delay_ns = 0.0;

            if (!cascade_steps.empty()) {
                const auto& top_step = cascade_steps.back();

                if (top_step.tau_ns > 0.0) {
                    initial_delay_ns = SampleDelay(top_step.tau_ns);
                }
            }

            // Initialize cumulative time with the initial delay, if any
            // This represents how much time has passed since the creation of the nucleus
            // before the first gamma is emitted
            cumulative_time_ns = initial_delay_ns;

            // Loop through the cascade steps from highest to lowest nuclear levels
            // and accumulate delays for each gamma emission
            for (auto it = cascade_steps.rbegin(); it != cascade_steps.rend(); ++it) {
                auto& step = *it;

                // Add the sampled delay of THIS level to the cumulative time for the next step going down
                cumulative_time_ns += step.sampled_delay_ns;
                // Store the total cumulative time at which this gamma is emitted
                step.cumulative_time_ns = cumulative_time_ns;
            }

            // Re-invert to restore physical emission order
            // Originally, we reversed it to reconstruct levels from bottom-up
            // Now we restore it to the order of emission (from highest to lowest levels)
            std::reverse(gamma_candidates.begin(), gamma_candidates.end());

            // Assign the correct emission time to each gamma particle
            // This ensures that each gamma in the tree has the proper time offset
            // corresponding to delays in the nuclear cascade
            for (size_t i = 0; i < gamma_candidates.size(); ++i) {
                auto* g = gamma_candidates[i];
                // We need to access the cascade steps in reverse order,
                // because cascade_steps was filled from lowest to highest,
                // while gamma_candidates is now in physical emission order (highest to lowest)
                const auto& step = cascade_steps[cascade_steps.size() - 1 - i];

                double new_time = g->GetTime() + step.cumulative_time_ns;
                g->SetTime(new_time);

                log_info("Gamma energy %.3f keV assigned new time %.3f ns "
                         "(cumulative delay %.3f ns)",
                         g->GetEnergy() * 1e3,
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
            frames_with_cascade++;

            // Save the energy difference if it falls into the meta-stable
            double metastable_energy_diff = -1.0;
            double highest_level_energy = 0.0;

            for (const auto& step : cascade_steps) {

                // Check metastable (use index 4 directly)
                if (step.initial_level_index == 4 || step.final_level_index == 4) {
                    has_metastable = true;
                }

                // Save highest energy reached in the cascade
                if (step.initial_level_index >= 0) {
                    double E = levels_map_.at(step.initial_level_index).energy_keV;
                    if (E > highest_level_energy) {
                        highest_level_energy = E;
                    }
                }
            }

            if (has_metastable) {
                frames_with_metastable++;

                // Calculate difference between max E and metastable
                const double metastable_energy = levels_map_.at(4).energy_keV;
                metastable_energy_diff = highest_level_energy - metastable_energy;
                metastable_energy_diffs.push_back(metastable_energy_diff);

                log_info("Frame with metastable level. Energy diff from highest level = %.3f keV",
                         metastable_energy_diff);
            }
            // Saves in the frame if has metastable level n the cascade
            frame->Put("HasMetaStable", I3BoolPtr(new I3Bool(has_metastable)));
            // And the difference in Energy from max to meta stable
            frame->Put("MarleyMetastableDeltaE",
                       boost::make_shared<I3Double>(metastable_energy_diff));
            // ==End of counters=====


            frame->Put("MarleyGammaCascadeInfo", cascade_info);
            log_info("Stored MarleyGammaCascadeInfo in frame.");
            log_info("=== Finished reconstructing gamma cascade ===");

            //Save other important variables for analysis

            //Reverse cascade_steps so we save gammas in physical production order
            std::reverse(cascade_steps.begin(), cascade_steps.end());
            auto cumulative_times_vec = boost::make_shared<I3VectorDouble>();
            auto gamma_energies_vec = boost::make_shared<I3VectorDouble>();

            for (const auto& step : cascade_steps){
                cumulative_times_vec->push_back(step.cumulative_time_ns);
                gamma_energies_vec->push_back(step.gamma_energy_keV);
            }

            frame->Put("MarleyGammaCumulativeTimes", cumulative_times_vec);
            frame->Put("MarleyGammaEnergies", gamma_energies_vec);

        } //end of 'if' for K40

    }//end of main 'for' mcTree

    frame->Put("HasK40", I3BoolPtr(new I3Bool(has_k40)));

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
    return -mean_lifetime_ns * std::log(1.0 - u);
}


void MarleySimulator::FillSimulationFrame(I3FramePtr frame) {
    I3Int32Ptr obj = boost::make_shared<I3Int32>(1);
    frame->Put("MarleyConfiguration", obj);
}

void MarleySimulator::Finish(){
    log_info("====== MARLEY SIMULATOR SUMMARY ======");
    log_info("Frames with Marley Events: %zu", frames_in_total);
    log_info("Frames with MarleyGammaCascadeInfo: %zu", frames_with_cascade);
    log_info("Frames passing through metastable level [n=4] (1.64364 MeV): %zu", frames_with_metastable);

    //I want to check if the E of these prompt gamma is around 2.74 MeV as with the solar neutrinos
    if (!metastable_energy_diffs.empty()) {
        double sum = 0.0;
        for (auto x : metastable_energy_diffs) sum += x;
        double mean_diff = sum / metastable_energy_diffs.size();
        log_info("Average energy difference to metastable level: %.3f keV", mean_diff);
    }

    log_info("=======================================");
}
