// standard library stuff

#include <cmath> // math functions
#include <regex> // searches for patterns
#include <cctype>
#include <cstdio>
#include <limits>
#include <random> // for the seed
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <filesystem>

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


I3_MODULE(MarleySimulator);

void MarleySimulator::GetNucleonContent(int code, int & strange_count, int & neutron_count, int & proton_count, int & nucleon_count) {
    int prefix = 0;
    int excitation = 0;

    char buf[CHAR_BUF_SIZE];
    snprintf(buf, CHAR_BUF_SIZE, "%d", code);
    int nread = sscanf(buf, "%2d%1d%3d%3d%1d", &prefix, &strange_count, &proton_count, &nucleon_count, &excitation);
    if (nread != 5) {
        throw std::runtime_error("Failed to convert nuclear pdg to 10LZZZAAAI "
                "prefix "+std::to_string(prefix)+", L "+std::to_string(strange_count)+", Z "+std::to_string(proton_count)+", A "+std::to_string(nucleon_count)+", I "+std::to_string(excitation));
    }
    neutron_count = nucleon_count - proton_count - strange_count;
}

std::map<I3Particle::ParticleType, std::vector<std::tuple<double, double>>> const MarleySimulator::delayed_levels_ = {
    {
        I3Particle::ParticleType::K40Nucleus,
        {
            {1.09  * I3Units::ns, 2542.79 * I3Units::keV},
            {336.0 * I3Units::ns, 1643.64 * I3Units::keV},
            {4.25  * I3Units::ns, 29.8299 * I3Units::keV},
        }
    }
};

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

        AddParameter("Nuclei", "Nuclei with delayed gamma rays", I3Vector<I3Particle::ParticleType>({I3Particle::ParticleType::K40Nucleus}));

        AddParameter("Debug", "Enable debugging outputs", bool(false));
}

void MarleySimulator::Configure() {
    GetParameter("MarleySearchPath", marley_search_path_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
    GetParameter("InputMCTreeName", input_mc_tree_name_);

    int random_seed;
    GetParameter("RandomSeed",random_seed);
    log_debug("Using random seed for Sample Delay in gamma times = %d", random_seed);
    rng_ = boost::make_shared<I3GSLRandomService>(random_seed);
    unsigned int marley_seed = (unsigned int)rng_->Uniform(std::numeric_limits<unsigned int>::max());

    GetParameter("EnableGammaTimeOffset", enable_gamma_time_offset_);
    log_debug("EnableGammaTimeOffset = %s", enable_gamma_time_offset_ ? "True" : "False");

    GetParameter("SaveLevelsFile", levels_filename_);
    save_levels_file_ = levels_filename_ != std::string("");

    I3Vector<I3Particle::ParticleType> nuclei;
    GetParameter("Nuclei", nuclei);

    GetParameter("Debug", debug_);

    if(marley_search_path_.empty()) {
        // Check the environment variable for a value if no argument is provided
        const char* marley_path = std::getenv("MARLEY_SEARCH_PATH");
        if(!marley_path) {
            log_fatal("MARLEY_SEARCH_PATH is not set and marley_search_path_ is empty! Cannot load K.dat.");
        }
        marley_search_path_ = marley_path;
    } else {
        // If an argument is provided then store the value for potential downstream use
        setenv("MARLEY_SEARCH_PATH", marley_search_path_.c_str(), 0);
    }

    if(FileManager_::get_search_path().empty()) {
        // An empty search path indicated that MARLEY has not been initialized

        // Set the MARLEY environment variable to bypass internal MARLEY check
        setenv("MARLEY", "", 0);

        // Set up MARLEY and override the search path
        marley::FileManager::Instance();
        FileManager_::set_search_path(marley_search_path_);
    }

    // String that contains the contents of the config file (from examples/config/annotated.js, we only will need the reaction, in this case ve40arCC, but there are more.)
    std::stringstream config_ss;
    config_ss << R"({reactions: ["ve40ArCC_Bhattacharya2009.react"], seed: )";
    config_ss << marley_seed;
    config_ss << R"(})";
    std::string marley_config = config_ss.str();
    // Now convert the string into a JSON object of MARLEY (Using include/marley/JSON.hh)
    marley::JSON marley_json = marley::JSON::load(marley_config);
    // Pass the json object to the marley::JSONConfig constructor to create a config object (as in examples/marg4/src/MarleyPrimaryGeneratorAction.cc)
    marley::JSONConfig config(marley_json);
    // Call the config.create_generator function to get a marley::Generator object
    marley_generator_ = config.create_generator();

    for(std::pair<I3Particle::ParticleType const, std::vector<std::tuple<double, double>>> const & p : delayed_levels_) {
        I3Particle::ParticleType type = p.first;
        std::string type_string = i3particle_type_string(static_cast<int32_t>(type));
        std::string element_name = type_string.substr(0, type_string.size() - std::string("Nucleus").size());
        std::string fname = element_name;
        fname.erase(std::remove_if(fname.begin(), fname.end(), [](char c)->bool{return std::isdigit(c);}), fname.end());
        fname += ".dat";

        // Loads all the info of levels and transitions for K40
        // This works for the version with MARLEY-SIREN-CCMAnalysis shared framework
        // The data files are in local/share/marley/structure
        std::string element_file = marley::FileManager::Instance().find_file(fname, marley_search_path_);
        log_debug("Loading %s transitions from: %s", element_name.c_str(), element_file.c_str());
        if(not std::filesystem::exists(element_file)) {
            log_fatal("%s file \"%s\" does not exist!", element_name.c_str(), element_file.c_str());
        }
        this->LoadK40Transitions(type, element_file);
    }
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

    // The script crashes if the tree does not exist
    if(!inputMCTree) {
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

    if(debug_)
        ++frames_in_total; // Count how many frames have Marley events

    // Pass the parameters to marley_generator_.create_event()
    int pdg_a = neutrino->GetPdgEncoding();
    double KEa = neutrino->GetEnergy() / I3Units::MeV;  // Energy in MeV (it was in GeV in the input Tree)
    int pdg_atom = I3Particle::ParticleType::Ar40Nucleus;
    std::array<double, 3> dir_vec = { neutrino->GetDir().GetX(),
                                      neutrino->GetDir().GetY(),
                                      neutrino->GetDir().GetZ() };

    // Call the function create_event from marley generator
    marley::Event ev = marley_generator_.create_event(pdg_a, KEa, pdg_atom, dir_vec);

    // Create a new I3MCTree to add the particles from Marley event
    I3MCTreePtr outputMCTree = boost::make_shared<I3MCTree>(*inputMCTree);
    outputMCTree->erase_children(neutrino->GetID());

    // Put particles from Marley into output tree
    //
    // Loop over each of the final particles (fp) in the MARLEY event
    // and save their properties into the output I3MCTree
    for(const auto& fp : ev.get_final_particles()) {
        I3Particle p;
        p.SetType(static_cast<I3Particle::ParticleType>(fp->pdg_code()));;
        double total_energy = fp->total_energy();
        double mass = fp->mass();
        double kinetic_energy = total_energy - mass;
        p.SetEnergy(kinetic_energy * I3Units::MeV); // kinetic_energy is in MeV. The energies in the tree are in GeV

        // Set the position of the interaction using the initial neutrino position, direction and length
        p.SetPos(neutrino->GetStopPos());

        // Calculate the Direction using the momentum
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

    // Pass the outputMCTree to the function for gammas
    AdjustGammaTimes(outputMCTree, frame);
    // Now save to a frame
    frame->Put(output_mc_tree_name_, outputMCTree);
    PushFrame(frame);
}

std::vector<LevelInfo>::iterator MarleySimulator::ClosestLevel(std::vector<LevelInfo> & levels, double energy) {
    // high points to the first element that is *not* ordered before the energy
    // i.e. high >= energy
    std::vector<LevelInfo>::iterator high = std::lower_bound(levels.begin(), levels.end(), energy);
    if(high == levels.begin()) {
        log_debug("Hit begin");
        log_debug("Energy: %f keV, Level: %d", energy / I3Units::keV, (int)std::distance(levels.begin(), levels.begin()));
        return levels.begin();
    }
    if(high == levels.end()) {
        log_debug("Hit end");
        auto it = levels.begin() + (std::max(levels.size(), size_t(1)) - 1);
        log_debug("Energy: %f keV, Level: %d", energy / I3Units::keV, (int)std::distance(levels.begin(), it));
        return it;
    }
    std::vector<LevelInfo>::iterator low = high-1;
    double diff_high = std::abs(high->energy - energy);
    double diff_low = std::abs(energy - low->energy);
    if(diff_low < diff_high) {
        log_debug("Using low");
        log_debug("Energy: %f keV, Level: %d", energy / I3Units::keV, (int)std::distance(levels.begin(), low));
        return low;
    } else {
        log_debug("Using high");
        log_debug("Energy: %f keV, Level: %d", energy / I3Units::keV, (int)std::distance(levels.begin(), high));
        return high;
    }
}

// This method gets the transitions from the K40.dat file
// That file is an extract of the marley/data/structure/K.dat only for the K40
// and is located in the same path as the main python script
// In this method we save all the info of the levels in levels_map_
void MarleySimulator::LoadK40Transitions(I3Particle::ParticleType type, std::string const & filename) {
    int code = static_cast<int32_t>(type);
    int strange_count;
    int neutron_count;
    int proton_count;
    int nucleon_count;
    GetNucleonContent(code, strange_count, neutron_count, proton_count, nucleon_count);

    std::string type_string = i3particle_type_string(static_cast<int32_t>(type));
    std::string element_name = type_string.substr(0, type_string.size() - std::string("Nucleus").size());

    std::ifstream infile(filename);
    std::string line;

    // Excitation energy of level, 2*Spin, Parity, Total gammas (transitions)
    //  Raw string()
    //  ^ Start of the line
    //  \s+ One or more blank spaces
    // ([\d\.Ee+-]+) Digit, point, scientific notation E, sign
    //  \d+ One or more digits
    // [-+]? Optional sign
    // [+-] Mandatory sign (Parity)
    std::regex level_regex(R"(^\s*([\d\.Ee+-]+)\s+([-+]?\d+)\s+([+-])\s+(\d+))");

    // Gamma energy, Relative intensity (Branching ratio), Final level (to which it descends).
    std::regex gamma_regex(R"(^\s*([\d\.Ee+-]+)\s+([\d\.Ee+-]+)\s+(\d+))");

    std::vector<LevelInfo> levels_temp;

    int current_level_index = -1;

    // We need to look for the block of K40 in the K.dat file
    bool inside_K40_block = false;

    if(!infile.is_open()) {
        log_error("Could not open file %s", filename.c_str());
        return;
    }

    // Start reading the file line by line
    while (std::getline(infile, line)) {
        std::smatch match;

        // Match header lines (Z A num_levels)
        std::istringstream line_stream(line);
        int Z, A, num_levels;

        if(line_stream >> Z >> A >> num_levels) {
            if(Z == proton_count && A == nucleon_count) {
                inside_K40_block = true;
                log_debug("Found start of the K40 block: Z=%d A=%d N=%d", Z, A, num_levels);
            } else if(inside_K40_block) {
                inside_K40_block = false;
                log_debug("Exiting block");
            }
            continue;
        }
        // Ignore everything outside the K40 block
        if(!inside_K40_block) continue;

        // Now match level line
        if(std::regex_match(line, match, level_regex)) {
            try{
                // Get Energy of the level, 2J, parity, and number of gammas (or transitions)
                double energy_MeV = std::stod(match[1]);
                double spin_times_two = std::stod(match[2]);
                std::string parity_sign = match[3];
                int num_transitions = std::stoi(match[4]);

                log_debug("Found level: %.5f MeV with %d gammas", energy_MeV, num_transitions);

                // Then save the info:
                LevelInfo lvl;
                lvl.level_index = static_cast<int>(levels_temp.size());
                lvl.energy = energy_MeV * I3Units::MeV;
                lvl.spin = spin_times_two / 2.0;
                lvl.parity = (parity_sign == "+") ? +1 : -1;

                lvl.T12_ns = 0.0; // half life time
                lvl.tau_ns = 0.0; // mean life time

                current_level_index = lvl.level_index;

                levels_temp.push_back(lvl);
            } catch (const std::exception& e) {
                log_warn("Failed to parse level line: '%s'. Exception: %s", line.c_str(), e.what());
                continue;
            }
        } else if(inside_K40_block && std::regex_match(line, match, gamma_regex) && current_level_index >= 0) {
        // Now we look at the lines for gammas
            try {
                double E_gamma_MeV = std::stod(match[1]);  // Energy of the gamma
                double branching_ratio = std::stod(match[2]); // Relative intensity
                int final_level_idx = std::stoi(match[3]);  // index of the level of de-excitation (to which it descends)

                LevelInfo::Transition t;
                t.gamma_energy = E_gamma_MeV * I3Units::MeV;
                t.branching_ratio = branching_ratio;
                t.final_level_index = final_level_idx;

                levels_temp[current_level_index].transitions.push_back(t);
            } catch(const std::exception& e) {
                log_warn("Failed to parse gamma line: '%s'. Error: %s", line.c_str(), e.what());
                continue;
            }
        }
    }

    // Copy levels_temp into levels_map_
    std::vector<LevelInfo> & levels = levels_map_.insert({type, levels_temp}).first->second;
    log_debug("Loaded %zu %s levels.", levels.size(), element_name.c_str());
    // There are only 3 levels that will add a significant T1/2
    // All the other levels are of the order of pico or femto seconds
    // Manually set T1/2 where known from ENSDF

    for(std::tuple<double, double> delayed_level : delayed_levels_.at(type)) {
        double half_life = std::get<0>(delayed_level);
        double energy = std::get<1>(delayed_level);
        std::vector<LevelInfo>::iterator level = ClosestLevel(levels, energy);
        level->T12_ns = half_life; // from www.nndc.bnl.gov/
        level->tau_ns = half_life / std::log(2.0);
        log_debug("Set T1/2 for level %d (%.3f keV) = %.3f ns. Tau = (%.3f ns) ",
                 level->level_index, level->energy / I3Units::keV, level->T12_ns, level->tau_ns);
    }

    double metastable_energy = std::max_element(levels.begin(), levels.end(),
            [](LevelInfo const & a, LevelInfo const & b) -> bool { return a.tau_ns < b.tau_ns; }
    )->energy;

    metastable_index_.insert({type, ClosestLevel(levels, metastable_energy)->level_index});

    // Check that everything is ok for the first levels
    if(debug_) {
        size_t N = std::min(levels_map_.size(), size_t(5));
        for(size_t i = 0; i < N; ++i) {
            const LevelInfo& lvl = levels[i];
            std::string parity_init = (lvl.parity > 0) ? "+" : "-";

            log_debug("Level [%d]: %.3f keV, spin %.1f%s, T1/2 = %.3f ns",
                     lvl.level_index,
                     lvl.energy / I3Units::keV,
                     lvl.spin,
                     parity_init.c_str(),
                     lvl.T12_ns);

            for(const auto& t : lvl.transitions) {
                const LevelInfo& lvl_final = levels[t.final_level_index];
                std::string parity_final = (lvl_final.parity > 0) ? "+" : "-";

                log_debug("   Transition: [%d] %.3f keV J=%.1f%s → [%d] %.3f keV J=%.1f%s | Gamma = %.3f keV | BR = %.5f",
                         lvl.level_index,
                         lvl.energy / I3Units::keV,
                         lvl.spin,
                         parity_init.c_str(),
                         t.final_level_index,
                         lvl_final.energy / I3Units::keV,
                         lvl_final.spin,
                         parity_final.c_str(),
                         t.gamma_energy / I3Units::keV,
                         t.branching_ratio);
            }
        }
    }


    // Print totals
    size_t total_transitions = 0;
    for(LevelInfo const & level : levels) {
        total_transitions += level.transitions.size();
    }

    log_debug("Finished reading %s.dat", element_name.c_str());
    log_debug("Total number of levels parsed: %zu", levels_map_.size());
    log_debug("Total number of transitions parsed: %zu", total_transitions);

    // Finally we can write all the levels and transitions into a file for reference
    // Write all levels and transitions to a file for reference
    if(save_levels_file_) {
        std::ofstream outfile(levels_filename_);
        if(!outfile.is_open()) {
            log_error("Could not open output file %s", levels_filename_.c_str());
            return;
        }

        outfile << "# " << element_name << " Nuclear Levels and Transitions\n";
        outfile << "# =================================\n\n";

        outfile << "# ENERGY LEVELS (in keV)\n";
        for(LevelInfo const & lvl : levels) {
            outfile << "Level[" << lvl.level_index << "] "
                    << lvl.energy / I3Units::keV << " keV, "
                    << "J=" << lvl.spin
                    << ((lvl.parity > 0) ? "+" : "-");

            if(lvl.T12_ns > 0.0) {
                outfile << ", T1/2=" << lvl.T12_ns << " ns";
            }

            outfile << "\n";
        }

        outfile << "\n# TRANSITIONS\n";
        for(LevelInfo const & lvl : levels) {
            for(const auto& t : lvl.transitions) {
                const LevelInfo& lvl_final = levels.at(t.final_level_index);
                outfile << "Transition: "
                        << "[" << lvl.level_index << "] "
                        << lvl.energy / I3Units::keV << " keV → "
                        << "[" << lvl_final.level_index << "] "
                        << lvl_final.energy / I3Units::keV << " keV, "
                        << "Gamma = " << t.gamma_energy / I3Units::keV << " keV, "
                        << "BR = " << t.branching_ratio
                        << "\n";
            }
        }

        outfile.close();
        log_debug("%s levels and transitions written to %s", element_name.c_str(), levels_filename_.c_str());
        log_debug("Loaded %zu levels and %zu transitions.", levels.size(), total_transitions);
    }
}


// Look for events with K40 in the I3MCTree and save save pointers to the gamma particles
// Reconstruct the cascade for each event
// Save the cascade info in a key
// Sample a time delay for those de excitations coming from levels with a tau > ns
// Add new time delay in the I3MCTree
void MarleySimulator::AdjustGammaTimes(I3MCTreePtr mcTree, I3FramePtr frame) {
    log_trace("AdjustGammaTimes called.");

    if(!enable_gamma_time_offset_) {
        log_debug("Skipping time offsets for gammas because EnableGammaTimeOffset=False.");
        return;
    }

    std::function<bool(I3Particle const &)> is_nucleus_leaf = [&](I3Particle const & p) -> bool {
        return levels_map_.find(p.GetType()) != levels_map_.end() and I3MCTreeUtils::GetDaughtersPtr(mcTree, p).size() == 0;
    };

    std::vector<typename I3MCTree::fast_const_iterator> const K40_vec = I3MCTreeUtils::GetFilterPtr(mcTree, is_nucleus_leaf);

    std::function<bool(I3Particle const &)> is_gamma = [](I3Particle const & p) -> bool {
        return p.GetType() == I3Particle::ParticleType::Gamma;
    };

    //std::vector<typename I3MCTree::fast_const_iterator> const K40_vec = I3MCTreeUtils::GetFilterPtr(mcTree, is_K40);
    bool has_K40 = K40_vec.size() > 0;
    if(K40_vec.size() > 1) {
        log_fatal("Found more than one K40 in the MARLEY tree");
        return;
    }

    if(not has_K40) {
        frame->Put("HasK40", boost::make_shared<I3Bool>(has_K40));
        return;
    }

    I3Particle const & K40 = *K40_vec.front();
    I3Particle::ParticleType type = K40.GetType();
    I3Particle const & parent = *I3MCTreeUtils::GetParentPtr(mcTree, K40);
    std::vector<I3Particle*> const siblings = I3MCTreeUtils::GetDaughtersPtr(mcTree, parent);

    std::function<bool(I3Particle const *)> is_gamma_ptr = [=](I3Particle const * p) -> bool {
        return p->GetType() == I3Particle::ParticleType::Gamma;
    };

    std::vector<I3Particle*> gamma_rays; gamma_rays.reserve(siblings.size() - 1);
    std::copy_if(siblings.begin(), siblings.end(), std::back_inserter(gamma_rays), is_gamma_ptr);

    // Now we add the logic for the time delays
    // First we need to reconstruct the cascade

    log_trace("=== Begin reconstructing gamma cascade ===");

    // We start from the ground state (E=0)
    double running_energy = 0.0;

    // Vector for saving info of the cascade and then pass it to a key in the frame
    std::vector<NuclearCascadeStep> cascade_steps(gamma_rays.size());

    bool has_metastable = false;

    std::vector<LevelInfo> & levels = levels_map_.at(type);

    // Iterate in reverse order since we're working our way up from the ground state
    for(int i = gamma_rays.size()-1; i >= 0; --i) {
        I3Particle* g = gamma_rays.at(i);
        double gamma_energy = g->GetEnergy();

        // Find the matching level in the nuclear levels map
        // Remember the map wad made from the .dat file
        std::vector<LevelInfo>::iterator final_level = ClosestLevel(levels, running_energy);
        running_energy += gamma_energy;
        std::vector<LevelInfo>::iterator initial_level = ClosestLevel(levels, running_energy);

        // Update the runnning energy for the next step
        running_energy = initial_level->energy;

        // Sample delay with SampleDelay function if level has lifetime > ns
        double delay_ns = 0.0;
        if(initial_level->tau_ns > 0.0) {
            delay_ns = SampleDelay(initial_level->tau_ns);
        }

        has_metastable |= initial_level->level_index == metastable_index_.at(type);

        // And save the step in the cascade
        NuclearCascadeStep step;
        step.initial_level_index = initial_level->level_index;
        step.initial_level_energy = initial_level->energy;
        step.initial_level_spin = initial_level->spin;
        step.initial_level_parity = initial_level->parity;

        step.final_level_index = final_level->level_index;
        step.final_level_energy = final_level->energy;
        step.final_level_spin = final_level->spin;
        step.final_level_parity = final_level->parity;

        step.gamma_energy = gamma_energy;
        step.T12_ns = initial_level->T12_ns;
        step.tau_ns = initial_level->tau_ns;
        step.sampled_delay_ns = delay_ns;

        cascade_steps.at(i) = step;

        log_trace("Cascade step: [%d] %.3f keV J=%.1f%s → [%d] %.3f keV J=%.1f%s | gamma = %.3f keV | Sampled delay = %.3f ns",
                step.initial_level_index,
                step.initial_level_energy / I3Units::keV,
                step.initial_level_spin,
                (step.initial_level_parity > 0) ? "+" : "-",
                step.final_level_index,
                step.final_level_energy / I3Units::keV,
                step.final_level_spin,
                (step.final_level_parity > 0) ? "+" : "-",
                step.gamma_energy / I3Units::keV,
                step.sampled_delay_ns
        );

    } // end of 'for' for gammas
    // END OF LEVELS RECONSTRUCTION

    // Assign the new time to the gamma
    double cumulative_time_ns = 0.0;
    for(size_t i=0; i<gamma_rays.size(); ++i) {
        I3Particle* g = gamma_rays.at(i);
        NuclearCascadeStep & step = cascade_steps.at(i);
        cumulative_time_ns += step.sampled_delay_ns;
        step.cumulative_time_ns = cumulative_time_ns;
        double new_time = g->GetTime() + cumulative_time_ns;
        g->SetTime(new_time);
        log_trace("Gamma energy %.3f keV assigned new time %.3f ns "
                "(cumulative delay %.3f ns)",
                g->GetEnergy() / I3Units::keV,
                new_time,
                step.cumulative_time_ns
        );
    }

    if(debug_) {
        // Save the cascade into the frame
        I3VectorStringPtr cascade_info = boost::make_shared<I3VectorString>();
        for(const auto& step : cascade_steps) {
            std::ostringstream ss;

            ss << "Step:\n";
            ss << "  From Level [" << (step.initial_level_index >= 0 ? std::to_string(step.initial_level_index) : "UNKNOWN") << "] "
                << step.initial_level_energy / I3Units::keV << " keV, "
                << "J=" << step.initial_level_spin
                << (step.initial_level_parity > 0 ? "+" : "-") << "\n";

            ss << "  To Level [" << (step.final_level_index >= 0 ? std::to_string(step.final_level_index) : "UNKNOWN") << "] "
                << step.final_level_energy / I3Units::keV << " keV, "
                << "J=" << step.final_level_spin
                << (step.final_level_parity > 0 ? "+" : "-") << "\n";

            ss << "  Gamma Energy = " << step.gamma_energy / I3Units::keV << " keV\n";
            ss << "  T1/2 = " << step.T12_ns << " ns\n";
            ss << "  Tau = " << step.tau_ns << " ns\n";
            ss << "  Sampled Delay = " << step.sampled_delay_ns << " ns\n";
            ss << "  Cumulative time = " << step.cumulative_time_ns << " ns\n";

            cascade_info->push_back(ss.str());
        }
        frame->Put("MarleyGammaCascadeInfo", cascade_info);
        log_trace("Stored MarleyGammaCascadeInfo in frame.");
        // ===Some counters ===
        ++frames_with_cascade;

        if(has_metastable) {
            // Save the energy difference if it falls into the meta-stable
            double highest_level_energy = cascade_steps.front().initial_level_energy;
            double metastable_energy_diff = -1.0;
            ++frames_with_metastable;

            // Calculate difference between max E and metastable
            double metastable_energy = levels.at(metastable_index_.at(type)).energy;
            metastable_energy_diff = highest_level_energy - metastable_energy;
            metastable_energy_diffs.push_back(metastable_energy_diff);

            log_debug("Frame with metastable level. Energy diff from highest level = %.3f keV",
                    metastable_energy_diff / I3Units::keV);

            // And the difference in Energy from max to meta stable
            frame->Put("MarleyMetastableDeltaE",
                    boost::make_shared<I3Double>(metastable_energy_diff));
        }
    }
    // Saves in the frame if has metastable level n the cascade
    frame->Put("HasMetaStable", boost::make_shared<I3Bool>(has_metastable));

    log_trace("=== Finished reconstructing gamma cascade ===");

    // Save other important variables for analysis
    auto cumulative_times_vec = boost::make_shared<I3VectorDouble>();
    auto gamma_energies_vec = boost::make_shared<I3VectorDouble>();

    for(const auto& step : cascade_steps) {
        cumulative_times_vec->push_back(step.cumulative_time_ns);
        gamma_energies_vec->push_back(step.gamma_energy);
    }

    frame->Put("MarleyGammaCumulativeTimes", cumulative_times_vec);
    frame->Put("MarleyGammaEnergies", gamma_energies_vec);

    frame->Put("HasK40", boost::make_shared<I3Bool>(has_K40));

} // End of AdjustGammaTimes

// This method samples the delay for a given lifetime
// using an exponential distribution.
// t= -tau * ln (1-u)
// tau in ns
// u= random number between [0,1]
double MarleySimulator::SampleDelay(double mean_lifetime_ns) {
    if(mean_lifetime_ns <= 0.0)
        return 0.0;

    double u = rng_->Uniform(0.0, 1.0);
    return -mean_lifetime_ns * std::log1p(-u);
}


void MarleySimulator::FillSimulationFrame(I3FramePtr frame) {
    I3Int32Ptr obj = boost::make_shared<I3Int32>(1);
    frame->Put("MarleyConfiguration", obj);
}

void MarleySimulator::Finish() {
    if(debug_) {
        log_debug("====== MARLEY SIMULATOR SUMMARY ======");
        log_debug("Frames with Marley Events: %zu", frames_in_total);
        log_debug("Frames with MarleyGammaCascadeInfo: %zu", frames_with_cascade);
        log_debug("Frames passing through metastable level [n=4] (1.64364 MeV): %zu", frames_with_metastable);

        // I want to check if the E of these prompt gamma is around 2.74 MeV as with the solar neutrinos
        if(!metastable_energy_diffs.empty()) {
            double sum = 0.0;
            for(auto x : metastable_energy_diffs) sum += x;
            double mean_diff = sum / metastable_energy_diffs.size();
            log_debug("Average energy difference to metastable level: %.3f keV", mean_diff / I3Units::keV);
        }

        log_debug("=======================================");
    }
}
