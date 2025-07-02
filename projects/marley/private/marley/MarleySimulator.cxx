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

#include "phys-services/I3RandomService.h"

#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/MarleySimulator.h"

#include <fstream>
#include <sstream>
#include <regex> //searches for patterns

I3_MODULE(MarleySimulator);

MarleySimulator::MarleySimulator(const I3Context& context) : I3ConditionalModule(context),
    output_mc_tree_name_("I3MCTree"), input_mc_tree_name_("SIRENMarleyInjectionTree") {
    AddParameter("MarleySearchPath", "Search path for Marley environment file", marley_search_path_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
    AddParameter("InputMCTreeName", "Name of the MCTree in the frame.", input_mc_tree_name_);
}

void MarleySimulator::Configure() {
    GetParameter("MarleySearchPath", marley_search_path_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
    GetParameter("InputMCTreeName", input_mc_tree_name_);

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

    this->LoadK40Transitions("K40.dat"); //Test.. I would like to use the original K.dat from marley that has all the isotopes but later
    //Saves the info of energy levels and transitions obtained fromLoadK40Transitions
    this->PrintLevelsAndTransitions("K40_levels_and_transitions.txt");
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
    // Now save to a frame
    frame->Put(output_mc_tree_name_, outputMCTree);
    PushFrame(frame);
}

//This method gets the transitions from the K40.dat file
//That file is an extract of the marley/data/structure/K.dat only for the K40
//and is located in the same path as the main python script
//In this method we obtain two important vectors:
//energy_levels_ -> energies of all the levels
//k40_transitions_ -> list of all the gamma transitions
void MarleySimulator::LoadK40Transitions(const std::string& filename) {
    std::ifstream infile(filename);
    std::string line;
    std::vector<double> level_energies;
    int level_index = 0;

    //Header: Z, A, Num of levels
    std::regex header_regex(R"(^\d+\s+\d+\s+\d+$)");
    //Excitation energy of level, 2*Spin, Parity, Total gammas
    std::regex level_regex(R"(^\s*([\d\.Ee+-]+)\s+[-+]?\d+\s+[+-]\s+(\d+))");

    if (!infile.is_open()) {
        log_error("Could not open file %s", filename.c_str());
        return;
    }

    //Start reading the file line by line
    while (std::getline(infile, line)) {
        log_info("Line raw content: '%s'", line.c_str());

        // Skip the Z A num_levels if header line (header_regex)
        if (std::regex_match(line, header_regex)) {
            log_info("Skipping header line.");
            continue;
        }

        std::smatch match; //que es esto
        //If it is not header check if it is a line with exitation energy level (level_regex)
        if (std::regex_match(line, match, level_regex)) {
            //Get Initial energy and number of gammas
            double E_initial = std::stod(match[1]);
            int num_gammas = std::stoi(match[2]);
            log_info("Found level: %.5f MeV with %d gammas", E_initial, num_gammas);

            //Then save the energy in a vector
            level_energies.push_back(E_initial);
            level_index++;

            //Now we loop over the number of gammas that we found previously
            for (int i = 0; i < num_gammas && std::getline(infile, line); ++i) {
                std::istringstream gamma_iss(line); //que es esto
                double E_gamma, RI; //Energy of the gamma, Relative intensity
                int tmp_index; //index of the level of de-excitation (to which it descends)

                if (gamma_iss >> E_gamma >> RI >> tmp_index) {
                    //The final level index indicates to what level the transition goes down:
                    //If the index is valid, gets the final energy, if not = 0.0
                    size_t Lf_index = static_cast<size_t>(tmp_index);
                    double E_final = (Lf_index < level_energies.size())
                        ? level_energies[Lf_index]
                        : 0.0;

                    //Save the transition in this structure:
                    k40_transitions_.push_back({E_initial, E_final, E_gamma, RI});

                    log_info("   Transition: %.5f MeV -> %.5f MeV | Gamma Energy = %.5f MeV | RI = %.5f",
                        E_initial, E_final, E_gamma, RI);
                } else {
                    log_warn("Could not parse gamma line: '%s'", line.c_str());
                }
            }
        }
    }
    //saves the energies in a vector
    energy_levels_ = level_energies;


    // Summary printout
    //stdcout::Transitions for K40 loaded in MarleySimulator
    log_info("Finished reading K40.dat");
    log_info(" Total number of levels parsed: %zu", energy_levels_.size());
    log_info("Total number of transitions parsed: %zu", k40_transitions_.size());

    //Some checks...
    // Print the first few transitions
    size_t N = std::min(k40_transitions_.size(), size_t(10));
    for (size_t i = 0; i < N; ++i) {
        const auto& tr = k40_transitions_[i];
        log_info("Transition %zu: %.5f MeV → %.5f MeV | Gamma Energy= %.5f MeV | RI = %.5f",
                 i, tr.initial_energy, tr.final_energy, tr.gamma_energy, tr.relative_intensity);
    }

}

//This is for testing but it could be useful
void MarleySimulator::PrintLevelsAndTransitions(const std::string& outfilename) {
    std::ofstream outfile(outfilename);
    if (!outfile.is_open()) {
        log_error("Could not open output file %s", outfilename.c_str());
        return;
    }

    outfile << "# ENERGY LEVELS (in MeV)\n";
    for (size_t i = 0; i < energy_levels_.size(); ++i) {
        outfile << "Level[" << i << "] = " << energy_levels_[i] << " MeV\n";
    }

    outfile << "\n# TRANSITIONS\n";
    for (size_t i = 0; i < k40_transitions_.size(); ++i) {
        const auto& tr = k40_transitions_[i];
        outfile << "Transition[" << i << "]: "
                << tr.initial_energy << " MeV -> "
                << tr.final_energy << " MeV, "
                << "gamma = " << tr.gamma_energy << " MeV, "
                << "RI = " << tr.relative_intensity << "\n";
    }

    outfile.close();
    log_info("K40 levels and transitions written to %s", outfilename.c_str());
}




void MarleySimulator::FillSimulationFrame(I3FramePtr frame) {
    I3Int32Ptr obj = boost::make_shared<I3Int32>(1);
    frame->Put("MarleyConfiguration", obj);
}

