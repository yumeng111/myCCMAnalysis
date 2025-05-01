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

class MarleySimulator : public I3ConditionalModule {
public:
    MarleySimulator(const I3Context& context);
    ~MarleySimulator() = default;

    void Configure();
    void Simulation(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void FillSimulationFrame(I3FramePtr frame);

private:
    bool seen_s_frame_ = false;

    std::string output_mc_tree_name_;
    std::string input_mc_tree_name_;
    std::string marley_search_path_;
    // Add marley_generator_as member of the class
    marley::Generator marley_generator_;
};

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

void MarleySimulator::FillSimulationFrame(I3FramePtr frame) {
    I3Int32Ptr obj = boost::make_shared<I3Int32>(1);
    frame->Put("MarleyConfiguration", obj);
}
