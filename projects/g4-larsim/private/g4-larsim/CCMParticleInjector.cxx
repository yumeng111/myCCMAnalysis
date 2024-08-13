#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"

#include "g4-larsim/CCMParticleInjector.h"

#include <vector>
#include <string>
#include <algorithm>

CCMParticleInjector::CCMParticleInjector(const I3Context& context) : I3ConditionalModule(context),
    mcPrimaryName_("MCPrimaries"), output_mc_tree_name_("I3MCTree") {
    AddParameter("PrimariesName", "Name of the output frame key for the primary particles.", mcPrimaryName_);
    AddParameter("OutputMCTreeName", "Name of the MCTree in the frame.", output_mc_tree_name_);
}

void CCMParticleInjector::Configure() {
    GetParameter("PrimariesName", mcPrimaryName_);
    GetParameter("OutputMCTreeName", output_mc_tree_name_);
}

void CCMParticleInjector::Simulation(I3FramePtr frame) {
    seen_s_frame_ = true;
    FillSimulationFrame(frame);
    PushFrame(frame);
}

void CCMParticleInjector::DAQ(I3FramePtr frame) {
    if(not seen_s_frame_) {
        I3FramePtr sim_frame = boost::make_shared<I3Frame>(I3Frame::Simulation);
        FillSimulationFrame(sim_frame);
        PushFrame(sim_frame);
        seen_s_frame_ = true;
    }

    I3MCTreePtr mcTree = GetMCTree();

    std::vector<I3Particle> primaries = I3MCTreeUtils::GetPrimaries(*mcTree);
    I3VectorI3ParticlePtr primary(new I3Vector<I3Particle>(primaries.begin(), primaries.end()));

    // now save to a frame
    frame->Put(mcPrimaryName_, primary);
    frame->Put(output_mc_tree_name_, mcTree);
    PushFrame(frame);
}

void CCMParticleInjector::FillSimulationFrame(I3FramePtr frame) {
    I3FrameObjectPtr obj = this->GetSimulationConfiguration();
    frame->Put("InjectorConfiguration", obj);
}
