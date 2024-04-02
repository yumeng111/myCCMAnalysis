#include <icetray/I3Context.h>
#include <icetray/I3Frame.h>
#include <icetray/I3Bool.h>
#include <g4-larsim/CCMSimulator.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3EventHeader.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3ScintRecoPulseSeriesMap.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Double.h>

#include <stdexcept>

I3_MODULE(CCMSimulator);

// Name of the IncEventID boolean in the frame (see I3MCTimeGeneratorService)
const std::string CCMSimulator::INC_ID_NAME = "MCTimeIncEventID";

CCMSimulator::CCMSimulator(const I3Context& context): I3Module(context) {
    injector_ = CCMParticleInjectorPtr();
    response_ = CCMDetectorResponsePtr();

    responseServiceName_ = "CCM200Response";
    AddParameter("ResponseServiceName", "Name of the detector response service.", responseServiceName_);

    injectorServiceName_ = "CCMSimpleInjector"; 
    AddParameter("InjectorServiceName", "Name of the injector service.", injectorServiceName_);

    mcPrimaryName_ = "MCPrimary";
    AddParameter("PrimaryName", "Name of the primary particle in the frame.", mcPrimaryName_);

    PMTHitSeriesName_ = "PMTMCHitsMap";
    AddParameter("PMTHitSeriesName", "Name of the resulting PMT hit map in the frame.", PMTHitSeriesName_);
    
    LArHitSeriesName_ = "LArMCHitsSeries";
    AddParameter("LArHitSeriesName", "Name of the resulting LAr hit series in the frame.", LArHitSeriesName_);

    AddOutBox("OutBox");
}


CCMSimulator::~CCMSimulator() { }

 
void CCMSimulator::Configure() {
    log_info("Configuring the CCMSimulator:");

    GetParameter("InjectorServiceName", injectorServiceName_);
    log_info("+ Injector Service: %s", injectorServiceName_.c_str());
    injector_ = GetContext().Get<CCMParticleInjectorPtr>(injectorServiceName_);
    if(!injector_) log_fatal("No injector service \"%s\" in context", injectorServiceName_.c_str());

    GetParameter("ResponseServiceName", responseServiceName_);
    log_info("+ Response Service: %s", responseServiceName_.c_str());
    response_ = GetContext().Get<CCMDetectorResponsePtr>(responseServiceName_);
    if(!response_) log_fatal("No response service \"%s\" in context", responseServiceName_.c_str());

    GetParameter("PrimaryName", mcPrimaryName_);
    log_info("+ MC primary: %s", mcPrimaryName_.c_str());

    GetParameter("PMTHitSeriesName", PMTHitSeriesName_);
    log_info("+ PMT hit series : %s", PMTHitSeriesName_.c_str());
    
    GetParameter("LArHitSeriesName", LArHitSeriesName_);
    log_info("+ LAr hit series : %s", LArHitSeriesName_.c_str());

    // initialize injector and response services
    injector_->Configure();
    response_->Configure();
    response_->Initialize();

    // set our controls over SD
    PMTSDStatus_ = response_->GetPMTSDStatus();
    LArSDStatus_ = response_->GetLArSDStatus();
}


// we don't have DetectorStatus frames at this moment...
void CCMSimulator::DetectorStatus(I3FramePtr frame) {
    // Configure and initialize response and injector services
    //response_->Configure();
    //response_->Initialize();

    //injector_->SetResponseService(response_);
    //injector_->Configure();

    PushFrame(frame);
}

void CCMSimulator::DAQ(I3FramePtr frame) {
    log_debug("   Simulating CCM");
    
    // let's grab the mcPrimary from the injector
    mcTree_ = injector_->GetMCTree();
    
    std::vector<I3Particle*> primary_particles = I3MCTreeUtils::GetPrimariesPtr(mcTree_);
    
    for (size_t p = 0; p < primary_particles.size(); p++){
        // Tell the response service of a new event
        response_->BeginEvent(*primary_particles[p]);

        // Tell the response service that the event has ended
        response_->EndEvent();

        // Now grab the map between CCMPMTKey and std::vector<CCMMCPE> to save to frame
        // also grab hits in fiducial argon for voxelization
        // note -- if SD not enabled, these will just be empty objects
        CCMMCPEMap = response_->GetHitsMap(); 
        CCMMCPEList = response_->GetVoxelHits();
    }
    
    if (!saved_simulation_setup_) {
        // now save simulation set up info (just MC tree for now) into S frame
        I3FramePtr simframe(new I3Frame(I3Frame::Simulation));
        // put mcTree into frame
        simframe->Put(mcPrimaryName_, mcTree_);
        // push simframe
        PushFrame(simframe);
        
        saved_simulation_setup_ = true; 
    }
    
    // Put the hits to the DAQ frame 
    if (PMTSDStatus_){
        frame->Put(PMTHitSeriesName_, CCMMCPEMap);
    }
    if (LArSDStatus_){
        frame->Put(LArHitSeriesName_, CCMMCPEList);
    }

    // push daq frame
    PushFrame(frame);
}




