
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3EventHeader.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/physics/I3ScintRecoPulseSeriesMap.h"

#include "g4-larsim/CCMSimulator.h"

#include "icetray/I3Bool.h"
#include "icetray/I3Frame.h"
#include "icetray/I3Context.h"

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
    
    LArMCTreeName_ = "LArMCTree";
    AddParameter("LArMCTreeName", "Name of the MC tree containing energy depositions in LAr", LArMCTreeName_);

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
    
    GetParameter("LArMCTreeName", LArMCTreeName_);
    log_info("+ LAr MC Tree : %s", LArMCTreeName_.c_str());

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
    
    std::cout << "simulating " << primary_particles.size() << " events" << std::endl;
    for (size_t p = 0; p < primary_particles.size(); p++){
        //std::cout << "injecting " << primary_particles[p]->GetType() << " and pos = " << primary_particles[p]->GetPos() << std::endl;

        // Tell the response service of a new event
        response_->BeginEvent(*primary_particles[p]);

        // Tell the response service that the event has ended
        // this will also populate the map between CCMPMTKey and std::vector<CCMMCPE> to save to frame
        // also grab hits in fiducial argon for voxelization
        // note -- if SD not enabled, these will just be empty objects
        I3MCTreePtr LArEnergyDep = boost::make_shared<I3MCTree>();
        boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap = boost::make_shared<CCMMCPESeriesMap> ();
        response_->EndEvent(LArEnergyDep, CCMMCPEMap);

        // now save to put into frames and push at the end 
        AllEventsLArEnergyDep.push_back(boost::make_shared<I3MCTree>(*LArEnergyDep));
        AllEventsCCMMCPEMap.push_back(boost::make_shared<CCMMCPESeriesMap>(*CCMMCPEMap));

    }

    // terminate geant4    
    response_->TerminateRun();
}


void CCMSimulator::Finish(){
    
    // now save simulation set up info (just MC tree for now) into S frame
    I3FramePtr simframe(new I3Frame(I3Frame::Simulation));
    // put mcTree into frame
    simframe->Put(mcPrimaryName_, mcTree_);
    // push simframe
    PushFrame(simframe);

    // now go through our deque,
    // save information for each event to frames,
    // pop from deque, and push frame
    while (AllEventsLArEnergyDep.size() > 0){
        I3FramePtr tempframe(new I3Frame(I3Frame::DAQ));
        
        tempframe->Put(PMTHitSeriesName_, AllEventsCCMMCPEMap.at(0));
        // remove CCMMCPEMap from cache
        AllEventsCCMMCPEMap.pop_front();
        
        tempframe->Put(LArMCTreeName_, AllEventsLArEnergyDep.at(0));
        // remove LArMCTree from cache
        AllEventsLArEnergyDep.pop_front();
        
        PushFrame(tempframe);
    }
    
}



