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
    injector_ = CCMSimpleInjectorPtr();
    response_ = CCM200ResponsePtr();

    responseServiceName_ = I3DefaultName<CCM200Response>::value();
    AddParameter("ResponseServiceName",
         "Name of the detector response service.", 
         responseServiceName_);

    injectorServiceName_ = I3DefaultName<CCMSimpleInjector>::value();
    AddParameter("InjectorServiceName",
         "Name of the injector service.", 
         injectorServiceName_);

    mcPrimaryName_ = "MCPrimary";
    AddParameter("PrimaryName",
         "Name of the primary particle in the frame.",
         mcPrimaryName_);

    hitSeriesName_ = "CCMMCHits";
    AddParameter("CCMHitSeriesName",
         "Name of the resulting detector hit series in the frame.",
         hitSeriesName_);

    writeEvtHeader_ = true;
    AddParameter("WriteEventHeader",
         "Create an event header to store the run and event numbers.",
         writeEvtHeader_);

    createSframe_ = true;
    AddParameter("CreateSFrame", "Create an SFrame for CCM simulator.", createSframe_);

    AddOutBox("OutBox");

    sampleCount_ = 0;
}


CCMSimulator::~CCMSimulator() { }

 
void CCMSimulator::Configure() {
  log_info("Configuring the CCMSimulator:");

  GetParameter("InjectorServiceName", injectorServiceName_);
  log_info("+ Injector Service: %s", injectorServiceName_.c_str());
  injector_ = GetContext().Get<CCMSimpleInjectorPtr>(injectorServiceName_);
  if(!injector_) log_fatal("No injector service \"%s\" in context", injectorServiceName_.c_str());
  
  GetParameter("ResponseServiceName", responseServiceName_);
  log_info("+ Response Service: %s", responseServiceName_.c_str());
  response_ = GetContext().Get<CCM200ResponseServicePtr>(responseServiceName_);
  if(!response_) log_fatal("No response service \"%s\" in context", responseServiceName_.c_str());
  
  GetParameter("PrimaryName", mcPrimaryName_);
  log_info("+ MC primary: %s", mcPrimaryName_.c_str());
  
  GetParameter("CCMHitSeriesName", hitSeriesName_);
  log_info("+ CCM hit series : %s", hitSeriesName_.c_str());
  
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
    
    // let's initialize our response and injector services 
    response_->Configure();
    response_->Initialize();
    injector_->Configure();
    injector_->Initialize();
    
    // let's grab the mcPrimary from the frame
    I3MCTreePtr mcTree = frame->Get<I3MCTreePtr>(mcPrimaryName_); 

    // Tell the response service of a new event
    response_->BeginEvent(mcTree);

    // Tell the response service that the event has ended
    response_->EndEvent();

    // Now grab the map between CCMPMTKey and std::vector<CCMMCPE> to save to frame
    CCMMCPEMap response_->GetHitsMap(); 

    // now save everything into S frame
    
    I3FramePtr simframe(new I3Frame(I3Frame::Simulation));
    // Put primary to the frame (if it is valid)
    if((!mcPrimaryName_.empty()) && (mcPrimary->GetType()!=I3Particle::unknown)) {
        simframe->Put(mcPrimaryName_, mcPrimary);
    }

    // Put the hits to the frame 
    if (!hitSeriesName_.empty()) {
        simframe->Put(hitSeriesName_, CCMMCPEMap);
    }

    // push simframe
    PushFrame(simframe);
    
    // push daq frame
    PushFrame(frame);
}




