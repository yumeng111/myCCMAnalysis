#ifndef CCM200RESPONSE_H
#define CCM200RESPONSE_H
// standard library stuff
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/g4classes/G4Interface.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3ServiceBase.h"

#include "phys-services/I3RandomService.h"

#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>
#include <algorithm>

class CCM200Response : public CCMDetectorResponse {
private:
    CCM200Response operator= (const CCM200Response& rhs);

    std::string visMacroFile_;

    bool PMTSDStatus_; // turn PMT SD on/off
    bool LArSDStatus_; // turn fiducial LAr SD on/off
    bool SodiumSourceRun_; // turn source run on/off 
    double SodiumSourceLocation_; // location of end of source rod 

    G4Interface* g4Interface_;
        
    SET_LOGGER("CCM200Response");

public:

    CCM200Response(const I3Context& context);
    virtual ~CCM200Response() override;

    virtual void Configure() override;
    virtual void Initialize() override;
    virtual void BeginEvent(const I3Particle& primary) override;
    virtual void EndEvent(size_t event_idx) override;
    virtual void TerminateRun() override;
    virtual boost::shared_ptr<CCMMCPESeriesMap> GetHitsMap() override {return CCMMCPEMap;};
    virtual I3MCTreePtr GetLArEnergyDep() override { return LArEnergyDep; }
    virtual bool GetPMTSDStatus() override { return PMTSDStatus_; }
    virtual bool GetLArSDStatus() override { return LArSDStatus_; }

    boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap = boost::make_shared<CCMMCPESeriesMap> ();
    I3MCTreePtr LArEnergyDep = boost::make_shared<I3MCTree>();

};

#endif // CCM200RESPONSE_H


