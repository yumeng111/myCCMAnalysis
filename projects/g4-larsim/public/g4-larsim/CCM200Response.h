#ifndef CCM200RESPONSE_H
#define CCM200RESPONSE_H
// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

#include <g4-larsim/g4classes/G4Interface.h>
#include "g4-larsim/CCMDetectorResponse.h"

#include <icetray/I3Frame.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <simclasses/CCMMCPE.h>
#include <icetray/I3ServiceBase.h>

class CCM200Response : public CCMDetectorResponse {
private:
    CCM200Response operator= (const CCM200Response& rhs);

    std::string mc_tree_name_;
    std::string output_hits_map_name_;
    std::string visMacroFile_;

    bool PMTSDStatus_; // turn PMT SD on/off
    bool LArSDStatus_; // turn fiducial LAr SD on/off

    G4Interface* g4Interface_;
        
    SET_LOGGER("CCM200Response");

public:

    CCM200Response(const I3Context& context);
    virtual ~CCM200Response() override;

    virtual void Configure() override;
    virtual void Initialize() override;
    virtual void BeginEvent(const I3Particle& primary) override;
    virtual void EndEvent() override;
    virtual boost::shared_ptr<CCMMCPESeriesMap> GetHitsMap() override {return CCMMCPEMap;};
    virtual boost::shared_ptr<CCMMCPESeries> GetVoxelHits() override { return CCMMCPEList; }
    virtual bool GetPMTSDStatus() override { return PMTSDStatus_; }
    virtual bool GetLArSDStatus() override { return LArSDStatus_; }

    boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap = boost::make_shared<CCMMCPESeriesMap> ();
    boost::shared_ptr<CCMMCPESeries> CCMMCPEList = boost::make_shared<CCMMCPESeries> ();

};

#endif // CCM200RESPONSE_H


