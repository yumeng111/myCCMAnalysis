#ifndef CCM200RESPONSE_H
#define CCM200RESPONSE_H
// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

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
#include "g4-larsim/CCMDetectorResponse.h"

class CCM200Response : public CCMDetectorResponse {
private:
    CCM200Response operator= (const CCM200Response& rhs);

    std::string geometry_name_;
    std::string mc_tree_name_;
    std::string pulse_series_output_name_;
    std::string visMacroFile_;

    G4Interface* g4Interface_;
        
    // return CCMMCPEMap
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<CCMMCPE>>> GetHitsMap(){ return CCMMCPEMap; }

    SET_LOGGER("CCM200Response");

public:

    CCM200Response(const I3Context& context);
    virtual ~CCM200Response() override;

    void Configure();
    virtual void Initialize() override;
    virtual void BeginEvent(const I3Particle& primary) override;
    virtual void EndEvent() override;
    virtual boost::shared_ptr<I3Map<CCMPMTKey, std::vector<CCMMCPE>>> GetHitsMap() override;
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<CCMMCPE>>> CCMMCPEMap = boost::make_shared<I3Map<CCMPMTKey, std::vector<CCMMCPE>>> ();
};

#endif // CCM200RESPONSE_H


