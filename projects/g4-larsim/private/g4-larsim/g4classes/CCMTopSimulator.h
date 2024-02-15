#ifndef _TOPSIMULATOR_I3TOPSIMULATOR_H_
#define _TOPSIMULATOR_I3TOPSIMULATOR_H_


#include <icetray/I3Module.h>
#include <topsimulator/interface/I3InjectorService.h>
#include <topsimulator/interface/I3IceTopResponseService.h>
#include <topsimulator/ExtendedI3Particle.h>
//#include <dataclasses/TankKey.h>
//#include <dataclasses/ScintKey.h>
#include <dataclasses/LightKey.h>
/**
 * \brief The CCMTopimulator module handles the whole simulation and 
 * stores the results in the frame.
 */

class CCMTopSimulator : public I3Module 
{
public:
  CCMTopSimulator(const I3Context& context);

  ~CCMTopSimulator();

  void Configure();
  
  void DetectorStatus(I3FramePtr frame);
  
  void DAQ(I3FramePtr frame);

private:
  
  /// This function generates an event header and put it to the frame
  void WriteEventHeader(I3FramePtr frame, int32_t runID, int32_t evtID);

  //TankKey GetTankKey(std::string key) const;
  //ScintKey GetScintKey(std::string key) const;
  LightKey GetLightKey(std::string key) const;
  //AirShowerComponent GetAirShowerComponent(const I3Particle& p) const;
  //std::map<std::string, int32_t> GetAirShowerComponentNameMap() const;

  std::string injectorServiceName_;
  std::string responseServiceName_;
  std::string mcPrimaryName_;
  std::string hitSeriesName_;

  std::string peSeriesName_;
  std::string scintpeSeriesName_;

  std::string cherSeriesName_;

  std::string scintSeriesName_;
  
  std::string testPulsesName_;
  std::string icMCTreeName_;
  std::string itMCTreeName_;
  double hitBinWidth_;
  double muEnergyCut_;
  bool writeEvtHeader_;
  bool createSframe_;
  int32_t compressPEs_;
  bool useInjectorComponents_;
  //bool useScintillator_;
  //std::set<TankKey> tankKeys_;
  //std::set<ScintKey> scintKeys_;
  std::set<LightKey> lightKeys_;
  int32_t sampleCount_;

  I3InjectorServicePtr injector_;
  I3IceTopResponseServicePtr response_;
  
  static const std::string INC_ID_NAME;
  
  SET_LOGGER("CCMTopSimulator");
};

#endif
