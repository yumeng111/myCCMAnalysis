#ifndef _CCMSIMULATOR_CCMSIMULATOR_H_
#define _CCMSIMULATOR_CCMSIMULATOR_H_


#include <icetray/I3Module.h>
#include <topsimulator/interface/I3InjectorService.h>
#include <topsimulator/interface/I3IceTopResponseService.h>
#include <topsimulator/ExtendedI3Particle.h>
#include <simclasses/CCMMCPE.h>
/**
 * \brief The CCMTopimulator module handles the whole simulation and 
 * stores the results in the frame.
 */

class CCMSimulator : public I3Module 
{
public:
  CCMSimulator(const I3Context& context);

  ~CCMSimulator();

  void Configure();
  
  void DetectorStatus(I3FramePtr frame);
  
  void DAQ(I3FramePtr frame);

private:
  
  /// This function generates an event header and put it to the frame
  void WriteEventHeader(I3FramePtr frame, int32_t runID, int32_t evtID);

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
  int32_t sampleCount_;

  I3InjectorServicePtr injector_;
  I3IceTopResponseServicePtr response_;
  
  static const std::string INC_ID_NAME;
  
  SET_LOGGER("CCMSimulator");
};

#endif
