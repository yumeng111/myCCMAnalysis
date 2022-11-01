#ifndef CCMTimeShift_h
#define CCMTimeShift_h

#include <random>

#include "CCMAnalysis/CCMFramework/CCMModule.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

class TFile;
class TH1D;

enum CCMTimeDists_t {
  kCCMUniform = 0,
  kCCMBeamPi0 = 1,
  kCCMBeamPipm = 2,
  kCCMBeamMupm = 3,
  kCCMTimeTotal = 4
};

class CCMTimeShift : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMTimeShift(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMTimeShift(const CCMTimeShift& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMTimeShift();

    /*!
     *  \brief This is where the action takes place
     *  \return CCMResult_t the result of running this module
     */
    CCMResult_t ProcessTrigger();

    /*!
     *  \brief Returns true of the job has ended
     *  \return CCMResult_t result if the Job had ended
     */
    CCMResult_t EndOfJob();

    /*!
     *  \brief Configures things that are hardware specific 
     *  \param c holds hardware and event specific settings
     */
    void Configure(const CCMConfig& c);

    void ConnectPulses(std::shared_ptr<Pulses> pulses) { fPulses = pulses; }

    double GetTimeFromUniform(double start, double end);
    double GetTimeFromBeamPi0(double start, double end, double beamWidth);
    double GetTimeFromBeamPipm(double start, double end, double beamWidth);
    double GetTimeFromBeamMupm(double start, double end, double beamWidth);

  private:

    //private methods

  private:

    //private data members
    std::shared_ptr<Pulses> fPulses;

    std::random_device fRD;
    std::mt19937_64 fMT;

    CCMTimeDists_t fTimeDist;
    double fStartNS;
    double fEndNS;
    double fBeamWidthNS;

    double fTotalEvents;
    double fOutOfWindowBeamPi0;
    double fOutOfWindowBeamPipm;
    double fOutOfWindowBeamMupm;

};

#endif // CCMTimeShift_h


