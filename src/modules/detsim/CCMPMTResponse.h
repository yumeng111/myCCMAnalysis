#ifndef CCMPMTResponse_h
#define CCMPMTResponse_h

#include "Utility.h"
#include "CCMModule.h"

#include <random>
#include <map>
#include <utility>

class TFile;
class TH1D;

class CCMPMTResponse : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMPMTResponse(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMPMTResponse(const CCMPMTResponse& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMPMTResponse();

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

    void ConnectMCTruth(std::shared_ptr<MCTruth> mcTruth) { fMCTruth = mcTruth; }
    void ConnectPulses(std::shared_ptr<Pulses> pulses) { fPulses = pulses; }

   bool PMTQE(double energy, double angle); 
   void FillSPEWeights();
   void GetADCValueAndLength(size_t key, double & adc, double & length);
   void AddToWaveform(int key, double time, double adc, double length);

  private:

    //private methods

  private:

    //private data members
    std::shared_ptr<MCTruth> fMCTruth;
    std::shared_ptr<Pulses> fPulses;

    std::random_device fRD;
    std::mt19937_64 fMT;
    std::uniform_real_distribution<double> fUniform;

    std::map<int,std::pair<double,double>> fSPEWeights;
    std::map<int,std::vector<double>> fWaveforms;
    std::map<int,std::shared_ptr<TH1D>> fWaveformsHist;

    TFile * fFile;

    double fHighEnergy;// = 4.0;
    double fLowEnergy;// = 2.0;
    double fQE;// = 0.2;

    bool fSquareWF;
    bool fTriangleWF;

    double fTotalHits;
    double fAfterQE;
    double fAfterADC;
    double fTotalPulses;
    double fAvgPulseLength;
    double fAvgPulseIntegral;

};

#endif // CCMPMTResponse_h

