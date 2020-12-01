#ifndef CCMPMTResponse_h
#define CCMPMTResponse_h

#include "Utility.h"
#include "CCMModule.h"

#include <vector>
#include <iterator>
#include <algorithm>
#include <random>
#include <map>
#include <utility>

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
   double GetADCValue(int row, int col, double initialCharge);

  private:

    //private methods

  private:

    //private data members
    std::shared_ptr<MCTruth> fMCTruth;
    std::shared_ptr<Pulses> fPulses;

    std::random_device fRD;
    std::mt19937_64 fMT;
    std::uniform_real_distribution<double> fUniform;

    std::map<std::pair<int,int>,double> fSPEWeights;
};

#endif // CCMPMTResponse_h

