#ifndef CCMQuenchingFactor_h
#define CCMQuenchingFactor_h

#include "Utility.h"
#include "CCMModule.h"

#include <vector>
#include <iterator>
#include <algorithm>
#include <random>
#include <map>
#include <utility>

class CCMQuenchingFactor : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMQuenchingFactor(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMQuenchingFactor(const CCMQuenchingFactor& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMQuenchingFactor();

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

    bool TestQuench(double quench);

    void ResetQuenchingFactor(std::shared_ptr<MCTruth> mcTruth, double quench);

  private:

    //private methods

  private:

    //private data members
    std::shared_ptr<MCTruth> fMCTruth;

    std::random_device fRD;
    std::mt19937_64 fMT;
    std::uniform_real_distribution<double> fUniform;

};

#endif // CCMQuenchingFactor_h

//Process an event, then process a trigger, then assign the QFvalue of the MCTruth according to user input. Have a method to rethrow according to QF. 
