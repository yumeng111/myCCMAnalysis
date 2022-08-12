#ifndef CCMStrobeOverlay_h
#define CCMStrobeOverlay_h

#include <vector>
#include <iterator>
#include <algorithm>

#include "CCMAnalysis/modules/framework/CCMModule.h"
#include "CCMAnalysis/utils/Utility.h"

class CCMRootIO;

class CCMStrobeOverlay : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMStrobeOverlay(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMStrobeOverlay(const CCMStrobeOverlay& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMStrobeOverlay();

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

    void ConnectAccumWaveform(std::shared_ptr<AccumWaveform> accumWF) { fAccumWaveform = accumWF; }

  private:

    //private methods

  private:

    //private data members
    std::shared_ptr<AccumWaveform> fAccumWaveform;

    std::shared_ptr<CCMRootIO> fStrobeData;

    long int fTotalTriggers;
    long int fNumOverLap;
};

#endif // CCMStrobeOverlay_h

