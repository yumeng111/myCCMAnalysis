#ifndef CCMFindEvents_h
#define CCMFindEvents_h

#include <vector>
#include <iterator>
#include <algorithm>

#include "CCMAnalysis/modules/framework/CCMModule.h"
#include "CCMAnalysis/utils/Utility.h"

class CCMFindEvents : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMFindEvents(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMFindEvents(const CCMFindEvents& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMFindEvents();

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

    void ConnectAccumWaveform(std::shared_ptr<AccumWaveform> accum) { fAccumWaveform = accum; }
    void ConnectEvents(std::shared_ptr<Events> evt) { fEvents = evt; }

    void SetThreshold(double threshold) { fThreshold = threshold; }
    void SetEventFinderID(CCMEventFinderID_t evtFinder) { fEventFinderID = evtFinder; }
    void SetAccumWaveformMethodID(CCMAccumWaveformMethod_t method) {fAccumWaveformMethodID = method; }
    void SetResetEvents(bool flag) { fResetEvents = flag; }
    void SetFixedLength(int length) { fFixedLength = length; }

  private:

    //private methods

    // TODO Try to implement this correctly
    //int ExtrapolateStartTime(int startBin);

    void SaveEvent(int startBin, int endBin);

  private:

    //private data members
    double fThreshold;

    std::shared_ptr<AccumWaveform> fAccumWaveform;
    std::shared_ptr<Events> fEvents;

    CCMEventFinderID_t fEventFinderID;
    CCMAccumWaveformMethod_t fAccumWaveformMethodID;

    unsigned long int fNumTriggers;
    unsigned long int fNumEvents;

    int fResetEvents;
    int fFixedLength;
};

#endif // CCMFindEvents_h

