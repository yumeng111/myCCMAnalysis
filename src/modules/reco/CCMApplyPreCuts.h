#ifndef CCMApplyPreCuts_h
#define CCMApplyPreCuts_h

#include "Utility.h"
#include "CCMModule.h"

class Events;
class SimplifiedEvent;
class TH1D;


class CCMApplyPreCuts : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMApplyPreCuts(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMApplyPreCuts(const CCMApplyPreCuts& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMApplyPreCuts();

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

    void ConnectEvents(std::shared_ptr<Events> evt) { fEvents = evt; }

  private:

    //private methods
    void RecalculateStartTime(const SimplifiedEvent & events, double & st, double & charge, double & hits, double & length);
    void SetupHists();

  private:

    //private data members
    std::string fOutFileName;

    std::shared_ptr<Events> fEvents;

    std::vector<TH1D*> fTimeHist;
    TH1D* fPreEventHist;
    TH1D* fBeamWindow1Hist;
    TH1D* fBeamWindow2Hist;
    TH1D* fBeamWindow3Hist;
    TH1D* fBeamWindow4Hist;
    TH1D* fBeamWindow5Hist;
    TH1D* fBeamWindow6Hist;
    TH1D* fBeamWindow7Hist;
    TH1D* fBeamWindow8Hist;
    
};

#endif // CCMApplyPreCuts_h

