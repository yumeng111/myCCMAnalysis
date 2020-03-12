#ifndef CCMApplyPreCuts_h
#define CCMApplyPreCuts_h

#include "Utility.h"
#include "CCMModule.h"

class SimplifiedEvent;


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
    CCMResult_t ProcessEvent();

    /*!
     *  \brief Returns true of the job has ended
     *  \return CCMResult_t result if the Job had ended
     */
    CCMResult_t EndOfJob() { return kCCMSuccess;}

    /*!
     *  \brief Configures things that are hardware specific 
     *  \param c holds hardware and event specific settings
     */
    void Configure(const CCMConfig& c);

  private:

    //private methods
    void RecalculateStartTime(SimplifiedEvent * events, double & st, double & charge, double & hits, double & length);

  private:

    //private data members
    std::string fOutFileName;
    std::string fInFileName;
    
};

#endif // CCMApplyPreCuts_h

