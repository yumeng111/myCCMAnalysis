#ifndef CCMFindEvents_h
#define CCMFindEvents_h

#include "Utility.h"
#include "CCMModule.h"



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

  private:

    //private data members
    std::string fTriggerType;
    std::string fHVOffList;
    std::string fCalibrationFile;
    std::string fOutFileName;
    std::string fInFileName;
    
};

#endif // CCMFindEvents_h

