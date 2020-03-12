#ifndef CCMConvertBinary2ROOT_h
#define CCMConvertBinary2ROOT_h

#include "Utility.h"
#include "CCMModule.h"
#include "data_structures.hh"

class Pulses;
class RawData;

#include <vector>
#include <map>

class CCMConvertBinary2ROOT : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMConvertBinary2ROOT(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMConvertBinary2ROOT(const CCMConvertBinary2ROOT& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMConvertBinary2ROOT();

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
    std::string fInFileName;
    std::string fOutFileName;

    bool fWriteDBEntry;
    std::string fDBHost;
    std::string fDBUser;
    std::string fDBPwd;

    event_t fReadData;
    std::vector<float> fTriggerTime;
    Pulses * fPulses;
    RawData * fRawData;
    
};

#endif // CCMConvertBinary2ROOT_h

