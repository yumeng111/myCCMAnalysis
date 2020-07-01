#ifndef CCMFindPulses_h
#define CCMFindPulses_h

#include "Utility.h"
#include "CCMModule.h"
#include "data_structures.hh"

class Pulses;
class RawData;

#include <vector>
#include <map>

class CCMFindPulses : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMFindPulses(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMFindPulses(const CCMFindPulses& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMFindPulses();

    /*!
     *  \brief This is where the action takes place
     *  \return CCMResult_t the result of running this module
     */
    CCMResult_t ProcessTrigger();

    /*!
     *  \brief What to do if we start a new run
     *  \return CCMResult_t the result of if everything happened
     */
    CCMResult_t NewRun(uint32_t run, uint32_t subRun);

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

    void ConnectBinaryRawData(std::shared_ptr<RawData> rawData)  { fReadData = rawData; }
    void ConnectRawData(std::shared_ptr<RawData> rawData)  { fRawData = rawData; }
    void ConnectPulses(std::shared_ptr<Pulses> pulses) {fPulses = pulses; }

  private:

    //private methods
    void WriteDB();

  private:

    //private data members
    std::string fTriggerType;
    int fTruncateWaveform;
    int fFromRootFile;
    int fResetPulses;


    int fWriteDBEntry;
    std::string fDBHost;
    std::string fDBUser;
    std::string fDBPwd;

    std::vector<float> fTriggerTime;
    std::shared_ptr<RawData> fReadData;
    std::shared_ptr<Pulses> fPulses;
    std::shared_ptr<RawData> fRawData;

    long fNEventsTotal;
    long fNEventsSkipped;
    unsigned short fHighestTemp;
    std::time_t fFirstTriggerTime;
    std::time_t fLastTriggerTime;
    
};

#endif // CCMFindPulses_h

