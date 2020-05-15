#ifndef CCMRateCalc_h
#define CCMRateCalc_h

#include "Utility.h"
#include "CCMModule.h"

#include <array>


class CCMRateCalc : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMRateCalc(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMRateCalc(const CCMRateCalc& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMRateCalc();

    /*!
     *  \brief This is where the action takes place
     *  \return CCMResult_t the result of running this module
     */
    CCMResult_t ProcessEvent();

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
    void ConnectRawData(std::shared_ptr<RawData> rawData)  { fRawData = rawData; }
    void ConnectPulses(std::shared_ptr<Pulses> pulses) {fPulses = pulses; }

  private:

    //private methods
    void AddPulses();
    void AddEvents(bool shiftTime = true);

  private:

    //private data members
    std::shared_ptr<RawData> fRawData;
    std::shared_ptr<Pulses> fPulses;
    std::shared_ptr<Events> fEvents;

    std::string fTriggerType;
    bool fWriteDBEntry;
    std::string fDBHost;
    std::string fDBUser;
    std::string fDBPwd;

    ///spe count for each of the 160 pmts in the first 10 boards. use as fSPECount[key] in the pulse loop.       
    std::array<int,160> fSPECount;
    /// pmt type (0=1inveto,1=uncoated,2=coated) for the 160 pmts. use as fPMTType[key] in pulse loop.              
    std::array<int,160> fPMTType;

    std::time_t fFirstTriggerTime;
    std::time_t fLastTriggerTime;

    size_t fTotalTriggers;
    size_t fPreBeamTriggers;
    float fInBeamIntegral;
    
};

#endif // CCMRateCalc_h

