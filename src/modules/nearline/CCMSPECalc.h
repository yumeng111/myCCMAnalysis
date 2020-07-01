#ifndef CCMSPECalc_h
#define CCMSPECalc_h

#include "Utility.h"
#include "CCMModule.h"

#include <memory>


class RawData;
class Pulses;
class NearlineSPEDiag;

class CCMSPECalc : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMSPECalc(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMSPECalc(const CCMSPECalc& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMSPECalc();

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

    void ConnectRawData(std::shared_ptr<RawData> rawData)  { fRawData = rawData; }
    void ConnectPulses(std::shared_ptr<Pulses> pulses) {fPulses = pulses; }

  private:

    //private methods

  private:

    //private data members
    std::shared_ptr<RawData> fRawData;
    std::shared_ptr<Pulses> fPulses;

    std::string fTriggerType;
    std::string fCalibrationFile;
    int fRedoLEDCalib;
    double fDAQWindowStart;
    double fDAQWindowEnd;
    std::string fSaveParameters;

    std::unique_ptr<NearlineSPEDiag> fSPEFinder;
    
};

#endif // CCMSPECalc_h

