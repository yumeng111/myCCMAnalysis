#ifndef CCMNearlineDiag_h
#define CCMNearlineDiag_h

#include "Utility.h"
#include "CCMModule.h"



class CCMNearlineDiag : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMNearlineDiag(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMNearlineDiag(const CCMNearlineDiag& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMNearlineDiag();

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
    void CalculateRates(const std::vector<std::string> & fileNames);
    void LEDCalib(const std::vector<std::string> & fileNames);

  private:

    //private data members
    std::string fTriggerType;
    std::string fHVOffList;
    std::string fCalibrationFile;
    bool fDoLEDCalib;
    bool fRedoLEDCalib;
    bool fCalculateRates;
    bool fWriteDBEntry;
    std::string fDBHost;
    std::string fDBUser;
    std::string fDBPwd;
    std::string fInFileName;
    double fDAQWindowStart;
    double fDAQWindowEnd;
    std::string fSaveParameters;
    
};

#endif // CCMNearlineDiag_h

