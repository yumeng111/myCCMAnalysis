#ifndef CCMNa22Cuts_h
#define CCMNa22Cuts_h

#include "Utility.h"
#include "CCMModule.h"

class SimplifiedEvent;


class CCMNa22Cuts : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMNa22Cuts(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMNa22Cuts(const CCMNa22Cuts& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMNa22Cuts();

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
    std::vector<std::string> fInFileNames;
    int fDoVetoCut;
    int fDoPositionCut;
    int fDoEnergyCut;
    int fDoPrevCut;
    int fShiftTime;
    int fDoLengthCut;
    int fDoTimeCut;
    int fDoNicenessCut;
    int fReverseVeto;

    int fNumVetoCut;
    double fRadiusCutValueLow;
    double fZCutValueLow;
    double fRadiusCutValueHigh;
    double fZCutValueHigh;
    double fEnergyCutValueLow;
    double fEnergyCutValueHigh;
    double fPrevCutTime;
    double fLengthCutValueLow;
    double fLengthCutValueHigh;
    double fTimeCutValueLow;
    double fTimeCutValueHigh;
    double fNicenessCutValueLow;
    double fNicenessCutValueHigh;
    
};

#endif // CCMNa22Cuts_h

