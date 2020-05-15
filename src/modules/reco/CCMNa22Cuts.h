#ifndef CCMNa22Cuts_h
#define CCMNa22Cuts_h

#include "Utility.h"
#include "CCMModule.h"

class TTree;
class TFile;


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
    CCMResult_t EndOfJob();

    /*!
     *  \brief Configures things that are hardware specific 
     *  \param c holds hardware and event specific settings
     */
    void Configure(const CCMConfig& c);

    void ConnectEvents(std::shared_ptr<Events> evt) { fEvents = evt; }

  private:

    //private methods
    void SetupOutFile();

  private:

    //private data members

    std::shared_ptr<Events> fEvents;

    std::string fOutFileName;

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

    // variables for the output tree
    TFile * fOutfile;
    TTree * fTree;
    double fEnergy;
    double fLength;
    double fHits;
    int fNumPMTs;
    double fTime;
    int fVetoTop;
    int fVetoBottom;
    int fVetoCLeft;
    int fVetoCRight;
    int fVetoCFront;
    int fVetoCBack;
    double fWeight;
    double fX;
    double fY;
    double fZ;
    double fLargestPMTFraction;
    unsigned int fEpochSec;

    
};

#endif // CCMNa22Cuts_h

