#ifndef CCMNa22Cuts_h
#define CCMNa22Cuts_h

#include "Utility.h"
#include "CCMModule.h"

class TTree;
class TFile;
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
    void ConnectOutFileName(std::string name) { fOutFileName = name; SetupOutFile(); }

    void SetDoVetoCut(bool flag) {fDoVetoCut = flag;}
    void SetDoPositionCut(bool flag) {fDoPositionCut = flag;}
    void SetDoEnergyCut(bool flag) {fDoEnergyCut = flag;}
    void SetDoPrevCut(bool flag) {fDoPrevCut = flag;}
    void SetDoLengthCut(bool flag) {fDoLengthCut = flag;}
    void SetDoTimeCut(bool flag) {fDoTimeCut = flag;}
    void SetDoNicenessCut(bool flag) {fDoNicenessCut = flag;}
    void SetReverseVetoCut(bool flag) {fReverseVeto = flag;}
    void SetDoWaveformMaxCut(bool flag) {fWaveformMaxCut = flag;}

    void SetNumVetoCut(int value) {fNumVetoCut = value;}
    void SetPrevCutTime(double value) {fPrevCutTime = value;}
    void SetZCutValueLow(double value) {fZCutValueLow = value;}
    void SetZCutValueHigh(double value) {fZCutValueHigh = value;}
    void SetTimeCutValueLow(double value) {fTimeCutValueLow = value;}
    void SetTimeCutValueHigh(double value) {fTimeCutValueHigh = value;}
    void SetEnergyCutValueLow(double value) {fEnergyCutValueLow = value;}
    void SetEnergyCutValueHigh(double value) {fEnergyCutValueHigh = value;}
    void SetLengthCutValueLow(double value) {fLengthCutValueLow = value;}
    void SetLengthCutValueHigh(double value) {fLengthCutValueHigh = value;}
    void SetRadiusCutValueLow(double value) {fRadiusCutValueLow = value;}
    void SetRadiusCutValueHigh(double value) {fRadiusCutValueHigh = value;}
    void SetNicenessCutValueLow(double value) {fNicenessCutValueLow = value;}
    void SetNicenessCutValueHigh(double value) {fNicenessCutValueHigh = value;}

  private:

    //private methods
    void SetupOutFile();
    void RecalculatePosition(const SimplifiedEvent & simplifiedEvent, 
        const double fitLength, double & x, double & y, double & z);
    int WaveformMaxPosition(const SimplifiedEvent & simplifiedEvent);

  private:

    //private data members

    std::shared_ptr<Events> fEvents;

    std::string fOutFileName;

    int fDoVetoCut;
    int fDoPositionCut;
    int fDoEnergyCut;
    int fDoPrevCut;
    int fDoLengthCut;
    int fDoTimeCut;
    int fDoNicenessCut;
    int fReverseVeto;
    int fWaveformMaxCut;

    bool fPassedVetoCut;
    bool fPassedPositionCut;
    bool fPassedEnergyCut;
    bool fPassedPrevCut;
    bool fPassedLengthCut;
    bool fPassedTimeCut;
    bool fPassedNicenessCut;
    bool fPassedWaveformMaxCut;
    bool fPassedNumHitsCut;

    int fNumVetoCut;
    double fRadiusCutValueLow;
    double fRadiusCutValueHigh;
    double fZCutValueLow;
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
    int fVetoPromptTop;
    int fVetoPromptBottom;
    int fVetoPromptCLeft;
    int fVetoPromptCRight;
    int fVetoPromptCFront;
    int fVetoPromptCBack;
    double fWeight;
    double fX;
    double fY;
    double fZ;
    double fLargestPMTFraction;
    unsigned int fEpochSec;
};

#endif // CCMNa22Cuts_h

