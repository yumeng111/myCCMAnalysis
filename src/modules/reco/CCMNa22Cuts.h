#ifndef CCMNa22Cuts_h
#define CCMNa22Cuts_h

#include "Utility.h"
#include "CCMModule.h"

class TTree;
class TFile;
class TProfile;
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

    void SetDoBCMCut(bool flag) {fDoBCMCut = flag;}
    void SetDoVetoCut(bool flag) {fDoVetoCut = flag;}
    void SetDoPositionCut(bool flag) {fDoPositionCut = flag;}
    void SetDoEnergyCut(bool flag) {fDoEnergyCut = flag;}
    void SetDoPrevCut(bool flag) {fDoPrevCut = flag;}
    void SetDoLengthCut(bool flag) {fDoLengthCut = flag;}
    void SetDoTimeCut(bool flag) {fDoTimeCut = flag;}
    void SetDoNicenessCut(bool flag) {fDoNicenessCut = flag;}
    void SetReverseVetoCut(bool flag) {fReverseVeto = flag;}
    void SetDoWaveformMaxCut(bool flag) {fWaveformMaxCut = flag;}

    void SetBCMTimeLowCut( double value) { fBCMTimeLowCut = value; }
    void SetBCMTimeHighCut( double value) { fBCMTimeHighCut = value; }
    void SetBCMLengthLowCut( double value) { fBCMLengthLowCut = value; }
    void SetBCMLengthHighCut( double value) { fBCMLengthHighCut = value; }
    void SetBCMIntegralLowCut( double value) { fBCMIntegralLowCut = value; }
    void SetBCMIntegralHighCut( double value) { fBCMIntegralHighCut = value; }
    void SetNumVetoCut(int value) {fNumVetoCut = value;}
    void SetPrevCutTime(double value) {fPrevCutTime = value;}
    void SetPrevCutHoldOff(double value) {fPrevCutHoldOff = value;}
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

    void SetEventFinderID(CCMEventFinderID_t evtFinder) { fEventFinderID = evtFinder; }
    void SetAccumWaveformMethodID(CCMAccumWaveformMethod_t method) {fAccumWaveformMethodID = method; }
    void SetTreeName(std::string treeName) { fTreeName = treeName; }
    void SetRemoveOtherEventsFlag(bool flag) { fRemoveOtherEvents = flag;}
    void SetRemovePrimaryEventsFlag(bool flag) { fRemovePrimaryEvents = flag;}

    bool PassedVetoCut(int vetoTotal);
    bool PassedPositionCut(double x, double y, double z);
    bool PassedEnergyCut(double energy);
    bool PassedPrevCut(const double kStartTime, const long kStartingIndex, const double kEnergy);
    bool PassedLengthCut(double length);
    bool PassedTimeCut(double time);
    bool PassedNicenessCut(double energy, int hits, double largestFrac);
    bool PassedWaveformMaxCut(double length, double maxWFTime);
    bool PassedNumHitsCut(double hits);
    bool PassedBCMCut(double bcmTime, double bcmLength, double bcmIntegral);

  private:

    //private methods
    void SetupOutFile();
    void RecalculatePosition(const SimplifiedEvent & simplifiedEvent, 
        const double fitLength, double & x, double & y, double & z);
    void RecalculateStartTime(const SimplifiedEvent & events, 
        double & st, double & charge, double & hits, double & length);
    //int WaveformMaxPosition(const SimplifiedEvent & simplifiedEvent);

  private:

    //private data members

    std::shared_ptr<Events> fEvents;

    CCMEventFinderID_t fEventFinderID;
    CCMAccumWaveformMethod_t fAccumWaveformMethodID;

    std::string fOutFileName;
    std::string fTreeName;

    int fRemovePrimaryEvents;
    int fRemoveOtherEvents;
    int fDoBCMCut;
    int fDoVetoCut;
    int fDoPositionCut;
    int fDoEnergyCut;
    int fDoPrevCut;
    int fDoLengthCut;
    int fDoTimeCut;
    int fDoNicenessCut;
    int fReverseVeto;
    int fWaveformMaxCut;

    int fPassedVetoCut;
    int fPassedPositionCut;
    int fPassedEnergyCut;
    int fPassedPrevCut;
    int fPassedLengthCut;
    int fPassedTimeCut;
    int fPassedNicenessCut;
    int fPassedWaveformMaxCut;
    int fPassedNumHitsCut;
    int fPassedBCMCut;

    unsigned long int fNumFailedVetoCut;
    unsigned long int fNumFailedPositionCut;
    unsigned long int fNumFailedEnergyCut;
    unsigned long int fNumFailedPrevCut;
    unsigned long int fNumFailedLengthCut;
    unsigned long int fNumFailedTimeCut;
    unsigned long int fNumFailedNicenessCut;
    unsigned long int fNumFailedWaveformMaxCut;
    unsigned long int fNumFailedNumHitsCut;
    unsigned long int fNumFailedBCMCut;

    int fNumVetoCut;
    double fBCMTimeLowCut;
    double fBCMTimeHighCut;
    double fBCMLengthLowCut;
    double fBCMLengthHighCut;
    double fBCMIntegralLowCut;
    double fBCMIntegralHighCut;
    double fRadiusCutValueLow;
    double fRadiusCutValueHigh;
    double fZCutValueLow;
    double fZCutValueHigh;
    double fEnergyCutValueLow;
    double fEnergyCutValueHigh;
    double fPrevCutTime;
    double fPrevCutHoldOff;
    double fPrevCutEnergyFrac;
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
    unsigned int fEpochNSSec;
    double fBCMTime;
    double fBCMLength;
    double fBCMIntegral;
    double fMaxPrevEnergy1600;
    double fMaxPrevEnergy3200;
    double fMaxPrevEnergy4800;
    double fMaxPrevEnergy;
    double fMaxPrevEnergyTime1600;
    double fMaxPrevEnergyTime3200;
    double fMaxPrevEnergyTime4800;
    double fMaxPrevEnergyTime;
    double fMaxPrevEnergyLength1600;
    double fMaxPrevEnergyLength3200;
    double fMaxPrevEnergyLength4800;
    double fMaxPrevEnergyLength;

    std::shared_ptr<TFile> fEnergyLengthFile;
    TProfile* fEnergyLengthProf;

    unsigned long int fNumInitEvents;
    unsigned long int fNumFinalEvents;
};

#endif // CCMNa22Cuts_h

