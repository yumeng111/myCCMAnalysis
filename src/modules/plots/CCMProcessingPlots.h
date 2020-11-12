#ifndef CCMProcessingPlots_h
#define CCMProcessingPlots_h

#include "Utility.h"
#include "CCMModule.h"

#include <map>
#include <ctime>

#include "THnSparse.h"

class TFile;
class TH1D;
class TTree;


class CCMProcessingPlots : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMProcessingPlots(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMProcessingPlots(const CCMProcessingPlots& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMProcessingPlots();

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

    void ConnectOutFileName(std::string name) { fOutFileName = name; SetupOutFile(); }
    void ConnectInFileName(std::string name) {fInFileName = name; SetupInputFile(); }
    void ConnectEvents(std::shared_ptr<Events> events) {fEvents = events; }

  private:

    //private methods
    void SetupOutFile();
    void SetupInputFile();

  private:

    //private data members
    std::shared_ptr<Events> fEvents;

    std::string fOutFileName;
    std::string fInFileName;
    std::string fPrevInFileName;

    // variables for the input and output files
    TFile * fOutfile;
    TFile * fInfile;
    TTree * fTree;

    // these histograms are taken from the inputfile
    std::shared_ptr<TH1D> fCCMPulsesTimeHist;
    std::shared_ptr<TH1D> fCCMPulsesTimeIntHist;

    std::shared_ptr<TH1D> fCCMVetoTimeHist;
    std::shared_ptr<TH1D> fCCMVetoTimeIntHist;

    std::shared_ptr<TH1D> fFP3TimeHist;
    std::shared_ptr<TH1D> fFP3TimeIntHist;

    // these histograms are created by looping over the events
    std::shared_ptr<TH1D> fCCMEventsTimeHist;
    std::shared_ptr<TH1D> fCCMEventsTimeIntHist;

    std::shared_ptr<THnSparseF> fBCM3DHist;

    std::shared_ptr<TH1D> fPromptTopVeto;
    std::shared_ptr<TH1D> fPromptBottomVeto;
    std::shared_ptr<TH1D> fPromptCFrontVeto;
    std::shared_ptr<TH1D> fPromptCBackVeto;
    std::shared_ptr<TH1D> fPromptCLeftVeto;
    std::shared_ptr<TH1D> fPromptCRightVeto;
    std::shared_ptr<TH1D> fPromptTotalVeto;

    // these values are copied from the CCMBeamInfo class
    std::pair<std::time_t,double> fPrevCurrentInfo;
    std::pair<std::time_t,double> fNextCurrentInfo;

    // these maps are saved in the fOutfile for plotting time plots with
    // input from other time frames
    unsigned int fCurrentTime;
    double fNumberOfTriggersInSum;
    std::vector<double> fBCMIntegralTime;
    std::vector<double> fBCMTimeTime;
    std::vector<double> fBCMWidthTime;
    std::vector<double> fPromptTopVetoTime;
    std::vector<double> fPromptBottomVetoTime;
    std::vector<double> fPromptCLeftVetoTime;
    std::vector<double> fPromptCRightVetoTime;
    std::vector<double> fPromptCFrontVetoTime;
    std::vector<double> fPromptCBackVetoTime;
    std::vector<double> fPromptTotalVetoTime;
    //double fSumBCMIntegralTime;
    //double fSumBCMTimeTime;
    //double fSumBCMWidthTime;
    //double fSumPromptTopVetoTime;
    //double fSumPromptBottomVetoTime;
    //double fSumPromptCLeftVetoTime;
    //double fSumPromptCRightVetoTime;
    //double fSumPromptCFrontVetoTime;
    //double fSumPromptCBackVetoTime;
    //double fSumPromptTotalVetoTime;


};

#endif // CCMProcessingPlots_h

