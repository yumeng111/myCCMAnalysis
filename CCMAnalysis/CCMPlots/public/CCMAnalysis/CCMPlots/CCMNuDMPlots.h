#ifndef CCMNuDMPlots_h
#define CCMNuDMPlots_h

#include <map>
#include <ctime>
#include <unordered_map>

#include "CCMAnalysis/CCMFramework/CCMModule.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

#include "THnSparse.h"

class TFile;
class TH1D;
class TTree;
class TH1;
class TH2;


class CCMNuDMPlots : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMNuDMPlots(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMNuDMPlots(const CCMNuDMPlots& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMNuDMPlots();

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

    void DefineHists();
    void Get(TTree * tree);
    void Fill(std::string name, double par0, double par1 = 1.0, double par2 = 1.0);


  private:

    //private methods
    void SetupOutFile();
    void SetupInputFile();

  private:

    //private data members
    std::string fOutFileName;
    std::string fInFileName;
    std::string fPrevInFileName;

    // variables for the input and output files
    TFile * fOutfile;
    TFile * fInfile;

    static const size_t fgkNumVetoRegions = 6;
    static const size_t fgkNumEnergyRegions = 4;
    static const size_t fgkNumNicenessRegions = 6;

    size_t fNumTimeRegions;
    std::string fTitle;
    std::string fName;

    typedef std::unordered_map<std::string,std::shared_ptr<TH1>> MapStringHist1D;
    typedef std::unordered_map<std::string,std::shared_ptr<TH2>> MapStringHist2D;

    MapStringHist1D fHists1D;
    MapStringHist2D fHists2D;
    std::vector<std::string> fCutNames;

    std::vector<double> fTimeWindows;

    std::string fTreeName;
    bool fIncludeBCM;

    double fTotalPOT;
    double fTotalTriggers;
    constexpr static const double fgkBCMIntToPOT = 1351271580.158518;
    constexpr static const double fgkBCMIntToPOTRMSPerc = 0.007805;

};

#endif // CCMNuDMPlots_h


