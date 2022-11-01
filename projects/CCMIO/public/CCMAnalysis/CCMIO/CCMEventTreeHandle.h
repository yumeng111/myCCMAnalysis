#ifndef CCMEEVENTTREEHANDLE_H
#define CCMEEVENTTREEHANDLE_H
/*----------------------------------------------------------
 *
 *   CCMEventTreeHandle
 *
 *     Based on NiffteEventHandle on NIFFTE
 *
 *       Provide an interface between the SBRootIO object and
 *         the SBRootEvent (ROOT streamable) object.  
 *           Allows for partial IO of various major compoents of the event
 *
 *             Adapter: R. T. Thornton (IU)
 *               Date: Dec-2011
 *
 *-----------------------------------------------------------*/

#include <list>
#include <memory>
#include <vector>
#include <stdint.h>

#include "CCMAnalysis/CCMUtils/Utility.h"

class TFile;
class TTree;
class TBranch;
class RawData;
class Pulses;
class MCTruth;
class Events;
class AccumWaveform;

class CCMEventTreeHandle
{

  public:
    CCMEventTreeHandle();
    virtual ~CCMEventTreeHandle();

    void Clear();
    int  ClearLoadFlags();
    void SetLoaded(int branchId) const;

    bool ReadTree();

    uint32_t Index() { return fIndex; }

    int SetupInputFile(TFile & f);
    int SetupOutputFile(TFile & f);
    uint32_t Advance(uint32_t n = 1);
    uint32_t Rewind(uint32_t n = 1);

    void SetReadBranches(const std::vector<std::string> & branchList) { fReadBranches = branchList; }
    void SetSaveBranches(const std::vector<std::string> & branchList) { fSaveBranches = branchList; }

    void SetCurrentEvents(const Events& event);
    void SetCurrentRawData(const RawData& rawData);
    void SetCurrentPulses(const Pulses& pulses);
    void SetCurrentMCTruth(const MCTruth& mcTruth);
    void SetCurrentAccumWaveform(const AccumWaveform& waveform);

    Events& CurrentEvents();
    RawData& CurrentRawData();
    Pulses& CurrentPulses();
    MCTruth& CurrentMCTruth();
    AccumWaveform& CurrentAccumWaveform();

    uint32_t NumOfEntries() const { return fNumOfEntries; }

    void Close();
    void Report();
    int Write();

    //const TTree * CurrentEventTree() const { return fEventsTree; }

  protected:

    //If you set split level to 0, the tree will be written in one
    //branch, which is not very efficient
    //static const int  fgkBucketSize = 320000; //Bytes
    static const int  fgkSplitLevel = 99;


    //Define ID codes for the different branches in a file
    bool IsLoaded(CCMEventBranchID_t branchID) const;
    int Load(int branchID) const;

  protected:

    uint32_t fIndex; ///< Location in current input file
    uint32_t fNumOfEntries; ///< number of entries in the tree

    TFile * fInputFile;  ///< Pointer to input source.  Not owned.
    TFile * fOutputFile; ///< Pointer to output file.  Not owned.
    TBranch* fBranch[kNEventBranch]; ///< pointer to each branch

    TTree * fEventsTree; ///< The input event tree (also rawData tree backwards compatibility)
    TTree * fPulsesTree; ///< The input pulses tree (backwards compatibility)
    TTree * fOutEventTree; ///< The output event tree

    int fObjectCnt;
    mutable bool fIsLoaded[kNEventBranch];  ///< Which parts of the event are loaded?

    RawData * fRawData;
    Pulses * fPulses;
    Events * fEvents;
    AccumWaveform * fAccumWaveform;
    MCTruth * fMCTruth;

    std::vector<std::string> fReadBranches;
    std::vector<std::string> fSaveBranches;
};

//////////////////////////////////////////////////////////////
// Template functions must be inlines
//////////////////////////////////////////////////////////////

#endif // CCMEEVENTTREEHANDLE_H

