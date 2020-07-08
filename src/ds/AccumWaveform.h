/*!**********************************************
 * \file AccumWaveform.h
 * \brief Header file for the #AccumWaveform class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef AccumWaveform_h
#define AccumWaveform_h

#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <algorithm>
#include <iterator>

#include "Utility.h"
#include "MsgLog.h"

#include "TObject.h"

/*!**********************************************
 * \class AccumWaveform
 * \brief Container of the accumulated waveforms built for a given trigger
 *
 * Container for the accumulated waveforms built for a given trigger,
 * The class is saved in the data file.
 ***********************************************/
class AccumWaveform : public TObject
{
  public:
    AccumWaveform();
    AccumWaveform(const AccumWaveform & p);
    ~AccumWaveform();

    void Reset();
    void ClearAccumWaveform();

    /// \brief Set the event (trigger) number of the DAQ window currently looking at
    /// \param[in] value The event number
    void SetEventNumber(unsigned int value) { fEventNumber = value; }

    /// \brief Set the computer time of the DAQ window in s
    /// \param[in] value The computer time in s
    void SetComputerSecIntoEpoch(unsigned int value) { fComputerSecIntoEpoch = value; }

    /// \brief Set the computer time of the DAQ window in ns
    /// \param[in] value The computer time in ns
    void SetComputerNSIntoSec(unsigned int value) { fComputerNSIntoSec = value; }

    /// \return #fEventNumber
    unsigned int GetEventNumber() const { return fEventNumber; }
    
    /// \return #fComputerSecIntoEpoch
    unsigned int GetComputerSecIntoEpoch() const { return fComputerSecIntoEpoch; }
    
    /// \return #fComputerNSIntoSec
    unsigned int GetComputerNSIntoSec() const { return fComputerNSIntoSec; }

    /// \brief Sets the time offset of the BCM signal (0 if not BEAM trigger)
    void SetBeamOffset(int time) { fBeamTime = time; }
    
    /// \brief Get the time offset of the BCM signal
    int GetBeamOffset() { return fBeamTime; }

    
    /// \return #fEventNumber
    unsigned int GetEventNumber() { return fEventNumber; }
    
    /// \return #fComputerSecIntoEpoch
    unsigned int GetComputerSecIntoEpoch() { return fComputerSecIntoEpoch; }
    
    /// \return #fComputerNSIntoSec
    unsigned int GetComputerNSIntoSec() { return fComputerNSIntoSec; }

    float GetTriggerTime() { return fTriggerTime; }

    float GetTriggerTime() const { return fTriggerTime; }

    void SetIndex(size_t index, float weight, CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);
    void FillIndex(size_t index, float weight, CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);
    float Integrate(size_t start, size_t end, CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);
    float GetIndex(size_t index, CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);
    int FindFirstNoneEmptyBin(size_t start, size_t end, CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);

    template <class T>
    void CopyVec(typename std::vector<T> & outVec, size_t start, size_t end, 
        CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);

    template <class T>
      void Max(size_t & loc, T & value, size_t start, size_t end, 
          CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);
    template <class T>
      void Min(size_t & loc, T & value, size_t start, size_t end, 
          CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);

    AccumWaveform & operator=(const AccumWaveform & rhs);

    void DumpInfo();
    void AddMethod(CCMAccumWaveformMethod_t method);

  private:
    std::array<float,Utility::fgkNumBins> * Get(CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);

  private:
    /// The event number
    unsigned int fEventNumber;
    /// The time of the window in s
    unsigned int fComputerSecIntoEpoch;
    /// The time of the window in ns
    unsigned int fComputerNSIntoSec;

    int fBeamTime;

    /// The trigger time of the trigger
    float fTriggerTime;

    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fPulsesTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fIntegralTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fIntegralDer;

    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoBottomTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoTopTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoCRightTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoCLeftTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoCFrontTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoCBackTime;
    std::map<CCMAccumWaveformMethod_t,std::array<float,Utility::fgkNumBins>> fVetoTotalTime;

    std::map<CCMAccumWaveformMethod_t,std::vector<std::array<float,Utility::fgkNumBins>>> fPMTWaveform;
    std::map<CCMAccumWaveformMethod_t,std::vector<std::array<float,Utility::fgkNumBins>>> fPMTWaveformCount;

  ClassDef(AccumWaveform,1)

};

//-------------------------------------------------------------------------------------------------
template<class T>
void AccumWaveform::Max(size_t & loc, T & value, size_t start, size_t end, 
    CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  auto it = std::max_element(waveform->begin()+start, waveform->begin()+end);
  loc = std::distance(waveform->begin(),it);
  value = *it;
  return;
}

//-------------------------------------------------------------------------------------------------
template<class T>
void AccumWaveform::Min(size_t & loc, T & value, size_t start, size_t end,
    CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  auto it = std::min_element(waveform->begin()+start, waveform->begin()+end);
  loc = std::distance(waveform->begin(),it);
  value = *it;
  return;
}

//-------------------------------------------------------------------------------------------------
template<class T>
void AccumWaveform::CopyVec(typename std::vector<T> & outVec, size_t start, size_t end,
    CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID)
{
  if (outVec.size() != end-start+1) {
    outVec.resize(end-start+1);
  }
  auto waveform = Get(method,waveformType,pmtID);
  std::copy(waveform->begin()+start, waveform->begin()+end,outVec.begin());
  return;
}

#endif // #ifndef AccumWaveform_h

