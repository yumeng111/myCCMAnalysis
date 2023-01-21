/*!**********************************************
 * \file AccumWaveform.h
 * \brief Header file for the #AccumWaveform class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef AccumWaveform_h
#define AccumWaveform_h

#include <map>
#include <vector>
#include <utility>
#include <iterator>
#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
#include <icetray/serialization.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3DefaultName.h>
#endif // __MAKECINT__ || __ROOTCLING__

#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"

#include "TObject.h"

static const unsigned legacy_accum_waveform_version_ = 3;

/*!**********************************************
 * \class AccumWaveform
 * \brief Container of the accumulated waveforms built for a given trigger
 *
 * Container for the accumulated waveforms built for a given trigger,
 * The class is saved in the data file.
 ***********************************************/
class AccumWaveform : public TObject
#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
    , public I3FrameObject
#endif // __MAKECINT__ || __ROOTCLING__
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

    /// \brief Sets the integral of the BCM signal (0 if not BEAM trigger)
    void SetBeamIntegral(float integral) { fBeamIntegral = integral; }
    
    /// \brief Get the integral of the BCM signal
    float GetBeamIntegral() { return fBeamIntegral; }

    /// \brief Sets the length of the BCM signal (0 if not BEAM trigger)
    void SetBeamLength(float length) { fBeamLength = length; }
    
    /// \brief Get the length of the BCM signal
    float GetBeamLength() { return fBeamLength; }
    
    /// \return #fEventNumber
    unsigned int GetEventNumber() { return fEventNumber; }
    
    /// \return #fComputerSecIntoEpoch
    unsigned int GetComputerSecIntoEpoch() { return fComputerSecIntoEpoch; }
    
    /// \return #fComputerNSIntoSec
    unsigned int GetComputerNSIntoSec() { return fComputerNSIntoSec; }

    void SetTriggerTime(float time) { fTriggerTime = time; }
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

    void Max(size_t & loc, float & value, size_t start, size_t end, 
        CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);
    void Min(size_t & loc, float & value, size_t start, size_t end, 
        CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID = 0);

    AccumWaveform & operator=(const AccumWaveform & rhs);

    AccumWaveform & operator+=(const AccumWaveform & rhs);
    void AddMaps(MapDAQWF1D & lhs, const MapDAQWF1D & rhs);
    void AddMaps(MapDAQWF2D & lhs, const MapDAQWF2D & rhs);

    void DumpInfo();
    void AddMethod(CCMAccumWaveformMethod_t method);

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive & ar, unsigned version);
#endif // __MAKECINT__ || __ROOTCLING__

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
    float fBeamIntegral;
    float fBeamLength;

    /// The trigger time of the trigger
    float fTriggerTime;

    MapDAQWF1D fPulsesTime;
    MapDAQWF1D fIntegralTime;
    MapDAQWF1D fIntegralDer;

    MapDAQWF1D fVetoBottomTime;
    MapDAQWF1D fVetoTopTime;
    MapDAQWF1D fVetoCRightTime;
    MapDAQWF1D fVetoCLeftTime;
    MapDAQWF1D fVetoCFrontTime;
    MapDAQWF1D fVetoCBackTime;
    MapDAQWF1D fVetoTotalTime;

    MapDAQWF2D fPMTWaveform;
    MapDAQWF2D fPMTWaveformCount;

  ClassDef(AccumWaveform,legacy_accum_waveform_version_)

};

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

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
I3_DEFAULT_NAME(AccumWaveform, LegacyAccumWaveform);
I3_POINTER_TYPEDEFS(AccumWaveform);
I3_CLASS_VERSION(AccumWaveform, legacy_accum_waveform_version_);
#endif // __MAKECINT__ || __ROOTCLING__

#endif // #ifndef AccumWaveform_h

