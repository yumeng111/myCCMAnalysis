/*!**********************************************
 * \file SimplifiedEvent.h
 * \brief Header file for the #SimplifiedEvent class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef SimplifiedEvent_h
#define SimplifiedEvent_h

#include "CCMAnalysis/CCMDataStructures/SinglePulse.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
#include <icetray/serialization.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3DefaultName.h>
#endif // __MAKECINT__ || __ROOTCLING__

#include <map>
#include <vector>
#include <utility>
#include <iostream>

static const unsigned legacy_simplified_event_version_ = 12;

#include "TObject.h"

//class BoardInfo;
//class ChannelInfo;


// from Utility header we get the definition of the
// different event builder types
//typedef enum {
//  kCCMDynamicLengthEventID = 0, 
//  kCCMFixedLengthEventID = 1, 
//} CCMEventFinderID_t;

// from Utility header we get the definition
// of the different accumulated waveform builder
// types
//typedef enum {
//  kCCMAccumWaveformTriangleID = 0,
//  kCCMAccumWaveformStartID = 1,
//  kCCMAccumWaveformTrianglePulseCutID = 2,
//  kCCMAccumWaveformStartPulseCutID = 3,
//  kCCMAccumWaveformTotalID = 4
//} CCMAccumWaveformMethod_t;

/*!**********************************************
 * \class SimplifiedEvent
 * \brief All the information needed after finding an event
 *
 * Not all parameters are filled with every event definition
 * but this is a good list of all the various parameters
 * one might want to calculate
 ***********************************************/
class SimplifiedEvent : public TObject
#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
    , public I3FrameObject
#endif // __MAKECINT__ || __ROOTCLING__
{
  public:
    SimplifiedEvent();
    SimplifiedEvent(const SimplifiedEvent & rhs);
    ~SimplifiedEvent();

    void Reset();

    ////////////////////////////////////
    // Set Functions
    ////////////////////////////////////
    void SetEventFinderMethod  ( CCMEventFinderID_t value  ) { fEventFinderMethod = static_cast<int>(value); }
    void SetAccumWaveformMethod  ( CCMAccumWaveformMethod_t value  ) { fAccumWaveformMethod = static_cast<int>(value); }
    void SetThreshold          ( float value  ) { fThreshold = value; }
    void SetLargestPMTFraction ( float value ) { fLargestPMTFraction = value; }

    void SetStartTime               (float value) { fStartTime = value; }
    void SetLength                  (float value) { fLength = value; }
    void SetStartTimeCoated         (float value) { fCoatedStartTime = value; }
    void SetStartTimeUncoated       (float value) { fUncoatedStartTime = value; }
    void SetStartTimeVetoBottom     (float value) { fVetoBottomStartTime = value; }
    void SetStartTimeVetoTop        (float value) { fVetoTopStartTime = value; }
    void SetStartTimeVetoBack       (float value) { fVetoBackStartTime = value; }
    void SetStartTimeVetoRight      (float value) { fVetoRightStartTime = value; }
    void SetStartTimeVetoLeft       (float value) { fVetoLeftStartTime = value; }
    void SetStartTimeVetoFront      (float value) { fVetoFrontStartTime = value; }

    // Integral variables
    void SetIntegralCoated         (float value, bool prompt = true) { fCoatedInt[prompt] = value; }
    void SetIntegralUncoated       (float value, bool prompt = true) { fUncoatedInt[prompt] = value; }

    void SetXPosition (float value) { fHitXPosition = value; }
    void SetYPosition (float value) { fHitYPosition = value; }
    void SetZPosition (float value) { fHitZPosition = value; }

    // Multiplicity variables
    void SetNumCoated     (float value, bool prompt = true) { fNCoated[prompt] = value; }
    void SetNumUncoated   (float value, bool prompt = true) { fNUncoated[prompt] = value; }
    void SetNumVetoBottom (float value, bool prompt = true) { fNVetoBottom[prompt] = value; }
    void SetNumVetoTop    (float value, bool prompt = true) { fNVetoTop[prompt] = value; }
    void SetNumVetoBack   (float value, bool prompt = true) { fNVetoBack[prompt] = value; }
    void SetNumVetoLeft   (float value, bool prompt = true) { fNVetoLeft[prompt] = value; }
    void SetNumVetoRight  (float value, bool prompt = true) { fNVetoRight[prompt] = value; }
    void SetNumVetoFront  (float value, bool prompt = true) { fNVetoFront[prompt] = value; }

    ////////////////////////////////////
    // Add Functions
    ////////////////////////////////////
    void AddWaveforms(const std::vector<float> & count, const std::vector<float> & integral);
    void AddPMTWaveforms(const int key, const std::vector<int> & count, const std::vector<float> & integral);

    const std::vector<float> & GetWaveformCount() { return fAccWaveformCount; }
    const std::vector<float> & GetWaveformInt() { return fAccWaveformInt; }
    const std::vector<float> & GetWaveformCount() const { return fAccWaveformCount; }
    const std::vector<float> & GetWaveformInt() const { return fAccWaveformInt; }

    void ResetPMTWaveformItr();
    void ResetPMTWaveformItr() const;
    bool NextPMTWaveform();
    bool NextPMTWaveform() const;

    void GetPMTWaveform(int & key, std::vector<float> & vecInt, std::vector<int> & vecCount);
    void GetPMTWaveform(int & key, std::vector<float> & vecInt, std::vector<int> & vecCount) const;

    // Integral variables
    void AddIntegralCoated           (float value, bool prompt = true) { fCoatedInt[prompt] += value; }
    void AddIntegralUncoated         (float value, bool prompt = true) { fUncoatedInt[prompt] += value; }

    // Multiplicity variables
    void AddNumCoated     (int value = 1, bool prompt = true) { fNCoated[prompt] += value; }
    void AddNumUncoated   (int value = 1, bool prompt = true) { fNUncoated[prompt] += value; }
    void AddNumVetoBottom (int value = 1, bool prompt = true) { fNVetoBottom[prompt] += value; }
    void AddNumVetoTop    (int value = 1, bool prompt = true) { fNVetoTop[prompt] += value; }
    void AddNumVetoBack   (int value = 1, bool prompt = true) { fNVetoBack[prompt] += value; }
    void AddNumVetoLeft   (int value = 1, bool prompt = true) { fNVetoLeft[prompt] += value; }
    void AddNumVetoRight  (int value = 1, bool prompt = true) { fNVetoRight[prompt] += value; }
    void AddNumVetoFront  (int value = 1, bool prompt = true) { fNVetoFront[prompt] += value; }

    ////////////////////////////////////
    // Get Functions
    ////////////////////////////////////
    CCMEventFinderID_t GetEventFinderMethod  () { return static_cast<CCMEventFinderID_t>(fEventFinderMethod); }
    CCMAccumWaveformMethod_t GetAccumWaveformMethod  () { return static_cast<CCMAccumWaveformMethod_t>(fAccumWaveformMethod); }

    float GetThreshold           ( ) { return fThreshold; }
    float GetLargestPMTFraction ( ) { return fLargestPMTFraction; }


    // Time variables
    float GetStartTime               () const { return fStartTime; }
    float GetStartTime               () { return fStartTime; }
    float GetLength                  () { return fLength; }
    float GetStartTimeCoated         () { return fCoatedStartTime; }
    float GetStartTimeUncoated       () { return fUncoatedStartTime; }
    float GetStartTimeVetoBottom     () { return fVetoBottomStartTime; }
    float GetStartTimeVetoTop        () { return fVetoTopStartTime; }
    float GetStartTimeVetoBack       () { return fVetoBackStartTime; }
    float GetStartTimeVetoRight      () { return fVetoRightStartTime; }
    float GetStartTimeVetoLeft       () { return fVetoLeftStartTime; }
    float GetStartTimeVetoFront      () { return fVetoFrontStartTime; }
    float GetStartTimeTank();
    float GetStartTimeVeto();

    // Integral variables
    float GetIntegralCoated         (bool prompt = true) { return fCoatedInt[prompt]; }
    float GetIntegralUncoated       (bool prompt = true) { return fUncoatedInt[prompt]; }
    float GetIntegralTank(bool prompt = true);

    float GetXPosition () { return fHitXPosition; }
    float GetYPosition () { return fHitYPosition; }
    float GetZPosition () { return fHitZPosition; }

    // Multiplicity variables
    float    GetNumCoated     (bool prompt = true) { return fNCoated[prompt]; }
    float    GetNumUncoated   (bool prompt = true) { return fNUncoated[prompt]; }
    float    GetNumVetoBottom (bool prompt = true) { return fNVetoBottom[prompt]; }
    float    GetNumVetoTop    (bool prompt = true) { return fNVetoTop[prompt]; }
    float    GetNumVetoBack   (bool prompt = true) { return fNVetoBack[prompt]; }
    float    GetNumVetoLeft   (bool prompt = true) { return fNVetoLeft[prompt]; }
    float    GetNumVetoRight  (bool prompt = true) { return fNVetoRight[prompt]; }
    float    GetNumVetoFront  (bool prompt = true) { return fNVetoFront[prompt]; }
    float    GetNumTank(bool prompt = true);
    float    GetNumVeto(bool prompt = true);
    float    GetNumVetoSide(bool prompt = true);

    void  SetPMTHits (std::vector<std::pair<int,float>> & vec) { fPMTHits = vec.size(); fPMTHitsVec = vec; }
    void  AddPMTHits (std::pair<int,float> amount) { ++fPMTHits; fPMTHitsVec.emplace_back(amount); }
    int   GetPMTHits () { return fPMTHits; }
    const std::vector<std::pair<int,float>> & GetPMTHitsVec () { return fPMTHitsVec; }

    void  SetPMTHits20 (std::vector<std::pair<int,float>> & vec) { fPMTHits20 = vec.size(); fPMTHits20Vec = vec; }
    void  AddPMTHits20 (std::pair<int,float> amount) { ++fPMTHits20; fPMTHits20Vec.emplace_back(amount); }
    int   GetPMTHits20 () { return fPMTHits20; }
    const std::vector<std::pair<int,float>> & GetPMTHits20Vec () { return fPMTHits20Vec; }

    void  SetPMTHits40 (std::vector<std::pair<int,float>> & vec) { fPMTHits40 = vec.size(); fPMTHits40Vec = vec; }
    void  AddPMTHits40 (std::pair<int,float> amount) { ++fPMTHits40; fPMTHits40Vec.emplace_back(amount); }
    int   GetPMTHits40 () { return fPMTHits40; }
    const std::vector<std::pair<int,float>> & GetPMTHits40Vec () { return fPMTHits40Vec; }

    void  SetPMTHits60 (std::vector<std::pair<int,float>> & vec) { fPMTHits60 = vec.size(); fPMTHits60Vec = vec; }
    void  AddPMTHits60 (std::pair<int,float> amount) { ++fPMTHits60; fPMTHits60Vec.emplace_back(amount); }
    int   GetPMTHits60 () { return fPMTHits60; }
    const std::vector<std::pair<int,float>> & GetPMTHits60Vec () { return fPMTHits60Vec; }

    void  SetPMTHitsStart (std::vector<std::pair<int,float>> & vec) { fPMTHitsStartVec = vec; }
    void  AddPMTHitsStart (std::pair<int,float> amount) { fPMTHitsStartVec.emplace_back(amount); }
    const std::vector<std::pair<int,float>> & GetPMTHitsStartVec () { return fPMTHitsStartVec; }

    float GetMaxAccumWaveformTime() { return fMaxWaveformTime; }
    float GetMaxAccumWaveformValue() { return fMaxWaveformValue; }
    void SetMaxAccumWaveformTime(float time) { fMaxWaveformTime = time; }
    void SetMaxAccumWaveformValue(float value) { fMaxWaveformValue = value; }

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive & ar, unsigned version);
#endif // __MAKECINT__ || __ROOTCLING__

  protected:
    int fEventFinderMethod;        ///< the event finder method converted from #CCMEventFinderID_t
    int fAccumWaveformMethod;        ///< the method used to build the accumulated waveform #CCMAccumWaveformMethod_t
    float fThreshold;                    ///< threshold used for finding the event
    int fPMTHits;                    ///< number of PMTs with charge in the prompt region
    int fPMTHits20;                    ///< number of PMTs with charge in the first 20ns
    int fPMTHits40;                    ///< number of PMTs with charge in the first 40ns
    int fPMTHits60;                    ///< number of PMTs with charge in the first 60ns
    std::vector<std::pair<int,float>> fPMTHitsVec; ///< amount of charge in each PMT for the prompt region
    std::vector<std::pair<int,float>> fPMTHits20Vec; ///< amount of charge in each PMT for the first 20 ns
    std::vector<std::pair<int,float>> fPMTHits40Vec; ///< amount of charge in each PMT for the first 40 ns
    std::vector<std::pair<int,float>> fPMTHits60Vec; ///< amount of charge in each PMT for the first 60 ns
    std::vector<std::pair<int,float>> fPMTHitsStartVec; ///< time when first light was seen by each PMT
    float fLargestPMTFraction;          ///< fraction of PMT with the most charge


    // Time variables
    float fStartTime;               ///< time of event in mus
    float fLength;                  ///< length of event
    float fCoatedStartTime;         ///< time of first coated tank event
    float fUncoatedStartTime;       ///< time of first uncoated event
    float fVetoBottomStartTime;     ///< time of first veto on the bottom
    float fVetoTopStartTime;        ///< time of first veto on the top
    float fVetoBackStartTime;       ///< time of first veto on the back
    float fVetoRightStartTime;      ///< time of first veto on the right
    float fVetoLeftStartTime;       ///< time of first veto on the left
    float fVetoFrontStartTime;      ///< time of first veto on the front

    float fMaxWaveformTime;         ///< time of where the accumulated waveform is at its maximum
    float fMaxWaveformValue;        ///< value corresponding to when the accumulated waveform is at its maximum

    // Integral variables
    // The integral is calculated for the event only if the PMT was above threshold
    // independent of the value of the integral or when the PMT crossed threshold
    float fCoatedInt[2];         ///< of coated tubes
    float fUncoatedInt[2];       ///< of uncoated tubes

    float fHitXPosition; ///< of tank tubes (coated + uncoated)
    float fHitYPosition; ///< of tank tubes (coated + uncoated)
    float fHitZPosition; ///< of tank tubes (coated + uncoated)

    // Multiplicity variables
    float fNCoated[2];     ///< for coated tank tubes
    float fNUncoated[2];   ///< for uncoated tank tubes
    float fNVetoBottom[2]; ///< for 8" bottom veto tubes
    float fNVetoTop[2];    ///< for 1" top veto tubes
    float fNVetoBack[2];   ///< for 1" cylinder veto tubes (columes 22 through 3)
    float fNVetoLeft[2];   ///< for 1" cylinder veto tubes (columes 4 through 9)
    float fNVetoRight[2];  ///< for 1" cylinder veto tubes (columes 16 through 21)
    float fNVetoFront[2];  ///< for 1" cylinder veto tubes (columes 10 through 15)

    /// Vector of the single pulses that made up the event
    std::vector<float> fAccWaveformCount; //
    std::vector<float> fAccWaveformInt; //
    std::map<int,std::pair<std::vector<float>,std::vector<int>>> fPMTWaveform; //!
    mutable std::map<int,std::pair<std::vector<float>,std::vector<int>>>::const_iterator fPMTWaveformItr; //! do not save to file

  ClassDef(SimplifiedEvent,legacy_simplified_event_version_)  //SimplifiedEvent class
};

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
I3_DEFAULT_NAME(SimplifiedEvent, LegacySimplifiedEvent);
I3_POINTER_TYPEDEFS(SimplifiedEvent);
I3_CLASS_VERSION(SimplifiedEvent, legacy_simplified_event_version_);
#endif // __MAKECINT__ || __ROOTCLING__

#endif
