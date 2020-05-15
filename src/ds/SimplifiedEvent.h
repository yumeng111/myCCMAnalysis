/*!**********************************************
 * \file SimplifiedEvent.h
 * \brief Header file for the #SimplifiedEvent class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef SimplifiedEvent_h
#define SimplifiedEvent_h

#include "SinglePulse.h"
#include "TObject.h"

#include <iostream>
#include <vector>
#include <map>
#include <utility>

//class BoardInfo;
//class ChannelInfo;

/*!**********************************************
 * \class SimplifiedEvent
 * \brief All the information needed after finding an event
 *
 * Not all parameters are filled with every event definition
 * but this is a good list of all the various parameters
 * one might want to calculate
 ***********************************************/
class SimplifiedEvent : public TObject
{
  public :
    SimplifiedEvent();
    SimplifiedEvent(const SimplifiedEvent & rhs);
    ~SimplifiedEvent();

    void Reset();

    ////////////////////////////////////
    // Set Functions
    ////////////////////////////////////
    void SetThreshold          ( float value  ) { fThreshold = value; }
    void SetLargestPMTFraction ( double value ) { fLargestPMTFraction = value; }

    void SetStartTime               (double value) { fStartTime = value; }
    void SetLength                  (double value) { fLength = value; }
    void SetStartTimeCoated         (double value) { fCoatedStartTime = value; }
    void SetStartTimeUncoated       (double value) { fUncoatedStartTime = value; }
    void SetStartTimeVetoBottom     (double value) { fVetoBottomStartTime = value; }
    void SetStartTimeVetoTop        (double value) { fVetoTopStartTime = value; }
    void SetStartTimeVetoBack       (double value) { fVetoBackStartTime = value; }
    void SetStartTimeVetoRight      (double value) { fVetoRightStartTime = value; }
    void SetStartTimeVetoLeft       (double value) { fVetoLeftStartTime = value; }
    void SetStartTimeVetoFront      (double value) { fVetoFrontStartTime = value; }

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
    float GetThreshold           ( ) { return fThreshold; }
    double GetLargestPMTFraction ( ) { return fLargestPMTFraction; }


    // Time variables
    double GetStartTime               () const { return fStartTime; }
    double GetStartTime               () { return fStartTime; }
    double GetLength                  () { return fLength; }
    double GetStartTimeCoated         () { return fCoatedStartTime; }
    double GetStartTimeUncoated       () { return fUncoatedStartTime; }
    double GetStartTimeVetoBottom     () { return fVetoBottomStartTime; }
    double GetStartTimeVetoTop        () { return fVetoTopStartTime; }
    double GetStartTimeVetoBack       () { return fVetoBackStartTime; }
    double GetStartTimeVetoRight      () { return fVetoRightStartTime; }
    double GetStartTimeVetoLeft       () { return fVetoLeftStartTime; }
    double GetStartTimeVetoFront      () { return fVetoFrontStartTime; }
    double GetStartTimeTank();
    double GetStartTimeVeto();

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

    void  SetPMTHits (int value, std::vector<float> & vec) { fPMTHits = value; fPMTHitsVec = vec; }
    void  AddPMTHits (int value = 1, float percent = 0) { fPMTHits += value; fPMTHitsVec.push_back(percent); }
    int   GetPMTHits () { return fPMTHits; }
    const std::vector<float> & GetPMTHitsVec () { return fPMTHitsVec; }

  protected:
    static float fgkStartDAQWindow; //!
    static float fgkSampleWidth; //!

    float fThreshold;                    //!
    int fPMTHits;                    ///< true if there is another triggerable event
    std::vector<float> fPMTHitsVec; ///< time of all the other events
    double fLargestPMTFraction;          ///< fraction of PMT with the most charge


    // Time variables
    double fStartTime;               ///< time of event in mus
    double fLength;                  ///< length of event
    double fCoatedStartTime;         ///< time of first coated tank event
    double fUncoatedStartTime;       ///< time of first uncoated event
    double fVetoBottomStartTime;     ///< time of first veto on the bottom
    double fVetoTopStartTime;        ///< time of first veto on the top
    double fVetoBackStartTime;       ///< time of first veto on the back
    double fVetoRightStartTime;      ///< time of first veto on the right
    double fVetoLeftStartTime;       ///< time of first veto on the left
    double fVetoFrontStartTime;      ///< time of first veto on the front

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
    std::map<int,std::pair<std::vector<float>,std::vector<int>>> fPMTWaveform; //
    mutable std::map<int,std::pair<std::vector<float>,std::vector<int>>>::const_iterator fPMTWaveformItr; //! do not save to file

  ClassDef(SimplifiedEvent,10)  //SimplifiedEvent class
};

#endif
