#ifndef Events_v2_h
#define Events_v2_h

#include <array>
#include <vector>
#include <iostream>
#include <sys/types.h>

#include "CCMAnalysis/CCMDataStructures/Events.h"
#include "CCMAnalysis/CCMDataStructures/SimplifiedEvent.h"

#include "TObject.h"

/*!**********************************************
 * \class Events_v2
 * \brief Container of the #SimplifiedEvent found in a given DAQ window.
 *
 * Container for a vector of #SimplifiedEvent events found in the detector
 * and the methods to find such events. The class is saved in the
 * data file.
 ***********************************************/
class Events_v2 : public TObject
{
  public:
    Events_v2();
    Events_v2(const Events_v2 & p);
    ~Events_v2();

    void Reset();
    void ClearEvents();

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

    
    /// \return #fEventNumber
    unsigned int GetEventNumber() { return fEventNumber; }
    
    /// \return #fComputerSecIntoEpoch
    unsigned int GetComputerSecIntoEpoch() { return fComputerSecIntoEpoch; }
    
    /// \return #fComputerNSIntoSec
    unsigned int GetComputerNSIntoSec() { return fComputerNSIntoSec; }
    
    /// \brief Sets the time of the BCM signal (0 if not BEAM trigger)
    void SetBeamTime(int time) { fBeamTime= time; }
    /// \brief Sets the integral of the BCM signal (0 if not BEAM trigger)
    void SetBeamIntegral(float integral) { fBeamIntegral = integral; }
    /// \brief Sets the length of the BCM signal (0 if not BEAM trigger)
    void SetBeamLength(float length) { fBeamLength = length; }
    
    /// \brief Get the time of the BCM signal
    int GetBeamTime() { return fBeamTime; }
    /// \brief Get the time of the BCM signal
    int GetBeamTime() const { return fBeamTime; }

    /// \brief Get the integral of the BCM signal
    float GetBeamIntegral() { return fBeamIntegral; }
    /// \brief Get the integral of the BCM signal
    float GetBeamIntegral() const { return fBeamIntegral; }

    /// \brief Get the length of the BCM signal
    float GetBeamLength() { return fBeamLength; }
    /// \brief Get the length of the BCM signal
    float GetBeamLength() const { return fBeamLength; }

    size_t GetNumEvents() { return fNumEvents; }
    float GetTriggerTime() { return fTriggerTime; }
    void SetTriggerTime(float time) { fTriggerTime = time; }

    size_t GetNumEvents() const { return fNumEvents; }
    float GetTriggerTime() const { return fTriggerTime; }

    void AddSimplifiedEvent(const SimplifiedEvent & event) { fEvents.push_back(event); ++fNumEvents;}
    void RemoveSimplifiedEvent(size_t pos) { fEvents.erase(fEvents.begin()+pos); --fNumEvents;}
    void RemoveSimplifiedEvent(std::vector<SimplifiedEvent>::iterator it) { fEvents.erase(it); --fNumEvents;}
    const SimplifiedEvent & GetSimplifiedEvent(size_t pos) { return fEvents[pos]; }
    const SimplifiedEvent & GetSimplifiedEvent(size_t pos) const { return fEvents[pos]; }
    
    void UpdatePosition(size_t index, double x, double y, double z) {
      fEvents[index].SetXPosition(x);
      fEvents[index].SetYPosition(y);
      fEvents[index].SetZPosition(z);
      return;
    }

    Events_v2 & operator=(const Events_v2 & rhs);

  private:
    /// The time of the window with sample number == 0 (not saved)
    /// in microseconds
    constexpr static const double fgkTimeWindowStart = -9.92; //!
    /// The sample time width (not saved)
    /// in microseconds
    constexpr static const double fgkSampleTimeWidth = 2e-3; //!

    /// The event number
    unsigned int fEventNumber;
    /// The time of the window in s
    unsigned int fComputerSecIntoEpoch;
    /// The time of the window in ns
    unsigned int fComputerNSIntoSec;

    /// The number of events found
    size_t fNumEvents;

    /// The trigger time of the trigger
    float fTriggerTime;

    /// The BCM time in the trigger
    int fBeamTime;
    /// The BCM integral in the trigger
    float fBeamIntegral;
    /// The BCM length in the trigger
    float fBeamLength;

    /// The vector of events found
    std::vector<SimplifiedEvent> fEvents;

    /// Pointer to the current #SimplifiedEvent (not saved)
    SimplifiedEvent * fCurrEvent; //!
};

#endif // #ifndef Events_v2_h
