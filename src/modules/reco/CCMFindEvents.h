#ifndef CCMFindEvents_h
#define CCMFindEvents_h

#include "Utility.h"
#include "CCMModule.h"

#include <vector>
#include <iterator>
#include <algorithm>

class CCMFindEvents : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMFindEvents(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMFindEvents(const CCMFindEvents& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMFindEvents();

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
    void ConnectRawData(std::shared_ptr<RawData> rawData)  { fRawData = rawData; }
    void ConnectPulses(std::shared_ptr<Pulses> pulses) {fPulses = pulses; }

    void SetBuildEventFlag(bool flag) { fBuildEventFlag = flag; }
    void SetThreshold(double threshold) { fThreshold = threshold; }
    void SetTriggerType(std::string triggerType) { fTriggerType = triggerType; }
    void SetPulseRepresentation(int rep) { fPulseRep = rep; }
    void SetPulseTimeCut(bool flag) { fPulseTimeCut = flag; }
    void SetPulseTimeLowValue(double value) { fPulseTimeLowValue = value; }
    void SetPulseTimeHighValue(double value) { fPulseTimeHighValue = value; }

    std::vector<float> & GetWaveformInt() { return fIntegralTime; }

  private:

    //private methods
    /*!
     *  \brief Define the vectors used in finding events
     */
    void DefineVectors();

    /*!
     *  \brief Loop through the pulses to create the accumulated waveforms
     *
     *  The pulses for the veto tubes always is represented by a hit at the start time
     *  of the pulse. The charge of the veto tubes is not used.
     *
     *  For the tank tubes two possition representations are allowed, which representation
     *  is done is dependent on the value of #PulseRep.
     *
     *  <b>Triangular Representation</b>
     *  The pulse is represented as a triangle where the integral of the triangle is 
     *  equal to the charge of the pulse. The formula for the triangle is
     *  
     *  \f[
     *  I_{B} = \left\{\begin{array}{ll}
     *  \frac{2I}{L}\frac{B - B_{f}}{B_{m} - B_{f}} & B < B_{m}\\
     *  \frac{2I}{L}\frac{B_{l} - B}{B_{l} - B_{m}} & B > B_{m}\\
     *  \frac{2I}{L} & B = B_{m}
     *  \end{array}\right.
     *  \f]
     *
     *  where \f$B\f$ is the current bin, \f$B_{x}\f$ is either the first \f$f\f$, 
     *  middle \f$m\f$, or last \f$l\f$ bin location. \f$I\f$ is the integral of the pulse, \f$L\f$ is the length of the pulse, and \f$I_{B}\f$ is
     *  the amount of charge that goes into bin \f$B\f$.
     *
     *  <b>Start Time Representation</b>
     *  The pusle is represented as a single instance in time where all of the charge of the pulse
     *  goes into the start bin of the pulse
     */
    void BuildAccumulatedWaveform();

    /*!
     * \brief Loop over the accumulated waveforms and find events
     *
     * Loop over the accumulated waveforms to find the events based off the threshold
     * that is passed in the configuration file
     */
    void FindEvents();

    /*!
     *  \brief Shift the time of the event based on the BCM and FP3 time offsets
     *  \param[in] start The start bin of the pulse
     *  \return The time of the event in ns (input is bin count) 
     *
     *  Shift the time of the pulse to account for the jitter of the BCM. #Utility::fgkFP3Offset 
     *  is used to account for the time difference between the PMTs in CCM and the 
     *  EJ301 detector in FP3 
     *
     *  If #fBeamTime == 0 then the trigger was STROBE or LED so shift the time of the event 
     *  based on the DAQ window true start time (#Utility::fgkWindowStartTime)
     */
    double ShiftTime(int start);

    void ResetVectors();

  private:

    //private data members
    std::string fTriggerType;
    double fThreshold;

    std::shared_ptr<Events> fEvents;
    std::shared_ptr<RawData> fRawData;
    std::shared_ptr<Pulses> fPulses;

    unsigned long int fNumTriggers;
    
    std::vector<float> fPulsesTime;
    std::vector<float> fIntegralTime;
    std::vector<float> fIntegralDer;

    std::vector<int> fVetoBottomTime;
    std::vector<int> fVetoTopTime;
    std::vector<int> fVetoCRightTime;
    std::vector<int> fVetoCLeftTime;
    std::vector<int> fVetoCFrontTime;
    std::vector<int> fVetoCBackTime;
    std::vector<int> fVetoTotalTime;

    std::vector<std::vector<float>> fPMTWaveform;
    std::vector<std::vector<int>> fPMTWaveformCount;

    int fBuildEventFlag;
    int fPulseRep;
    int fPulseTimeCut;
    double fPulseTimeLowValue;
    double fPulseTimeHighValue;

    int fBeamTime;
};

#endif // CCMFindEvents_h

