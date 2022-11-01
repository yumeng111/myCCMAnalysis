#ifndef CCMBuildAccumWaveform_h
#define CCMBuildAccumWaveform_h

#include <vector>
#include <iterator>
#include <algorithm>

#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMFramework/CCMModule.h"

class CCMBuildAccumWaveform : public CCMModule
{
  public:
    /*!
     *  \brief The constructor
     *  \param version FIXME
     */
    CCMBuildAccumWaveform(const char* version);

    /*!
     *  \brief The copy constructor
     *  \param clufdr the object being copied
     */
    CCMBuildAccumWaveform(const CCMBuildAccumWaveform& clufdr);

    /*!
     *  \brief The destructor
     */
    ~CCMBuildAccumWaveform();

    /*!
     *  \brief This is where the action takes place
     *  \return CCMResult_t the result of running this module
     *
     *  The pulses for the veto tubes always is represented by a hit at the start time
     *  of the pulse. The charge of the veto tubes is not used.
     *
     *  For the tank tubes two position representations are allowed.
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
     *  middle \f$m\f$, or last \f$l\f$ bin location. \f$I\f$ is the integral of the pulse, 
     *  \f$L\f$ is the length of the pulse, and \f$I_{B}\f$ is the amount of charge that goes into bin \f$B\f$.
     *
     *  <b>Start Time Representation</b>
     *  The pusle is represented as a single instance in time where all of the charge of the pulse
     *  goes into the start bin of the pulse
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

    void ConnectAccumWaveform(std::shared_ptr<AccumWaveform> evt) { fAccumWaveform = evt; }
    void ConnectRawData(std::shared_ptr<RawData> rawData)  { fRawData = rawData; }
    void ConnectPulses(std::shared_ptr<Pulses> pulses) {fPulses = pulses; }

    void SetTriggerType(std::string triggerType) { fTriggerType = triggerType; }
    void SetPulseTimeLowValue(double value) { fPulseTimeLowValue = value; }
    void SetPulseTimeHighValue(double value) { fPulseTimeHighValue = value; }
    void AddBuildMethod(CCMAccumWaveformMethod_t method) { fBuildMethods.emplace_back(method); }

    void Dump();

  private:

    //private methods

  private:

    //private data members
    std::string fTriggerType;

    std::shared_ptr<AccumWaveform> fAccumWaveform;
    std::shared_ptr<RawData> fRawData;
    std::shared_ptr<Pulses> fPulses;

    unsigned long int fNumTriggers;
    
    double fPulseTimeLowValue;
    double fPulseTimeHighValue;

    std::vector<CCMAccumWaveformMethod_t> fBuildMethods;

    int fBeamTime;
};

#endif // CCMBuildAccumWaveform_h

