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
    CCMResult_t ProcessEvent();

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


  private:

    //private methods
    void DefineVectors();
    template<typename T>
    int FindFirstNoneEmptyBin(typename std::vector<T>::iterator begin, 
        typename std::vector<T>::iterator start, typename std::vector<T>::iterator end);

  private:

    //private data members
    std::string fTriggerType;
    double fThreshold;

    std::shared_ptr<Events> fEvents;
    std::shared_ptr<RawData> fRawData;
    std::shared_ptr<Pulses> fPulses;

    unsigned long int fNumTriggers;
    
    // Set the number of bins and bin width
    // This is hard coded and should be taken 
    // from either data_structures.hh or a data base
    // since they (in principle) could change
    constexpr static const int fgkNumBins = 8000;
    constexpr static const double fgkBinWidth = 2e-3;
    constexpr static const double fgkNumPMTs = 160;

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
    
};

template <class T>
int CCMFindEvents::FindFirstNoneEmptyBin(typename std::vector<T>::iterator begin, 
    typename std::vector<T>::iterator start, typename std::vector<T>::iterator end)
{
  auto time = std::find_if_not(start,end,[](const T & i){return i == static_cast<T>(0);});
  if (std::distance(time,end) == 0) {
    return -1;
  }
  return std::distance(begin,time);
}

#endif // CCMFindEvents_h

