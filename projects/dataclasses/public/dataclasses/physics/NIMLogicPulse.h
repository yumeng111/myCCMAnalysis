
#ifndef NIMLogicPulse_H_INCLUDED
#define NIMLogicPulse_H_INCLUDED

#include <icetray/I3Logging.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/Utility.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/I3Map.h>

/**
 * @brief The basic NIM pulse class
 *
 */
static const unsigned nimlogicpulse_version_ = 1;

class NIMLogicPulse {
private:
    double     time_;             // Time (ns) at which the pulse was observed in the DAQ window
    double     length_;           // Duration of pulse
public:
    /**
     * Default constructor.
     */
    NIMLogicPulse()
        : time_(0.0), length_(0.0) {}

    /**
     * Destructor.
     */
    ~NIMLogicPulse();

    std::ostream& Print(std::ostream&) const;

    /**
     * Retrieves time at which the NIM logic pulse was observed.
     *
     * @return NIM time.
     */
    double GetNIMPulseTime() const {return time_;}
    /**
     * Sets time at which the NIM logic pulse was observed
     *
     * @param time NIM time.
     */
    void SetNIMPulseTime(double time) {time_ = time;}

    /**
     * Retrieves duration of the NIM logic pulse
     *
     * @return NIM duration.
     */
    double GetNIMPulseLength() const {return length_;}
    /**
     * Sets duration of NIMed readout window.
     *
     * @param length NIM duration.
     */
    void SetNIMPulseLength(double length) {length_ = length;}

    bool operator==(const NIMLogicPulse& rhs) const;
    bool operator!=(const NIMLogicPulse& rhs) const;

private:

    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);

    // logging
    SET_LOGGER("NIMLogicPulse");
};

std::ostream& operator<<(std::ostream& oss, const NIMLogicPulse& t);

I3_CLASS_VERSION(NIMLogicPulse, nimlogicpulse_version_);
/**
 * pointer type to insulate users from memory management
 */
I3_POINTER_TYPEDEFS(NIMLogicPulse);


typedef I3Vector<NIMLogicPulse> NIMLogicPulseSeries;
typedef I3Map<CCMTriggerKey, NIMLogicPulseSeries> NIMLogicPulseSeriesMap;
I3_POINTER_TYPEDEFS(NIMLogicPulseSeries);
I3_POINTER_TYPEDEFS(NIMLogicPulseSeriesMap);

#endif

