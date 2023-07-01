/**
 * copyright  (C) 2004
 * the icecube collaboration
 * $Id$
 *
 * @file CCMWaveform.h
 * @version $Revision$
 * @date $Date$
 */
#ifndef CCMWAVEFORM_H_INCLUDED
#define CCMWAVEFORM_H_INCLUDED

#include <utility>
#include <vector>
#include <dataclasses/I3Map.h>
#include <dataclasses/ScintKey.h>
#include <icetray/CCMPMTKey.h>

/**
 * List the names of enumeration members defined in this file
 * here. These can be used for e.g. pybindings, which require
 * the names of the enumeration members to be known. This list
 * should be updated whenever members or new enums are added to
 * the class.
 */
#define CCMWAVEFORM_H_CCMWaveform_Source (V1730)
#define CCMWAVEFORM_H_CCMWaveform_Status (VIRGINAL)(COMBINED)(SATURATED)(UNDERSHOT)

enum CCMSource {
    Unknown = 0,
    V1730 = 1,
};

enum CCMStatus {
    VIRGINAL = 0, // Waveform looks good with nothing weird going on
    COMBINED = (1 << 1), // Waveform section that comes from a partial trigger
    SATURATED = (1 << 2), // Waveform section that is saturated
    UNDERSHOT = (1 << 3) // Waveform section that undershoots the pedestal
};

class CCMStatusCompound {
    private:
        std::pair<unsigned long long int, unsigned long long int>
            interval_;
        CCMStatus status_;

    public:
        CCMStatusCompound() : interval_(std::make_pair(0, 0)), status_(SATURATED) {}

        ~CCMStatusCompound();

        std::ostream& Print(std::ostream&) const;

        const std::pair<unsigned long long int, unsigned long long int>&
            GetInterval() const { return interval_; }

        std::pair<unsigned long long int, unsigned long long int>&
            GetInterval() { return interval_; }

        CCMStatus GetStatus() const { return status_; }

        void SetStatus(CCMStatus status) { status_ = status; }

        bool operator==(const CCMStatusCompound& rhs) const
        {
            return status_ == rhs.status_
                && interval_ == rhs.interval_;
        }
    private:
        friend class icecube::serialization::access;
        template<class Archive> void save(Archive& ar, unsigned version) const;
        template<class Archive> void load(Archive& ar, unsigned version);
        I3_SERIALIZATION_SPLIT_MEMBER();
};

static const unsigned ccmwaveform_version_ = 1;
template<typename T>
class CCMWaveform {
public:
    typedef T WaveformType;
    typedef CCMStatusCompound StatusCompound;
    typedef CCMStatus Status;
    typedef CCMSource Source;

    /**
     * Returns a summary of a given waveform/status information.
     * 
     * @return ADULTERATED/SHADY if any included status compound is ADULTERATED/SHADY,
     * or VIRGINAL.
     */
    static unsigned GetStatus(const std::vector<StatusCompound>& waveformInfo);

    unsigned GetStatus() const;

private:
    double startTime_;
    double binWidth_;
    std::vector<WaveformType> waveform_;
    std::vector<StatusCompound> waveformInfo_;
    Source source_;

public:
    CCMWaveform() :  startTime_(0), binWidth_(0), source_(V1730) {}

    ~CCMWaveform();

    std::ostream& Print(std::ostream&) const;

    double GetStartTime() const {return startTime_;}

    void SetStartTime(double startTime) {startTime_ = startTime;}

    double GetBinWidth() const {return binWidth_;}

    void SetBinWidth(double binWidth) {binWidth_ = binWidth;}

    const std::vector<WaveformType>& GetWaveform() const {return waveform_;}

    std::vector<WaveformType>& GetWaveform() {return waveform_;}

    void SetWaveform(const std::vector<WaveformType>& waveform) {waveform_ = waveform;}

    /**
     * Returns a status information for this waveform.
     * 
     * @return A collection of status compounds.
     * A status compound consists of an interval and a status information.
     * It describes the status of all waveform bins with indices
     * [GetInterval().first, GetInterval().second). If there are waveform bins not
     * described by a status compound, these bins are assumed to have a status equal
     * to VIRGINAL, e. g. GetWaveformInformation().empty() equal true means, all
     * waveform bins are good.
     */
    const std::vector<StatusCompound>& GetWaveformInformation() const
    {return waveformInfo_;}

    /**
     * Returns a status information for this waveform.
     * 
     * @return A collection of status compounds.
     * A status compound consists of an interval and a status information.
     * It describes the status of all waveform bins with indices
     * [GetInterval().first, GetInterval().second). If there are waveform bins not
     * described by a status compound, these bins are assumed to have a status equal
     * to VIRGINAL, e. g. GetWaveformInformation().empty() equal true means, all
     * waveform bins are good.
     */
    std::vector<StatusCompound>& GetWaveformInformation() {return waveformInfo_;}

    /**
     * Set the status information for this waveform.
     * 
     * Note: This method was added since there is a set method for the waveforms, too.
     * One might want to use the non-const get methods instead.
     * 
     * @param info The status information.
     */
    void SetWaveformInformation(const std::vector<StatusCompound>& info)
    {waveformInfo_ = info;}

    Source GetSource() const { return source_; }
    void SetSource(Source source) { source_ = source; }

private:
    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

template<typename T>
bool operator==(CCMWaveform<T> const & lhs, CCMWaveform<T> const & rhs);
template<typename T>
std::ostream& operator<<(std::ostream& oss, CCMWaveform<T> const & wf);
std::ostream& operator<<(std::ostream& oss, CCMStatusCompound const & wf);

typedef CCMWaveform<uint16_t> CCMWaveformUInt16;
typedef CCMWaveform<double> CCMWaveformDouble;

I3_CLASS_VERSION(CCMStatusCompound, 1);
I3_CLASS_VERSION(CCMWaveformUInt16, ccmwaveform_version_);
I3_CLASS_VERSION(CCMWaveformDouble, ccmwaveform_version_);
I3_POINTER_TYPEDEFS(CCMWaveformUInt16);
I3_POINTER_TYPEDEFS(CCMWaveformDouble);

typedef I3Vector<CCMWaveformUInt16> CCMWaveformUInt16Series;
typedef I3Vector<CCMWaveformDouble> CCMWaveformDoubleSeries;
I3_POINTER_TYPEDEFS(CCMWaveformUInt16Series);
I3_POINTER_TYPEDEFS(CCMWaveformDoubleSeries);

#include "CCMWaveform.tcc"

#endif // CCMWAVEFORM_H_INCLUDED
