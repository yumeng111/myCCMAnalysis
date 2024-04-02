#ifndef CCMTIME_H_INCLUDED
#define CCMTIME_H_INCLUDED

#include <string>
#include <time.h>

#include <icetray/I3FrameObject.h>
#include <dataclasses/Utility.h>

static const unsigned ccmtime_version_ = 0;

class CCMTime : public I3FrameObject {
public:
    CCMTime();

    /**
     * @brief creates the object with the given times as the DAQ time
     */
    CCMTime(int64_t tv_sec, double tv_nsec);

    ~CCMTime();

    std::ostream& Print(std::ostream&) const override;

    void SetUnixTime(time_t unixTime, double ns=0);

    int64_t GetSeconds() const;
    double GetNanoSeconds() const;

    void SetSeconds(int64_t seconds);
    void SetNanoSeconds(double nano_seconds);

    /**
     * @brief Gets the time in Unix convention
     * @returns The number of seconds since the Epoch (00:00:00 UTC, January 1, 1970)
     */
    time_t GetUnixTime() const;

    /**
     * @brief Gets a string representing the time in UTC
     */
    std::string GetUTCString(bool sub_second=true) const;

    /**
     * equality operator.  
     * @return true if the times are the same
     * @param rhs the CCMTime to compare this one to.
     */
    bool operator==(const CCMTime& rhs) const;

    /**
     * inequality operator
     * @return false if the times are different
     * @param rhs the CCMTime to compare this one to.
     */
    bool operator!=(const CCMTime& rhs) const
    {
        if(rhs == *this)
            return false;
        return true;
    }

    /**
     * comparison operator.
     * Compares first the year and then the DAQ time
     * @return true if the lhs should be ordered before the rhs
     */
    bool operator<(const CCMTime& rhs) const;

    /**
     * comparison operator.
     * Compares first the year and then the DAQ time
     * @return true if the lhs should be ordered after the rhs
     */
    bool operator>(const CCMTime& rhs) const;

    bool operator<=(const CCMTime& rhs) const;

    bool operator>=(const CCMTime& rhs) const;

    /**
     *Adds a double (please use I3Units of time) to CCMTime
     *Takes into account rounding and leap years
     */
    CCMTime operator+(const double) const;
    CCMTime operator-(const double) const;

    friend double operator-(const CCMTime t1,const CCMTime t2);

private:

    int64_t tv_sec_;
    double tv_nsec_;

    friend class icecube::serialization::access;

    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

/**
 * Returns the difference between two CCMTimes in nanoseconds
 */
double operator-(const CCMTime t1,const CCMTime t2);

std::ostream& operator<<(std::ostream& oss, const CCMTime& d);

I3_POINTER_TYPEDEFS(CCMTime);
I3_CLASS_VERSION(CCMTime, ccmtime_version_);

#endif //CCMTIME_H_INCLUDED
