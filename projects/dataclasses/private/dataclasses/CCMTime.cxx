#include <string>
#include <dataclasses/jday.h>

#include <icetray/serialization.h>
#include <dataclasses/CCMTime.h>
#include <icetray/I3Units.h>
#include <boost/algorithm/string.hpp>
#include <boost/multiprecision/cpp_int.hpp>


CCMTime::CCMTime() {
    tv_sec_ = 0;
    tv_nsec_ = 0;
    tv_nsec_frac_ = 0;
}

CCMTime::~CCMTime() {}

CCMTime::CCMTime(int64_t tv_sec, int64_t tv_nsec, double nsec_frac) :
    tv_sec_(tv_sec), tv_nsec_(tv_nsec), tv_nsec_frac_(nsec_frac) {}

void CCMTime::SetUnixTime(time_t unixTime, double ns) {
    if(unixTime < 0) log_fatal("invalid Unix time");
    tv_sec_ = unixTime;
    tv_nsec_frac_ = ns % 1.0;
    tv_nsec_ = ns - tv_nsec_frac_;
}

time_t CCMTime::GetUnixTime() const {
    return tv_sec_;
}

std::string CCMTime::GetUTCString(std::string format) const {
    time_t t=GetUnixTime();
    if( IsLeapSecond())
    {
        t-=1;
        boost::replace_all(format, "%S", "60");
    }
    struct tm *tm=gmtime(&t);
    char datestring[256];
    strftime (datestring, sizeof(datestring), format.c_str(), tm);
    return datestring;
}

bool CCMTime::operator<(const CCMTime& rhs) const {
    if(tv_sec_ < rhs.tv_sec_)
        return true;
    if(tv_sec_ > rhs.tv_sec_)
        return false;
    if(tv_nsec_ < rhs.tv_nsec_)
        return true;
    if(tv_nsec_ > rhs.tv_nsec_)
        return false;
    if(tv_nsec_frac_ < rhs.tv_nsec_frac_)
        return true;
    return false;
}

bool CCMTime::operator>(const CCMTime& rhs) const {
    if(tv_sec_ > rhs.tv_sec_)
        return true;
    if(tv_sec_ < rhs.tv_sec_)
        return false;
    if(tv_nsec_ > rhs.tv_nsec_)
        return true;
    if(tv_nsec_ < rhs.tv_nsec_)
        return false;
    if(tv_nsec_frac_ > rhs.tv_nsec_frac_)
        return true;
    return false;
}

bool CCMTime::operator==(const CCMTime& rhs) const {
    return tv_sec_ == rhs.tv_sec_ && tv_nsec_ == rhs.tv_nsec_ && tv_nsec_frac_ == rhs.tv_nsec_frac_;
}

bool CCMTime::operator<=(const CCMTime& rhs) const
{
    if(*this < rhs || *this == rhs)
        return true;
    return false;
}

bool CCMTime::operator>=(const CCMTime& rhs) const
{
    if(*this > rhs || *this == rhs)
        return true;
    return false;
}

CCMTime CCMTime::operator+(const double second_term) const {
    CCMTime result;
    result.tv_sec_ = tv_sec_ + (int64_t)(second_term / 1e9);
    result.tv_nsec_ = tv_nsec + (int64_t)(second_term % 1e9);
    result.tv_nsec_frac_ = tv_nsec_frac_ + second_term % 1;
    if(result.tv_nsec_frac_ >= 1) {
        result.tv_nsec_frac_ -= 1;
        result.tv_nsec_ += 1;
    } else if(result.tv_nsec_frac_ < 0) {
        result.tv_nsec_frac_ += 1;
        result.tv_nsec_ -= 1;
    }
    if(result.tv_nsec_ >= 1e9) {
        result.tv_nsec_ -= 1e9;
        result.tv_sec_ += 1;
    } else if(result.tv_nsec_ < 0) {
        result.tv_nsec_ += 1e9;
        result.tv_sec_ -= 1;
    }
    return result;
}

CCMTime CCMTime::operator-(const double second_term) const {
    return *this+(-second_term);
}

double operator-(const CCMTime t1,const CCMTime t2) {
    double result = (t1.tv_sec_ - t2.tv_sec_) * 1e9 + (t1.tv_nsec_ - t2.tv_nsec_) + t1.tv_nsec_frac_ - t2.tv_nsec_frac_;
    return result;
}

std::ostream& CCMTime::Print(std::ostream& oss) const {
    int64_t daqt = tv_nsec_ + tv_nsec_frac_;
    oss << GetUTCString("%Y-%m-%d %H:%M:%S.");
    oss << std::setw(3) << std::setfill('0') << (daqt/1000000)%1000  << ',';
    oss << std::setw(3) << std::setfill('0') << (daqt/1000)%1000 << ',';
    oss << std::setw(3) << std::setfill('0') << (daqt)%1000 << " UTC";

    return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMTime& t){
    return(t.Print(oss));
}

template <class Archive>
    void
CCMTime::serialize(Archive& ar, unsigned version)
{
    if (version>ccmtime_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMTime class.",version,ccmtime_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("Sec", tv_sec_);
    ar & make_nvp("NSec", tv_nsec_);
    ar & make_nvp("NSecFraction", tv_nsec_frac_);
}

I3_SERIALIZABLE(CCMTime);
