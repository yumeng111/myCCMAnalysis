#include <icetray/serialization.h>
#include <dataclasses/physics/CCMEventHeader.h>

#include <limits>


CCMEventHeader::CCMEventHeader()
    : runID_(std::numeric_limits<unsigned int>::max()),
    subRunID_(std::numeric_limits<unsigned int>::max()),
    eventID_(std::numeric_limits<unsigned int>::max()),
    subEventID_(0)
{
}

CCMEventHeader::~CCMEventHeader() {}

template <class Archive>
void CCMEventHeader::serialize(Archive& ar, unsigned version) {
    if(version > ccmeventheader_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMEventHeader class.", version, ccmeventheader_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("RunID", runID_);
    ar & make_nvp("SubRunID", subRunID_);
    ar & make_nvp("EventID", eventID_);
    ar & make_nvp("SubEventID", subEventID_);
    ar & make_nvp("SubEventStream", subEventStream_);
    ar & make_nvp("StartTime", startTime_);
    ar & make_nvp("EndTime", endTime_);
}

std::ostream& CCMEventHeader::Print(std::ostream& oss) const {
    oss << "[ CCMEventHeader:\n"
        << "        StartTime: " << GetStartTime() << '\n'
        << "         EndTime : " << GetEndTime() << '\n'
        << "           RunID : " << GetRunID() << '\n'
        << "        SubrunID : " << GetSubRunID() << '\n'
        << "         EventID : " << GetEventID() << '\n'
        << "      SubEventID : " << GetSubEventID() << '\n'
        << "  SubEventStream : " << GetSubEventStream() << '\n'
        << "]" ;
    return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMEventHeader& eh) {
    return(eh.Print(oss));
}

I3_SERIALIZABLE(CCMEventHeader);
