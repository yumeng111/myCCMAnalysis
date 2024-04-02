#ifndef CCMEVENTHEADER_H_INCLUDED
#define CCMEVENTHEADER_H_INCLUDED

// includes

#include <string>

#include <dataclasses/CCMTime.h>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>

static const unsigned ccmeventheader_version_ = 0;

/**
 * @brief The header for data on the Event stream.
 *
 * Supposed to be a header for the event that
 * you can store in a file if you don't want to store the event.
 * It's supposed to have enough data to reconstruct the full event
 * from the database
 */
class CCMEventHeader : public I3FrameObject {
private:
    /**
     * @brief Unique identifier assigned sequentially to every run the DAQ has ever recorded.
     */
    unsigned runID_;

    /**
     * @brief
     */
    unsigned subRunID_;

    /**
     * @brief Unique identifier assigned by the DAQ to every event (Q-frame) in the given run.
     */
    unsigned eventID_;

    /**
     * @brief Unique identifier assigned to every P-frame by the splitter.
     * 
     * This quantity is not assigned by the DAQ. When Q-frames are split into multiple
     * P-frames by a splitter module this is the unique number assigned to each P-frame.
     */
    unsigned subEventID_;

    /** 
     * @brief Name of the P-frame stream assigned by the splitter module
     */
    std::string subEventStream_;

    /**
     * @brief The 
     * 
     */
    CCMTime startTime_;
    CCMTime endTime_;

public:
    CCMEventHeader();

    ~CCMEventHeader();

    std::ostream& Print(std::ostream&) const override;

    bool operator==(const CCMEventHeader& other) const {
        return (runID_ == other.runID_ &&
                subRunID_ == other.subRunID_ &&
                eventID_ == other.eventID_ &&
                subEventID_ == other.subEventID_ &&
                subEventStream_ == other.subEventStream_ &&
                startTime_ == other.startTime_ &&
                endTime_ == other.endTime_);
    }
    bool operator!=(const CCMEventHeader& other) const {
        return !operator==(other);
    }

    CCMTime GetStartTime() const {
        return startTime_;
    }

    void SetStartTime(CCMTime time) {
        startTime_ = time;
    }

    CCMTime GetEndTime() const {
        return endTime_;
    }

    void SetEndTime(CCMTime time) {
        endTime_ = time;
    }

    /**
     * @return the run id for the event
     */
    unsigned GetRunID() const { return runID_; }

    /**
     * @param runid the new run id for the event
     */
    void SetRunID(unsigned runid) { runID_ = runid; }

    /**
     * @return the subrun id for the event
     */
    unsigned GetSubRunID() const { return subRunID_; }

    /**
     * @param subrunid the new subrun id for the event
     */
    void SetSubRunID(unsigned subrunid) { subRunID_ = subrunid; }

    /**
     * @return the event id for this event
     */
    unsigned GetEventID() const { return eventID_; }

    /**
     * @param eventid the new event id for the event
     */
    void SetEventID(unsigned eventid) { eventID_ = eventid; }

    /**
     * @return the subevent id for this subevent
     */
    unsigned GetSubEventID() const { return subEventID_; }

    /**
     * @param eventid the new subevent id for the subevent
     */
    void SetSubEventID(unsigned eventid) { subEventID_ = eventid; }

    /**
     * @return the subevent stream description
     */
    std::string GetSubEventStream() const { return subEventStream_; }

    /**
     * @param desc the new subevent stream description
     */
    void SetSubEventStream(const std::string &desc) { subEventStream_ = desc; }

    /**
     * @return the name of the stream this header is for.... "Physics"
     */
    const std::string GetDataStream(){ return "Physics";}

private:
    friend class icecube::serialization::access;

    template <class Archive> void serialize(Archive & ar, unsigned version);
};

std::ostream& operator<<(std::ostream& oss, const CCMEventHeader& eh);

I3_CLASS_VERSION(CCMEventHeader, ccmeventheader_version_);
I3_POINTER_TYPEDEFS(CCMEventHeader);
I3_DEFAULT_NAME(CCMEventHeader);

#endif //CCMEVENTHEADER_H_INCLUDED

