/**
 * copyright  (C) 2004
 * the icecube collaboration
 * @version $Id$
 * @file CCMDetectorStatus.h
 * @date $Date$
 */

#ifndef CCMDetectorStatus_H_INCLUDED
#define CCMDetectorStatus_H_INCLUDED

#include <map>

#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <dataclasses/I3Time.h>
#include <icetray/OMKey.h>
#include <dataclasses/TriggerKey.h>
#include <dataclasses/Utility.h>
#include <dataclasses/status/CCMPMTStatus.h>
#include <dataclasses/status/I3TriggerStatus.h>

/**
 * @brief This is the state of the aspects of the detector that people have 
 * direct control over.  Contains the "per run" settings.
 *
 * Stuff that is a 'knob' on the detector.  This is a
 * top-level object in the frame related to this 'Detector Status' information.
 * Contains: 
 * - map of per DOM configurations (also the list of active DOMs),
 * - map of per AOM configurations (including AOMs read out by the TWR DAQ),
 * - map of active icecube triggers (and their configurations),
 * - map of active amanda triggers (and their configurations) and
 * - map of active domhubs (and their settings) ... eventually.
 */
static const unsigned ccmdetectorstatus_version_ = 1;

class CCMDetectorStatus : public I3FrameObject {
public:
    I3Time startTime;
    I3Time endTime;

    CCMPMTStatusMap pmtStatus;

    CCMDetectorStatus() {}

    ~CCMDetectorStatus();

    bool operator==(const CCMDetectorStatus& rhs)
    {
        return (startTime == rhs.startTime &&
                endTime == rhs.endTime &&
                pmtStatus == rhs.pmtStatus);
    }
    bool operator!=(const CCMDetectorStatus& rhs)
    {
        return !operator==(rhs);
    }

private:
    friend class icecube::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    I3_SERIALIZATION_SPLIT_MEMBER();
};

I3_CLASS_VERSION(CCMDetectorStatus, ccmdetectorstatus_version_);
I3_DEFAULT_NAME(CCMDetectorStatus);
I3_POINTER_TYPEDEFS(CCMDetectorStatus);

#endif /*CCMDetectorStatus_H_INCLUDED*/

