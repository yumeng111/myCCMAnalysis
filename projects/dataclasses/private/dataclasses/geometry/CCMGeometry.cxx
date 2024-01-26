#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/geometry/CCMGeometry.h>

CCMGeometry::~CCMGeometry() {}

template <class Archive>
void CCMGeometry::save(Archive& ar, unsigned version) {
    if (version > ccmgeometry_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMGeometry class.", version, ccmgeometry_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTGeo", pmt_geo);
    ar & make_nvp("PMTChannelFormat", pmt_channel_map);
    ar & make_nvp("TriggerChannelMap", trigger_channel_map);
    ar & make_nvp("TriggerCopyMap", trigger_copy_map);
    ar & make_nvp("TriggerToTriggerCopyMap", trigger_to_trigger_copy_map);
    ar & make_nvp("StartTime", startTime);
    ar & make_nvp("EndTime", endTime);
}

template <class Archive>
void CCMGeometry::load(Archive& ar, unsigned version) {
    if (version > ccmgeometry_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMGeometry class.", version, ccmgeometry_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTGeo", pmt_geo);
    ar & make_nvp("PMTChannelFormat", pmt_channel_map);
    ar & make_nvp("TriggerChannelMap", trigger_channel_map);
    ar & make_nvp("TriggerCopyMap", trigger_copy_map);
    ar & make_nvp("TriggerToTriggerCopyMap", trigger_to_trigger_copy_map);
    ar & make_nvp("StartTime", startTime);
    ar & make_nvp("EndTime", endTime);

    if(version == 0) {
        for(std::pair<CCMPMTKey const, CCMOMGeo> & p : pmt_geo) {
            I3Position & pos = p.second.position;
            pos.SetX(pos.GetX() * I3Units::cm);
            pos.SetY(pos.GetY() * I3Units::cm);
            pos.SetZ(pos.GetZ() * I3Units::cm);
        }
    }
}

const CCMGeometry& CCMGeometry::operator=(const CCMGeometry& geometry) {
    pmt_geo = geometry.pmt_geo;
    pmt_channel_map = geometry.pmt_channel_map;
    trigger_channel_map = geometry.trigger_channel_map;
    trigger_copy_map = geometry.trigger_copy_map;
    trigger_to_trigger_copy_map = geometry.trigger_to_trigger_copy_map;
    startTime = geometry.startTime;
    endTime = geometry.endTime;

    return *this;
}

I3_SPLIT_SERIALIZABLE(CCMGeometry);
