// this module counts the number of cosmic triggers in a given dataset

#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>


class CountCosmicFrames: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    bool IsCosmicFrame(I3FramePtr frame);
public:
    CountCosmicFrames(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(CountCosmicFrames);

CountCosmicFrames::CountCosmicFrames(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void CountCosmicFrames::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void CountCosmicFrames::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    bool found_bcm_key = false;
    for(auto const & key_om : geo.pmt_geo) {
        if(key_om.second.omtype == CCMOMGeo::OMType::BeamCurrentMonitor) {
            bcm_key = key_om.first;
            found_bcm_key = true;
        }
    }
    if(not found_bcm_key) {
        log_fatal("CCMGeometry does not contain a channel corresponding to a BeamCurrentMonitor");
    }
    bcm_channel = geo.pmt_channel_map.at(bcm_key);
    cosmic_trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::CosmicTrigger, 1);
    PushFrame(frame);
}

bool CountCosmicFrames::IsCosmicFrame(I3FramePtr frame) {
    if(not frame->Has(nim_pulses_name_)) {
        log_warn(("No key named " + nim_pulses_name_ + " present in frame").c_str());
        return false;
    }
    boost::shared_ptr<NIMLogicPulseSeriesMap const> nim_pulses = frame->Get<boost::shared_ptr<NIMLogicPulseSeriesMap const>>(nim_pulses_name_);
    if(not nim_pulses) {
        log_warn(("No NIMLogicPulseSeriesMap named " + nim_pulses_name_ + " present in frame").c_str());
        return false;
    }
    NIMLogicPulseSeriesMap::const_iterator it = nim_pulses->find(cosmic_trigger_key);
    if(it == nim_pulses->end()) {
        log_warn(("NIMLogicPulseSeriesMap named " + nim_pulses_name_ + " does not contain the CosmicTrigger key").c_str());
        return false;
    }
    if(it->second.size() < 1) {
        // No NIM pulse on the cosmic trigger means we skip this frame
        return false;
    }

    // We found a NIM pulse on the cosmic trigger,
    // therefore this is a cosmic frame
    return true;
}

void CountCosmicFrames::DAQ(I3FramePtr frame) {
    
    // place to save tag on cosmics
    // 0 = does not contain a cosmic trigger
    // 1 = contains a cosmic trigger
    boost::shared_ptr<I3Vector<size_t>> cosmic_tag (new I3Vector<size_t>(1));

    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }

    if(not IsCosmicFrame(frame)) {
        size_t cosmic_status = 0;
        cosmic_tag->operator[](1) = cosmic_status;		
        frame->Put("CosmicTag", cosmic_tag);
	PushFrame(frame);
        return;
    }

    if(IsCosmicFrame(frame)) {
        size_t cosmic_status = 1;
        cosmic_tag->operator[](1) = cosmic_status;
        frame->Put("CosmicTag", cosmic_tag);
        std::cout << "found a cosmic!" << std::endl;
        PushFrame(frame);
        return;
    }

}
void CountCosmicFrames::Finish() {
}

