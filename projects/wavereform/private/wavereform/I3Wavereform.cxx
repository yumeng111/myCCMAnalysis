/**
 * @file 
 * @brief 
 *
 * (c) 2011 the IceCube Collaboration
 *
 * $Id:$
 * @version $Revision$
 * @date $Date$
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 */

#include "wavereform/I3Wavereform.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/physics/I3RecoPulse.h"
#include "dataclasses/physics/I3Waveform.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

I3_MODULE(I3Wavereform);

I3Wavereform::I3Wavereform(const I3Context &ctx) : I3ConditionalModule(ctx)
{
	waveform_name_ = "CalibratedATWD";
	AddParameter("Waveforms", "Name of calibrated waveforms in the frame",
	    waveform_name_);
	
	pulse_name_ = "RecoPulses";
	AddParameter("Pulses", "Name of the pulse series in the frame",
	    pulse_name_);
	
	chi_name_ = "ATWDChi_RecoPulses";
	AddParameter("Chi", "Where to store the chi^2 between waveform and pulses",
	    chi_name_);
	
	chi_threshold_ = 2e3;
	AddParameter("ChiThreshold",
	    "Flag any OM with a waveform/pulse chi^2 greater than this.",
	    chi_threshold_);
	
	flag_name_ = "";
	AddParameter("Flag", "Where to store the list of bad OMs", flag_name_);
	
	use_domsimulator_hacks_ = false;
	AddParameter("UseDOMsimulatorTemplates", "Use pulse shapes from DOMsimulator, "
	    "rather than those from reality.", use_domsimulator_hacks_);
	
	AddOutBox("OutBox");
}

void
I3Wavereform::Configure()
{
	GetParameter("Waveforms", waveform_name_);
	GetParameter("Pulses", pulse_name_);
	GetParameter("Chi", chi_name_);
	GetParameter("ChiThreshold", chi_threshold_);
	GetParameter("Flag", flag_name_);
}

void
I3Wavereform::Calibration(I3FramePtr frame)
{
	calibration_ = frame->Get<I3CalibrationConstPtr>();
	PushFrame(frame);
}

void
I3Wavereform::DetectorStatus(I3FramePtr frame)
{
	status_ = frame->Get<I3DetectorStatusConstPtr>();
	PushFrame(frame);
}

static double
chisquared(const std::vector<I3Wavereform::Refolded> &bins)
{
	double chi = 0.0;
	BOOST_FOREACH(const I3Wavereform::Refolded &bin, bins) {
		double diff = (bin.waveform-bin.refolded)/bin.step;
		chi += diff*diff;
	}
	
	return chi;
}

void
I3Wavereform::DAQ(I3FramePtr frame)
{
	I3RecoPulseSeriesMapConstPtr pulse_map =
	    frame->Get<I3RecoPulseSeriesMapConstPtr>(pulse_name_);
	I3WaveformSeriesMapConstPtr waveform_map = 
	    frame->Get<I3WaveformSeriesMapConstPtr>(waveform_name_);
	
	if (!(pulse_map && waveform_map)) {
		PushFrame(frame);
		return;
	} else if (!calibration_) {
		log_fatal("I haven't seen a calibration frame yet!");
		return;
	} else if (!status_) {
		log_fatal("I haven't seen a detector status frame yet!");
		return;
	}
	
	I3MapKeyVectorDoublePtr chi_map = boost::make_shared<I3MapKeyVectorDouble>();
	I3VectorOMKeyPtr screwy_doms = boost::make_shared<I3VectorOMKey>();
	
	I3MapKeyVectorDouble::iterator inserter = chi_map->begin();
	I3WaveformSeriesMap::const_iterator waveform_map_it = waveform_map->begin();
	
	for ( ; waveform_map_it != waveform_map->end(); waveform_map_it++) {
		const OMKey &om = waveform_map_it->first;
		const I3WaveformSeries &waveforms = waveform_map_it->second;
		
		I3RecoPulseSeriesMap::const_iterator series_it =
		    pulse_map->find(om);
		if (series_it == pulse_map->end())
			continue;
		
		std::map<OMKey, I3DOMCalibration>::const_iterator calib =
		    calibration_->domCal.find(om);
		std::map<OMKey, I3DOMStatus>::const_iterator status =
		    status_->domStatus.find(om);
		
		if (calib == calibration_->domCal.end() ||
		    status == status_->domStatus.end())
			continue;
		
		std::vector<double> chis;
		chis.reserve(waveforms.size());
		
		bool borked = false;
		I3WaveformSeries::const_iterator waveform_it = waveforms.begin();
		for ( ; waveform_it != waveforms.end(); waveform_it++) {
			chis.push_back(chisquared(GetRefolded(series_it->second,
			    *waveform_it, calib->second, status->second)));
			
			if (chis.back() > chi_threshold_)
				borked = true;
		}
			
		inserter = chi_map->insert(inserter, std::make_pair(om, chis));
		
		if (borked)
			screwy_doms->push_back(waveform_map_it->first);
	}
	
	frame->Put(chi_name_, chi_map);
	
	if ((flag_name_.size() > 0) && (screwy_doms->size() > 0))
		frame->Put(flag_name_, screwy_doms);
		
	PushFrame(frame);
}
