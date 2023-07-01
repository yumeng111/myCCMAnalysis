#include <cctype>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <cholmod.h>

#include <boost/make_shared.hpp>

#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/I3ConditionalModule.h>

#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

#include <dataclasses/CCMPMTFunctions.h>
#include <dataclasses/I3TimeWindow.h>
#include <dataclasses/calibration/CCMCalibration.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/status/CCMDetectorStatus.h>
#include <dataclasses/status/CCMPMTStatus.h>

#include "nnls.h"
#include "rnnls.h"

/* Simple class to hold together ATWD and FADC
 * templates along with a validity flag and waveform features */
struct CCMWaveformTemplate {

    std::vector<double> digitizer_template;
    double digitizerStart;
    double digitizerStop;
    bool filled;
};

class CCMWavedeform : public I3ConditionalModule {
    public:
        bool geo_seen;
        std::string geometry_name_;
        std::string nim_pulses_name_;
        I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
        CCMWavedeform(const I3Context &);
        virtual ~CCMWavedeform();

        void Geometry(I3FramePtr frame);
        void Configure();
        void Calibration(I3FramePtr frame);
        void DAQ(I3FramePtr frame);
    private:

        CCMRecoPulseSeriesPtr GetPulses(
                CCMWaveformDouble const & wf,
                const CCMWaveformTemplate& wfTemplate,
                const CCMPMTCalibration& calibration,
                const double spe_charge);

        std::string waveforms_name_;
        std::string waveform_range_name_;
        std::string output_name_;

        double spes_per_bin_;
        double tolerance_;
        double noise_threshold_;
        double basis_threshold_;

        int template_bins_;
        double template_bin_spacing_;

        bool apply_spe_corr_;
        bool reduce_;

        CCMWaveformTemplate template_;

        void FillTemplate(CCMWaveformTemplate& wfTemplate,
                const CCMPMTCalibration& calibration);

        cholmod_common c;
};

I3_MODULE(CCMWavedeform);

// End time (ns) of pulse template
inline double PulseWidth() {
    return 50;
}

// Beginning time (ns) of pulse template
inline double PulseMin() {
    return -2;
}

CCMWavedeform::CCMWavedeform(const I3Context& context) : I3ConditionalModule(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
    AddParameter("SPEsPerBin",
            "Number of basis functions to unfold per waveform bin", 4.);
    AddParameter("Tolerance",
            "Stopping tolerance, in units of bin mV^2/PE", 9.);
    AddParameter("NoiseThreshold",
            "Consider bins with amplitude below this number of counts as noise", 2.);
    AddParameter("BasisThreshold",
            "Require a bin with amplitude at least this number of counts "
            "within the FWHM of the template waveform in order to include "
            "a given start time in the basis set", 3.);
    AddParameter("Waveforms", "Name of input waveforms",
            "CalibratedWaveforms");
    AddParameter("WaveformTimeRange", "Name of maximum time range of "
            "calibrated waveforms for this event", "CalibratedWaveformRange");
    AddParameter("Output", "Name of output pulse series",
            "WavedeformPulses");
    AddParameter("ApplySPECorrections", "Whether to apply DOM-by-DOM"
            " corrections to the pulse charge scaling if available", true);
    AddParameter("Reduce", "Find the optimal NNLS solution, then eliminate"
            " basis members until tolerance is reached", true);

    cholmod_l_start(&c);
}

void CCMWavedeform::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("Waveforms", waveforms_name_);
    GetParameter("WaveformTimeRange", waveform_range_name_);
    GetParameter("Output", output_name_);
    GetParameter("SPEsPerBin", spes_per_bin_);
    GetParameter("Tolerance", tolerance_);
    GetParameter("NoiseThreshold", noise_threshold_);
    GetParameter("BasisThreshold", basis_threshold_);
    GetParameter("ApplySPECorrections", apply_spe_corr_);
    GetParameter("Reduce", reduce_);

    double range;

    /* Determine the number of bins and the bin resolution when creating
     * the ATWD template waveform cache.  Since the ATWD bin size is ~3.3 ns,
     * we will need to sample our cached template at a resolution of
     * ~3.3 / spes_per_bin_.  Build the cached waveform templates with a
     * factor 10 better resolution to minimize binning artifacts when
     * sampling it.
     */
    // CCM has a bin width of 2ns
    template_bin_spacing_ = 2.0 / spes_per_bin_ / 10;
    range = PulseWidth() - PulseMin();

    template_bins_ = (int)ceil(range / template_bin_spacing_);

    template_.filled = false;
    template_.digitizer_template.resize(template_bins_);
}

CCMWavedeform::~CCMWavedeform() {
    cholmod_l_finish(&c);
}

void CCMWavedeform::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void CCMWavedeform::Calibration(I3FramePtr frame) {

    /* Void the waveform templates since they could possibly
     * change frame-by-frame.
     */
    template_.filled = false;
    PushFrame(frame, "OutBox");
}

void CCMWavedeform::DAQ(I3FramePtr frame) {
    if (!frame->Has(waveforms_name_)) {
        PushFrame(frame);
        return;
    }

    const CCMCalibration& calibration = frame->Get<CCMCalibration>();
    const CCMDetectorStatus& status = frame->Get<CCMDetectorStatus>();
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    boost::shared_ptr<const CCMWaveformDoubleSeries> droopy_waveforms = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>("DroopCorrectedWf");
    //const std::vector<CCMWaveformUInt16> & waveforms = frame->Get<std::vector<CCMWaveformUInt16>>(waveforms_name_);
    I3Map<CCMPMTKey, BaselineEstimate> const & baselines = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");
    boost::shared_ptr<CCMRecoPulseSeriesMap> output(new CCMRecoPulseSeriesMap);

    // Unfold each set of waveforms into a pulse series
    // TODO Either loop over channels and check for PMTs or find all the PMTs and then loop and grab the right channels
    // loop over each pmt
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baselines){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second;
        double baseline = value.baseline; 
        uint32_t channel = pmt_channel_map_[key];

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);
        CCMWaveformDouble const & droopy_waveform = droopy_waveforms->at(channel);

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);
        std::map<CCMPMTKey, CCMPMTStatus>::const_iterator stat = status.pmtStatus.find(key);

        // Get/create the appropriate template
        if (!template_.filled) {
            FillTemplate(template_, calib->second);
        }

        (*output)[key] = *GetPulses(droopy_waveform, template_, calib->second, SPEMean(stat->second, calib->second));

    }


    frame->Put(output_name_, output);

    // Add a time window for the earliest possible pulse in the event
    I3TimeWindowConstPtr waveform_range = frame->Get<I3TimeWindowConstPtr>(waveform_range_name_);
    if (!waveform_range)
        log_fatal("Waveform time range \"%s\" not found in frame",
                waveform_range_name_.c_str());

    // XXX: Assume FADC has the widest binning
    I3TimeWindowPtr pulse_range(new I3TimeWindow(
                waveform_range->GetStart() - 2.0*I3Units::ns,
                waveform_range->GetStop()
                ));
    frame->Put(output_name_ + "TimeRange", pulse_range);

    PushFrame(frame, "OutBox");
}

/*
 *  GetPulses: Solve for the best fit pulse series to a set
 *  of DOM waveforms using NNLS.  We do the following:
 *  1.  Load waveform data (FADC + ATWD, if available) into a linear vector
 *  2.  Eliminate data bins below the set noise threshold
 *  3.  Determine the least-squares weights for each data bin.  Weights
 *          are reduced or set to zero for data we trust less or known
 *          bad data.
 *  4.  Create the set of possible pulses given the time range of nonzero data.
 *          These pulses are defined entirely by their start times.
 *  5.  Remove data bins with zero weight and those outside the support
 *          of all possible pulses.
 *  6.  Build the response matrix A, such Ax = y, where x is the set of
 *          amplitudes corresponding to each pulse and y is the data.
 *  7.  Solve the above for x using NNLS, yielding the pulse amplitudes.
 */
CCMRecoPulseSeriesPtr CCMWavedeform::GetPulses(CCMWaveformDouble const & wf,  const CCMWaveformTemplate& wfTemplate, const CCMPMTCalibration& calibration, const double spe_charge) {

    boost::shared_ptr<CCMRecoPulseSeries> output(new CCMRecoPulseSeries);
    cholmod_triplet *basis_trip;
    cholmod_sparse *basis;
    cholmod_dense *data, *unfolded;
    unsigned j, k;
    int nbins = 0;

    // Determine the total number of WF bins
    nbins = wf.GetWaveform().size();

    // If we have no data, nothing to do
    if (nbins == 0 || !std::isfinite(spe_charge) || spe_charge == 0)
        return output;

    // Original code defined `std::vector<int> sources;` to specify a different source per bin
    // Each of our sensors should only have a single digitizer associated with it,
    //  so every bin will have the same source
    int source; // Bitmask of CCMRecoPulse::PulseFlags
    std::vector<double> redges(nbins+1);
    std::vector<double> weights(nbins);
    std::vector<bool> passBasisThresh(nbins,false);
    // Channel 0 unless otherwise noted
    data = cholmod_l_zeros(nbins, 1, CHOLMOD_REAL, &c);

    // Fill data vector
    // Set pulse flags to use for this waveform
    if (wf.GetSource() == CCMSource::V1730) 
        source = CCMSource::V1730;
    else {
        log_error("Unknown waveform source (%d), assuming V1730 " 
                "pulse template", wf.GetSource());
        source = CCMSource::V1730;
    }

    double base_weight = 1;

    // If the waveform is shorter than a pulse width,
    // increase the per-bin weights so the aggregate weight
    // is closer to what it should be
    if (wf.GetWaveform().size() * wf.GetBinWidth() < PulseWidth()){
        base_weight *= PulseWidth() / (wf.GetWaveform().size() * wf.GetBinWidth());
    }

    double noise = noise_threshold_;
    double basisThreshmV = basis_threshold_;

    // Read waveform
    for (k = 0; k < wf.GetWaveform().size(); k++) {
        redges[k] = wf.GetStartTime() + (1. + k) * wf.GetBinWidth();
        ((double *)(data->x))[k] = wf.GetWaveform()[k] / I3Units::mV;

        weights[k] = base_weight;

        // Remove waveform bins that were crazy for some reason
        if (!std::isfinite(((double *)(data->x))[k])) {
            ((double *)(data->x))[k] = 0;
            weights[k] = 0;
        }

        // Deweight and zero below noise-floor bins
        if (fabs(((double *)(data->x))[k]) < noise) {
            ((double *)(data->x))[k] = 0;
            weights[k] /= 4.;
        } else if (fabs(((double *)(data->x))[k]) > basisThreshmV) {
            passBasisThresh[k] = true;
        }
    }

    // Remove saturated bins
    for (k = 0; k < wf.GetWaveformInformation().size(); k++) {
        const CCMWaveformDouble::StatusCompound &status = wf.GetWaveformInformation()[k];

        if (status.GetStatus() & CCMStatus::SATURATED) {
            for (uint64_t a = status.GetInterval().first; a < status.GetInterval().second; ++a){
                weights[a] = 0;
            }
        }
    }


    // Apply weights
    for (int i = 0; i < nbins; i++) {
        ((double *)(data->x))[i] *= weights[i];
    }

    // Precalculate the max number of basis functions to avoid reallocation
    int maxspes = int(spes_per_bin_*(wf.GetWaveform().size()));

    std::vector<std::pair<double, double> > start_times;
    start_times.reserve(maxspes);
    double min_spe_spacing = DBL_MAX;

    // Start, end two bins early
    double present = wf.GetStartTime() - 2.*wf.GetBinWidth();
    double max = present + wf.GetWaveform().size()*wf.GetBinWidth();

    double spacing = wf.GetBinWidth() / spes_per_bin_;
    if (spacing < min_spe_spacing) {
        min_spe_spacing = spacing;
    }

    // Get the peak FWHM start, stop times for this waveform
    double fwhmStart = wfTemplate.digitizerStart;
    double fwhmStop = wfTemplate.digitizerStop;

    // Generate the set of pulse start times that have reasonable
    // non-zero data within the FWHM of the corresponding pulse
    for (k = 0; k < wf.GetWaveform().size(); ++k) {
        if (((double *)(data->x))[k] != 0. && passBasisThresh[k]) {

            // Move the present time forward if necessary
            double binTime = redges[k];
            // Don't jump if we're moving less than the basis spacing
            if (present < (binTime - fwhmStop - spacing)) {
                present = binTime - fwhmStop;
            }

            // Add start times to the set
            while (present < (binTime - fwhmStart) && present < max) {
                start_times.push_back(std::pair<double, double>(present, spacing));
                present += spacing;
            }
        }
    }

    int nspes = start_times.size();
    if (nspes == 0) {
        // We have no nonzero data left
        return output;
    }
    std::sort(start_times.begin(), start_times.end());

    // Deduplicate SPE start times to improve basis matrix conditioning.
    // Due to the superposition of multiple basis function rows, multiple
    // basis functions can become interleaved or even duplicate each other.
    // This can slow down the solution and lead to pulse splitting, and
    // is at best redundant.
    {
        int i, j;
        for (i = 0, j = 0; i < nspes-1; i++) {
            // NB: i is always >= j, so the below is always safe
            // min_spe_spacing is inflated slightly here to avoid
            // round-off error
            if (start_times[i+1].first - start_times[j].first >=
                    0.9 * min_spe_spacing)
                start_times[++j] = start_times[i+1];
        }
        nspes = j+1;
    }
    start_times.resize(nspes);

    // Recompute spacings based on deduplicated basis
    if (nspes > 1) {
        int i;
        // Compute rightdiff
        for (i = 0; i < nspes-1; ++i) {
            start_times[i].second =
                start_times[i+1].first - start_times[i].first;
        }
        // Compute leftdiff and use the smaller
        for (i = 1; i < nspes; ++i) {
            double leftdiff =
                start_times[i].first - start_times[i-1].first;
            if (leftdiff < start_times[i].second) {
                start_times[i].second = leftdiff;
            }
        }
    }

    // Don't use data bins that are not in the support of any basis function
    double start = DBL_MIN;
    double end = DBL_MIN;

    // commented out this check since we're no longer looping over the different wfs the continue statement doesnt make sense
    //if (wf.GetWaveform().size() == 0) {
    //	continue;
    //}

    k = 0;
    for (std::vector<std::pair<double, double> >::const_iterator it = start_times.begin();
            it != start_times.end(); ++it) {
        start = it->first + PulseMin();
        end = it->first + PulseWidth();

        // Evaluate bins up until we pass the end of the current time range
        for (; k < wf.GetWaveform().size() && redges[k] < end; ++k) {
            // Check if bin time is inside the max of
            // the last range and the min of the current range
            if (redges[k] < start) {
                weights[k] = 0.;
            }
        }

        if (k == wf.GetWaveform().size()) {
            break;
        }
    }

    while (k < wf.GetWaveform().size()) {
        // Get rid of bins with times later than the end of the support of the
        // last basis function
        weights[k] = 0.;
        ++k;
    }

    // Chop out the data bins that we don't use since we don't fit them
    j = 0;
    for (int i = 0; i < nbins; ++i) {
        if (weights[i] > 0.) {
            weights[j] = weights[i];
            redges[j] = redges[i];
            //sources[j] = sources[i];
            ((double *)(data->x))[j] = ((double *)(data->x))[i];
            ++j;
        }
    }
    nbins = j;
    data->nrow = nbins;

    // Compute a reasonable upper bound on the number of non-zero matrix
    // elements.
    long nzmax = 0;
    int first_spe = 0;
    double last_t = 0;
    for (int i = 0; i < nbins; i++) {
        // See below for comments about how this loop works
        if (redges[i] < last_t)
            first_spe = 0;
        last_t = redges[i];
        while (first_spe < nspes-1 && redges[i] -
                start_times[first_spe].first > PulseWidth())
            first_spe++;
        for (int j = first_spe; j < nspes; j++) {
            if (((redges[i] - start_times[j].first) -
                        PulseMin()) < -template_bin_spacing_)
                break;
            nzmax++;
        }
    }

    // Create model matrix
    basis_trip = cholmod_l_allocate_triplet(nbins, nspes, nzmax, 0,
            CHOLMOD_REAL, &c);
    basis_trip->nnz = 0;

    first_spe = 0;
    last_t = 0;
    std::vector<int> col_counts(nspes,0);
    std::vector<int> pflags(nspes,0);

    for (int i = 0; i < nbins; i++) {
        // SPEs are sorted in time, as are the bins (though in blocks)
        // Exploit this to minimize the amount of the matrix through
        // which we have to loop.
        if (redges[i] < last_t)
            first_spe = 0;
        last_t = redges[i];

        // The earliest pulse influencing this bin is PULSE_WIDTH in the past.
        // The template is defined up to (but not including) PULSE_WIDTH.
        while (first_spe < nspes && redges[i] -
                start_times[first_spe].first >= PulseWidth())
            first_spe++;
        if (first_spe == nspes) {
            continue;
        }

        // Precompute the multiplicative term in the basis
        double weighted_charge = (spe_charge/I3Units::mV) * weights[i];
        if (weighted_charge == 0) // Removed from fit, so ignore
            continue;

        // Precache which pulse template we're using
        double templ_bin_spacing_inv = 1./template_bin_spacing_;
        const std::vector<double>* pulse_templ = &(wfTemplate.digitizer_template);

        // The last pulse for this bin is 2 ns in the future
        for (int j = first_spe; j < nspes; j++) {
            int templ_bin = int(((redges[i] - start_times[j].first) -
                        PulseMin())*templ_bin_spacing_inv);
            if (templ_bin < 0)
                break;

            ((long *)(basis_trip->i))[basis_trip->nnz] = i;
            ((long *)(basis_trip->j))[basis_trip->nnz] = j;
            ((double *)(basis_trip->x))[basis_trip->nnz] =
                (*pulse_templ)[templ_bin] * weighted_charge;
            //pflags[j] |= sources[i];
            pflags[j] |= source;

            basis_trip->nnz++;
            assert(basis_trip->nnz <= basis_trip->nzmax);
            col_counts[j]++;
        }
    }

    //  Convert to column-ordered sparse matrix
    //  Note: This is handrolled instead of using
    //  cholmod_l_triplet_to_sparse() in order to exploit some of the
    //  structure of our specific triplet matrix, which lets this
    //  run in less than one third the time of cholmod_l_triplet_to_sparse.
    basis = cholmod_l_allocate_sparse(basis_trip->nrow, basis_trip->ncol,
            basis_trip->nnz, true, true, 0, CHOLMOD_REAL, &c);
    int accum = 0;
    for (int i = 0; i < nspes; i++) {
        ((long *)(basis->p))[i] = accum;
        accum += col_counts[i];
    }
    // Need to set data end pointer for the last column.  Otherwise
    // SuiteSparse will ignore the last column.
    ((long *)(basis->p))[nspes] = accum;
    std::vector<long> col_indices(nspes,0);

    for (unsigned i = 0; i < basis_trip->nnz; i++) {
        long col = ((long *)(basis_trip->j))[i];
        long index = ((long *)(basis->p))[col] + col_indices[col]++;
        ((long *)(basis->i))[index] = ((long *)(basis_trip->i))[i];
        ((double *)(basis->x))[index] = ((double *)(basis_trip->x))[i];
    }
    cholmod_l_free_triplet(&basis_trip, &c);

    // Solve for SPE heights
    if (reduce_) {
        unfolded = rnnls(basis, data, tolerance_, 1000, 0, &c);
    } else {
        unfolded = nnls_lawson_hanson(basis, data, tolerance_,
                0, 1000, nspes, 0, 1, 0, &c);
    }

    cholmod_l_free_sparse(&basis, &c);
    cholmod_l_free_dense(&data, &c);

    // Load the SPE corrections
    double speCorrection = 1.;
    if (apply_spe_corr_) {
        speCorrection = 1. / calibration.GetMeanPMTCharge();
    }

    // Convert to pulse series
    for (int i = 0; i < nspes; i++) {
        if (((double *)(unfolded->x))[i] == 0)
            continue;

        CCMRecoPulse pulse;

        pulse.SetTime(start_times[i].first);
        pulse.SetCharge((((double *)(unfolded->x))[i]) * speCorrection);
        pulse.SetWidth(start_times[i].second);
        output->push_back(pulse);
    }
    cholmod_l_free_dense(&unfolded, &c);
    return output;
}

namespace {
/* Simple routine to calculate FWHM start and stop given a waveform template */
void FillFWHM(double& start, double& stop,
        const std::vector<double>& data, double spacing, double min) {

    // Find the maximum
    double max = DBL_MIN;
    for (unsigned i = 0; i < data.size(); ++i) {
        if (data[i] > max) {
            max = data[i];
        }
    }

    double hm = max*0.5;
    unsigned hmbin = 0;
    for (unsigned i = 0; i < data.size(); ++i) {
        if (data[i] > hm) {
            hmbin = i;
            break;
        }
    }

    // Be conservative: choose the bin before the first rise above hm
    if (hmbin > 0) {
        --hmbin;
    }

    start = min + hmbin*spacing;

    for (int i = data.size() - 1; i >= 0; --i) {
        if (data[i] > hm) {
            hmbin = i;
            break;
        }
    }

    if (hmbin < data.size() - 1) {
        ++hmbin;
    }

    stop = min + hmbin*spacing;
    if (stop == start) {
        stop += spacing;
    }
}
} // namespace

void CCMWavedeform::FillTemplate(CCMWaveformTemplate& wfTemplate,
        const CCMPMTCalibration& calibration) {
    CCMPMTCalibration::DroopedSPETemplate chan_template =
        calibration.PMTPulseTemplate();
    wfTemplate.digitizer_template.resize(template_bins_);
    for (int i = 0; i < template_bins_; i++) {
        wfTemplate.digitizer_template[i] =
            chan_template(PulseMin() +
                    i*template_bin_spacing_);
    }

    FillFWHM(wfTemplate.digitizerStart,
            wfTemplate.digitizerStop,
            wfTemplate.digitizer_template,
            template_bin_spacing_, PulseMin());

    wfTemplate.filled = true;
}

