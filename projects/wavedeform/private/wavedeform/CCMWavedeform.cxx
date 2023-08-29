#include <cctype>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <thread>
#include <mutex>
#include <cholmod.h>

#include <boost/make_shared.hpp>

#include <icetray/I3Frame.h>
#include <icetray/ctpl.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/I3ConditionalModule.h>

#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/geometry/CCMOMGeo.h>
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
void CCMFillFWHM(double& start, double& stop, const std::vector<double>& data, double spacing, double min) {

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

void GetPulses(CCMWaveformDouble const & wf, CCMWaveformTemplate const & wfTemplate, CCMPMTCalibration const & calibration, double spe_charge, double start_time, double pulse_width, double template_bin_spacing_, double noise_threshold_, double basis_threshold_, double spes_per_bin_, bool reduce_, double tolerance_, bool apply_spe_corr_, cholmod_common & chol_common, CCMRecoPulseSeries & output, I3FramePtr frame) {
    output.clear();
    cholmod_triplet *basis_trip;
    cholmod_sparse *basis;
    cholmod_dense *data, *unfolded;
    unsigned j, k;
    int nbins = 0;

    // Determine the total number of WF bins
    nbins = wf.GetWaveform().size();
    // If we have no data, nothing to do
    if (nbins == 0 || !std::isfinite(spe_charge) || spe_charge == 0)
        return;

    // Original code defined `std::vector<int> sources;` to specify a different source per bin
    // Each of our sensors should only have a single digitizer associated with it,
    //  so every bin will have the same source
    int source; // Bitmask of CCMRecoPulse::PulseFlags
    std::vector<double> redges(nbins+1);
    std::vector<double> weights(nbins);
    std::vector<bool> passBasisThresh(nbins,false);
    // Channel 0 unless otherwise noted
    data = cholmod_l_zeros(nbins, 1, CHOLMOD_REAL, &chol_common);

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
    // HARDCODING WF BIN WIDTH TO AVOID NANS
    double wf_bin_width = 2; // 2 nsec

    // If the waveform is shorter than a pulse width,
    // increase the per-bin weights so the aggregate weight
    // is closer to what it should be
    if (wf.GetWaveform().size() * wf_bin_width < pulse_width){
        base_weight *= pulse_width / (wf.GetWaveform().size() * wf_bin_width);
    }

    double noise = noise_threshold_;
    double basisThreshmV = basis_threshold_;

    // Read waveform
    for (k = 0; k < wf.GetWaveform().size(); k++) {
        redges[k] = (1. + k) * wf_bin_width;
        ((double *)(data->x))[k] = wf.GetWaveform()[k];

        weights[k] = base_weight;

        // Remove waveform bins that were crazy for some reason
        if (!std::isfinite(((double *)(data->x))[k])) {
            ((double *)(data->x))[k] = 0;
            weights[k] = 0;
        }

        // Deweight and zero below noise-floor bins
        if (((double *)(data->x))[k] < noise) {
            ((double *)(data->x))[k] = 0;
            weights[k] /= 4.;
        } else if (fabs(((double *)(data->x))[k]) > basisThreshmV) {
            passBasisThresh[k] = true;
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
    double present = - 2. * wf_bin_width;
    double max = present + wf.GetWaveform().size() * wf_bin_width;

    double spacing = wf_bin_width / spes_per_bin_;
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
        return;
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

    k = 0;
    for (std::vector<std::pair<double, double> >::const_iterator it = start_times.begin();
            it != start_times.end(); ++it) {
        start = it->first + start_time;
        end = it->first + start_time + pulse_width;

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
        if (((double *)(data->x))[i] > 0.) {
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
                start_times[first_spe].first > pulse_width)
            first_spe++;
        for (int j = first_spe; j < nspes; j++) {
            if (((redges[i] - start_times[j].first) -
                        start_time) < -template_bin_spacing_)
                break;
            nzmax++;
        }
    }

    // Create model matrix
    basis_trip = cholmod_l_allocate_triplet(nbins, nspes, nzmax, 0,
            CHOLMOD_REAL, &chol_common);
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
                start_times[first_spe].first - start_time >= pulse_width)
            first_spe++;
        if (first_spe == nspes) {
            continue;
        }

        // Precompute the multiplicative term in the basis
        //double weighted_charge = (spe_charge/I3Units::mV) * weights[i];
        double weighted_charge = spe_charge * weights[i];
        if (weighted_charge == 0) // Removed from fit, so ignore
            continue;

        // Precache which pulse template we're using
        double templ_bin_spacing_inv = 1./template_bin_spacing_;
        std::vector<double> const & pulse_templ = wfTemplate.digitizer_template;

        // The last pulse for this bin is 2 ns in the future
        for (int j = first_spe; j < nspes; j++) {
            int templ_bin = int(((redges[i] - start_times[j].first) - start_time)*templ_bin_spacing_inv);
            if (templ_bin < 0)
                break;

            ((long *)(basis_trip->i))[basis_trip->nnz] = i;
            ((long *)(basis_trip->j))[basis_trip->nnz] = j;
            ((double *)(basis_trip->x))[basis_trip->nnz] = pulse_templ.at(templ_bin) * weighted_charge;

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
            basis_trip->nnz, true, true, 0, CHOLMOD_REAL, &chol_common);
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
    cholmod_l_free_triplet(&basis_trip, &chol_common);
    // Solve for SPE heights
    if (reduce_) {
        unfolded = rnnls(basis, data, tolerance_, 2000, 0, &chol_common);
    } else {
        unfolded = nnls_lawson_hanson(basis, data, tolerance_,
                0, 1000, nspes, 0, 1, 0, &chol_common);
    }

    cholmod_l_free_sparse(&basis, &chol_common);
    cholmod_l_free_dense(&data, &chol_common);

    // Load the SPE corrections
    double speCorrection = 1.;
    if (apply_spe_corr_) {
        //speCorrection = 1. / calibration.GetMeanPMTCharge();
        speCorrection = 1.;
    }

    // Convert to pulse series
    for (int i = 0; i < nspes; i++) {
        if (((double *)(unfolded->x))[i] == 0)
            continue;

        CCMRecoPulse pulse;

        pulse.SetTime(start_times[i].first);
        pulse.SetCharge((((double *)(unfolded->x))[i]) * speCorrection);
        pulse.SetWidth(start_times[i].second);
        output.push_back(pulse);
    }
    cholmod_l_free_dense(&unfolded, &chol_common);
}

void RunPulsesThread(
        std::tuple<size_t, size_t> thread_range,
        std::vector<CCMPMTKey> const & pmt_keys,
        I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map,
        CCMWaveformDoubleSeries const & waveforms,
        I3Map<CCMPMTKey, CCMWaveformTemplate> const & templates,
        CCMCalibration const & calibration,
        std::vector<cholmod_common> & cholmod_common_vec,
        bool reduce,
        std::vector<std::reference_wrapper<CCMRecoPulseSeries>> & output_pulses_references,
        double start_time,
        double pulse_width,
        double template_bin_spacing_,
        double noise_threshold_,
        double basis_threshold_,
        double spes_per_bin_,
        double tolerance_,
        bool apply_spe_corr_,
        I3FramePtr frame
        ) {
    for(size_t i=std::get<0>(thread_range); i<std::get<1>(thread_range); ++i) {
        CCMPMTKey pmt_key = pmt_keys.at(i);
        size_t channel = pmt_channel_map.at(pmt_key);
        CCMWaveformDouble const & waveform = waveforms.at(channel);

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(pmt_key);

        if(calib == calibration.pmtCal.cend()) {
            //std::cout << "oops! no calibration for " << pmt_key << std::endl;
            continue;
        }

        double start_time = 2 * -10;
        double end_time = 2 * 60;
        double pulse_width = end_time - start_time;

        double template_bin_spacing_ = 2.0 / spes_per_bin_ / 10;
        double range = end_time - start_time;

        double placeholder = 1.0;

        GetPulses(
                waveform,
                templates.at(pmt_key),
                calib->second,
                placeholder,
                start_time,
                pulse_width,
                template_bin_spacing_,
                noise_threshold_,
                basis_threshold_,
                spes_per_bin_,
                reduce,
                tolerance_,
                apply_spe_corr_,
                cholmod_common_vec.at(i),
                output_pulses_references[i].get(),
                frame
                );


    }
}

void FillTemplate(CCMWaveformTemplate& wfTemplate, const CCMPMTCalibration& calibration, double const & start_time, double const & pulse_width, int const & template_bins, double const & template_bin_spacing_) {
    CCMSPETemplate channel_template = calibration.GetSPETemplate();
    wfTemplate.digitizer_template.resize(template_bins);
    for (int i = 0; i < template_bins; i++) {
        wfTemplate.digitizer_template[i] =
            channel_template.Evaluate(start_time + i*template_bin_spacing_ ); // hard coding bin spacing...use calib one day 
    }

    CCMFillFWHM(wfTemplate.digitizerStart,
            wfTemplate.digitizerStop,
            wfTemplate.digitizer_template,
            template_bin_spacing_, start_time);

    wfTemplate.filled = true;
}


class CCMWavedeform : public I3ConditionalModule {
    public:
        size_t num_threads;
        bool geo_seen;
        std::string geometry_name_;
        std::string nim_pulses_name_;
        I3Vector<CCMPMTKey> allowed_pmt_keys_;
        I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
        CCMWavedeform(const I3Context &);

        std::vector<std::tuple<size_t, size_t>> thread_ranges_;
        CCMCalibration calibration;

        virtual ~CCMWavedeform();

        void Geometry(I3FramePtr frame);
        void Configure();
        void Calibration(I3FramePtr frame);
        void DAQ(I3FramePtr frame);
    private:
        std::string waveforms_name_;
        std::string waveform_range_name_;
        std::string output_name_;

        double spes_per_bin_;
        double tolerance_;
        double noise_threshold_;
        double basis_threshold_;

        std::map<CCMPMTKey, int> start_bins_;
        std::map<CCMPMTKey, int> end_bins_;
        double template_bin_spacing_;

        bool apply_spe_corr_;
        bool reduce_;

        I3Map<CCMPMTKey, CCMWaveformTemplate> template_;
        std::vector<CCMPMTKey> pmt_keys_;
        std::vector<cholmod_common> cholmod_common_vec_;

        double start_time;
        double pulse_width;
        double end_time;
        double range;
};

I3_MODULE(CCMWavedeform);

CCMWavedeform::CCMWavedeform(const I3Context& context) : I3ConditionalModule(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("NumThreads", "Number of worker threads to use for pulse fitting", (size_t)(0));
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
    AddParameter("SPEsPerBin", "Number of basis functions to unfold per waveform bin", 4.);
    AddParameter("Tolerance", "Stopping tolerance, in units of bin ADC^2/PE", 9.);
    AddParameter("NoiseThreshold","Consider bins with amplitude below this number of counts as noise", 5.0);
    AddParameter("BasisThreshold",
            "Require a bin with amplitude at least this number of counts "
            "within the FWHM of the template waveform in order to include "
            "a given start time in the basis set", 10.0);
    AddParameter("Waveforms", "Name of input waveforms",
            "CCMCalibratedWaveforms");
    AddParameter("WaveformTimeRange", "Name of maximum time range of "
            "calibrated waveforms for this event", "CalibratedWaveformRange");
    AddParameter("Output", "Name of output pulse series",
            "WavedeformPulses");
    AddParameter("ApplySPECorrections", "Whether to apply DOM-by-DOM"
            " corrections to the pulse charge scaling if available", false);
    AddParameter("Reduce", "Find the optimal NNLS solution, then eliminate"
            " basis members until tolerance is reached", true);
    AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>());
    //AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>({CCMPMTKey(3,8,0)}));

}

void CCMWavedeform::Configure() {
    GetParameter("NumThreads", num_threads);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
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
    GetParameter("PMTKeys", allowed_pmt_keys_);
}

CCMWavedeform::~CCMWavedeform() {
    for(size_t i=0; i<cholmod_common_vec_.size(); ++i) {
        cholmod_l_finish(&(cholmod_common_vec_[i]));
    }
}

void CCMWavedeform::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    template_.clear();
    for(auto pmt_channel_pair : pmt_channel_map_) {
        template_[pmt_channel_pair.first] = CCMWaveformTemplate();
        template_[pmt_channel_pair.first].filled = false;
    }

    I3Map<CCMPMTKey, CCMOMGeo> const & pmt_geo_ = geo.pmt_geo;
    std::set<CCMPMTKey> allowed_pmt_keys(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end());
    bool filter_pmts = allowed_pmt_keys.size() > 0;
    for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : pmt_geo_) {
        if(filter_pmts and allowed_pmt_keys.find(p.first) == allowed_pmt_keys.end()) {
            continue;
        }
        if(p.second.omtype == CCMOMGeo::OMType::CCM8inCoated or p.second.omtype == CCMOMGeo::OMType::CCM8inUncoated or p.second.omtype == CCMOMGeo::OMType::CCM1in) {
            pmt_keys_.push_back(p.first);
        }
    }

    size_t min_channels_per_thread = pmt_keys_.size() / num_threads;
    size_t max_channels_per_thread = min_channels_per_thread + 1;
    size_t num_threads_with_max = pmt_keys_.size() - num_threads * min_channels_per_thread;
    size_t channel_start = 0;
    for(size_t i=0; i<num_threads; ++i) {
        size_t n_channels = min_channels_per_thread;
        if(i < num_threads_with_max) {
            n_channels = max_channels_per_thread;
        }
        thread_ranges_.emplace_back(channel_start, channel_start + n_channels);
        channel_start += n_channels;
    }

    cholmod_common_vec_.resize(pmt_keys_.size());
    for(size_t i=0; i<pmt_keys_.size(); ++i) {
        cholmod_l_start(&(cholmod_common_vec_[i]));
    }

    PushFrame(frame);
}

void CCMWavedeform::Calibration(I3FramePtr frame) {
    calibration = frame->Get<CCMCalibration>("CCMPMTCalibration");

    start_time = 2 * -10;
    end_time = 2 * 60;
    pulse_width = end_time - start_time;

    template_bin_spacing_ = 2.0 / spes_per_bin_ / 10;
    range = end_time - start_time;

    for(size_t i=0; i<pmt_keys_.size(); ++i) {
        // let's get our wf and what not
        CCMPMTKey key = pmt_keys_[i];
        uint32_t channel = pmt_channel_map_[key];

        if(calibration.pmtCal.count(key) == 0) {
            //std::cout << "oops! no calibration for " << key << std::endl;
            continue;
        }

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);

        double placeholder = 1.0;

        if(not template_.at(key).filled) {
            int template_bins = (int)ceil(range / template_bin_spacing_);
            template_.at(key).digitizer_template.resize(template_bins);
            FillTemplate(template_.at(key), calib->second, start_time, pulse_width, template_bins, template_bin_spacing_);
        }
    }

    PushFrame(frame);
}

void CCMWavedeform::DAQ(I3FramePtr frame) {
    if (!frame->Has(waveforms_name_)) {
        PushFrame(frame);
        return;
    }

    //const CCMDetectorStatus& status = frame->Get<CCMDetectorStatus>();
    boost::shared_ptr<const CCMWaveformDoubleSeries> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>(waveforms_name_);
    I3Map<CCMPMTKey, BaselineEstimate> const & baselines = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    if(waveforms == nullptr) {
        log_warn("CCMWaveformDoubleSeries named %s not present in frame", "CCMCalibratedWaveforms");
    }

    size_t num_pulse_series = pmt_keys_.size();

    if(num_threads == 0) {
        num_threads = pmt_keys_.size();
    }

    std::vector<std::thread> threads;
    threads.reserve(thread_ranges_.size());

    // place to store pulses
    boost::shared_ptr<CCMRecoPulseSeriesMap> output(new CCMRecoPulseSeriesMap);
    for(size_t i = 0; i < pmt_keys_.size(); ++i) {
        output->operator[](pmt_keys_[i]) = CCMRecoPulseSeries();
    }

    std::vector<std::reference_wrapper<CCMRecoPulseSeries>> output_pulses_references;
    output_pulses_references.reserve(pmt_keys_.size());
    for(size_t i = 0; i < pmt_keys_.size(); ++i) {
        output_pulses_references.emplace_back(output->at(pmt_keys_[i]));
    }

    // loop over each channel in waveforms
    for(size_t i=0; i<num_threads; ++i) {
        threads.emplace_back(
            RunPulsesThread,
            thread_ranges_[i],
            std::cref(pmt_keys_),
            std::cref(pmt_channel_map_),
            std::cref(*waveforms),
            std::cref(template_),
            std::cref(calibration),
            std::ref(cholmod_common_vec_),
            reduce_,
            std::ref(output_pulses_references),
            start_time,
            pulse_width,
            template_bin_spacing_,
            noise_threshold_,
            basis_threshold_,
            spes_per_bin_,
            tolerance_,
            apply_spe_corr_,
            frame
        );
    }

    for(size_t i=0; i<num_threads; ++i) {
        threads[i].join();
    }

    frame->Put(output_name_, output);

    PushFrame(frame, "OutBox");
}



