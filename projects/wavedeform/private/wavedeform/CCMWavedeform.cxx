#include <cctype>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <thread>
#include <mutex>
#include <stack>
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

#include "timer.h"

/* Simple class to hold together ATWD and FADC
 * templates along with a validity flag and waveform features */
struct CCMWaveformTemplate {

    std::vector<double> digitizer_template;
    double digitizerStart;
    double digitizerStop;
    double start_time;
    double end_time;
    double pulse_width;
    double total_mass;
    bool filled;
};

void rebin_bayesian_blocks(std::vector<double> const & raw_bins,
    double const * raw_charges, std::vector<std::vector<size_t>> & bin_indices,
    double ncp_prior, double poisson_scale_factor) {
    size_t N_raw = raw_bins.size() - 1;
    std::vector<double> best(N_raw, -std::numeric_limits<double>::infinity());
    std::vector<int> last(N_raw, 0);
    
    // sum of raw_charges[:k+1]
    double total = 0;
    for(int k = 0; k < N_raw; k++) {
        total += raw_charges[k] * poisson_scale_factor;
        // sum of raw_charges[:j]
        double subtotal = 0.;
        for(int j = 0; j < k+1; j++) {
            // sum of raw_charges[j:k+1] (contents of proposedblock)
            double counts = total - subtotal;
            // width of proposed block
            double width = raw_bins[k+1]-raw_bins[j];
            // the fitness of the block is a saturated Poisson
            // likelihood, penalized by a constant factor for each
            // extra block
            double fitness = (counts > 0 ? counts*(std::log(counts) - std::log(width)) : 0) + (j > 0 ? best[j-1] : 0) - ncp_prior;
            if (fitness > best[k]) {
                best[k] = fitness;
                last[k] = j;
            }
            subtotal += raw_charges[j] * poisson_scale_factor;
        }
    }
    
    // Recover changepoints by iteratively peeling off the last block
    std::stack<int> change_points;
    for(int i = N_raw; i > 0; i = last[i-1])
        change_points.push(i);
    
    // Fill contents of blocks into output
    int p;
    for(p = 0; p < N_raw; change_points.pop()) {
        bin_indices.emplace_back();
        assert(!change_points.empty());
        while(p < change_points.top()) {
            bin_indices.back().push_back(p);
            p++;
        }
    }
    // Don't include the last bin edge index because the output is for specifying which bins to merge
    // bin_indices.back().push_back(p);
}

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

void GetPulses(CCMWaveformDouble const & wf, CCMWaveformTemplate const & wfTemplate, CCMPMTCalibration const & calibration, double spe_charge, double template_bin_spacing_, double noise_threshold_, double basis_threshold_, double spes_per_bin_, bool reduce_, double tolerance_, bool apply_spe_corr_, cholmod_common & chol_common, CCMRecoPulseSeries & output, std::vector<double> & output_data_times, std::vector<double> & output_rebin_data_times, I3FramePtr frame) {
    double GetPulsesInternal_s = 0.0;
    double GetPulsesInternal_u = 0.0;
    NNLSTimer GetPulsesInternal_timer("GetPulsesInternal", GetPulsesInternal_s, GetPulsesInternal_u, true);
    output.clear();
    cholmod_triplet *basis_trip;
    cholmod_sparse *basis;
    cholmod_dense *data, *unfolded;
    unsigned j, k;
    int nbins = 0;

    double dummy_s = 0.0;
    double dummy_u = 0.0;
    //double eigen_s = 0.0;
    //double eigen_u = 0.0;
    //NNLSTimer eigen_timer("Eigen", eigen_s, eigen_u, false);
    NNLSTimer pre_compute("Pre compute total", dummy_s, dummy_u, true);

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
    if (wf.GetWaveform().size() * wf_bin_width < wfTemplate.pulse_width){
        base_weight *= wfTemplate.pulse_width / (wf.GetWaveform().size() * wf_bin_width);
    }

    double noise = noise_threshold_;
    double basisThreshmV = basis_threshold_;

    std::vector<double> data_times;
    data_times.reserve(wf.GetWaveform().size());

    // Read waveform
    for (k = 0; k < wf.GetWaveform().size(); k++) {
        redges[k] = (1. + k) * wf_bin_width;
        ((double *)(data->x))[k] = wf.GetWaveform()[k];
        data_times.push_back(k * 2.0);

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

    double fine_spacing = wf_bin_width / spes_per_bin_;
    if (fine_spacing < min_spe_spacing) {
        min_spe_spacing = fine_spacing;
    }

    double coarse_spacing = wf_bin_width;
    double spacing = coarse_spacing;

    // Get the peak FWHM start, stop times for this waveform
    double fwhmStart = wfTemplate.digitizerStart;
    double fwhmStop = wfTemplate.digitizerStop;
    // Generate the set of pulse start times that have reasonable
    // non-zero data within the FWHM of the corresponding pulse
    bool prev_below = true;
    size_t fine_spacing_tau = 1;
    size_t n_since_fine = 0;
    for (k = 0; k < wf.GetWaveform().size(); ++k) {
        if (((double *)(data->x))[k] != 0. && passBasisThresh[k]) {
            n_since_fine += 1;
            if(prev_below) {
                n_since_fine = 0;
                spacing = fine_spacing;
            } else if(n_since_fine > fine_spacing_tau) {
                spacing = coarse_spacing;
            }
            double binTime = redges[k];
            // Don't jump if we're moving less than the basis spacing
            if (present + fwhmStop < (binTime - fine_spacing)) {
                present = binTime - fwhmStop;
            }

            // Add start times to the set
            while (present < (binTime - fwhmStart) && present < max) {
                start_times.push_back(std::pair<double, double>(present, spacing));
                present += spacing;
            }
            prev_below = false;
        } else {
            prev_below = true;
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

    /*
    // now let's save the start times after de-duplication
    boost::shared_ptr<I3Vector<double>> start_times_copy = boost::make_shared<I3Vector<double>>();

    for (size_t i = 0; i < start_times.size(); ++i) {
        start_times_copy->push_back(start_times[i].first);
    }
    frame->Put("StartTimes" , start_times_copy);

    */

    //Don't use data bins that are not in the support of any basis function
    double start = DBL_MIN;
    double end = DBL_MIN;

    k = 0;
    for (std::vector<std::pair<double, double> >::const_iterator it = start_times.begin();
            it != start_times.end(); ++it) {
        start = it->first + wfTemplate.start_time;
        end = it->first + wfTemplate.end_time;

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
            data_times[j] = data_times[i];
            //sources[j] = sources[i];
            ((double *)(data->x))[j] = ((double *)(data->x))[i];
            ++j;
        }
    }
    nbins = j;
    data->nrow = nbins;
    data_times.resize(nbins);
    data_times.push_back(data_times.back() + 2.0);

    output_data_times = data_times;
    std::vector<std::vector<size_t>> rebin_ranges;
    std::vector<double> rebinned_data;
    double ncp_prior = 4.0;
    rebin_bayesian_blocks(data_times, ((double *)(data->x)), rebin_ranges, ncp_prior, wfTemplate.total_mass);
    output_rebin_data_times.clear();
    for(size_t i=0; i<rebin_ranges.size(); ++i) {
        output_rebin_data_times.push_back(rebin_ranges[i][0]);
        rebinned_data.emplace_back(0);
        for(size_t rebin_idx : rebin_ranges[i]) {
            rebinned_data.back() += ((double *)(data->x))[rebin_idx];
        }
    }
    output_rebin_data_times.push_back(rebin_ranges.back().back());

    std::cout << "Original data bins: " << data_times.size() - 1 << std::endl;
    std::cout << "Rebinned data bins: " << output_rebin_data_times.size() - 1 << std::endl;
    std::cout << "Rebinned last index: " << rebin_ranges.back().back() << std::endl;
    std::cout << std::endl;

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
                start_times[first_spe].first - wfTemplate.start_time > wfTemplate.end_time - wfTemplate.pulse_width)
            first_spe++;
        for (int j = first_spe; j < nspes; j++) {
            if (((redges[i] - start_times[j].first) - wfTemplate.start_time) < -template_bin_spacing_)
                break;
            nzmax++;
        }
    }
    nzmax = (int(wfTemplate.digitizer_template.size() / spes_per_bin_) + 1) * nspes;

    //std::vector<double> SBB_seed_vector;
    size_t data_idx = 0;
    double template_avg = 0.0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        template_avg += wfTemplate.digitizer_template[i];
    }
    template_avg /= wfTemplate.digitizer_template.size();
    for(size_t i=0; i<start_times.size(); ++i) {
        double start_time = start_times[i].first;
        double spacing = start_times[i].second;
        while(data_idx < data_times.size() - 1 and data_times[data_idx] <= start_time) {
            data_idx += 1;
        }
        if(data_idx > 0 and std::abs(data_times[data_idx - 1] - start_time) < std::abs(data_times[data_idx] - start_time)) {
            data_idx -= 1;
        }
        double seed_value = ((double *)data->x)[data_idx] / (template_avg * (spes_per_bin_ * template_bin_spacing_) * spacing);
        //SBB_seed_vector.push_back(seed_value);
    }
    //nsNNLS::vector SBB_seed(SBB_seed_vector.size(), SBB_seed_vector.data());


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
                start_times[first_spe].first - wfTemplate.start_time >= wfTemplate.pulse_width)
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
            int templ_bin = int(((redges[i] - start_times[j].first) - wfTemplate.start_time)*templ_bin_spacing_inv);
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

    {
    int k = 0;
    double data_running_total = 0.0;
    std::vector<std::vector<size_t>> spe_original_support;
    std::vector<std::vector<size_t>> spe_bb_support;
    std::vector<double> rebinned_data;
    std::vector<size_t> first_spe_idxs;
    std::vector<std::vector<double>> rebinned_entries;
    // Merge entries in the matrix according to the data rebinning
    for(size_t bb_idx=0; bb_idx<rebin_ranges.size(); ++bb_idx) {
        rebinned_entries.emplace_back();
        for(size_t data_idx : rebin_ranges[bb_idx]) {
            while(((long *)(basis_trip->i))[k] < data_idx) {
                ++k;
            }
            assert(((long *)(basis_trip->i))[k] == data_idx);
            size_t spe_idx = ((long *)(basis_trip->j))[k];
            first_spe_idxs.push_back(spe_idx);
            while(((long *)(basis_trip->i))[k] == data_idx) {
                while(rebinned_entries.back().size() < spe_idx - first_spe_idxs.back() + 1) {
                    rebinned_entries.back().emplace_back(0);
                }
                while(spe_original_support.size() < spe_idx + 1) {
                    spe_original_support.emplace_back();
                    spe_bb_support.emplace_back();
                }
                spe_original_support[spe_idx].push_back(data_idx);
                if(spe_bb_support[spe_idx].size() > 0 and spe_bb_support[spe_idx].back() != bb_idx) {
                    spe_bb_support[spe_idx].push_back(bb_idx);
                }
                rebinned_entries.back()[spe_idx - first_spe_idxs.back()] += ((double *)(basis_trip->x))[k];
                ++k;
                spe_idx = ((long *)(basis_trip->j))[k];
            }
        }
    }

    std::vector<double> new_spe_start_times;
    // Store representative SPE index, and a vector of combined SPE times and ratios
    std::vector<std::tuple<size_t, std::vector<std::tuple<double, double>>>> merged_spes;

    // Determine if any SPEs have the same support for the rebinned data and merge them accordingly
    // We want to leave the leading edge SPEs alone
    // Only merge if we have 3 or more SPEs in a row with the same support, and merge all SPEs except the first one
    // We must choose a relative magnitude for the merged SPEs when combining them and when uncombining after fitting,
    //   we can use the dot product of the SPE row with the data to determine this relative magnitude
    size_t num_same_support = 1;
    size_t matching_support_idx = 0;
    for(size_t spe_idx=1; spe_idx<=spe_bb_support.size(); ++spe_idx) {
        if(spe_idx<spe_bb_support.size() and spe_bb_support[spe_idx] == spe_bb_support[matching_support_idx]) {
            num_same_support += 1;
        } else {
            // The support has now changed
            // Store the first spe
            merged_spes.emplace_back(new_spe_start_times.size(), std::vector<std::tuple<double, double>>(1, {start_times[matching_support_idx].first, 1.0}));
            new_spe_start_times.emplace_back(start_times[matching_support_idx].first);
            // Check for how many identical supports we had in a row
            if(num_same_support >= 3) {
                // Merge all but the first SPE with matching support
                double total_dot_product = 0.0;
                std::vector<double> dot_products;
                dot_products.reserve(num_same_support-1);
                for(size_t j=spe_idx-num_same_support+1; j<spe_idx; ++j) {
                    double dot_product = 0.0;
                    for(size_t i_sub=0; i_sub<spe_original_support[matching_support_idx].size(); ++i_sub) {
                        size_t i = i_sub + spe_original_support[matching_support_idx][0];
                        dot_product += rebinned_entries[j][i_sub] * ((double *)(data->x))[i];
                    }
                    dot_products.emplace_back(dot_product);
                    total_dot_product += dot_product;
                }
                size_t j_sub = 0;
                std::vector<std::tuple<double, double>> degenerate_spes;
                for(size_t j=spe_idx-num_same_support+1; j<spe_idx; ++j, ++j_sub) {
                    degenerate_spes.emplace_back(start_times[j].first, dot_products[j_sub] / total_dot_product);
                }
            }
            num_same_support = 1;
            matching_support_idx = spe_idx;
        }
    }
    }

    /*
    // let's save our basis trip vector
    boost::shared_ptr<I3Vector<long int>> basis_trip_i_copy = boost::make_shared<I3Vector<long int>>();
    boost::shared_ptr<I3Vector<long int>> basis_trip_j_copy = boost::make_shared<I3Vector<long int>>();
    boost::shared_ptr<I3Vector<double>> basis_trip_x_copy = boost::make_shared<I3Vector<double>>();

    for (unsigned i = 0; i < basis_trip->nnz; ++i){
        basis_trip_i_copy->push_back(((long int *)(basis_trip->i))[i]);
        basis_trip_j_copy->push_back(((long int *)(basis_trip->j))[i]);
        basis_trip_x_copy->push_back(((double *)(basis_trip->x))[i]);
    }

    frame->Put("BasisTripI", basis_trip_i_copy);
    frame->Put("BasisTripJ", basis_trip_j_copy);
    frame->Put("BasisTripX", basis_trip_x_copy);
    */

    //  Convert to column-ordered sparse matrix
    //  Note: This is handrolled instead of using
    //  cholmod_l_triplet_to_sparse() in order to exploit some of the
    //  structure of our specific triplet matrix, which lets this
    //  run in less than one third the time of cholmod_l_triplet_to_sparse.
    basis = cholmod_l_allocate_sparse(basis_trip->nrow, basis_trip->ncol,
            basis_trip->nnz, true, true, 0, CHOLMOD_REAL, &chol_common);


    //std::vector<size_t> SBB_ridx(basis_trip->nnz);
    //std::vector<size_t> SBB_cptr(basis_trip->ncol + 1);
    //std::vector<double> SBB_val(basis_trip->nnz);
    //nsNNLS::sparseMatrix SBB_basis(basis_trip->nrow, basis_trip->ncol, basis_trip->nnz, SBB_ridx.data(), SBB_cptr.data(), SBB_val.data());
    //eigen_timer.start();
    //Eigen::SparseMatrix<double> eigen_basis(basis_trip->nrow, basis_trip->ncol);
    //eigen_basis.reserve(basis_trip->nnz);
    //eigen_timer.end();
    int accum = 0;
    for (int i = 0; i < nspes; i++) {
        ((long *)(basis->p))[i] = accum;
        //eigen_timer.start();
        //eigen_basis.outerIndexPtr()[i] = accum;
        //eigen_timer.end();
        //SBB_cptr[i] = accum;
        accum += col_counts[i];
        //eigen_basis.innerNonZeroPtr()[i] = accum; // the matrix is in the compressed format so we do not need this
    }
    // Need to set data end pointer for the last column.  Otherwise
    // SuiteSparse will ignore the last column.
    ((long *)(basis->p))[nspes] = accum;
    //SBB_cptr[nspes] = accum;
    //eigen_timer.start();
    //eigen_basis.outerIndexPtr()[nspes] = accum;
    //eigen_timer.end();
    std::vector<long> col_indices(nspes,0);


    for (unsigned i = 0; i < basis_trip->nnz; i++) {
        long col = ((long *)(basis_trip->j))[i];
        long index = ((long *)(basis->p))[col] + col_indices[col]++;
        ((long *)(basis->i))[index] = ((long *)(basis_trip->i))[i];
        ((double *)(basis->x))[index] = ((double *)(basis_trip->x))[i];
        //SBB_ridx[i] = ((long *)(basis_trip->i))[i];
        //SBB_val[i] = ((double *)(basis_trip->x))[i];
        //eigen_timer.start();
        //eigen_basis.innerIndexPtr()[index] = ((long *)(basis_trip->i))[i];
        //eigen_basis.valuePtr()[index] = ((double *)(basis_trip->x))[i];
        //eigen_timer.end();
    }
    cholmod_l_free_triplet(&basis_trip, &chol_common);
    //eigen_timer.start();
    //Eigen::Map<fnnls::VectorX_<double>> data_view(((double *)(data->x)), data->nrow);
    //eigen_timer.end();
    //std::vector<double> SBB_data_vector(((double *)data->x), ((double *)data->x) + data->nrow);
    //nsNNLS::vector SBB_data(SBB_data_vector.size(), SBB_data_vector.data());
    pre_compute.end();
    //eigen_timer.print();


    // Solve for SPE heights
    if (reduce_) {
        {
            double s = 0.0;
            double u = 0.0;
            std::string name = "RNNLS";
            NNLSTimer nnls_timer(name, s, u, true);
            unfolded = rnnls(basis, data, tolerance_, 2000, 0, &chol_common);
        }
        /*{
            double s = 0.0;
            double u = 0.0;
            std::string name = "FNNLS";
            NNLSTimer nnls_timer(name, s, u, true);
            fnnls::fnnls_solver(eigen_basis, data_view, 2000, tolerance_);
        }
        nsNNLS::vector* result;
        {
            double s = 0.0;
            double u = 0.0;
            std::string name = "SBB_NNLS";
            NNLSTimer nnls_timer(name, s, u, true);
            
            nsNNLS::nnls solver(&SBB_basis, &SBB_data, &SBB_seed, 2000);
            solver.optimize();
            result = solver.getSolution();
        }
        delete result; */
    } else {
        unfolded = nnls_lawson_hanson(basis, data, tolerance_,
            0, 2000, nspes, 0, 1, 0, &chol_common);
    }
    double PostNNLS_s = 0.0;
    double PostNNLS_u = 0.0;
    NNLSTimer PostNNLS_timer("PostNNLS", PostNNLS_s, PostNNLS_u, true);

    cholmod_l_free_sparse(&basis, &chol_common);
    cholmod_l_free_dense(&data, &chol_common);

    // Load the SPE corrections
    double speCorrection = 1.;
    if (apply_spe_corr_) {
        //speCorrection = 1. / calibration.GetMeanPMTCharge();
        speCorrection = 1.;
    }

    // let's account for the electron transit time!!
    double electron_time = calibration.GetPMTDeltaT();

    // Convert to pulse series
    for (int i = 0; i < nspes; i++) {
        if (((double *)(unfolded->x))[i] == 0)
            continue;

        CCMRecoPulse pulse;

        pulse.SetTime(start_times[i].first - electron_time);
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
        std::vector<std::reference_wrapper<std::vector<double>>> & output_data_times_references,
        std::vector<std::reference_wrapper<std::vector<double>>> & output_rebin_data_times_references,
        double template_bin_spacing_,
        double noise_threshold_,
        double basis_threshold_,
        double spes_per_bin_,
        double tolerance_,
        bool apply_spe_corr_,
        I3FramePtr frame
        ) {
    double RunPulsesThread_s = 0.0;
    double RunPulsesThread_u = 0.0;
    NNLSTimer RunPulsesThread_timer("RunPulsesThread", RunPulsesThread_s, RunPulsesThread_u, true);
    for(size_t i=std::get<0>(thread_range); i<std::get<1>(thread_range); ++i) {
        CCMPMTKey pmt_key = pmt_keys.at(i);
        size_t channel = pmt_channel_map.at(pmt_key);
        CCMWaveformDouble const & waveform = waveforms.at(channel);

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(pmt_key);

        if(calib == calibration.pmtCal.cend()) {
            //std::cout << "oops! no calibration for " << pmt_key << std::endl;
            continue;
        }

        double template_bin_spacing_ = 2.0 / spes_per_bin_ ;

        double placeholder = 1.0;

        double GetPulses_s = 0.0;
        double GetPulses_u = 0.0;
        NNLSTimer GetPulses_timer("GetPulses", GetPulses_s, GetPulses_u, true);
        GetPulses_timer.start();
        GetPulses(
                waveform,
                templates.at(pmt_key),
                calib->second,
                placeholder,
                template_bin_spacing_,
                noise_threshold_,
                basis_threshold_,
                spes_per_bin_,
                reduce,
                tolerance_,
                apply_spe_corr_,
                cholmod_common_vec.at(i),
                output_pulses_references[i].get(),
                output_data_times_references[i].get(),
                output_rebin_data_times_references[i].get(),
                frame
                );
        GetPulses_timer.end();
    }
}

void FillTemplate(CCMWaveformTemplate& wfTemplate, const CCMPMTCalibration& calibration, double const & start_time, int const & template_bins, double const & template_bin_spacing_, I3FramePtr frame) {
    CCMSPETemplate channel_template = calibration.GetSPETemplate();
    wfTemplate.digitizer_template.resize(template_bins);
    double total = 0.0;
    double wf_max = DBL_MIN;
    for(int i = 0; i < template_bins; i++) {
        double value = channel_template.Evaluate(start_time + i*template_bin_spacing_); // hard coding bin spacing...use calib one day
        wfTemplate.digitizer_template[i] = value;
        total += value;
        wf_max = std::max(wf_max, value);
    }

    double tail_fraction = 0.005;
    double target_amount = tail_fraction * total;
    double height_cut = tail_fraction * wf_max;

    double running_total = 0.0;
    size_t begin_idx = 0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        running_total += wfTemplate.digitizer_template[i];
        //if(running_total >= target_amount) {
        if(wfTemplate.digitizer_template[i] >= height_cut) {
            begin_idx = i;
            break;
        }
    }
    if(begin_idx > 0)
        begin_idx -= 1;

    target_amount = tail_fraction * total;
    running_total = 0.0;
    size_t end_idx = wfTemplate.digitizer_template.size();
    for(int i=wfTemplate.digitizer_template.size()-1; i>=0; --i) {
        running_total += wfTemplate.digitizer_template[i];
        //if(running_total >= target_amount) {
        if(wfTemplate.digitizer_template[i] >= height_cut) {
            end_idx = i;
            break;
        }
    }
    if(end_idx < wfTemplate.digitizer_template.size() - 1)
        end_idx += 1;

    wfTemplate.digitizer_template.erase(wfTemplate.digitizer_template.begin(), wfTemplate.digitizer_template.begin() + begin_idx);
    wfTemplate.start_time = start_time + template_bin_spacing_ * begin_idx;

    end_idx -= begin_idx;

    wfTemplate.digitizer_template.erase(wfTemplate.digitizer_template.begin() + end_idx, wfTemplate.digitizer_template.end());
    wfTemplate.end_time = wfTemplate.start_time + wfTemplate.digitizer_template.size() * template_bin_spacing_;

    wfTemplate.pulse_width = wfTemplate.end_time - wfTemplate.start_time;

    total = 0.0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        total += wfTemplate.digitizer_template[i];
    }

    wfTemplate.total_mass = total;

    /*
    boost::shared_ptr<I3Vector<double>> wf_template_copy = boost::make_shared<I3Vector<double>>();
    for (int i = 0; i < wfTemplate.digitizer_template.size(); ++i){
        wf_template_copy->push_back(wfTemplate.digitizer_template.at(i));
    }

    frame->Put("WaveformTemplate", wf_template_copy);
    */

    CCMFillFWHM(wfTemplate.digitizerStart,
            wfTemplate.digitizerStop,
            wfTemplate.digitizer_template,
            template_bin_spacing_, wfTemplate.start_time);

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

        double rnnls_u = 0.0;
        double rnnls_s = 0.0;
        double fnnls_u = 0.0;
        double fnnls_s = 0.0;

        std::vector<std::tuple<size_t, size_t>> thread_ranges_;
        CCMCalibration calibration;

        virtual ~CCMWavedeform();

        void Geometry(I3FramePtr frame);
        void Configure();
        void Calibration(I3FramePtr frame);
        void DAQ(I3FramePtr frame);
    private:
        std::string waveforms_name_;
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
    calibration = frame->Get<CCMCalibration>("CCMCalibration");

    double start_time = 2 * -10;
    double end_time = 2 * 60;

    template_bin_spacing_ = 2.0 / spes_per_bin_;
    double range = end_time - start_time;

    for(size_t i=0; i<pmt_keys_.size(); ++i) {
        // let's get our wf and what not
        CCMPMTKey key = pmt_keys_[i];
        uint32_t channel = pmt_channel_map_[key];

        if(calibration.pmtCal.count(key) == 0) {
            // std::cout << "oops! no calibration for " << key << std::endl;
            continue;
        }

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);

        double placeholder = 1.0;

        if(not template_.at(key).filled) {
            int template_bins = (int)ceil(range / template_bin_spacing_);
            FillTemplate(template_.at(key), calib->second, start_time, template_bins, template_bin_spacing_, frame);
            std::cout << key << " Template(" << template_.at(key).start_time << ", " << template_.at(key).end_time << ") FWHM(" << template_.at(key).digitizerStart << ", " << template_.at(key).digitizerStop << ")" << std::endl;
        }
    }

    PushFrame(frame);
}

void CCMWavedeform::DAQ(I3FramePtr frame) {
    if (!frame->Has(waveforms_name_)) {
        PushFrame(frame);
        return;
    }
    double DAQ_s = 0.0;
    double DAQ_u = 0.0;
    NNLSTimer DAQ_timer("DAQ", DAQ_s, DAQ_u, true);

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
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> output_data_times(new I3Map<CCMPMTKey, std::vector<double>>());
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> output_rebin_data_times(new I3Map<CCMPMTKey, std::vector<double>>());
    for(size_t i = 0; i < pmt_keys_.size(); ++i) {
        output->operator[](pmt_keys_[i]) = CCMRecoPulseSeries();
        output_data_times->operator[](pmt_keys_[i]) = std::vector<double>();
        output_rebin_data_times->operator[](pmt_keys_[i]) = std::vector<double>();
    }

    std::vector<std::reference_wrapper<CCMRecoPulseSeries>> output_pulses_references;
    std::vector<std::reference_wrapper<std::vector<double>>> output_data_times_references;
    std::vector<std::reference_wrapper<std::vector<double>>> output_rebin_data_times_references;
    output_pulses_references.reserve(pmt_keys_.size());
    for(size_t i = 0; i < pmt_keys_.size(); ++i) {
        output_pulses_references.emplace_back(output->at(pmt_keys_[i]));
        output_data_times_references.emplace_back(output_data_times->at(pmt_keys_[i]));
        output_rebin_data_times_references.emplace_back(output_rebin_data_times->at(pmt_keys_[i]));
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
            std::ref(output_data_times_references),
            std::ref(output_rebin_data_times_references),
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

    frame->Put("OriginalDataBins", output_data_times);
    frame->Put("RebinnedDataBins", output_rebin_data_times);
    frame->Put(output_name_, output);

    PushFrame(frame, "OutBox");
}



