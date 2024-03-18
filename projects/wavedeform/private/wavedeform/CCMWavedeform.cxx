#include <cctype>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <thread>
#include <mutex>
#include <stack>
#include <numeric>

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

#include <dataclasses/I3Double.h>
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

void linear_fit(double const * begin, double const * end, double & linear_b, double & linear_m) {
    double N = std::distance(begin, end);
    double X = 0.0;
    double T = 0.0;
    double T2 = 0.0;
    double XT = 0.0;

    size_t t = 0;
    double const * it = begin;
    for(; it != end; ++it) {
        double x = *it;
        X += x;
        T += t;
        T2 += t*t;
        XT += t*x;

        ++t;
    }
    double b = X/N - (N*XT*T-X*T*T) / (N*N*T2 - N*T*T);
    double m = (N*XT - X*T) / (N*T2 - T*T);
    linear_b = b;
    linear_m = m;
}

void rebin_bayesian_blocks(std::vector<double> const & raw_bins,
    double const * raw_charges, std::vector<std::vector<size_t>> & bin_indices,
    double ncp_prior, double poisson_scale_factor) {
    int N_raw = int(raw_bins.size()) - 1;
    std::vector<double> best(N_raw, -std::numeric_limits<double>::infinity());
    std::vector<int> last(N_raw, 0);

    // sum of raw_charges[:k+1]
    double total = 0;
    for(int k = 0; k < N_raw; k++) {
        total += raw_charges[k] / poisson_scale_factor;
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
            subtotal += raw_charges[j] / poisson_scale_factor;
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

bool GetPulses(CCMWaveformDouble const & wf, size_t wf_begin, size_t wf_end, CCMWaveformTemplate const & wfTemplate, CCMPMTCalibration const & calibration, double spe_charge, double wf_bin_width, double noise_threshold, double basis_threshold, double spes_per_bin, bool reduce, double tolerance, cholmod_common & chol_common, CCMRecoPulseSeries & output, std::vector<double> & output_data_times, std::vector<double> & output_rebin_data_times, DurationTimer & timer, size_t & chunk_iterations, I3Frame * frame) {
    cholmod_triplet *basis_trip;
    cholmod_sparse *basis;
    cholmod_dense *data, *unfolded;
    unsigned j, k;
    int nbins = 0;

    size_t wf_size = wf_end - wf_begin;

    // Increase size by 1 to account for zero padding at the start of the waveform
    wf_size += 1;

    // Determine the total number of WF bins
    nbins = wf_size;
    // If we have no data, nothing to do
    if (nbins == 0 || !std::isfinite(spe_charge) || spe_charge == 0)
        return true;

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
    double template_bin_spacing = wf_bin_width / spes_per_bin;

    // If the waveform is shorter than a pulse width,
    // increase the per-bin weights so the aggregate weight
    // is closer to what it should be
    if (wf_size * wf_bin_width < wfTemplate.pulse_width){
        base_weight *= wfTemplate.pulse_width / (wf_size * wf_bin_width);
    }

    double noise = noise_threshold;
    double basisThreshmV = basis_threshold;

    std::vector<double> data_times;
    data_times.reserve(wf_size);

    redges[0] = 0;
    data_times.push_back(-2.0);
    weights[0] = base_weight / 1.0;
    passBasisThresh[0] = false;

    // Read waveform
    for (k = 1; k < wf_size; k++) {
        redges[k] = (1. + k - 1) * wf_bin_width;
        ((double *)(data->x))[k] = wf.GetWaveform()[k+wf_begin-1];
        data_times.push_back((k-1) * wf_bin_width);

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
        } else if (fabs(((double *)(data->x))[k]) > basis_threshold) {
            passBasisThresh[k] = true;
        } else {
            ((double *)(data->x))[k] = 0;
            weights[k] = 0;
        }
    }

    // Apply weights
    for (int i = 0; i < nbins; i++) {
        ((double *)(data->x))[i] *= weights[i];
    }

    // Precalculate the max number of basis functions to avoid reallocation
    int maxspes = int(spes_per_bin*(wf_size));

    std::vector<std::pair<double, double> > start_times;
    start_times.reserve(maxspes);
    double min_spe_spacing = DBL_MAX;
    double max_spe_spacing = DBL_MIN;

    // Start, end two bins early
    double present = - 2. * wf_bin_width;
    double max = present + wf_size * wf_bin_width;

    double fine_spacing = wf_bin_width / spes_per_bin;
    min_spe_spacing = std::min(fine_spacing, min_spe_spacing);
    max_spe_spacing = std::max(fine_spacing, max_spe_spacing);

    double coarse_spacing = wf_bin_width * 1.5;
    min_spe_spacing = std::min(coarse_spacing, min_spe_spacing);
    max_spe_spacing = std::max(coarse_spacing, max_spe_spacing);

    double ccoarse_spacing = wf_bin_width * 2.0;
    min_spe_spacing = std::min(ccoarse_spacing, min_spe_spacing);
    max_spe_spacing = std::max(ccoarse_spacing, max_spe_spacing);

    double spacing = coarse_spacing;

    // Get the peak FWHM start, stop times for this waveform
    double fwhmStart = wfTemplate.digitizerStart;
    double fwhmStop = wfTemplate.digitizerStop;
    // Generate the set of pulse start times that have reasonable
    // non-zero data within the FWHM of the corresponding pulse
    bool prev_below = true;
    size_t fine_spacing_tau = 20;
    size_t coarse_spacing_tau = 40;
    size_t n_since_fine = 0;
    for (k = 0; k < wf_size; ++k) {
        if (((double *)(data->x))[k] != 0. && passBasisThresh[k]) {
            if(prev_below) {
                n_since_fine = 0;
                spacing = fine_spacing;
            } else if(n_since_fine > coarse_spacing_tau) {
                spacing = ccoarse_spacing;
            } else if(n_since_fine > fine_spacing_tau) {
                spacing = coarse_spacing;
            }
            double binTime = redges[k];
            // Don't jump if we're moving less than the basis spacing
            if (present < (binTime - fwhmStop)) {
                present = binTime - fwhmStop;
            }

            // Add start times to the set
            while (present < (binTime - fwhmStart) && present < max && present >= (binTime - fwhmStop)) {
                start_times.push_back(std::pair<double, double>(present, spacing));
                present += spacing;
            }
            n_since_fine += 1;
            prev_below = false;
        } else {
            prev_below = true;
        }
    }

    int nspes = start_times.size();
    if (nspes == 0) {
        // We have no nonzero data left
        return true;
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
        start_times[0].second = start_times[1].first - start_times[0].first;
        start_times[nspes-1].second = start_times[nspes-1].first - start_times[nspes-2].first;
        // Compute rightdiff
        for (i = 1; i < nspes-1; ++i) {
            start_times[i].second =
                start_times[i+1].first - start_times[i].first;
        }
        // Compute leftdiff and average
        for (i = 1; i < nspes-1; ++i) {
            double leftdiff =
                start_times[i].first - start_times[i-1].first;
            start_times[i].second = (leftdiff + start_times[i].second) / 2.0;
        }
    }

    //Don't use data bins that are not in the support of any basis function
    double start = DBL_MIN;
    double end = DBL_MIN;

    k = 0;
    for (std::vector<std::pair<double, double> >::const_iterator it = start_times.begin();
            it != start_times.end(); ++it) {
        start = it->first + wfTemplate.start_time;
        end = it->first + wfTemplate.end_time;

        // Evaluate bins up until we pass the end of the current time range
        for (; k < wf_size && redges[k] < end; ++k) {
            // Check if bin time is inside the max of
            // the last range and the min of the current range
            if (redges[k] < start) {
                weights[k] = 0.;
            }
        }

        if (k == wf_size) {
            break;
        }
    }

    while (k < wf_size) {
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
            ((double *)(data->x))[j] = ((double *)(data->x))[i];
            ++j;
        }
    }
    nbins = j;
    data->nrow = nbins;
    data_times.resize(nbins);
    data_times.push_back(data_times.back() + wf_bin_width);
    output_data_times = data_times;


    // Attempt to rebin the data with bayesian blocks
    std::vector<std::vector<size_t>> bb_rebin_ranges;
    double ncp_prior = 100.0;
    rebin_bayesian_blocks(data_times, ((double *)(data->x)), bb_rebin_ranges, ncp_prior, wfTemplate.total_mass);
    output_rebin_data_times.clear();

    // Modify the rebinning so that only data below threshold is rebinned
    std::vector<std::vector<size_t>> rebin_ranges;
    rebin_ranges.reserve(bb_rebin_ranges.size());
    std::vector<size_t> current_binning;
    std::function<void(size_t)> new_bin = [&current_binning, &rebin_ranges] (size_t idx) {
        if(current_binning.size() == 0) {
            current_binning.push_back(idx);
            rebin_ranges.push_back(current_binning);
            current_binning.clear();
        } else {
            rebin_ranges.push_back(current_binning);
            current_binning.clear();
            current_binning.push_back(idx);
            rebin_ranges.push_back(current_binning);
            current_binning.clear();
        }
    };

    bool zero_bin = true;
    for(size_t i=0; i<bb_rebin_ranges.size(); ++i) {
        for(size_t j=0; j<bb_rebin_ranges[i].size(); ++j) {
            size_t idx = bb_rebin_ranges[i][j];
            bool current_zero = ((double *)(data->x))[idx] <= 0;
            if((zero_bin != current_zero) or (not current_zero)) {
                new_bin(idx);
            } else {
                current_binning.push_back(idx);
            }
            zero_bin = current_zero;
            if(not zero_bin and current_binning.size() > 2) {
                rebin_ranges.push_back(current_binning);
                current_binning.clear();
            }
        }
        if(current_binning.size() > 0 and not zero_bin) {
            rebin_ranges.push_back(current_binning);
            current_binning.clear();
            zero_bin = true;
        }
    }
    if(current_binning.size() > 0) {
        rebin_ranges.push_back(current_binning);
        current_binning.clear();
    }

    for(size_t i=0; i<rebin_ranges.size(); ++i) {
        output_rebin_data_times.push_back(data_times[rebin_ranges[i][0]]);
    }
    output_rebin_data_times.push_back(data_times.back());

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
        while(first_spe < nspes-1 && redges[i] > start_times[first_spe].first + wfTemplate.end_time + max_spe_spacing)
            first_spe++;
        for(int j = first_spe; j < nspes; j++) {
            if(redges[i] < start_times[j].first + wfTemplate.start_time - max_spe_spacing)
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
        while (first_spe < nspes && start_times[first_spe].first + wfTemplate.end_time <= redges[i]) {
            first_spe++;
        }
        if (first_spe == nspes) {
            continue;
        }

        // Precompute the multiplicative term in the basis
        //double weighted_charge = (spe_charge/I3Units::mV) * weights[i];
        double weighted_charge = spe_charge * weights[i];
        if (weighted_charge == 0) // Removed from fit, so ignore
            continue;

        // Precache which pulse template we're using
        double templ_bin_spacing_inv = 1./template_bin_spacing;
        std::vector<double> const & pulse_templ = wfTemplate.digitizer_template;
        int n_template_bins = pulse_templ.size();

        // The last pulse for this bin is 2 ns in the future
        for (int j = first_spe; j < nspes; j++) {
            int templ_bin = int(((redges[i] - start_times[j].first) - wfTemplate.start_time)*templ_bin_spacing_inv);
            if (templ_bin < 0)
                break;
            if (templ_bin >= n_template_bins)
                continue;

            ((long *)(basis_trip->i))[basis_trip->nnz] = i;
            ((long *)(basis_trip->j))[basis_trip->nnz] = j;
            ((double *)(basis_trip->x))[basis_trip->nnz] = pulse_templ.at(templ_bin) * weighted_charge;

            pflags[j] |= source;

            basis_trip->nnz++;
            assert(basis_trip->nnz <= basis_trip->nzmax);
            col_counts[j]++;
        }
    }

    k = 0;
    double data_running_total = 0.0;
    std::vector<std::vector<size_t>> spe_original_support(nspes);
    std::vector<std::vector<double>> spe_original_entries(nspes);
    std::vector<std::vector<size_t>> spe_bb_support(nspes);
    std::vector<std::vector<double>> spe_bb_entries(nspes);
    std::vector<double> rebinned_data(rebin_ranges.size(), 0.0);

    // Merge entries in the matrix according to the data rebinning
    for(size_t bb_idx=0; bb_idx<rebin_ranges.size(); ++bb_idx) {
        for(size_t data_idx : rebin_ranges[bb_idx]) {
            while(((long *)(basis_trip->i))[k] < (long int)(data_idx)) {
                ++k;
            }
            assert(((long *)(basis_trip->i))[k] == (long int)(data_idx));
            rebinned_data[bb_idx] += ((double *)(data->x))[data_idx];
            while(k < basis_trip->nnz and ((long *)(basis_trip->i))[k] == (long int)(data_idx)) {
                size_t spe_idx = ((long *)(basis_trip->j))[k];
                spe_original_support[spe_idx].push_back(data_idx);
                double x = ((double *)(basis_trip->x))[k];
                spe_original_entries[spe_idx].push_back(x);
                if(spe_bb_support[spe_idx].size() == 0 or spe_bb_support[spe_idx].back() != bb_idx) {
                    spe_bb_support[spe_idx].push_back(bb_idx);
                    spe_bb_entries[spe_idx].push_back(x);
                } else {
                    spe_bb_entries[spe_idx].back() += x;
                }
                ++k;
            }
        }
    }

    size_t reduced_nnz = 0;
    // Parent SPE, child SPE
    std::vector<std::vector<double>> merged_spe_magnitudes;
    std::vector<std::vector<double>> merged_spe_start_times;
    std::vector<std::vector<double>> merged_spe_widths;

    // Parent SPE, data index
    std::vector<std::vector<size_t>> merged_spe_support;
    std::vector<std::vector<double>> merged_spe_entries;
    // Store representative SPE index, and a vector of combined SPE times and ratios

    // Determine if any SPEs have the same support for the rebinned data and merge them accordingly
    // Only merge SPEs that are supported by a data bin that has been rebinned
    // We must choose a relative magnitude for the merged SPEs when combining them and when uncombining after fitting,
    //   we can use a linear fit of the supporting data to determine this relative magnitude
    size_t min_spe_idx = 0;
    size_t max_spe_idx = 0;
    size_t pad_left = 0;
    size_t pad_right = 0;
    for(size_t rebin_data_idx=0; rebin_data_idx<rebin_ranges.size() and min_spe_idx<spe_bb_support.size(); ++rebin_data_idx) {
        max_spe_idx += 1;

        while(rebin_data_idx < rebin_ranges.size()-1
                and max_spe_idx < spe_bb_support.size()
                and std::abs(start_times[max_spe_idx].first - data_times[rebin_ranges[rebin_data_idx].back()])
                    < std::abs(start_times[max_spe_idx].first - data_times[rebin_ranges[rebin_data_idx].back()+1])) {
            max_spe_idx += 1;
        }
        if(rebin_ranges[rebin_data_idx].size() == 1) {
            // Put the unmerged SPEs into the final basis
            for(size_t jj=0; jj<(max_spe_idx-min_spe_idx); ++jj) {
                size_t j = min_spe_idx + jj;
                merged_spe_magnitudes.push_back({1.0});
                merged_spe_start_times.push_back({start_times[j].first});
                merged_spe_widths.push_back({start_times[j].second});

                merged_spe_support.push_back(spe_bb_support[j]);
                merged_spe_entries.push_back(spe_bb_entries[j]);
                reduced_nnz += merged_spe_entries.back().size();
                for(size_t ii=0; ii<spe_bb_entries[j].size(); ++ii) {
                    if(std::isnan(spe_bb_entries[j][ii])) {
                        std::cout << "NaN spe entry 2" << std::endl;
                    }
                }
            }
            min_spe_idx = max_spe_idx;
        } else {
            bool skip_left = (rebin_data_idx > 0 and rebin_ranges[rebin_data_idx-1].size() == 1);
            bool skip_right = (rebin_data_idx < rebin_ranges.size() - 1 and rebin_ranges[rebin_data_idx+1].size() == 1);
            size_t merge_min_idx = std::min(min_spe_idx + pad_left, spe_bb_support.size()-1);
            size_t merge_max_idx = std::max(max_spe_idx - std::min(pad_right, max_spe_idx), merge_min_idx);
            size_t n_merge = merge_max_idx - merge_min_idx;
            if(skip_left) {
                for(size_t jj=0; jj<pad_left; ++jj) {
                    size_t j = merge_min_idx + jj;
                    merged_spe_magnitudes.push_back({1.0});
                    merged_spe_start_times.push_back({start_times[j].first});
                    merged_spe_widths.push_back({start_times[j].second});

                    merged_spe_support.push_back(spe_bb_support[j]);
                    merged_spe_entries.push_back(spe_bb_entries[j]);
                    reduced_nnz += merged_spe_entries.back().size();
                    for(size_t ii=0; ii<spe_bb_entries[j].size(); ++ii) {
                        if(std::isnan(spe_bb_entries[j][ii])) {
                            std::cout << "NaN spe entry 2" << std::endl;
                        }
                    }
                }
            }
            if(skip_right) {
                for(size_t jj=0; jj<pad_right; ++jj) {
                    size_t j = merge_max_idx + jj;
                    merged_spe_magnitudes.push_back({1.0});
                    merged_spe_start_times.push_back({start_times[j].first});
                    merged_spe_widths.push_back({start_times[j].second});

                    merged_spe_support.push_back(spe_bb_support[j]);
                    merged_spe_entries.push_back(spe_bb_entries[j]);
                    reduced_nnz += merged_spe_entries.back().size();
                    for(size_t ii=0; ii<spe_bb_entries[j].size(); ++ii) {
                        if(std::isnan(spe_bb_entries[j][ii])) {
                            std::cout << "NaN spe entry 2" << std::endl;
                        }
                    }
                }
            }
            // Use the dot product with the data to set the relative magnitudes
            double total_dot_product = 0.0;
            std::vector<double> spe_magnitudes;
            spe_magnitudes.reserve(n_merge);
            // Iterate over SPEs in the center
            size_t min_data_idx = spe_original_support[merge_min_idx].front();
            size_t max_data_idx = spe_original_support[merge_max_idx-1].back();
            while(redges[min_data_idx] < start_times[merge_min_idx].first + fwhmStart - wf_bin_width and min_data_idx < max_data_idx) {
                min_data_idx += 1;
            }
            while(redges[max_data_idx] > start_times[merge_max_idx-1].first + fwhmStop - wf_bin_width and max_data_idx > min_data_idx) {
                max_data_idx -= 1;
            }
            for(size_t i=0; i<2 and max_data_idx == min_data_idx; ++i) {
                if(min_data_idx > 0)
                    min_data_idx -= 1;
                if(max_data_idx < redges.size() - 1)
                    max_data_idx += 1;
            }
            double linear_b = 0;
            double linear_m = 0;
            if(max_data_idx > min_data_idx) {
                linear_fit((double *)(data->x) + min_data_idx, (double *)(data->x) + max_data_idx + 1, linear_b, linear_m);
                linear_m /= wf_bin_width; // convert slope to ns
            }

            for(size_t jj=0; jj<n_merge; ++jj) {
                size_t j = merge_min_idx + jj;
                double dot_product = std::max(0.0, linear_b + linear_m * (start_times[j].first - redges[min_data_idx]));
                dot_product /= start_times[j].second;
                spe_magnitudes.emplace_back(dot_product);
                total_dot_product += dot_product;
            }
            std::vector<double> spe_start_times;
            spe_start_times.reserve(n_merge);
            std::vector<double> spe_widths;
            spe_widths.reserve(n_merge);
            size_t min_support_idx = spe_bb_support[merge_min_idx].front();
            size_t max_support_idx = spe_bb_support[merge_max_idx-1].back();
            std::vector<double> spe_entries(max_support_idx - min_support_idx + 1, 0.0);
            std::vector<size_t> spe_support(max_support_idx - min_support_idx + 1);
            std::iota(spe_support.begin(), spe_support.end(), min_support_idx);

            for(size_t jj=0; jj<n_merge; ++jj) {
                if(total_dot_product > 0) {
                    spe_magnitudes[jj] /= total_dot_product;
                } else {
                    spe_magnitudes[jj] = 1.0 / (n_merge);
                }
                size_t j = merge_min_idx + jj;
                spe_start_times.emplace_back(start_times[j].first);
                spe_widths.emplace_back(start_times[j].second);
                for(size_t ii = 0; ii<spe_bb_support[j].size(); ++ii) {
                    spe_entries[ii + spe_bb_support[j].front() - min_support_idx] += spe_bb_entries[j][ii] * spe_magnitudes[jj];
                    if(std::isnan(spe_entries[ii + spe_bb_support[j].front() - min_support_idx])) {
                        std::cout << "NaN spe entry 1" << std::endl;
                    }
                }
            }

            merged_spe_magnitudes.push_back(spe_magnitudes);
            merged_spe_start_times.push_back(spe_start_times);
            merged_spe_widths.push_back(spe_widths);

            merged_spe_support.push_back(spe_support);
            merged_spe_entries.push_back(spe_entries);
            reduced_nnz += merged_spe_entries.back().size();
        }
        min_spe_idx = max_spe_idx;
    }

    basis = cholmod_l_allocate_sparse(rebinned_data.size(), merged_spe_entries.size(),
            reduced_nnz, true, true, 0, CHOLMOD_REAL, &chol_common);

    size_t nnz_idx = 0;
    for(size_t spe_idx=0; spe_idx<merged_spe_entries.size(); ++spe_idx) {
        ((long *)(basis->p))[spe_idx] = nnz_idx;
        for(size_t jj=0; jj<merged_spe_support[spe_idx].size(); ++jj) {
            ((long *)(basis->i))[nnz_idx] = merged_spe_support[spe_idx][jj];
            ((double *)(basis->x))[nnz_idx] = merged_spe_entries[spe_idx][jj];
            ++nnz_idx;
        }
    }

    for(size_t i=0; i<rebinned_data.size(); ++i) {
        ((double *)(data->x))[i] = rebinned_data[i];
    }
    data->nrow = rebin_ranges.size();

    cholmod_l_free_triplet(&basis_trip, &chol_common);

    // Solve for SPE heights
    if (reduce) {
        unfolded = rnnls(basis, data, tolerance, 1000, 0, &chol_common, timer, chunk_iterations);
    } else {
        unfolded = nnls_lawson_hanson(basis, data, tolerance,
            0, 1000, nspes, 0, 1, 0, &chol_common, timer, chunk_iterations);
    }

    cholmod_l_free_sparse(&basis, &chol_common);
    cholmod_l_free_dense(&data, &chol_common);

    if(unfolded == nullptr) {
        return false;
    }

    // Convert to pulse series
    for(size_t i = 0; i < merged_spe_start_times.size(); i++) {
        if (((double *)(unfolded->x))[i] == 0)
            continue;

        for(size_t j=0; j<merged_spe_start_times[i].size(); ++j) {
            CCMRecoPulse pulse;
            pulse.SetTime(merged_spe_start_times[i][j]);
            pulse.SetCharge((((double *)(unfolded->x))[i]) * merged_spe_magnitudes[i][j]);
            pulse.SetWidth(merged_spe_widths[i][j]);
            output.push_back(pulse);
        }
    }
    cholmod_l_free_dense(&unfolded, &chol_common);
    return true;
}

bool GetPulses(CCMWaveformDouble const & wf, CCMWaveformTemplate const & wfTemplate, CCMPMTCalibration const & calibration, double spe_charge, double wf_bin_width, double noise_threshold, double basis_threshold, double spes_per_bin, bool reduce, double tolerance, cholmod_common & chol_common, CCMRecoPulseSeries & output, std::vector<double> & output_data_times, std::vector<double> & output_rebin_data_times, DurationTimer & timer, I3Vector<unsigned int> & channel_iterations, I3Frame * frame) {
    output.clear();
    std::vector<double> const & w = wf.GetWaveform();
    if(w.size() == 0)
        return true;

    std::vector<std::pair<size_t, size_t>> regions;
    size_t start = 0;
    size_t end = 0;
    bool prev_in_region = w[0] > basis_threshold;
    for(size_t i=0; i<w.size(); ++i) {
        bool in_region = w[i] > basis_threshold;
        if(in_region) {
            if(prev_in_region) {
                //
            } else {
                start = i;
            }
        } else {
            if(prev_in_region) {
                end = i;
                if(regions.size() > 0 and end - regions.back().second < wfTemplate.digitizer_template.size()*2) {
                    regions.back() = std::pair<size_t, size_t>(regions.back().first, end);
                } else {
                    regions.emplace_back(start, end);
                }
            } else {
                //
            }
        }
        prev_in_region = in_region;
    }
    if(prev_in_region) {
        end = w.size();
        regions.emplace_back(start, end);
        prev_in_region = false;
    }

    size_t front_ext = wfTemplate.digitizer_template.size();
    size_t back_ext = wfTemplate.digitizer_template.size();

    end = 0;
    for(size_t i=0; i<regions.size(); ++i) {
        if(i+1 < regions.size())
            start = regions[i+1].first;
        else
            start = w.size();
        size_t start_idx = regions[i].first;
        size_t ext = std::min(start_idx, back_ext);
        start_idx -= ext;
        start_idx = std::max(end, start_idx);
        size_t end_idx = regions[i].second;
        ext = std::min(w.size() - end_idx, front_ext);
        end_idx += ext;
        end_idx = std::min(start, end_idx);
        end = end_idx;
        double start_time = start_idx * 2.0;
        double end_time = end_idx * 2.0;
        CCMRecoPulseSeries region_output;
        std::vector<double> region_output_data_times;
        std::vector<double> region_output_rebin_data_times;
        size_t chunk_iterations = 0;
        bool success = GetPulses(wf, start_idx, end_idx, wfTemplate, calibration, spe_charge, wf_bin_width, noise_threshold, basis_threshold, spes_per_bin, reduce, tolerance, chol_common, region_output, region_output_data_times, region_output_rebin_data_times, timer, chunk_iterations, frame);
        channel_iterations.push_back(chunk_iterations);
        if(not success) {
            return false;
        }
        for(size_t j=0; j<region_output.size(); ++j) {
            CCMRecoPulse const & pulse = region_output[j];
            CCMRecoPulse new_pulse;
            new_pulse.SetTime(pulse.GetTime() + start_time);
            new_pulse.SetCharge(pulse.GetCharge());
            new_pulse.SetWidth(pulse.GetWidth());
            output.push_back(new_pulse);
        }
        for(size_t j=0; j<region_output_data_times.size(); ++j)
            output_data_times.push_back(region_output_data_times[j] + start_time);
        for(size_t j=0; j<region_output_rebin_data_times.size(); ++j)
            output_rebin_data_times.push_back(region_output_rebin_data_times[j] + start_time);
    }
    return true;
}

void RunPulsesThread(
        std::atomic<bool> & running,
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
        double wf_bin_width,
        double noise_threshold,
        double basis_threshold,
        double spes_per_bin,
        double tolerance,
        std::atomic<bool> & success,
        double time_limit_seconds,
        double & elapsed_time,
        I3Vector<I3Vector<unsigned int>> & iterations,
        I3Frame * frame
        ) {
    running.store(true);
    DurationTimer timer(success, elapsed_time, time_limit_seconds);
    for(size_t i=std::get<0>(thread_range); (i<std::get<1>(thread_range) and success.load() == true); ++i) {
        CCMPMTKey pmt_key = pmt_keys.at(i);
        size_t channel = pmt_channel_map.at(pmt_key);
        CCMWaveformDouble const & waveform = waveforms.at(channel);

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(pmt_key);

        if(calib == calibration.pmtCal.cend()) {
            continue;
        }

        double placeholder = 1.0;

        I3Vector<unsigned int> channel_iterations;
        bool channel_success = GetPulses(
                waveform,
                templates.at(pmt_key),
                calib->second,
                placeholder,
                wf_bin_width,
                noise_threshold,
                basis_threshold,
                spes_per_bin,
                reduce,
                tolerance,
                cholmod_common_vec.at(i),
                output_pulses_references[i].get(),
                output_data_times_references[i].get(),
                output_rebin_data_times_references[i].get(),
                timer,
                channel_iterations,
                frame
                );
        iterations.at(i) = channel_iterations;
        if(not channel_success) {
            success.store(false);
            break;
        }
    }
    running.store(false);
}

void FillTemplate(CCMWaveformTemplate& wfTemplate, const CCMPMTCalibration& calibration, double const & start_time, int const & template_bins, double const & template_bin_spacing, I3FramePtr frame) {
    CCMSPETemplate channel_template = calibration.GetSPETemplate();
    wfTemplate.digitizer_template.resize(template_bins);
    double total = 0.0;
    double wf_max = DBL_MIN;
    for(int i = 0; i < template_bins; i++) {
        double value = channel_template.Evaluate(start_time + i*template_bin_spacing); // hard coding bin spacing...use calib one day
        wfTemplate.digitizer_template[i] = value;
        total += value;
        wf_max = std::max(wf_max, value);
    }

    double tail_fraction = 0.005;
    double height_cut = tail_fraction * wf_max;

    double running_total = 0.0;
    size_t begin_idx = 0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        running_total += wfTemplate.digitizer_template[i];
        if(wfTemplate.digitizer_template[i] >= height_cut) {
            begin_idx = i;
            break;
        }
    }
    if(begin_idx > 0)
        begin_idx -= 1;

    running_total = 0.0;
    size_t end_idx = wfTemplate.digitizer_template.size();
    for(int i=wfTemplate.digitizer_template.size()-1; i>=0; --i) {
        running_total += wfTemplate.digitizer_template[i];
        if(wfTemplate.digitizer_template[i] >= height_cut) {
            end_idx = i;
            break;
        }
    }
    if(end_idx < wfTemplate.digitizer_template.size() - 1)
        end_idx += 1;

    wfTemplate.digitizer_template.erase(wfTemplate.digitizer_template.begin(), wfTemplate.digitizer_template.begin() + begin_idx);
    wfTemplate.start_time = start_time + template_bin_spacing * begin_idx;

    end_idx -= begin_idx;

    wfTemplate.digitizer_template.erase(wfTemplate.digitizer_template.begin() + end_idx, wfTemplate.digitizer_template.end());
    wfTemplate.end_time = wfTemplate.start_time + wfTemplate.digitizer_template.size() * template_bin_spacing;

    wfTemplate.pulse_width = wfTemplate.end_time - wfTemplate.start_time;

    total = 0.0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        total += wfTemplate.digitizer_template[i];
    }

    wfTemplate.total_mass = total;

    CCMFillFWHM(wfTemplate.digitizerStart,
            wfTemplate.digitizerStop,
            wfTemplate.digitizer_template,
            template_bin_spacing, wfTemplate.start_time);

    wfTemplate.filled = true;
}

struct WavedeformJob {
    std::atomic<bool> running = false;
    std::thread thread;
    size_t thread_index = 0;
    I3FramePtr frame = nullptr;
    size_t frame_index = 0;
};

struct WavedeformResult {
    I3FramePtr frame = nullptr;
    bool done = false;
    size_t frame_index = 0;
};

struct FrameWorkspace {
    bool can_run_jobs = false;
    size_t total_jobs = 0;
    size_t jobs_queued = 0;
    size_t jobs_done = 0;

    std::atomic<bool> success = true;
    std::vector<double> elapsed_times;

    I3FramePtr frame;

    boost::shared_ptr<CCMWaveformDoubleSeries const> waveforms;
    boost::shared_ptr<I3Double const> total_adc_count;

    boost::shared_ptr<I3Vector<I3Vector<unsigned int>>> iterations;

    // place to store pulses
    boost::shared_ptr<CCMRecoPulseSeriesMap> output;
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> output_data_times;
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> output_rebin_data_times;

    std::vector<std::reference_wrapper<CCMRecoPulseSeries>> output_pulses_references;
    std::vector<std::reference_wrapper<std::vector<double>>> output_data_times_references;
    std::vector<std::reference_wrapper<std::vector<double>>> output_rebin_data_times_references;

    std::vector<std::tuple<size_t, size_t>> job_ranges;

    FrameWorkspace(I3FramePtr frame, std::vector<CCMPMTKey> const & pmt_keys_, std::string const & waveforms_name_) :
        frame(frame) {
        if (!frame->Has(waveforms_name_)) {
            can_run_jobs = false;
            return;
        }

        waveforms = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>(waveforms_name_);
        total_adc_count = frame->Get<boost::shared_ptr<const I3Double>>(waveforms_name_ + "TotalADC");

        // place to store pulses
        output = boost::make_shared<CCMRecoPulseSeriesMap>();
        output_data_times = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
        output_rebin_data_times = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();

        for(size_t i = 0; i < pmt_keys_.size(); ++i) {
            output->operator[](pmt_keys_[i]) = CCMRecoPulseSeries();
            output_data_times->operator[](pmt_keys_[i]) = std::vector<double>();
            output_rebin_data_times->operator[](pmt_keys_[i]) = std::vector<double>();
        }

        output_pulses_references.reserve(pmt_keys_.size());
        output_data_times_references.reserve(pmt_keys_.size());
        output_rebin_data_times_references.reserve(pmt_keys_.size());

        for(size_t i = 0; i < pmt_keys_.size(); ++i) {
            output_pulses_references.emplace_back(output->at(pmt_keys_[i]));
            output_data_times_references.emplace_back(output_data_times->at(pmt_keys_[i]));
            output_rebin_data_times_references.emplace_back(output_rebin_data_times->at(pmt_keys_[i]));
        }

        total_jobs = 1;
        if(total_adc_count) {
            total_jobs = (unsigned int)(total_adc_count->value / 5e5);
            total_jobs = std::max(size_t(1), total_jobs);
            total_jobs = std::min(total_jobs, pmt_keys_.size());
        }
        if(total_jobs == 1) {
            job_ranges.push_back(std::tuple<size_t, size_t>(0, pmt_keys_.size()));
        } else {
            size_t max_size = pmt_keys_.size() / total_jobs; 
            for(size_t i=0; i<total_jobs-1; ++i) {
                job_ranges.push_back(std::tuple<size_t, size_t>(max_size*i, max_size*(i+1)));
            }
            job_ranges.push_back(std::tuple<size_t, size_t>(max_size*(total_jobs-1), pmt_keys_.size()));
        }

        iterations = boost::make_shared<I3Vector<I3Vector<unsigned int>>>(pmt_keys_.size());

        elapsed_times.resize(total_jobs, 0.0);

        success.store(true);
    }

    void StartJob(WavedeformJob * job,
        std::string const & waveforms_name,
        std::string const & output_name,
        std::vector<CCMPMTKey> const & pmt_keys,
        I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map,
        I3Map<CCMPMTKey, CCMWaveformTemplate> const & templates,
        CCMCalibration const & calibration,
        std::vector<std::vector<cholmod_common>> & cholmod_common_vec,
        bool reduce,
        double wf_bin_width,
        double noise_threshold,
        double basis_threshold,
        double spes_per_bin,
        double tolerance,
        size_t num_threads,
        bool remove_waveforms,
        double time_limit_seconds,
        std::exception_ptr & teptr) {

        job->running.store(true);

        size_t i = jobs_queued;
        job->thread = std::thread(RunPulsesThread,
            std::ref(job->running),
            job_ranges[i],
            std::cref(pmt_keys),
            std::cref(pmt_channel_map),
            std::cref(*waveforms),
            std::cref(templates),
            std::cref(calibration),
            std::ref(cholmod_common_vec[job->thread_index]),
            reduce,
            std::ref(output_pulses_references),
            std::ref(output_data_times_references),
            std::ref(output_rebin_data_times_references),
            wf_bin_width,
            noise_threshold,
            basis_threshold,
            spes_per_bin,
            tolerance,
            std::ref(success),
            time_limit_seconds,
            std::ref(elapsed_times.at(i)),
            std::ref(*iterations),
            job->frame.get()
        );
        jobs_queued += 1;
    }

    void Finish(bool remove_waveforms, std::string const & waveforms_name_, std::string const & output_name_) {
        double elapsed_time = std::accumulate(elapsed_times.begin(), elapsed_times.end(), 0.0);
        frame->Put(output_name_ + "RuntimeInSeconds", boost::make_shared<I3Double>(elapsed_time));
        frame->Put(output_name_ + "Iterations", iterations);
        if(not success.load())
            return;
        if(remove_waveforms) {
            frame->Delete(waveforms_name_);
        }
        frame->Put(output_name_ + "OriginalDataBins", output_data_times);
        frame->Put(output_name_ + "RebinnedDataBins", output_rebin_data_times);
        frame->Put(output_name_, output);
    }
};

class CCMWavedeform : public I3ConditionalModule {
    std::exception_ptr teptr;
public:
    bool remove_waveforms;
    double time_limit_seconds;
    size_t num_threads;
    size_t max_cached_frames;
    std::string geometry_name_;
    std::string pmt_channel_map_name_;
    bool geo_seen;
    std::string nim_pulses_name_;
    I3Vector<CCMPMTKey> allowed_pmt_keys_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    CCMWavedeform(const I3Context &);

    std::vector<std::tuple<size_t, size_t>> thread_ranges_;
    CCMCalibration calibration;

    virtual ~CCMWavedeform();

    void Configure();
    void Process();
    void Finish();
    void Geometry(I3FramePtr frame);
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
    double wf_bin_width_;

    bool reduce_;

    I3Map<CCMPMTKey, CCMWaveformTemplate> templates_;
    std::vector<CCMPMTKey> pmt_keys_;
    std::vector<std::vector<cholmod_common>> cholmod_common_vec_;

    size_t frame_index = 0;
    size_t min_frame_idx = 0;
    std::deque<WavedeformJob *> free_jobs;
    std::deque<WavedeformJob *> running_jobs;
    std::deque<WavedeformResult> results;
    std::map<size_t, FrameWorkspace *> frame_workspaces;
};

I3_MODULE(CCMWavedeform);

CCMWavedeform::CCMWavedeform(const I3Context& context) : I3ConditionalModule(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("NumThreads", "Number of worker threads to use for pulse fitting", (size_t)(0));
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTChannelMapName", "Key for PMTChannelMap", std::string(""));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
    AddParameter("SPEsPerBin", "Number of basis functions to unfold per waveform bin", 1.);
    AddParameter("Tolerance", "Stopping tolerance, in units of bin ADC^2/PE", 49.);
    AddParameter("NoiseThreshold","Consider bins with amplitude below this number of counts as noise", 5.0);
    AddParameter("BasisThreshold",
            "Require a bin with amplitude at least this number of counts "
            "within the FWHM of the template waveform in order to include "
            "a given start time in the basis set", 7.0);
    AddParameter("Waveforms", "Name of input waveforms",
            "CCMCalibratedWaveforms");
    AddParameter("Output", "Name of output pulse series",
            "WavedeformPulses");
    AddParameter("Reduce", "Find the optimal NNLS solution, then eliminate"
            " basis members until tolerance is reached", true);
    AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>());
    AddParameter("MaxCachedFrames", "The maximum number of frames this module is allowed to have cached", (size_t)(3000));
    AddParameter("RemoveWaveforms", "Remove the input waveforms?", bool(false));
    AddParameter("TimeLimitSeconds", "Time limit for each thread", (double)(60.0 * 10));

}

void CCMWavedeform::Configure() {
    GetParameter("NumThreads", num_threads);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTChannelMapName", pmt_channel_map_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("Waveforms", waveforms_name_);
    GetParameter("Output", output_name_);
    GetParameter("SPEsPerBin", spes_per_bin_);
    GetParameter("Tolerance", tolerance_);
    GetParameter("NoiseThreshold", noise_threshold_);
    GetParameter("BasisThreshold", basis_threshold_);
    GetParameter("Reduce", reduce_);
    GetParameter("PMTKeys", allowed_pmt_keys_);
    GetParameter("MaxCachedFrames", max_cached_frames);
    GetParameter("RemoveWaveforms", remove_waveforms);
    GetParameter("TimeLimitSeconds", time_limit_seconds);
    teptr = nullptr;
}

CCMWavedeform::~CCMWavedeform() {
    for(size_t i=0; i<cholmod_common_vec_.size(); ++i) {
        for(size_t j=0; j<cholmod_common_vec_[i].size(); ++j) {
            cholmod_l_finish(&(cholmod_common_vec_[i][j]));
        }
    }
}

void CCMWavedeform::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);

    if(pmt_channel_map_name_ != "") {
        if(not frame->Has(pmt_channel_map_name_)) {
            log_fatal("Could not find PMTChannelMap object with the key named \"%s\" in the Geometry frame.", pmt_channel_map_name_.c_str());
        }
        boost::shared_ptr<const I3Map<CCMPMTKey, uint32_t>> pmt_channel_map = frame->Get<boost::shared_ptr<const I3Map<CCMPMTKey, uint32_t>>>(pmt_channel_map_name_);
        boost::shared_ptr<const I3Map<CCMPMTKey, int>> pmt_channel_map_alt = frame->Get<boost::shared_ptr<const I3Map<CCMPMTKey, int>>>(pmt_channel_map_name_);
        if(pmt_channel_map) {
            pmt_channel_map_ = *pmt_channel_map;
        } else if(pmt_channel_map_alt) {
            for(std::pair<CCMPMTKey const, int> const & p : *pmt_channel_map_alt) {
                pmt_channel_map_[p.first] = int32_t(p.second);
            }
        } else {
            log_fatal("Could not find PMTChannelMap object with the key named \"%s\" in the Geometry frame.", pmt_channel_map_name_.c_str());
        }
    } else {
        pmt_channel_map_ = geo.pmt_channel_map;
    }

    geo_seen = true;
    templates_.clear();
    for(auto pmt_channel_pair : pmt_channel_map_) {
        templates_[pmt_channel_pair.first] = CCMWaveformTemplate();
        templates_[pmt_channel_pair.first].filled = false;
    }

    I3Map<CCMPMTKey, CCMOMGeo> const & pmt_geo_ = geo.pmt_geo;
    std::set<CCMPMTKey> allowed_pmt_keys(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end());
    bool filter_pmts = allowed_pmt_keys.size() > 0;
    for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : pmt_geo_) {
        if(filter_pmts and allowed_pmt_keys.find(p.first) == allowed_pmt_keys.end()) {
            continue;
        }
        if(pmt_channel_map_.find(p.first) == pmt_channel_map_.end()) {
            continue;
        }
        if(p.second.omtype == CCMOMGeo::OMType::CCM8inCoated or p.second.omtype == CCMOMGeo::OMType::CCM8inUncoated or p.second.omtype == CCMOMGeo::OMType::CCM1in) {
            pmt_keys_.push_back(p.first);
        }
    }

    /*
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
    */
    thread_ranges_.emplace_back(0, pmt_keys_.size());

    cholmod_common_vec_.resize(num_threads);
    for(size_t i=0; i<cholmod_common_vec_.size(); ++i) {
        cholmod_common_vec_[i].resize(pmt_keys_.size());
        for(size_t j=0; j<cholmod_common_vec_[i].size(); ++j) {
            cholmod_l_start(&(cholmod_common_vec_[i][j]));
        }
    }

    PushFrame(frame);
}

void CCMWavedeform::Calibration(I3FramePtr frame) {
    calibration = frame->Get<CCMCalibration>("CCMCalibration");

    double start_time = 2 * -10;
    double end_time = 2 * 60;

    wf_bin_width_ = 2.0;
    double template_bin_spacing = wf_bin_width_ / spes_per_bin_;
    double range = end_time - start_time;

    for(size_t i=0; i<pmt_keys_.size(); ++i) {
        // let's get our wf and what not
        CCMPMTKey key = pmt_keys_[i];
        uint32_t channel = pmt_channel_map_[key];

        if(calibration.pmtCal.count(key) == 0) {
            continue;
        }

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);

        double placeholder = 1.0;

        if(not templates_.at(key).filled) {
            int template_bins = (int)ceil(range / template_bin_spacing);
            FillTemplate(templates_.at(key), calib->second, start_time, template_bins, template_bin_spacing, frame);
            std::stringstream ss;
            ss << key << " Template(" << templates_.at(key).start_time << ", " << templates_.at(key).end_time << ") FWHM(" << templates_.at(key).digitizerStart << ", " << templates_.at(key).digitizerStop << ")" << std::endl;
            log_info("%s", ss.str().c_str());
        }
    }

    PushFrame(frame);
}

void CCMWavedeform::Process() {
  if (inbox_)
    log_trace("%zu frames in inbox", inbox_->size());

  I3FramePtr frame = PopFrame();

  if(!frame or frame->GetStop() == I3Frame::DAQ) {
    DAQ(frame);
    return;
  }

  if(frame->GetStop() == I3Frame::Physics && ShouldDoPhysics(frame))
    {
      Physics(frame);
    }
  else if(frame->GetStop() == I3Frame::Geometry && ShouldDoGeometry(frame))
    Geometry(frame);
  else if(frame->GetStop() == I3Frame::Calibration && ShouldDoCalibration(frame))
    Calibration(frame);
  else if(frame->GetStop() == I3Frame::DetectorStatus && ShouldDoDetectorStatus(frame))
    DetectorStatus(frame);
  else if(frame->GetStop() == I3Frame::Simulation && ShouldDoSimulation(frame))
    Simulation(frame);
  else if(frame->GetStop() == I3Frame::DAQ && ShouldDoDAQ(frame)) {
    DAQ(frame);
  } else if(ShouldDoOtherStops(frame))
    OtherStops(frame);
}

void CCMWavedeform::DAQ(I3FramePtr frame) {
    FrameWorkspace * workspace = new FrameWorkspace(frame, pmt_keys_, waveforms_name_);
    frame_workspaces[frame_index] = workspace;
	while(true) {

		// Check if any jobs have finished
		for(int i=int(running_jobs.size())-1; i>=0; --i) {
			if (teptr) {
				try{
					std::rethrow_exception(teptr);
				}
				catch(const std::exception &ex)
				{
					std::cerr << "Thread exited with exception: " << ex.what() << "\n";
				}
			}
			if(not running_jobs[i]->running.load()) {
				WavedeformJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                frame_workspaces[job->frame_index]->jobs_done += 1;
                if(frame_workspaces[job->frame_index]->jobs_done == frame_workspaces[job->frame_index]->total_jobs) {
                    results[job->frame_index - min_frame_idx].done = true;
                }
            } else {
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                FrameWorkspace * workspace_ptr = frame_workspaces[results[i].frame_index];
                workspace_ptr->Finish(remove_waveforms, waveforms_name_, output_name_);
                PushFrame(results[i].frame);
                delete workspace_ptr;
                frame_workspaces.erase(results[i].frame_index);
                results[i].frame = nullptr;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_frame_idx += results_done;
        }

        if(not frame)
            break;

        // Attempt to queue up a new job for the frame
        WavedeformJob * job = nullptr;

        if(free_jobs.size() > 0) {
            job = free_jobs.front();
            job->running.store(false);
            free_jobs.pop_front();
        } else if(running_jobs.size() < num_threads) {
            job = new WavedeformJob();
            job->running.store(false);
            job->thread_index = running_jobs.size();
        }

        if(job != nullptr and results.size() < max_cached_frames) {
            job->running.store(true);
            running_jobs.push_back(job);
            job->frame = frame;
            job->frame_index = frame_index;
            if(workspace->jobs_queued == 0) {
                results.emplace_back();
                results.back().frame = frame;
                results.back().frame_index = frame_index;
                results.back().done = false;
            }
            // Starts a job and increments workspace->jobs_queued
            workspace->StartJob(job,
                waveforms_name_,
                output_name_,
                pmt_keys_,
                pmt_channel_map_,
                templates_,
                calibration,
                cholmod_common_vec_,
                reduce_,
                wf_bin_width_,
                noise_threshold_,
                basis_threshold_,
                spes_per_bin_,
                tolerance_,
                num_threads,
                remove_waveforms,
                time_limit_seconds,
                teptr);
            if(workspace->jobs_queued == workspace->total_jobs) {
                break;
            }
        } else if(job != nullptr) {
            sleep(1);
            free_jobs.push_back(job);
        }
    }
    frame_index += 1;
}

void CCMWavedeform::Finish() {
    while(running_jobs.size() > 0) {
        // Check if any jobs have finished
        for(int i=int(running_jobs.size())-1; i>=0; --i) {
            if (teptr) {
                try{
                    std::rethrow_exception(teptr);
                }
                catch(const std::exception &ex)
                {
                    std::cerr << "Thread exited with exception: " << ex.what() << "\n";
                }
            }
            if(not running_jobs[i]->running.load()) {
                WavedeformJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                frame_workspaces[job->frame_index]->jobs_done += 1;
                if(frame_workspaces[job->frame_index]->jobs_done == frame_workspaces[job->frame_index]->total_jobs) {
                    results[job->frame_index - min_frame_idx].done = true;
                }
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                FrameWorkspace * workspace_ptr = frame_workspaces[results[i].frame_index];
                workspace_ptr->Finish(remove_waveforms, waveforms_name_, output_name_);
                PushFrame(results[i].frame);
                delete workspace_ptr;
                frame_workspaces.erase(results[i].frame_index);
                results[i].frame = nullptr;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_frame_idx += results_done;
        }
    }
}



