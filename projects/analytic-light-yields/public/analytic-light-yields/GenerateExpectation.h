#ifndef GenerateExpectation_H
#define GenerateExpectation_H

#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

#include <tuple>
#include <vector>
#include <random>
#include <cmath>
#include <map>
#include <memory>

#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/physics/PhotonYieldSummary.h"
#include "dataclasses/physics/HESodiumEvent.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"
#include <analytic-light-yields/SodiumVertexDistribution.h>
#include <analytic-light-yields/YieldsPerPMT.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>

class GenerateExpectation {
    // things we need to do in order to generate expected hits/pmt :
    // 1) get soidum vertex distribution
    // 2) get yields + time offsets per pmt
    // 3) get light profile
    // then we put it all together and bin!
    // so for now, we are going to set sodium vertices ONCE and get yields + offsets ONCE
    // and update light profile with every call to GetExpectation
    std::shared_ptr<SodiumVertexDistribution> sodium_events_constructor = nullptr;
    std::shared_ptr<YieldsPerPMT> yields_and_offset_constructor = nullptr;
    std::shared_ptr<LArScintillationLightProfile> LAr_scintillation_light_constructor = nullptr;

    boost::shared_ptr<HESodiumEventSeries> event_vertices = boost::make_shared<HESodiumEventSeries> ();
    std::vector<boost::shared_ptr<PhotonYieldSummarySeriesMap>> yields_per_pmt_per_event;
    
    // let's also add small noise to every bin
    double noise_photons = 1.0;
    double noise_triggers = 5.0;
    double digitization_time = 16.0 * std::pow(10.0, 3.0); //16 usec in nsec
    double noise_rate = noise_photons / (noise_triggers * digitization_time); // units of photons/nsec
    double noise_rate_per_time_bin = 2.0 * noise_rate; // 2nsec binning

public:
    GenerateExpectation();
    void GetSodiumVertices(size_t const & n_events_to_simulate, double const & z_position);
    void GetYieldsAndOffsets(I3FramePtr geo_frame, double const & uv_absorption);
    I3Vector<double> LightProfile(double const & Rs, double const & Rt, double const & tau_s, double const & tau_t, double const & tau_rec, double const & tau_TPB, double const & time_offset,
                                       AnalyticLightYieldGenerator::LArLightProfileType const & light_profile_type);
    I3Vector<double> DLightProfile(double const & Rs, double const & Rt, double const & tau_s, double const & tau_t, double const & tau_TPB, double const & time_offset, std::string deriv_variable);
    std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>,boost::shared_ptr<I3MapPMTKeyVectorDouble>> GetExpectation(AnalyticLightYieldGenerator analytic_light_yield_setup,
                                                                                                      I3FramePtr geo_frame, double const & time_offset);
    //boost::shared_ptr<I3MapPMTKeyVectorDouble> GetExpectation(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame);
    std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>>
               GetExpectationWithDerivs(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, double const & time_offset);
};
#endif
