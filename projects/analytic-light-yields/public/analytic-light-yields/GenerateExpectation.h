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
    SodiumVertexDistribution* sodium_events_constructor = nullptr;
    YieldsPerPMT* yields_and_offset_constructor = nullptr;
    LArScintillationLightProfile* LAr_scintillation_light_constructor = nullptr;

    boost::shared_ptr<HESodiumEventSeries> event_vertices = boost::make_shared<HESodiumEventSeries> ();
    boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_per_pmt = boost::make_shared<PhotonYieldSummarySeriesMap> ();

public:
    GenerateExpectation();
    void GetSodiumVertices(size_t const & n_events_to_simulate, double const & z_position);
    void GetYieldsAndOffsets(I3FramePtr geo_frame, double const & uv_absorption);
    I3Vector<double> LightProfile(double const & Rs, double const & Rt, double const & tau_s, double const & tau_t, double const & tau_rec, double const & tau_TPB,
                                       AnalyticLightYieldGenerator::LArLightProfileType const & light_profile_type);
    boost::shared_ptr<I3MapPMTKeyVectorDouble> GetExpectation(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame);
};
#endif
