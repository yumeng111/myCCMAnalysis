#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <chrono>
#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <math.h>

//#include <dataio/I3File.h>
//#include <icetray/ctpl.h>
//#include <icetray/open.h>
//#include <icetray/I3Frame.h>
//#include <icetray/I3Units.h>
//#include <icetray/I3TrayInfo.h>
//#include <icetray/I3Module.h>
//#include <icetray/I3Logging.h>
//#include <icetray/I3PODHolder.h>
//#include <icetray/CCMPMTKey.h>
//#include <icetray/CCMTriggerKey.h>
//#include <dataclasses/I3Double.h>
//#include <dataclasses/geometry/CCMGeometry.h>
//#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
//#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include <analytic-light-yields/G4YieldsPerPMT.h>

G4YieldsPerPMT::G4YieldsPerPMT() {}

std::vector<double> G4YieldsPerPMT::GetPlottingInformation(CCMPMTKey key, double uv_absorption, double scaling){
    
    // grab frames
    std::deque<I3FramePtr> G4Events;
    std::string g4_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumCenterHEEvents.i3.zst";
    dataio::I3File g4_file(g4_fname, dataio::I3File::Mode::read);
    while (g4_file.more()){
        I3FramePtr g4_frame = g4_file.pop_frame();
        G4Events.push_back(g4_frame);
    }
    
    // now get all yields
    size_t n_threads = 0;
    double max_time = 50.0;
    std::map<CCMPMTKey, std::vector<double>> binned_yields;
    std::map<CCMPMTKey, std::vector<double>> binned_square_yields;

    GetAllYields<double>(n_threads, {key}, G4Events, G4Events, false, max_time, 0.0, 0.0, uv_absorption, 0.0, scaling, binned_yields, binned_square_yields);

    // now grab yields for this pmt
    return binned_yields.at(key);
}

void G4YieldsPerPMT::TimeComparison(std::vector<CCMPMTKey> keys_to_fit, double z_offset){
    
    // make sure we have filled G4Events_
    if (grabbed_g4_events == false){
        std::vector<std::string> g4_fnames = {"/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumCenterHEEvents.i3.zst", "/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumPlus50HEEvents.i3.zst"};
        std::vector<double> z_locs = {0.0, 50.0};

        for (size_t f = 0; f < g4_fnames.size(); f++){
            std::deque<I3FramePtr> this_file_events;

            dataio::I3File g4_file(g4_fnames.at(f), dataio::I3File::Mode::read);
            while (g4_file.more()){
                I3FramePtr g4_frame = g4_file.pop_frame();
                this_file_events.push_back(g4_frame);
            }

            G4Events_[z_locs.at(f)] = this_file_events;
        }

    }
    
    // some constants
    double c = 2.998 * std::pow(10.0, 8.0); // speed of light in m/s
    double c_cm_per_nsec = c * std::pow(10.0, -7.0); // speed of light in cm/nsec
    double uv_index_of_refraction = 1.358;
    double vis_index_of_refraction = 1.23;
    
    for (size_t event_it = 0; event_it < G4Events_[z_offset].size(); event_it ++){
        I3FramePtr frame = G4Events_[z_offset].at(event_it);

        // lets grab necessary things from each frame
        boost::shared_ptr<CCMMCPESeriesMap const> CCMMCPEMap = frame->Get<boost::shared_ptr<CCMMCPESeriesMap const>>("PMTMCHitsMap");
        boost::shared_ptr<PhotonSummarySeries const> photon_summary_series = frame->Get<boost::shared_ptr<PhotonSummarySeries const>>("PhotonSummarySeries");
        boost::shared_ptr<I3Map<int, size_t> const> photon_summary_series_map = frame->Get<boost::shared_ptr<I3Map<int, size_t> const>>("PhotonSummaryMap");

        // now loop over all pmts we are fitting
        for(size_t k = 0; k < keys_to_fit.size(); k++){
            
            // grab the list of CCMMCPEs for this key
            CCMMCPESeries this_key_ccmmcpe = CCMMCPEMap->at(keys_to_fit.at(k));
            
            // loop over each CCMMCPE 
            for (size_t m = 0; m < this_key_ccmmcpe.size(); m++){
                CCMMCPE this_ccmmcpe = this_key_ccmmcpe.at(m);
                double g4_time = this_ccmmcpe.global_time;

                // first check the wavelength!! only continue if visible
                double wavelength = this_ccmmcpe.wavelength * 1e9; // units of nm
                if (wavelength > 325.0){
                    // grab photon summary information related with this photon hit
                    size_t track_id = this_ccmmcpe.track_id;
                    std::map<int, size_t>::const_iterator it = photon_summary_series_map->find(track_id);
                    if (it != photon_summary_series_map->end()) {
                        PhotonSummary this_photon_summary = photon_summary_series->at(it->second);
                        double distance_travelled_uv = this_photon_summary.distance_uv / I3Units::cm;
                        double distance_travelled_vis = this_photon_summary.distance_visible / I3Units::cm;
                        
                        double travel_time_uv = distance_travelled_uv / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec
                        double travel_time_visible = distance_travelled_vis / (c_cm_per_nsec / vis_index_of_refraction); // units of nsec
                        double calculated_time = travel_time_uv + travel_time_visible;
                        
                        // now save
                        G4Times_.push_back(g4_time);
                        CalculatedTimes_.push_back(calculated_time);
                    }
                }
            }
        }
    }
}


