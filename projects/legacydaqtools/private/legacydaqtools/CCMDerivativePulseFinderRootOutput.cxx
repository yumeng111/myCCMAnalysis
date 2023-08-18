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
#include <algorithm>
#include <thread>
#include <cmath>
#include <functional>
#include <memory>
#include <cstdlib>

#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/BaselineEstimate.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <daqtools/WaveformDerivative.h>

#include "CCMAnalysis/CCMDataStructures/Pulses.h"
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

#include <TTree.h>
#include <TFile.h>


class CCMDerivativePulseFinderRootOutput: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string ccm_trigger_name_;
    std::string ccm_waveform_name_;
    std::string nim_pulses_name_;
    bool throw_without_input_;
    std::string output_name_;

    std::string output_file_name_;

    double pulse_initial_deriv_threshold;
    size_t pulse_max_width;
    size_t pulse_min_width;
    double pulse_min_height;
    double pulse_min_deriv_magnitude;
    double pulse_min_integral;

    TTree * fOutEventTree;
    std::shared_ptr<TFile> fOutputFile;
    Pulses * fPulses;
    size_t current_event_number_;
    std::map<size_t, int> board_idx_map;

    CCMGeometry geo;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    public:
    CCMDerivativePulseFinderRootOutput(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
    void Finish();
    static std::vector<std::tuple<CCMRecoPulse, double, double>> CheckForPulse(WaveformDerivative & deriv, size_t start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold);
    static void FindPulses(
            CCMRecoPulseSeries & output,
            std::vector<double> & max_derivative,
            std::vector<double> & max_amplitude,
            WaveformDerivative & deriv,
            double deriv_threshold,
            size_t max_pulse_width,
            size_t min_pulse_width,
            double min_pulse_height,
            double min_deriv_magnitude,
            double min_integral);
    static size_t GetBoardEventNumber(boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger> const> triggers);
    static timespec GetEventTime(boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger> const> triggers);
};

I3_MODULE(CCMDerivativePulseFinderRootOutput);

CCMDerivativePulseFinderRootOutput::CCMDerivativePulseFinderRootOutput(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false), current_event_number_(0) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("CCMTriggersName", "Key for CCMTriggers", std::string("CCMTriggers"));
        AddParameter("NIMPulsesName", "Key for NIMPulses", std::string("NIMPulses"));
        AddParameter("CCMCalibratedWaveformsName", "Key for the input CCMWaveformDoubleSeries", std::string("CCMCalibratedWaveforms"));
        AddParameter("ThrowWithoutInput", "Whether to throw an error when there is no input", false);

        AddParameter("InitialDerivativeThreshold", "Initial positive derivative threshold for a pulse", double(0.3));
        AddParameter("MinPulseWidth", "Minimum width for defining a pulse", size_t(5));
        AddParameter("MaxPulseWidth", "Maxiumum width for defining a pulse", size_t(100));
        AddParameter("MinPulseHeight", "Minimum height for defining a pulse", double(5.0));
        AddParameter("MinPulseDerivativeMagnitude", "Minimum derivative magnitude for defining a pulse", double(0.65));
        AddParameter("MinPulseIntegral", "Minimum integral for defining a pulse", double(25.0));
        AddParameter("OutputName", "Key for output CCMRecoPulseSeriesMap", std::string("DerivativePulses"));
        AddParameter("RootOutputFileName", "A path and filename for the output root file", std::string("test.root"));
    }

void CCMDerivativePulseFinderRootOutput::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMTriggerName", ccm_trigger_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("CCMCalibratedWaveformsName", ccm_waveform_name_);
    GetParameter("ThrowWithoutInput", throw_without_input_);
    GetParameter("InitialDerivativeThreshold", pulse_initial_deriv_threshold);
    GetParameter("MinPulseWidth", pulse_min_width);
    GetParameter("MaxPulseWidth", pulse_max_width);
    GetParameter("MinPulseHeight", pulse_min_height);
    GetParameter("MinPulseDerivativeMagnitude", pulse_min_deriv_magnitude);
    GetParameter("MinPulseIntegral", pulse_min_integral);
    GetParameter("OutputName", output_name_);
    GetParameter("RootOutputFileName", output_file_name_);

    // #### ED ####
    // Open the root file and set up data structures here
    //Currently copying necessary elements from the CCMEventTreeHandle::SetupOutputFile method. takes in a TFile &f object and creates an EventTree Object with a Pulses branch. 
    //Initiate these variables earlier:
    //TTree * fOutEventTree, TFile * fOutputFile = "${filename}_pulses.root", Pulses * fPulses
    // You can open the file here, and use the output_file_name_ variable
    static const char* gsEvents= "EventTree";

    //Out with the old                                                                                
    if (fOutEventTree != nullptr) {
        fOutEventTree = nullptr;
    }

    //Need to read in an output file name from somewhere, modifiable.
    //Set up an output file for writing events 
    fOutputFile = std::make_shared<TFile>(output_file_name_.c_str(),"RECREATE","CCM ROOT Event File",3);
    //fOutputFile = &f;

    fOutputFile->cd();
    fOutEventTree = new TTree(gsEvents, gsEvents);
    if(fOutEventTree == nullptr) {
        log_fatal("Problem setting up output tree!");
        //Report error; not sure if MsgError still in code. 
        //return 0;
    }

    int sLvl = 99;

    //Don't need to check if want to save Pulses -> set up module to save pulses. 
    //if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"pulses") != fSaveBranches.end() ||
    //std::find(fSaveBranches.begin(),fSaveBranches.end(),"all") != fSaveBranches.end()) {

    //Create Pulses branch and add it to the created tree. 
    if (fPulses == nullptr) {
        fPulses = new Pulses();
    }
    fOutEventTree->Branch("pulses", &fPulses, 320000, sLvl);
}

/*
   We might need this Advance method in order to move through the root file. 

//_____________________________________________________________                                     
uint32_t CCMEventTreeHandle::Advance(uint32_t n)
{
//advance n places in the input stream.  Return the difference                                    
//between the new position and the old.                                                           

if(fEventsTree == nullptr || n < 1) {
return 0;
}

uint32_t indexSave = fIndex;

fIndex += n;
if(fIndex >= fNumOfEntries) {
fIndex = fNumOfEntries - 1;
}

//Since we've moved on, mark all branches as unloaded                                             
this->ClearLoadFlags();

return fIndex-indexSave;
}
*/


void CCMDerivativePulseFinderRootOutput::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    if(not frame->Has("CCMDAQConfig")) {
        log_fatal("Could not find CCMDAQConfig object with the key named \"CCMDAQConfig\" in the Geometry frame.");
    }
    geo = frame->Get<CCMGeometry const>(geometry_name_);
    CCMAnalysis::Binary::CCMDAQConfig const & config = frame->Get<CCMAnalysis::Binary::CCMDAQConfig const>("CCMDAQConfig");

    for(size_t board_idx=0; board_idx<config.digitizer_boards.size(); ++board_idx) {
        std::string board_name = config.digitizer_boards[board_idx].physical_board_id;
        std::string::size_type underscore_idx = board_name.rfind("_");
        if(underscore_idx == std::string::npos) {
            log_fatal("Could not find underscore in physical board id");
        }
        std::string board_id_str = board_name.substr(underscore_idx + 1, std::string::npos);
        int physical_board_id = std::atoi(board_id_str.c_str());
        board_idx_map[board_idx] = physical_board_id - 1;
    }


    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void CCMDerivativePulseFinderRootOutput::DAQ(I3FramePtr frame) {

    if(not frame->Has(ccm_waveform_name_)) {
        if(throw_without_input_) {
            log_fatal("Input waveform named %s not present", ccm_waveform_name_.c_str());
        } else {
            log_warn("Input waveform named %s not present", ccm_waveform_name_.c_str());
        }
        return;
    }
    if(not frame->Has(ccm_trigger_name_)) {
        if(throw_without_input_) {
            log_fatal("Input CCMTriggers name %s not present", ccm_trigger_name_.c_str());
        } else {
            log_warn("Input CCMTriggers named %s not present", ccm_trigger_name_.c_str());
        }
        return;
    }

    // let's read in our waveform
    boost::shared_ptr<CCMWaveformDoubleSeries const> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>(ccm_waveform_name_);
    boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger> const> triggers = frame->Get<boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger> const>>(ccm_trigger_name_);
    boost::shared_ptr<I3Map<CCMTriggerKey, NIMLogicPulse> const> nim_pulses = frame->Get<boost::shared_ptr<I3Map<CCMTriggerKey, NIMLogicPulse> const>>(nim_pulses_name_);

    std::vector<CCMPMTKey> pmt_keys;
    pmt_keys.reserve(pmt_channel_map_.size());
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        pmt_keys.push_back(p.first);
    }

    boost::shared_ptr<CCMRecoPulseSeriesMap> output(new CCMRecoPulseSeriesMap);

    // #### ED ####
    // Push the data to the file
    // Reset the Pulses object so it is empty for the next trigger
    // clear pulses and raw data for the new trigger window
    //Austin Note: I think this may need to be before the pulse loop. We would then call an Advance on the root file in order to move to the next trigger. 
    fPulses->Reset();

    size_t event_number = std::max(current_event_number_ + 1, GetBoardEventNumber(triggers));
    current_event_number_ = event_number;

    fPulses->SetEventNumber(event_number);
    struct timespec event_time = GetEventTime(triggers);
    fPulses->SetComputerSecIntoEpoch(event_time.tv_sec);
    fPulses->SetComputerNSIntoSec(event_time.tv_nsec);

    // loop over each channel in waveforms
    for(size_t i=0; i<pmt_keys.size(); ++i) {
        CCMPMTKey pmt_key = pmt_keys[i];
        uint32_t channel = pmt_channel_map_[pmt_key];
        CCMTriggerKey trigger_key =  geo.trigger_copy_map[pmt_key];
        double nim_pulse_time = nim_pulses->at(trigger_key).GetNIMPulseTime();
        //uint32_t trigger_channel = geo.trigger_channel_map[trigger_key];
        WaveformDerivative deriv(waveforms->at(channel).GetWaveform().begin(), waveforms->at(channel).GetWaveform().end(), 2.0);
        std::vector<double> max_derivative;
        std::vector<double> max_amplitude;
        FindPulses(output->operator[](pmt_key),
                max_derivative,
                max_amplitude,
                deriv,
                pulse_initial_deriv_threshold,
                pulse_max_width,
                pulse_min_width,
                pulse_min_height,
                pulse_min_deriv_magnitude,
                pulse_min_integral);
        // #### ED ####
        // Assume a SinglePulse object is available here
        //initiate earlier again: SinglePulse * fCurrPulse
        SinglePulse fCurrPulse;
        int key = (channel % 16) + board_idx_map[channel/16]*16;
        fCurrPulse.SetKey(key);

        float pmt_offset;
        if(geo.pmt_geo[pmt_key].omtype == CCMOMGeo::OMType::CCM1in)
            pmt_offset = 0.0;
        else
            pmt_offset = 21.0;

        for(size_t pulse_idx=0; pulse_idx<output->at(pmt_key).size(); ++pulse_idx) {
            CCMRecoPulse const & pulse = output->at(pmt_key)[pulse_idx];
            fCurrPulse.Reset();
            fCurrPulse.SetTriggerOffset(nim_pulse_time / 2.0);//Channel 15 offset for the PMT/digitizer board
            fCurrPulse.SetPMTOffset(pmt_offset);//Offset timing for Veto vs. 8" PMT. 21.0 for Veto, 0 for 8" according to CCMFindPulses.cxx
            fCurrPulse.SetADCToPE(1.0);//initial ADC to PE value - set to 1 for all PMTs in the FindPulses module. 
            fCurrPulse.SetTime(pulse.GetTime() / 2.0);//Time of the Pulse
            fCurrPulse.SetLength(pulse.GetWidth() / 2.0);//Length of the Pulse
            fCurrPulse.SetMaxDerValue(max_derivative[pulse_idx]);//Location of the maximum derivative in the pulse. 
            fCurrPulse.SetAmplitude(max_amplitude[pulse_idx]);//Max amplitude of the pulse.
            fCurrPulse.SetIntegral(pulse.GetCharge());//Total integral of the Pulse
            fCurrPulse.SetBaseline(0.0);//Baseline value for the waveform. 
            /*for (int sampleLoc2 = sampleLoc; sampleLoc2 < sampleLoc+5 && sampleLoc2 != length-4; ++sampleLoc2) {
                fCurrPulse.AddSample(fSmoothArray[key][sampleLoc2]/smoothWidth);
                //This saves the waveform of the pulse as an array of samples (2ns bins). Should be skippable. 
            }
            if (sampleLoc+5 < length-4) {
                fCurrPulse.SetWaveformEnd(sampleLoc);//Location of the start of the pulse
            } else {
                fCurrPulse.SetWaveformEnd(length);//Location of the end of the pulse. 
            }*/
            // Assume a Pulses object is available here
            // Add the SinglePulse object to the Pulses object
            fPulses->AddSinglePulse(fCurrPulse);
            //Counter object: questionable necessity. Not saved anywhere I can see. (Does get outputted maybe?)
        }
    }

    // #### ED ####
    //Somehow call the advance to move forward in the output file. 

    frame->Put(output_name_, output);
    PushFrame(frame);
}

void CCMDerivativePulseFinderRootOutput::FindPulses(
        CCMRecoPulseSeries & output,
        std::vector<double> & max_derivative,
        std::vector<double> & max_amplitude,
        WaveformDerivative & deriv,
        double deriv_threshold,
        size_t max_pulse_width,
        size_t min_pulse_width,
        double min_pulse_height,
        double min_deriv_magnitude,
        double min_integral) {
    // Determine the maximum starting sample for a pulse
    // We need at least 5 samples after the starting pulse
    size_t N = deriv.Size();
    // Define the maximum sample to consider for the start of a pulse
    // A pulse needs at least 5 samples, so subtract off 5 from the max sample we can look at
    if(N > 5)
        N -= 5;
    size_t pulse_first_index = 0;
    size_t pulse_last_index = 0;

    // Reset the smoother position to the start of the waveform
    deriv.Reset();
    for(size_t i=0; i<N; ++i) {
        // Check the condition for the beginning of a pulse
        if(deriv.Derivative() > deriv_threshold) {
            // Get the last index of the found pulse, 0 if no pulse is found
            std::vector<std::tuple<CCMRecoPulse, double, double>> pulses = CheckForPulse(deriv, i, max_pulse_width, min_pulse_width, min_pulse_height, min_deriv_magnitude, min_integral);
            for(size_t j=0; j<pulses.size(); ++j) {
                output.push_back(std::get<0>(pulses[j]));
                pulse_first_index = std::get<0>(pulses[j]).GetTime() / 2;
                pulse_last_index = pulse_first_index + std::get<0>(pulses[j]).GetWidth() / 2;
                i = pulse_last_index;
            }
            deriv.Reset(pulse_last_index);
        }
        deriv.Next();
    }
}

size_t CCMDerivativePulseFinderRootOutput::GetBoardEventNumber(boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger> const> triggers) {
    uint32_t min_event_num;
    bool first = true;
    for(CCMAnalysis::Binary::CCMTrigger const & trigger : (*triggers)) {
        size_t board_idx = 0;
        for(uint32_t event_num : trigger.board_event_numbers) {
            bool has_data = false;
            for(size_t channel_idx=0; channel_idx<16; ++channel_idx) {
                if(trigger.channel_masks[board_idx * 16 + channel_idx]) {
                    has_data = true;
                    break;
                }
            }
            if(not has_data)
                continue;
            if(first) {
                min_event_num = event_num;
                first = false;
                continue;
            }
            min_event_num = std::min(min_event_num, event_num);
            board_idx += 1;
        }
    }
    return size_t(min_event_num);
}

timespec CCMDerivativePulseFinderRootOutput::GetEventTime(boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger> const> triggers) {
    using namespace CCMAnalysis::Binary;
    struct timespec min_time;
    bool first = true;
    for(CCMAnalysis::Binary::CCMTrigger const & trigger : (*triggers)) {
        size_t board_idx = 0;
        for(struct timespec event_time : trigger.board_computer_times) {
            bool has_data = false;
            for(size_t channel_idx=0; channel_idx<16; ++channel_idx) {
                if(trigger.channel_masks[board_idx * 16 + channel_idx]) {
                    has_data = true;
                    break;
                }
            }
            if(not has_data)
                continue;
            if(first) {
                min_time = event_time;
                first = false;
                continue;
            }
            if(event_time < min_time)
                min_time = event_time;
            board_idx += 1;
        }
    }
    return min_time;
}

std::vector<std::tuple<CCMRecoPulse, double, double>> CCMDerivativePulseFinderRootOutput::CheckForPulse(WaveformDerivative & deriv, size_t start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold) {
    deriv.Reset(start_idx);
    // 0 before start of pulse checking [checking for positive derivative]
    // 1 derivative is positive (in rising edge of pulse) [checking for negative derivative]
    // 2 derivative is negative (in falling edge of pulse) [checking for positive derivative]
    // 3 derivative is positive (recovering from droop) [done]
    int state = 1;
    double max_value = deriv.Value();
    double max_abs_derivative = deriv.Derivative();
    double integral = 0;
    bool found_pulse = false;
    size_t final_idx = start_idx;
    size_t N = std::min(deriv.Size(), start_idx + max_samples);
    for(size_t i=start_idx; i<N; ++i) {
        max_value = std::max(max_value, deriv.Value());
        max_abs_derivative = std::max(max_abs_derivative, std::abs(deriv.Derivative()));
        integral += deriv.Value();
        if(state % 2) {
            if(deriv.Derivative() < 0) {
                state += 1;
            }
        } else {
            if(state == 2 and deriv.Value() > 0.1) {
                continue;
            }
            if(deriv.Derivative() > 0) {
                state += 1;
            }
        }
        if(state >= 3) {
            found_pulse = true;
            final_idx = i;
            break;
        }
        deriv.Next();
    }
    if(found_pulse
            and (max_value >= value_threshold)
            and (max_abs_derivative >= derivative_threshold)
            and (integral >= integral_threshold)
            and (final_idx - start_idx >= min_length)
      ) {
        CCMRecoPulse pulse;
        pulse.SetTime(double(start_idx) * 2.0);
        pulse.SetCharge(integral);
        pulse.SetWidth((double(final_idx) - double(start_idx)) * 2.0);
        return {std::make_tuple(pulse, max_abs_derivative, max_value)};
    } else {
        deriv.Reset(start_idx);
        return {};
    }
}

void CCMDerivativePulseFinderRootOutput::Finish() {
    if (fPulses != nullptr) {
        delete fPulses;
    }
    // #### ED ####
    // Close the root file here
    if (fOutputFile != nullptr && fOutEventTree != nullptr) {
        //Get to the output file, write anything still in the Event Tree, then clear the Event Tree. 
        fOutputFile->cd();
        fOutEventTree->Write();
        fOutEventTree = nullptr;
	//Then write the file and close it. 
	fOutputFile->Flush();
	fOutputFile->Write();
	fOutputFile = nullptr;
    }

}


/*
//_____________________________________________________________                                     
int CCMEventTreeHandle::Write()
{
if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"none") != fSaveBranches.end()) {
return 1;
}

if (fOutputFile != nullptr && fOutEventTree != nullptr) {
//Make sure we load any branches which haven't been loaded yet so                               
//they get written                                                                              

this->Load(kNEventBranch);

fOutEventTree->Fill();
}
return 1;
}

//_____________________________________________________________                                     
void CCMEventTreeHandle::Close()
{
//MsgDebug(3,"Calling Close method");                                                             

if (fOutputFile != nullptr && fOutEventTree != nullptr) {
fOutputFile->cd();
fOutEventTree->Write();
fOutEventTree = nullptr;
}

}
*/
