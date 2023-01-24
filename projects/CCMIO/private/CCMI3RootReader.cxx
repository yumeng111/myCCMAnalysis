#include <icetray/IcetrayFwd.h>
#include <icetray/I3DefaultName.h>

#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <fstream>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Context.h>
#include <icetray/I3Logging.h>
#include <icetray/I3FrameObject.h>
#include <icetray/serialization.h>

#include "dataio/I3FileStager.h"
#include "dataio/I3FrameSequence.h"

#include "CCMAnalysis/CCMIO/CCMRootHandle.h"

#include "CCMAnalysis/CCMDataStructures/AccumWaveform.h"
#include "CCMAnalysis/CCMDataStructures/Events.h"
#include "CCMAnalysis/CCMDataStructures/MCTruth.h"
#include "CCMAnalysis/CCMDataStructures/Pulses.h"
#include "CCMAnalysis/CCMDataStructures/RawData.h"

#include "TFile.h"

class CCMI3RootReader : public I3Module {
    unsigned nframes_;
    std::vector<std::string> filenames_;
    std::vector<std::string> load_;
    std::string tree_name_;

    boost::shared_ptr<TFile> file_ptr_;
    CCMRootHandlePtr file_handle_;
    std::string current_filename_;

    std::vector<std::string>::iterator filenames_iter_;

    void OpenNextFile();
    void LoadItems(std::map<I3Frame::Stream, I3FramePtr> & frame_map);
    template<typename T>
    void LoadItem(std::map<I3Frame::Stream, I3FramePtr> & frame_map, std::string const & branch_name, I3Frame::Stream stream) {
        boost::shared_ptr<T> ptr = file_handle_->Get<T>(branch_name);
        if(ptr == nullptr)
            return;
        if(frame_map.count(stream) == 0)
            frame_map[stream] = boost::shared_ptr<I3Frame>(new I3Frame(stream));
        frame_map[stream]->Put(I3DefaultName<T>::value(), ptr, stream);
    }

public:
    CCMI3RootReader(const I3Context&);
    void Configure();
    void Process();

    SET_LOGGER("CCMI3RootReader");
};

I3_MODULE(CCMI3RootReader);

namespace detail {
    std::string tolower(std::string const & s) {
        std::string r(s);
        std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c){return std::tolower(c);});
        return r;
    }
}

CCMI3RootReader::CCMI3RootReader(const I3Context& context) : I3Module(context), nframes_(0) {
  AddParameter("Filename",
	       "Filename to read.  Use either this or Filenamelist, not both.",
	       std::string());

  AddParameter("FilenameList",
	       "List of files to read, *IN SORTED ORDER*",
	       std::vector<std::string>());

  AddParameter("KeysToLoad",
          "Load frame objects with these keys",
          std::vector<std::string>());

  AddParameter("TreeName",
          "Tree name to load objects from",
          std::string());

  AddOutBox("OutBox");
}

void CCMI3RootReader::Configure() {
    std::string fname;

    GetParameter("FileNameList", filenames_);
    GetParameter("FileName", fname);
    if(!filenames_.empty() && !fname.empty())
        log_fatal("Both Filename and FileNameList were specified.");

    if (filenames_.empty() && fname.empty())
        log_fatal("Neither 'Filename' nor 'FilenameList' specified");

    if (filenames_.empty())
        filenames_.push_back(fname);

    GetParameter("KeysToLoad", load_);

    GetParameter("TreeName", tree_name_);

    std::vector<std::string> real_keys = {
        "rawData",
        "pulses",
        "events",
        "accumWaveform",
        "mcTruth",
    };

    std::map<std::string, std::string> key_map;
    for(std::string const & s : real_keys)
        key_map.insert({detail::tolower(s), s});

    if(load_.size() == 0) {
        load_ = real_keys;
    } else {
        std::vector<std::string> keys;
        for(std::string const & s : load_) {
            std::map<std::string, std::string>::iterator it = key_map.find(detail::tolower(s));
            if(it == key_map.end())
                continue;
            keys.push_back(it->second);
        }
        load_ = keys;
        if(load_.size() == 0) {
            throw std::runtime_error("KeysToLoad specified, but no keys match those in the valid keys list.");
        }
    }
    filenames_iter_ = filenames_.begin();
    OpenNextFile();
}

void CCMI3RootReader::LoadItems(std::map<I3Frame::Stream, I3FramePtr> & frame_map) {
    for(std::string const & s : load_) {
        switch(s[0]) {
            case 'a': // accumWaveform
                LoadItem<AccumWaveform>(frame_map, s, I3Frame::Physics);
                break;
            case 'e': // events
                LoadItem<Events>(frame_map, s, I3Frame::Physics);
                break;
            case 'm': // mcTruth
                LoadItem<MCTruth>(frame_map, s, I3Frame::DAQ);
                break;
            case 'p': // pulses
                LoadItem<Pulses>(frame_map, s, I3Frame::DAQ);
                break;
            case 'r': // rawData
                LoadItem<RawData>(frame_map, s, I3Frame::DAQ);
                break;
            default:
                break;
        }
    }
}

void CCMI3RootReader::Process() {
    while(file_handle_ == nullptr or not file_handle_->More()) {
        if(filenames_iter_ == filenames_.end()) {
            RequestSuspension();
            file_handle_ = nullptr;
            return;
        } else
            OpenNextFile();
    }

    std::map<I3Frame::Stream, I3FramePtr> frame_map;
    std::vector<I3Frame::Stream> frame_keys = {I3Frame::DAQ, I3Frame::Physics};
    try {
        LoadItems(frame_map);
    } catch (const std::exception &e) {
        log_fatal("Error reading %s at frame %d: %s!",
                current_filename_.c_str(), nframes_, e.what());
        return;
    }

    for(I3Frame::Stream stream : frame_keys) {
        if(frame_map.count(stream) == 0)
            continue;
        nframes_++;
        PushFrame(frame_map[stream], "OutBox");
    }
}

void CCMI3RootReader::OpenNextFile() {

  current_filename_ = *filenames_iter_;
  nframes_ = 0;
  filenames_iter_++;

  if(file_ptr_ != nullptr)
      file_ptr_->Close();
  file_ptr_ = boost::make_shared<TFile>(current_filename_.c_str(), "READ");
  if(tree_name_ != "") {
      if(load_.size() > 0) {
          file_handle_ = boost::make_shared<CCMRootHandle>(file_ptr_.get(), tree_name_, load_);
      } else {
          file_handle_ = boost::make_shared<CCMRootHandle>(file_ptr_.get(), tree_name_);
      }
  } else {
      if(load_.size() > 0) {
          file_handle_ = boost::make_shared<CCMRootHandle>(file_ptr_.get(), load_);
      } else {
          file_handle_ = boost::make_shared<CCMRootHandle>(file_ptr_.get());
      }
  }
  log_trace("Constructing with filename %s, %zu keys",
	    current_filename_.c_str(), load_);

  log_info("Opened file %s", current_filename_.c_str());
}
