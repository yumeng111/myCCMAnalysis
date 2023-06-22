/**
 *  $Id$
 *
 *  Copyright (C) 2007
 *  Troy D. Straszheim  <troy@icecube.umd.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 *  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 *  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 *
 *  SPDX-License-Identifier: BSD-2-Clause
 *
 */
#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <fstream>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>

#include "dataio/I3FileStager.h"

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

class CCMBinaryReader : public I3Module {

    unsigned nframes_;
    std::vector<std::string> filenames_;
    I3FileStagerPtr file_stager_;
    I3::dataio::shared_filehandle current_filename_;

    bool drop_duplicate_configs_;

    boost::iostreams::filtering_istream ifs_;
    std::ifstream file_stream_;

    boost::shared_ptr<CCMAnalysis::Binary::CCMDAQConfig> last_config_;
    size_t num_frames_in_file_;
    bool new_config_;
    boost::shared_ptr<CCMAnalysis::Binary::CCMTriggerReadout> readout_;

    std::vector<std::string>::iterator filenames_iter_;

    void OpenNextFile();

    public:

    CCMBinaryReader(const I3Context&);

    void Configure();
    void Process();

    ~CCMBinaryReader();
    SET_LOGGER("CCMBinaryReader");

};

namespace io = boost::iostreams;
namespace dataio = I3::dataio;

I3_MODULE(CCMBinaryReader);

CCMBinaryReader::CCMBinaryReader(const I3Context& context) : I3Module(context),
    nframes_(0), num_frames_in_file_(0), new_config_(false) {
    std::string fname;

    AddParameter("Filename",
            "Filename to read.  Use either this or Filenamelist, not both.",
            fname);

    AddParameter("FilenameList",
            "List of files to read, *IN SORTED ORDER*",
            std::vector<std::string>());

    AddParameter("DropDuplicateConfigs",
            "Drop frames that duplicate an existing configuration. "
            "Saves space when processing files from the same run.",
            true);

    AddOutBox("OutBox");
}

    void
CCMBinaryReader::Configure() {
    std::string fname;

    GetParameter("FileNameList", filenames_);
    GetParameter("FileName", fname);
    if(!filenames_.empty() && !fname.empty())
        log_fatal("Both Filename and FileNameList were specified.");

    if (filenames_.empty() && fname.empty())
        log_fatal("Neither 'Filename' nor 'FilenameList' specified");

    if (filenames_.empty())
        filenames_.push_back(fname);

    GetParameter("DropDuplicateConfigs", drop_duplicate_configs_);

    last_config_ = boost::make_shared<CCMAnalysis::Binary::CCMDAQConfig>();
    readout_ = boost::make_shared<CCMAnalysis::Binary::CCMTriggerReadout>();

    file_stager_ = context_.Get<I3FileStagerPtr>();
    if (!file_stager_)
        file_stager_ = I3TrivialFileStager::create();
    BOOST_FOREACH(const std::string &filename, filenames_)
        file_stager_->WillReadLater(filename);

    filenames_iter_ = filenames_.begin();
    OpenNextFile();
}

    void
CCMBinaryReader::Process() {
    while(file_stream_.peek() == EOF or num_frames_in_file_ == nframes_) {
        if(filenames_iter_ == filenames_.end()) {
            RequestSuspension();
            current_filename_.reset();
            return;
        }
        else
            OpenNextFile();
    }

    I3FramePtr frame(new I3Frame);
    try {
        nframes_++;
        if(new_config_) {
            frame->SetStop(I3Frame::Geometry);
            frame->Put("CCMDAQConfig", last_config_);
            new_config_ = false;
        } else {
            CCMAnalysis::Binary::read_binary(file_stream_, *readout_);
            frame->SetStop(I3Frame::DAQ);
            frame->Put("CCMTriggerReadout", readout_);
        }
    } catch (const std::exception &e) {
        log_fatal("Error reading %s at frame %d: %s!",
                current_filename_->c_str(), nframes_, e.what());
        return;
    }

    PushFrame(frame, "OutBox");
}

    void
CCMBinaryReader::OpenNextFile() {
    current_filename_.reset();
    current_filename_ = file_stager_->GetReadablePath(*filenames_iter_);
    nframes_ = 0;
    num_frames_in_file_ = 0;
    new_config_ = false;
    filenames_iter_++;

    I3::dataio::open(ifs_, *current_filename_);
    file_stream_ = std::ifstream(*current_filename_, std::ios::binary);
    log_trace("Constructing with filename %s",
            current_filename_->c_str());

    log_info("Opened file %s", current_filename_->c_str());
    boost::shared_ptr<CCMAnalysis::Binary::CCMDAQConfig> config = boost::make_shared<CCMAnalysis::Binary::CCMDAQConfig>();
    CCMAnalysis::Binary::read_config(file_stream_, *config);
    CCMAnalysis::Binary::read_size(file_stream_, num_frames_in_file_);
    if(drop_duplicate_configs_ and *config == *last_config_) {
        new_config_ = false;
        return;
    } else {
        num_frames_in_file_ += 1;
        new_config_ = true;
    }
}

CCMBinaryReader::~CCMBinaryReader() {
}

