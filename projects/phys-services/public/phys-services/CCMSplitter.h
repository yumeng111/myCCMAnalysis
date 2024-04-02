/**
 *  Copyright (C) 2011
 *  Nathan Whitehorn	<nwhitehorn@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */
#ifndef ICETRAY_CCMSPLITTERMODULE_H_INCLUDED
#define ICETRAY_CCMSPLITTERMODULE_H_INCLUDED

#include <icetray/I3Configuration.h>
#include <icetray/I3Frame.h>
#include <string>

/**
 * @brief This class is meant to be a mix-in base for modules that split
 * DAQ frames into multiple physics frames.
 *
 */

class CCMSplitter {
public:
    CCMSplitter(const I3Configuration& config);
    ~CCMSplitter();

    I3FramePtr GetNextSubEvent(I3FramePtr daq);

private:
    I3FramePtr last_daq;
    const I3Configuration& config_;
    int last_subevent_id;
    bool print_warning_;

protected:
    /**
     * The name of the SubEvent stream.  This should be set by
     * derived modules, which is why it's protected.
     * 
     * Modules where this is used as a mix-in simply need to add
     * the following lines to the constructor and Configure method, 
     * respectively:
     * 
     *   AddParameter("SubEventStreamName", "The name of the SubEvent stream.", "SplitSubStreamDefault");
     *   GetParameter("SubEventStreamName", sub_event_stream_name_);
     */
    std::string sub_event_stream_name_;
};

#endif // ICETRAY_CCMSPLITTERMODULE_H_INCLUDED
