/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Eike Middell <eike.middell@desy.de> $LastChangedBy$
 */

#include "tableio/converter/CCMEventHeaderConverter.h"

I3TableRowDescriptionPtr CCMEventHeaderConverter::CreateDescription(const CCMEventHeader& header) {
    
    I3TableRowDescriptionPtr desc(new I3TableRowDescription());
    desc->AddField<int64_t>("time_start_utc_nsec", "1e-9 s", "start of event - nanoseconds from the previous second");
    desc->AddField<int64_t> ("time_start_utc_sec",     "s",    "start of event - seconds from the Unix epoch");
    desc->AddField<int64_t>("time_end_utc_nsec",   "1e-9 s", "end of event in daq time - nanoseconds from the previous second");
    desc->AddField<int64_t> ("time_end_utc_sec",       "s",    "end of event - seconds from the Unix epoch");
    return desc;
}
    
size_t CCMEventHeaderConverter::FillRows(const CCMEventHeader& header,
                                        I3TableRowPtr rows) {
    CCMTime start = header.GetStartTime();
    CCMTime end = header.GetEndTime();

    rows->Set<int64_t>("time_start_utc_nsec", start.GetSeconds());
    rows->Set<int64_t>("time_start_utc_sec",  start.GetNanoSeconds());
    
    rows->Set<int64_t>("time_end_utc_nsec", end.GetSeconds());
    rows->Set<int64_t>("time_end_utc_sec",  end.GetNanoSeconds());

    return 1;
}
