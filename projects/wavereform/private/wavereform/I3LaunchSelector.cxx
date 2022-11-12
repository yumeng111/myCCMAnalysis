
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <icetray/I3ConditionalModule.h>
#include <dataclasses/physics/I3DOMLaunch.h>
#include <recclasses/I3CompactKeyList.h>

class I3LaunchSelector : public I3ConditionalModule {
public:
	I3LaunchSelector(const I3Context&);
	void Configure();
	void DAQ(I3FramePtr);
private:
	std::vector<std::string> flag_names_;
	std::string raw_launch_name_, selected_launch_name_;
	bool invert_;
};

I3_MODULE(I3LaunchSelector);

I3LaunchSelector::I3LaunchSelector(const I3Context &ctx) : I3ConditionalModule(ctx)
{
	raw_launch_name_ = "InIceRawData";
	AddParameter("Launches", "Name of DOMLaunch map in the frame.", raw_launch_name_);
	
	selected_launch_name_ = "InIceRawDataErrata";
	AddParameter("Output", "Name of DOMLaunch map to put in the frame. "
	    "The list of keys is also written to the frame as output+'Keys', "
	    "e.g. InIceRawDataErrataKeys", selected_launch_name_);
	
	AddParameter("Flags", "List of vector<OMKey>s in the frame. "
	    "Any OMKey that appears in one of these lists will be copied to the output map.",
	    flag_names_);

	invert_ = false;
	AddParameter("Invert", "Invert the selection: any key that does *not* appear in "
	    "Flags will be copied to the output map", invert_);
	
	AddOutBox("OutBox");
}

void
I3LaunchSelector::Configure()
{
	GetParameter("Launches", raw_launch_name_);
	GetParameter("Output", selected_launch_name_);
	GetParameter("Flags", flag_names_);
	GetParameter("Invert", invert_);
}

void
I3LaunchSelector::DAQ(I3FramePtr frame)
{
	/* OR the flag lists together */
	std::set<OMKey> flags;
	BOOST_FOREACH(const std::string &name, flag_names_) {
		boost::shared_ptr<const std::vector<OMKey> > flag_list =
		    frame->Get<boost::shared_ptr<const std::vector<OMKey> > >(name);
		if (!flag_list)
			continue;
		BOOST_FOREACH(const OMKey &key, *flag_list) {
			flags.insert(key);
		}
	}
	
	I3DOMLaunchSeriesMapConstPtr launches =
	    frame->Get<I3DOMLaunchSeriesMapConstPtr>(raw_launch_name_);
	
	if (!launches) {
		PushFrame(frame);
		return;
	}
	
	I3DOMLaunchSeriesMapPtr selection;
	I3CompactKeyListPtr keys = boost::make_shared<I3CompactKeyList>();
	
	if (invert_) {
		selection = boost::make_shared<I3DOMLaunchSeriesMap>(*launches);
		/* Remove the flagged keys, and fill a list to put in the frame. */
		BOOST_FOREACH(const OMKey &key, flags) {
			I3DOMLaunchSeriesMap::iterator target = selection->find(key);
			if (target != selection->end())
				selection->erase(target);
		}
		BOOST_FOREACH(const I3DOMLaunchSeriesMap::value_type &pair, *selection)
			keys->push_back(pair.first);
	} else {
		selection = boost::make_shared<I3DOMLaunchSeriesMap>();
		/* Copy the flagged keys, and fill a list to put in the frame. */
		BOOST_FOREACH(const OMKey &key, flags) {
			I3DOMLaunchSeriesMap::const_iterator target = launches->find(key);
			if (target != launches->end()) {
				selection->insert(*target);
				keys->push_back(key);
			}
		}
	}
	
	if (selection->size() > 0) {
		frame->Put(selected_launch_name_, selection);
		frame->Put(selected_launch_name_ + "Keys", keys);
	}
		
	PushFrame(frame);
}
