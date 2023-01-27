
#include <boost/foreach.hpp>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <icetray/python/dataclass_suite.hpp>
#include <icetray/python/function.hpp>

namespace bp = boost::python;

CCMRecoPulseSeriesMapPtr
underhanded_apply(const CCMRecoPulseSeriesMapMask& mask, const I3Frame &frame)
{
	return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(mask.Apply(frame));
}

bp::list
getbits(const CCMRecoPulseSeriesMapMask &mask)
{
	bp::list usermask;
	BOOST_FOREACH(const boost::dynamic_bitset<uint8_t> &bits, mask.GetBits()) {
		bp::list elements;
		for (size_t i = 0; i < bits.size(); i++)
			elements.append(bits.test(i));
		usermask.append(elements);
	}
	
	return usermask;
}

namespace {
CCMRecoPulseSeriesMapMaskPtr
from_callable(const I3Frame &frame, const std::string &key, bp::object callable)
{
	typedef boost::function<bool (const CCMPMTKey&, size_t, const CCMRecoPulse&)> callback_t;
	boost::shared_ptr<callback_t> predicate = bp::detail::function_from_object<callback_t>(callable);
	return CCMRecoPulseSeriesMapMaskPtr(new CCMRecoPulseSeriesMapMask(frame, key, *predicate));
}
};

void register_CCMRecoPulseSeriesMapMask()
{
	void (CCMRecoPulseSeriesMapMask::*set_om_all)(const CCMPMTKey&, bool) = &CCMRecoPulseSeriesMapMask::Set;
	void (CCMRecoPulseSeriesMapMask::*set_om_by_idx)(const CCMPMTKey&, const unsigned, bool) = &CCMRecoPulseSeriesMapMask::Set;
	void (CCMRecoPulseSeriesMapMask::*set_om_by_value)(const CCMPMTKey&, const CCMRecoPulse&, bool) = &CCMRecoPulseSeriesMapMask::Set;
	
	typedef boost::function<bool (const CCMPMTKey&, size_t, const CCMRecoPulse&)> callback_t;
	
	bp::scope mask_scope =
	bp::class_<CCMRecoPulseSeriesMapMask, bp::bases<I3FrameObject>,
	    CCMRecoPulseSeriesMapMaskPtr>("CCMRecoPulseSeriesMapMask",
	    bp::init<const I3Frame&, const std::string &>(bp::args("frame", "key")))
		.def("__init__", bp::make_constructor(&from_callable))
		.def(bp::init<const I3Frame&, const std::string &, const CCMRecoPulseSeriesMap &>())
		.def(bp::init<const I3Frame&, const std::string &, callback_t>())
		.add_property("source", &CCMRecoPulseSeriesMapMask::GetSource)
		.add_property("bits", &getbits)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
		.def(bp::self &  bp::self)
		.def(bp::self &= bp::self)
		.def(bp::self |  bp::self)
		.def(bp::self |= bp::self)
		.def(bp::self ^  bp::self)
		.def(bp::self ^= bp::self)
#ifdef __clang__
#pragma GCC diagnostic pop
#endif              
		.def("remove", &CCMRecoPulseSeriesMapMask::Remove)
		.def("has_ancestor", &CCMRecoPulseSeriesMapMask::HasAncestor)
		.def("repoint", &CCMRecoPulseSeriesMapMask::Repoint)
		.def("apply", &underhanded_apply, "Apply the mask to an I3Frame, returning an CCMRecoPulseSeries.")
		.def("any", &CCMRecoPulseSeriesMapMask::GetAnySet, "Are any of the bits set in the mask?")
		.def("all", &CCMRecoPulseSeriesMapMask::GetAllSet, "Are all of the bits set in the mask?")
		.def("sum", &CCMRecoPulseSeriesMapMask::GetSum, "Get the number of set bits in the mask.")
		.def("set", set_om_all, "Set/unset all bits for this CCMPMTKey.")
		.def("set", set_om_by_idx, "Set/unset the bit at this index.")
		.def("set", set_om_by_value, "Set/unset the bit corresponding to this RecoPulse.")
		.def("set_none", &CCMRecoPulseSeriesMapMask::SetNone, "Unset all bits in the mask.")
		.def(bp::dataclass_suite<CCMRecoPulseSeriesMapMask>())
	;
	
	bp::def_function<callback_t>("Predicate");
	
	register_pointer_conversions<CCMRecoPulseSeriesMapMask>();
}
