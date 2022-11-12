
#include <icetray/load_project.h>
#include <icetray/I3Units.h>
#include <wavereform/I3Wavereform.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <dataclasses/calibration/I3Calibration.h>
#include <dataclasses/status/I3DetectorStatus.h>

#include <numpy/numpyconfig.h>
#ifdef NPY_1_7_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include <numpy/ndarrayobject.h>

namespace bp = boost::python;

static PyObject *hack_import_array() {import_array(); return NULL;}

static bp::object
RefoldPulses(const std::vector<I3RecoPulse> &pulses, I3Waveform::Source source,
    unsigned channel, const I3DOMCalibration &calib, const I3DOMStatus &status,
    bp::object times, bool is_simulation = false)
{
	PyArrayObject *time_arr = (PyArrayObject*)PyArray_FromObject(times.ptr(), NPY_DOUBLE, 1, 1);
	if (!time_arr) {
		PyErr_SetString(PyExc_ValueError, "Can't turn 'times' into an array!");
		throw bp::error_already_set();
	}
	npy_intp n_bins = PyArray_DIM(time_arr, 0);
	
	PyArrayObject *bin_arr = (PyArrayObject*)PyArray_SimpleNew(1, &n_bins, NPY_DOUBLE);
	if (!bin_arr) {
		PyErr_SetString(PyExc_MemoryError, "Can't allocate bin array!");
		throw bp::error_already_set();
	}
	
	I3Wavereform::RefoldPulses(pulses, source, channel, calib, status,
	    (double*)PyArray_DATA(time_arr), (double*)PyArray_DATA(bin_arr),
	    n_bins, is_simulation);
	
	Py_XDECREF(time_arr);
	
	return bp::object(bp::handle<>((PyObject*)bin_arr));
}

BOOST_PYTHON_FUNCTION_OVERLOADS(RefoldPulsesOverloads, RefoldPulses, 6, 7);

static double
chisquared(const std::vector<I3Wavereform::Refolded> &bins)
{
	double chi = 0.0;
	BOOST_FOREACH(const I3Wavereform::Refolded &bin, bins) {
		double diff = (bin.waveform-bin.refolded)/bin.step;
		chi += diff*diff;
	}
	
	return chi;
}

static bp::tuple
GetChiSquared(const std::vector<I3RecoPulse> &pulses, const I3Waveform &wf,
    const I3DOMCalibration &calib, const I3DOMStatus &status)
{
	std::vector<I3Wavereform::Refolded> refolded
	    = I3Wavereform::GetRefolded(pulses, wf, calib, status);
	return bp::make_tuple(chisquared(refolded), refolded.size());
}

I3_PYTHON_MODULE(wavereform)
{
	load_project("wavereform", false);
	hack_import_array();
	
	bp::def("refold_pulses", RefoldPulses, RefoldPulsesOverloads(
	    bp::args("pulses", "source", "channel", "calibration", "status",
	    "times", "is_simulation")));
	bp::def("chi_squared", &GetChiSquared);
	bp::def("get_max_charge", &I3Wavereform::GetMaxCharge,
	    (bp::arg("pulses"), bp::arg("twindow")=6.4*I3Units::microsecond));
}
