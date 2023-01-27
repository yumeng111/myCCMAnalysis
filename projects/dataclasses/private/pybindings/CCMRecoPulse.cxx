//
//   Copyright (c) 2004, 2005, 2006, 2007   Troy D. Straszheim  
//   
//   $Id$
//
//   This file is part of IceTray.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//   
//   THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
//   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//   ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
//   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
//   OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
//   OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
//   SUCH DAMAGE.
//   
//   SPDX-License-Identifier: BSD-2-Clause
//   
//

#include <vector>

#include <dataclasses/physics/CCMRecoPulse.h>
#include <icetray/I3Frame.h>
#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

static CCMRecoPulseSeriesMapPtr
from_frame(I3Frame &frame, const std::string &name)
{
	if (!frame.Has(name)) {
		PyErr_SetString(PyExc_KeyError, name.c_str());
		throw_error_already_set();
	}
	CCMRecoPulseSeriesMapConstPtr rpsm =
	    frame.Get<CCMRecoPulseSeriesMapConstPtr>(name);
	if (!rpsm) {
		PyErr_SetString(PyExc_TypeError, name.c_str());
		throw_error_already_set();
	}
	return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(rpsm);
}

namespace {
static PyBufferProcs rps_bufferprocs;
static PyBufferProcs rpsm_bufferprocs;

static int
CCMRecoPulseSeries_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;
	view->internal = NULL;
	view->suboffsets = NULL;
	view->buf = NULL;

	handle<> self(boost::python::borrowed(obj));
	object selfobj(self);
	std::vector<CCMRecoPulse> &ts =
	    extract<std::vector<CCMRecoPulse> &>(selfobj)();
	if (ts.size() == 0) {
		PyErr_SetString(PyExc_BufferError, "Pulse series is empty.");
		view->obj = NULL;
		return -1;
	}
	if ((flags & PyBUF_WRITABLE) == PyBUF_WRITABLE) {
		PyErr_SetString(PyExc_BufferError, "Cannot provide writable "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}
	if ((flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS) {
		PyErr_SetString(PyExc_BufferError, "Cannot provide FORTRAN "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}

	view->obj = obj;
	view->len = 3*ts.size()*sizeof(double);
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;
	view->ndim = 2;

	view->shape = new Py_ssize_t[2];
	view->shape[0] = ts.size();
	view->shape[1] = 3;

	// To honor a contiguous buffer request, make a copy of the
	// data. This violates the spirit of the buffer protocol
	// slightly, but both simplifies the API and allows a
	// faster memory copy than iterating over the list in Python
	// would.

	view->buf = malloc(view->len);
	view->readonly = 1;

	view->strides = new Py_ssize_t[2];
	view->strides[0] = 3*view->itemsize;
	view->strides[1] = view->itemsize;

	int j = 0;
	for (auto i : ts) {
		((double *)view->buf)[3*j + 0] = i.GetTime();
		((double *)view->buf)[3*j + 1] = i.GetCharge();
		((double *)view->buf)[3*j + 2] = i.GetWidth();
		j++;
	}

	view->suboffsets = NULL;
	view->internal = view->buf;

        Py_INCREF(obj);

        return 0;
}

static void
CCMRecoPulseSeries_relbuffer(PyObject *obj, Py_buffer *view)
{
	if (view->strides != NULL)
		delete [] view->strides;
	if (view->shape != NULL)
		delete [] view->shape;
	if (view->suboffsets != NULL)
		delete [] view->suboffsets;
	if (view->buf != NULL)
		free(view->buf);
}

static int
CCMRecoPulseSeriesMap_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;
	view->internal = NULL;
	view->suboffsets = NULL;
	view->buf = NULL;

	handle<> self(boost::python::borrowed(obj));
	object selfobj(self);
	CCMRecoPulseSeriesMap &ts = extract<CCMRecoPulseSeriesMap &>(selfobj)();
	if (ts.size() == 0) {
		PyErr_SetString(PyExc_BufferError, "Pulse series is empty.");
		view->obj = NULL;
		return -1;
	}
	if ((flags & PyBUF_WRITABLE) == PyBUF_WRITABLE) {
		PyErr_SetString(PyExc_BufferError, "Cannot provide writable "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}
	if ((flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS) {
		PyErr_SetString(PyExc_BufferError, "Cannot provide FORTRAN "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}

	view->len = 0;
	for (auto i : ts)
		view->len += 6*i.second.size()*sizeof(double);

	view->obj = obj;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;
	view->ndim = 2;

	view->shape = new Py_ssize_t[2];
	view->shape[1] = 6;
	view->shape[0] = view->len/view->shape[1]/sizeof(double);

	// To honor a contiguous buffer request, make a copy of the
	// data. This violates the spirit of the buffer protocol
	// slightly, but both simplifies the API and allows a
	// faster memory copy than iterating over the list in Python
	// would.

	view->buf = malloc(view->len);
	view->readonly = 1;

	view->strides = new Py_ssize_t[2];
	view->strides[0] = 6*view->itemsize;
	view->strides[1] = view->itemsize;

	int j = 0;
	for (auto i : ts) {
		for (auto k : i.second) {
			((double *)view->buf)[6*j + 0] = i.first.GetRegion();
			((double *)view->buf)[6*j + 1] = i.first.GetSensor();
			((double *)view->buf)[6*j + 2] = i.first.GetSubsensor();
			((double *)view->buf)[6*j + 3] = k.GetTime();
			((double *)view->buf)[6*j + 4] = k.GetCharge();
			((double *)view->buf)[6*j + 5] = k.GetWidth();
			j++;
		}
	}

	view->suboffsets = NULL;
	view->internal = view->buf;

        Py_INCREF(obj);

        return 0;
}

static list
CCMRecoPulseSeriesMap_pmtoffsets(const CCMRecoPulseSeriesMap &rpsm)
{
	list out;

	size_t j = 0;
	for (auto i : rpsm) {
		out.append(j);
		j += i.second.size();
	}

	return out;
}
}


void register_CCMRecoPulse()
{
  object rps = class_<std::vector<CCMRecoPulse> >("vector_CCMRecoPulse")
    .def(dataclass_suite<std::vector<CCMRecoPulse> >())
    ;

  // Add buffer protocol interface
  PyTypeObject *rpsclass = (PyTypeObject *)rps.ptr();
  rps_bufferprocs.bf_getbuffer = CCMRecoPulseSeries_getbuffer;
  rps_bufferprocs.bf_releasebuffer = CCMRecoPulseSeries_relbuffer;
  rpsclass->tp_as_buffer = &rps_bufferprocs;

  object rpsm = class_<CCMRecoPulseSeriesMap, bases<I3FrameObject>, CCMRecoPulseSeriesMapPtr>("CCMRecoPulseSeriesMap")
    .def(dataclass_suite<CCMRecoPulseSeriesMap>())
    .def("from_frame", &from_frame, args("frame", "key"),
        "Get an CCMRecoPulseSeriesMap from the frame, performing any necessary "
        "format conversions behind the scenes.")
    .staticmethod("from_frame")
    .def("pmt_array_offsets", &CCMRecoPulseSeriesMap_pmtoffsets,
        "Provide a list of offsets into a numpy.asarray()-ed "
        "CCMRecoPulseSeriesMap corresponding to the beginning of the pulses for "
        "each OMKey. The format of the numpy.asarray() version of an "
        "CCMRecoPulseSeriesMap is one row per pulse, PMTs grouped together, "
        "with columns (String, OM, PMT, Time, Charge, Width).")
    ;
  register_pointer_conversions<CCMRecoPulseSeriesMap>();

  // Add buffer protocol interface
  PyTypeObject *rpsmclass = (PyTypeObject *)rpsm.ptr();
  rpsm_bufferprocs.bf_getbuffer = CCMRecoPulseSeriesMap_getbuffer;
  rpsm_bufferprocs.bf_releasebuffer = CCMRecoPulseSeries_relbuffer;
  rpsmclass->tp_as_buffer = &rpsm_bufferprocs;

  scope outer = 
  class_<CCMRecoPulse, boost::shared_ptr<CCMRecoPulse> >("CCMRecoPulse")
    #define PROPS (Time)(Charge)(Width)
    BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, CCMRecoPulse, PROPS)
    #undef PROPS
    .def(dataclass_suite<CCMRecoPulse>())
    ;
}
