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

#include <daqtools/PhotonPropagation.h>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/indexed_property.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>
#include <icetray/python/stream_to_string.hpp>

using namespace boost::python;

void register_PhotonPropagation() {
{

    class_<PhotonPropagation, boost::noncopyable>("PhotonPropagation")
    //class_<boost::shared_ptr<PhotonPropagation>> ("PhotonPropagation")
        .def(init<>())
        .def("SetData", &PhotonPropagation::SetData)
        .def("SetDataSampleSize", &PhotonPropagation::SetDataSampleSize)
        .def("SetNThreads", &PhotonPropagation::SetNThreads)
        .def("GetNFaceChunks", &PhotonPropagation::GetNFaceChunks)
        .def("GetNSideChunks", &PhotonPropagation::GetNSideChunks)
        .def("SetZOffset", &PhotonPropagation::SetZOffset)
        .def("GetNSimulatedEvents", &PhotonPropagation::GetNSimulatedEvents)
        .def("GetEventVertices", &PhotonPropagation::GetEventVertices)
        .def("GetPMTInformation", &PhotonPropagation::GetPMTInformation)
        .def("GetSecondaryLocs", &PhotonPropagation::GetSecondaryLocs)
        .def("GetTopCoatedPMTLocs", &PhotonPropagation::GetTopCoatedPMTLocs)
        .def("GetTopUncoatedPMTLocs", &PhotonPropagation::GetTopUncoatedPMTLocs)
        .def("GetBottomCoatedPMTLocs", &PhotonPropagation::GetBottomCoatedPMTLocs)
        .def("GetBottomUncoatedPMTLocs", &PhotonPropagation::GetBottomUncoatedPMTLocs)
        .def("GetSideCoatedPMTLocs", &PhotonPropagation::GetSideCoatedPMTLocs)
        .def("GetSideUncoatedPMTLocs", &PhotonPropagation::GetSideUncoatedPMTLocs)
        .def("GetTopLocWidth", &PhotonPropagation::GetTopLocWidth)
        .def("GetTopLocHeight", &PhotonPropagation::GetTopLocHeight)
        .def("GetTopLocXY", &PhotonPropagation::GetTopLocXY)
        .def("GetTopLocInsideCoatedPMTXY", &PhotonPropagation::GetTopLocInsideCoatedPMTXY)
        .def("GetTopLocOutsideCoatedPMTXY", &PhotonPropagation::GetTopLocOutsideCoatedPMTXY)
        .def("GetTopLocInsideDetectorRadiusXY", &PhotonPropagation::GetTopLocInsideDetectorRadiusXY)
        .def("GetTopLocOutsideDetectorRadiusXY", &PhotonPropagation::GetTopLocOutsideDetectorRadiusXY)
        .def("GetAllTopLocValidDots", &PhotonPropagation::GetAllTopLocValidDots)
        .def("GetAllTopLocInvalidDots", &PhotonPropagation::GetAllTopLocInvalidDots)
        .def("GetAllTopLocCoatedPMTDots", &PhotonPropagation::GetAllTopLocCoatedPMTDots)
        .def("GetAllBottomLocValidDots", &PhotonPropagation::GetAllBottomLocValidDots)
        .def("GetAllBottomLocInvalidDots", &PhotonPropagation::GetAllBottomLocInvalidDots)
        .def("GetAllBottomLocCoatedPMTDots", &PhotonPropagation::GetAllBottomLocCoatedPMTDots)
        .def("GetAllSideLocValidDots", &PhotonPropagation::GetAllSideLocValidDots)
        .def("GetAllSideLocInvalidDots", &PhotonPropagation::GetAllSideLocInvalidDots)
        .def("GetAllSideLocCoatedPMTDots", &PhotonPropagation::GetAllSideLocCoatedPMTDots)
        .def("GetSimulation", &PhotonPropagation::GetSimulation);

}
}
