#
#  $Id$
#  
#  Copyright (C) 2007   Troy D. Straszheim  <troy@icecube.umd.edu>
#  and the IceCube Collaboration <http://www.icecube.wisc.edu>
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions
#  are met:
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#  
#  THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
#  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
#  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
#  SUCH DAMAGE.
#  
#  SPDX-License-Identifier: BSD-2-Clause
#  
#  

# Find (serial or parallel) HDF5
TOOLDEF (hdf5
  include
  hdf5.h
  lib
  NONE
  hdf5 hdf5_hl
  )

if(NOT HDF5_FOUND)
  unset(HDF5_CONFIG_ERROR)
  unset(HDF5_INCLUDE_DIR)

  # serial HDF5 fallback
  TOOLDEF (hdf5
    /usr/include/hdf5/serial
    hdf5.h
    lib/x86_64-linux-gnu
    NONE
    hdf5_serial hdf5_serial_hl
    )

  if(HDF5_FOUND)
    set(HDF5_FOUND TRUE CACHE BOOL "Tool HDF5 found successfully")
    set(HDF5_LIBRARIES "${HDF5_LIBRARIES}" CACHE PATH "Libraries for tool HDF5")
  endif()
endif()

# Check whether hdf5 is a parallel build
#    Parallel HDF5 defines H5_HAVE_PARALLEL in H5pubconf.h and H5public.h includes mpi.h.
if(HDF5_FOUND)
  # Try to locate H5pubconf.h next to hdf5.h
  set(_H5CONF_CANDIDATES
      "${HDF5_INCLUDE_DIR}/H5pubconf.h"
      "${HDF5_INCLUDE_DIR}/../H5pubconf.h"
  )
  set(_HDF5_PARALLEL FALSE)
  foreach(_cand IN LISTS _H5CONF_CANDIDATES)
    if(EXISTS "${_cand}")
      file(STRINGS "${_cand}" _H5CONF_PAR_LINE REGEX "^#define[ \t]+H5_HAVE_PARALLEL[ \t]+1")
      if(_H5CONF_PAR_LINE)
        set(_HDF5_PARALLEL TRUE)
        break()
      endif()
    endif()
  endforeach()

  # If parallel, bring in MPI and append its include/lib paths to HDF5_* variables
  if(_HDF5_PARALLEL)
    find_package(MPI REQUIRED COMPONENTS C CXX)

    # Make sure these variables exist even if TOOLDEF left them as CACHE PATH, etc.
    if(NOT DEFINED HDF5_INCLUDE_DIR)
      set(HDF5_INCLUDE_DIR "")
    endif()
    if(NOT DEFINED HDF5_LIBRARIES)
      set(HDF5_LIBRARIES "")
    endif()

    # Append MPI includes/libs so any target linking HDF5 also gets MPI
    list(APPEND HDF5_INCLUDE_DIR ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
    list(APPEND HDF5_LIBRARIES ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

    # expose a boolean for debugging
    set(HDF5_IS_PARALLEL TRUE CACHE BOOL "HDF5 was detected as a parallel build")
  else()
    set(HDF5_IS_PARALLEL FALSE CACHE BOOL "HDF5 was detected as a serial build")
  endif()
endif()

