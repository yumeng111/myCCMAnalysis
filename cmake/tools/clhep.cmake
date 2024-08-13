#
#  cmake/tools/clhep.cmake
#
#  Copyright (C) 2010
#  the IceCube Collaboration <http://www.icecube.wisc.edu>
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

# - Try to find CLHEP
# Once done this will define
#
#  CLHEP_FOUND - system has CLHEP
#  CLHEP_INCLUDE_DIR - the CLHEP include directory
#  CLHEP_LIBRARIES - The libraries needed to use CLHEP
#  CLHEP_DEFINITIONS - Compiler switches required for using CLHEP
#

MESSAGE(STATUS "    Looking for CLHEP")

FIND_PATH(CLHEP_INCLUDE_DIR NAMES CLHEP
  PATHS ${CLHEP_BASE_DIR}/include
  NO_DEFAULT_PATH
)

IF(NOT ${CLHEP_INCLUDE_DIR} MATCHES ".*NOTFOUND$")

  FIND_PATH(CLHEP_LIBRARY_DIR NAMES libCLHEP.so
    PATHS ${CLHEP_BASE_DIR}/lib
    NO_DEFAULT_PATH
    )
  IF(NOT ${CLHEP_LIBRARY_DIR} MATCHES ".*NOTFOUND$")

    IF(CLHEP_BASE_DIR AND CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY_DIR)
      SET(CLHEP_FOUND TRUE)
    ENDIF(CLHEP_BASE_DIR AND CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY_DIR)

    IF(CLHEP_FOUND)
      MESSAGE(STATUS "    + found at ${CLHEP_BASE_DIR}")
      INCLUDE_DIRECTORIES(${CLHEP_INCLUDE_DIR})
      LINK_DIRECTORIES(${CLHEP_LIBRARY_DIR})
      LINK_LIBRARIES("CLHEP")
    ENDIF(CLHEP_FOUND)
  ENDIF(NOT ${CLHEP_LIBRARY_DIR} MATCHES ".*NOTFOUND$")
ENDIF(NOT ${CLHEP_INCLUDE_DIR} MATCHES ".*NOTFOUND$")

IF (CLHEP_FOUND)
  MESSAGE(STATUS "    Looking for CLHEP -- found")
ELSE (CLHEP_FOUND)
  MESSAGE(STATUS "    Looking for CLHEP -- not found")
ENDIF (CLHEP_FOUND)

tooldef (clhep
    include
    CLHEP/ClhepVersion.h
    lib
    NONE  # bin is n/a, placeholder
    CLHEP
    )
