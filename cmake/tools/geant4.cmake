#
#  cmake/tools/geant4.cmake
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

colormsg(HICYAN "")
colormsg(HICYAN "geant4")

find_package(Geant4 REQUIRED OPTIONAL_COMPONENTS ui_all vis_all)
if(GEANT4_FOUND)
    message(STATUS "Geant4 found")

    set(GEANT4_LIBRARIES ${Geant4_LIBRARIES})
    set(GEANT4_INCLUDE_DIR ${Geant4_INCLUDE_DIRS})
    set(GEANT4_LIB_DIR ${Geant4_ROOT_DIR}/lib)
    set(GEANT4_BIN_DIR ${Geant4_ROOT_DIR}/bin)

    execute_process (COMMAND
        cat ${GEANT4_INCLUDE_DIR}/G4Version.hh
        OUTPUT_VARIABLE VERSION_OUTPUT
    )

    # weird cmake way of extracting only a subexpression
    string (REGEX REPLACE
        ".*#define +G4VERSION_NUMBER +([0-9]+).*"
        "\\1"
        GEANT4_VERSION ${VERSION_OUTPUT}
    )

    found_ok("     version: ${GEANT4_VERSION}")
    found_ok("    includes: ${GEANT4_INCLUDE_DIR}")
    found_ok("        libs: ${GEANT4_LIBRARIES}")
    found_ok("c++ standard: ${Geant4_CXX_STANDARD}")
    found_ok(" definitions: ${Geant4_DEFINITIONS}")
    found_ok("    datasets: ${Geant4_DATASETS}")

    unset(CLHEP_FOUND CACHE)
    unset(CLHEP_INCLUDE_DIR CACHE)
    unset(CLHEP_FOUND)
    unset(CLHEP_INCLUDE_DIR)

    tooldef (clhep
      ${GEANT4_INCLUDE_DIR}
      CLHEP/Units/SystemOfUnits.h
      ${GEANT4_LIB_DIR}
      NONE  # bin is n/a, placeholder
      )

    # no extra libraries for the geant4-internal version of clhep
    unset(CLHEP_LIBRARIES CACHE)
    unset(CLHEP_LIBRARIES)

else()
    message(STATUS "Geant4 not found")
endif()
