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

# we don't use the tooldef() macro, so we have to fudge pretty-printing
colormsg("")
colormsg(HICYAN "ncurses (debuggable)")

# First, try to get WIDE ncurses
set(CURSES_NEED_NCURSES TRUE)
set(CURSES_NEED_WIDE    TRUE)
set(CURSES_USE_NCURSES  TRUE)
find_package(Curses)

# Helper to tell if the result is actually "wide"
function(_curses_result_is_wide out)
  set(_iswide FALSE)
  foreach(_lib IN LISTS CURSES_LIBRARIES)
    if(_lib MATCHES "ncursesw(|\\.so|\\.a|[-\\.]lib)($|[^a-zA-Z])")
      set(_iswide TRUE)
      break()
    endif()
  endforeach()
  set(${out} ${_iswide} PARENT_SCOPE)
endfunction()

if(CURSES_FOUND)
  _curses_result_is_wide(_got_wide)
else()
  set(_got_wide FALSE)
endif()

if(NOT _got_wide)
  # -> purge Curses cache vars and retry narrow
  if(COMMAND colormsg)
    colormsg(YELLOW "[ncurses] Wide stack not available; retrying narrow find_package(Curses).")
  else()
    message(STATUS "[ncurses] Wide stack not available; retrying narrow find_package(Curses).")
  endif()

  unset(CURSES_FOUND           CACHE)
  unset(CURSES_INCLUDE_DIR     CACHE)
  unset(CURSES_LIBRARY         CACHE)
  unset(CURSES_LIBRARIES       CACHE)
  unset(CURSES_FORM_LIBRARY    CACHE)
  unset(CURSES_MENU_LIBRARY    CACHE)
  unset(CURSES_PANEL_LIBRARY   CACHE)

  # Ask for ncurses, but not wide this time
  set(CURSES_NEED_NCURSES TRUE)
  set(CURSES_NEED_WIDE    FALSE)
  set(CURSES_USE_NCURSES  TRUE)
  find_package(Curses REQUIRED)
endif()

# Helpers
function(_dir_of out path)
  get_filename_component(_d "${path}" DIRECTORY)
  set(${out} "${_d}" PARENT_SCOPE)
endfunction()

function(_fatal_if_missing what var)
  if(NOT ${var})
    message(FATAL_ERROR "ncurses.cmake: missing ${what}")
  endif()
endfunction()

function(_require_same_prefix prefix libs_out)
  # Ensure all libs live in the same directory; return that list
  set(_first "")
  set(_accum "")
  foreach(_lib IN LISTS ARGN)
    if(NOT _lib)
      continue()
    endif()
    get_filename_component(_d "${_lib}" DIRECTORY)
    if(_first STREQUAL "")
      set(_first "${_d}")
      if(NCURSES_DEBUG)
        message(STATUS "[ncurses] prefix candidate: ${_first}")
      endif()
    elseif(NOT _d STREQUAL _first)
      message(FATAL_ERROR
        "ncurses.cmake: mixed prefixes detected:\n  ${_first}\n  ${_d}\n"
        "All curses libs must come from one prefix to avoid ABI mismatches.")
    endif()
    list(APPEND _accum "${_lib}")
  endforeach()
  set(${prefix} "${_first}" PARENT_SCOPE)
  set(${libs_out} "${_accum}" PARENT_SCOPE)
endfunction()

# Optional external hint
if(DEFINED ENV{NCURSES_ROOT})
  set(_CURSES_HINTS "$ENV{NCURSES_ROOT}/lib" "$ENV{NCURSES_ROOT}")
else()
  set(_CURSES_HINTS "")
endif()

# Locate candidates (prefer wide)
find_library(NCURSESW_LIBRARY NAMES ncursesw HINTS ${_CURSES_HINTS})
find_library(NCURSES_LIBRARY  NAMES ncurses  HINTS ${_CURSES_HINTS})

# If we found ncursesw, search companions in the same dir first.
set(_WIDE_MODE FALSE)
if(NCURSESW_LIBRARY)
  _dir_of(_ncw_dir "${NCURSESW_LIBRARY}")
  set(_WIDE_MODE TRUE)
  find_library(FORMW_LIBRARY  NAMES formw  HINTS "${_ncw_dir}" ${_CURSES_HINTS})
  find_library(PANELW_LIBRARY NAMES panelw HINTS "${_ncw_dir}" ${_CURSES_HINTS})
  find_library(MENUW_LIBRARY  NAMES menuw  HINTS "${_ncw_dir}" ${_CURSES_HINTS})
  find_library(TINFOW_LIBRARY  NAMES tinfow  HINTS "${_ncw_dir}" ${_CURSES_HINTS})
else()
    # Wide base (ncursesw) not found; will attempt narrow stack
endif()

# Decide stack and validate prefix
set(_libs "")
set(_prefix "")

if(_WIDE_MODE AND FORMW_LIBRARY AND PANELW_LIBRARY AND MENUW_LIBRARY)
  _require_same_prefix(_prefix _libs
  "${NCURSESW_LIBRARY}" "${FORMW_LIBRARY}" "${PANELW_LIBRARY}" "${MENUW_LIBRARY}" "${TINFOW_LIBRARY}")
else()
  if(NOT NCURSES_LIBRARY)
    message(FATAL_ERROR "ncurses.cmake: no ncursesw found and no narrow ncurses available.")
  endif()
  _dir_of(_nc_dir "${NCURSES_LIBRARY}")
  find_library(FORM_LIBRARY   NAMES form   HINTS "${_nc_dir}" ${_CURSES_HINTS})
  find_library(PANEL_LIBRARY  NAMES panel  HINTS "${_nc_dir}" ${_CURSES_HINTS})
  find_library(MENU_LIBRARY   NAMES menu   HINTS "${_nc_dir}" ${_CURSES_HINTS})
  find_library(TINFO_LIBRARY  NAMES tinfo  HINTS "${_nc_dir}" ${_CURSES_HINTS})

  _fatal_if_missing("form (narrow)"  FORM_LIBRARY)
  _fatal_if_missing("panel (narrow)" PANEL_LIBRARY)
  _fatal_if_missing("menu (narrow)"  MENU_LIBRARY)

  _require_same_prefix(_prefix _libs
    "${NCURSES_LIBRARY}" "${FORM_LIBRARY}" "${PANEL_LIBRARY}" "${MENU_LIBRARY}" "${TINFO_LIBRARY}")
endif()

# Build final ordered list and export
set(_CURSES_ALL_LIBS "${_libs}")
list(REMOVE_ITEM _CURSES_ALL_LIBS "")

# Cache in the "tool" style your build uses
set(NCURSES_FOUND TRUE CACHE BOOL "NCurses tool found" FORCE)
set(NCURSES_INCLUDE_DIR "${CURSES_INCLUDE_DIR}" CACHE PATH "Curses include dir" FORCE)
set(NCURSES_LIBRARIES   "${_CURSES_ALL_LIBS}"   CACHE PATH "Curses libs (+form/panel/menu/tinfo) — single prefix" FORCE)
set(NCURSES_PREFIX      "${_prefix}"            CACHE PATH "Prefix where curses libs were found" FORCE)

# Final summary line
if(COMMAND colormsg)
    colormsg(CYAN "+-- using curses from: ${NCURSES_PREFIX}")
else()
    message(STATUS "+-- using curses from: ${NCURSES_PREFIX}")
endif()

