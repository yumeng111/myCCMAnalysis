set(_CDK_HINT_PREFIXES ${CMAKE_PREFIX_PATH})

# Prefer wide library name first.
set(_CDK_LIBNAMES cdkw cdk)

# Find library path (absolute) and choose name/dir.
set(_CDK_LIB_PATH "")
set(_CDK_LIB_NAME "")
foreach(_name IN LISTS _CDK_LIBNAMES)
  find_library(_try_LIB
    NAMES ${_name}
    HINTS ${_CDK_HINT_PREFIXES}
    PATH_SUFFIXES lib lib64 lib/cdk
    NO_DEFAULT_PATH
  )
  if(NOT _try_LIB STREQUAL "_try_LIB-NOTFOUND")
    set(_CDK_LIB_PATH "${_try_LIB}")
    set(_CDK_LIB_NAME "${_name}")
    break()
  endif()
endforeach()

# Derive a candidate prefix from the found library (…/prefix/lib*/libcdk*.so)
set(_CDK_PREFIX "")
set(_CDK_LIB_DIR "")
if(_CDK_LIB_PATH)
  get_filename_component(_CDK_LIB_DIR "${_CDK_LIB_PATH}" DIRECTORY)
  # prefix should be parent of lib/lib64
  get_filename_component(_CDK_PREFIX "${_CDK_LIB_DIR}" DIRECTORY)
endif()

# Try headers under the same prefix first; prefer include/cdk/cdk.h, then include/cdk.h
set(_CDK_INC_DIR "")
set(_CDK_INC_FILE "")
if(_CDK_PREFIX)
  if(EXISTS "${_CDK_PREFIX}/include/cdk/cdk.h")
    set(_CDK_INC_DIR  "${_CDK_PREFIX}/include")
    set(_CDK_INC_FILE "cdk/cdk.h")
  elseif(EXISTS "${_CDK_PREFIX}/include/cdk.h")
    set(_CDK_INC_DIR  "${_CDK_PREFIX}/include")
    set(_CDK_INC_FILE "cdk.h")
  endif()
endif()

# If still unknown (library not found yet or headers not in same prefix), do a broader search.
if(NOT _CDK_INC_DIR)
  # Try hierarchical header first
  find_path(_try_INC_DIR
    NAMES cdk/cdk.h
    HINTS ${_CDK_HINT_PREFIXES}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
  )
  if(NOT _try_INC_DIR STREQUAL "_try_INC_DIR-NOTFOUND")
    set(_CDK_INC_DIR  "${_try_INC_DIR}")
    set(_CDK_INC_FILE "cdk/cdk.h")
  else()
    find_path(_try_INC_DIR2
      NAMES cdk.h
      HINTS ${_CDK_HINT_PREFIXES}
      PATH_SUFFIXES include
      NO_DEFAULT_PATH
    )
    if(NOT _try_INC_DIR2 STREQUAL "_try_INC_DIR2-NOTFOUND")
      set(_CDK_INC_DIR  "${_try_INC_DIR2}")
      set(_CDK_INC_FILE "cdk.h")
    endif()
  endif()
endif()

# Fallbacks if nothing was discovered: keep legacy defaults so tooldef will print a single failure report.
if(NOT _CDK_INC_DIR)
  set(_CDK_INC_DIR "include")
  # Try hierarchical header first; tooldef will report whichever fails.
  set(_CDK_INC_FILE "cdk/cdk.h")
endif()
if(NOT _CDK_LIB_DIR)
  set(_CDK_LIB_DIR "lib/cdk")  # legacy first
endif()
if(NOT _CDK_LIB_NAME)
  set(_CDK_LIB_NAME "cdk")     # legacy first; tooldef will say not found if only wide exists
endif()

# Single tooldef call to get user printouts

tooldef (cdk
  ${_CDK_INC_DIR}
  ${_CDK_INC_FILE}
  ${_CDK_LIB_DIR}
  NONE
  ${_CDK_LIB_NAME}
)

if (CDK_FOUND)
  # If headers live under include/cdk/, add that subdir too.
  if (EXISTS "${CDK_INCLUDE_DIR}/cdk")
    list(APPEND CDK_INCLUDE_DIR "${CDK_INCLUDE_DIR}/cdk")
    # de-duplicate
    list(REMOVE_DUPLICATES CDK_INCLUDE_DIR)
  endif()
endif()

