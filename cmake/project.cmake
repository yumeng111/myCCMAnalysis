# -*- tab-width: 8; indent-tabs-mode: t -*- ex: ts=8 noet: 
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
# rootcint() handles root dictionary generation
#
if(NOT ROOT_FOUND OR NOT USE_CINT)
  message(STATUS "ROOT CINT is OFF")
  macro(ROOTCINT)
  endmacro(ROOTCINT)
else()
  message(STATUS "ROOT CINT is ON")
  macro(ROOTCINT TARGET)
    parse_arguments(ARG
	"LINKDEF;SOURCES;INCLUDE_DIRECTORIES;USE_TOOLS;USE_PROJECTS;LIBNAME"
      ""
      ${ARGN}
      )

    get_directory_property(incdirs INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})
    foreach(dir ${incdirs})
      list(APPEND ROOTCINT_INCLUDE_FLAGS -I${dir})
    endforeach(dir ${ARG_INCLUDE_DIRECTORIES})

    foreach(TOOL_ ${ARG_USE_TOOLS})
      string(TOUPPER ${TOOL_} TOOL)
      foreach(PATH_ ${${TOOL}_INCLUDE_DIR})
    list(APPEND ROOTCINT_INCLUDE_FLAGS "-I${PATH_}")
      endforeach(PATH_ ${${TOOL}_INCLUDE_DIR})
    endforeach(TOOL_ ${ARG_USE_TOOLS})

    foreach(PROJECT ${ARG_USE_PROJECTS})
      list(APPEND ROOTCINT_INCLUDE_FLAGS -I${CMAKE_SOURCE_DIR}/projects/${PROJECT}/public)
    endforeach(PROJECT ${ARG_USE_PROJECTS})

    set(ROOTCINT_HEADERS "")
    set(ROOTINTERNAL_HEADERS "")
    foreach(header ${ARG_SOURCES})
      if(EXISTS ${ROOT_INCLUDE_DIR}/${header})
    # If this is a ROOT header, don't add to dependencies or
    # rootcint flags as root adds these automagically.
    list(APPEND ROOTINTERNAL_HEADERS ${header})
      elseif(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${header})
    # if it isn't in our project (and isn't a root header), it
    # isn't legal here, the other project needs to build the dict.
    # I guess.
    message("In ${CMAKE_CURRENT_SOURCE_DIR}:")
    message(FATAL_ERROR "Header '${header}' passed to rootcint does not exist")
      else()
    # okay, it exists in our project, add it to our commandline
    # and be dependent on it.
    list(APPEND ROOTCINT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/${header})
      endif()
    endforeach(header ${ARG_SOURCES})

    add_custom_command(
      OUTPUT ${TARGET}
      DEPENDS ${ARG_LINKDEF} ${ROOTCINT_HEADERS}
      COMMAND ${CMAKE_BINARY_DIR}/env-shell.sh
      # rootcint found and ROOTSYS set in env-shell.sh path
      ARGS rootcint -f ${TARGET} -c -DCCM_USE_ROOT -DCCM_USE_CINT ${ROOTCINT_INCLUDE_FLAGS} -p ${CCM_UBER_HEADER} ${ROOTCINT_HEADERS} ${ROOTINTERNAL_HEADERS} ${ARG_LINKDEF}
      COMMENT "Generating ${TARGET} with rootcint"
      VERBATIM
      )
  if(${CMAKE_VERSION} VERSION_LESS "3.19.0")
   add_custom_target(
       "ROOTCINT_${ARG_LIBNAME}"
       ALL
       DEPENDS ${TARGET})
  endif()
  get_filename_component(dict_parent_dir ${TARGET} DIRECTORY)
  install(FILES "${dict_parent_dir}/${ARG_LIBNAME}Dict_rdict.pcm"
      DESTINATION lib)
  ENDMACRO(ROOTCINT)
ENDIF()

macro(use_projects THIS_TARGET)
  parse_arguments(${THIS_TARGET}_USE_PROJECTS
    "PROJECTS"
    ""
    ${ARGN}
    )
  foreach(USED_PROJECT ${${THIS_TARGET}_USE_PROJECTS_PROJECTS})
    if(NOT IS_DIRECTORY ${CCM_SRC}/${USED_PROJECT})
      message(FATAL_ERROR "Attempt to use nonexistent project '${USED_PROJECT}'")
    endif(NOT IS_DIRECTORY ${CCM_SRC}/${USED_PROJECT})
    if(NOT EXISTS ${CCM_SRC}/${USED_PROJECT}/CMakeLists.txt)
      message(FATAL_ERROR "Attempt to use project '${USED_PROJECT}'. There is a directory but no CMakeLists.txt... don't know what to do.")
    endif(NOT EXISTS ${CCM_SRC}/${USED_PROJECT}/CMakeLists.txt)


    include_directories(${CCM_SRC}/${USED_PROJECT}/public)
    target_link_libraries(${THIS_TARGET} ${USED_PROJECT})
  endforeach(USED_PROJECT ${${THIS_TARGET}_USE_PROJECTS_PROJECTS})
endmacro(use_projects THIS_TARGET)


macro(ccm_add_library THIS_LIB_NAME)
  if (BUILD_${CCM_PROJECT})
    #
    # Grrr...  this *_ARGS variable has to be unique to the project, otherwise
    # if you have two instances of ccm_add_library the second will include the
    # parsed arg values from the first.
    #
    parse_arguments(${THIS_LIB_NAME}_ARGS
      "USE_TOOLS;USE_PROJECTS;ROOTCINT;INSTALL_DESTINATION;LINK_LIBRARIES;COMPILE_FLAGS"
      "MODULE;EXCLUDE_FROM_ALL;WITHOUT_CCM_HEADERS;NO_DOXYGEN;IWYU;PYBIND11;NOUNDERSCORE;INTERFACE"
      ${ARGN}
      )

    include_directories(
      ${PROJECT_SOURCE_DIR}/public
      ${PROJECT_SOURCE_DIR}/private
      )

    no_dotfile_glob(${THIS_LIB_NAME}_ARGS_SOURCES ${${THIS_LIB_NAME}_ARGS_DEFAULT_ARGS})

    if(${THIS_LIB_NAME}_ARGS_ROOTCINT)
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
        set(LINKDEF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
      else()
        set(LINKDEF_FILE "NOTFOUND")
      endif()
    endif(${THIS_LIB_NAME}_ARGS_ROOTCINT)

    if (LINKDEF_FILE AND ${THIS_LIB_NAME}_ARGS_ROOTCINT AND USE_CINT)

      set (DICTIONARY_SOURCEFILE
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${THIS_LIB_NAME}Dict.cxx)

      set (${THIS_LIB_NAME}_ARGS_SOURCES ${${THIS_LIB_NAME}_ARGS_SOURCES} ${DICTIONARY_SOURCEFILE})

      rootcint(${DICTIONARY_SOURCEFILE}
    LINKDEF ${LINKDEF_FILE}
    SOURCES ${${THIS_LIB_NAME}_ARGS_ROOTCINT}
    USE_PROJECTS ${${THIS_LIB_NAME}_ARGS_USE_PROJECTS} ${THIS_LIB_NAME}
    USE_TOOLS    ${${THIS_LIB_NAME}_ARGS_USE_TOOLS}
    LIBNAME ${THIS_LIB_NAME}
    )

    endif (LINKDEF_FILE AND ${THIS_LIB_NAME}_ARGS_ROOTCINT AND USE_CINT)

    set (ARGS)
    if (${THIS_LIB_NAME}_ARGS_EXCLUDE_FROM_ALL)
       set(ARGS EXCLUDE_FROM_ALL)
    endif (${THIS_LIB_NAME}_ARGS_EXCLUDE_FROM_ALL)
    if (${THIS_LIB_NAME}_ARGS_MODULE)
       set(ARGS ${ARGS} MODULE)
    endif (${THIS_LIB_NAME}_ARGS_MODULE)


    if(${THIS_LIB_NAME}_ARGS_INTERFACE)
  	if(${CMAKE_VERSION} VERSION_LESS "3.19.0")
            add_library(${THIS_LIB_NAME} INTERFACE ${ARGS})
        else(${CMAKE_VERSION} VERSION_LESS "3.19.0")
            add_library(${THIS_LIB_NAME} INTERFACE ${ARGS} ${${THIS_LIB_NAME}_ARGS_SOURCES})
  	endif(${CMAKE_VERSION} VERSION_LESS "3.19.0")
    else(${THIS_LIB_NAME}_ARGS_INTERFACE)
      add_library(${THIS_LIB_NAME} ${ARGS} ${${THIS_LIB_NAME}_ARGS_SOURCES})
    endif(${THIS_LIB_NAME}_ARGS_INTERFACE)
    set(LIB_${THIS_LIB_NAME}_EXISTS ON)
    add_dependencies(${THIS_LIB_NAME} env-check)

    if(${THIS_LIB_NAME}_ARGS_INTERFACE AND ${CMAKE_VERSION} VERSION_LESS "3.19.0")
    else(${THIS_LIB_NAME}_ARGS_INTERFACE AND ${CMAKE_VERSION} VERSION_LESS "3.19.0")
    set_target_properties(${THIS_LIB_NAME}
      PROPERTIES
      COMPILE_DEFINITIONS PROJECT=${PROJECT_NAME}
      NO_SYSTEM_FROM_IMPORTED TRUE
      )
    endif(${THIS_LIB_NAME}_ARGS_INTERFACE AND ${CMAKE_VERSION} VERSION_LESS "3.19.0")

    if(${THIS_LIB_NAME}_ARGS_IWYU AND USE_IWYU)
      set_target_properties(${THIS_LIB_NAME}
        PROPERTIES
        CXX_INCLUDE_WHAT_YOU_USE ${IWYU_PROGRAM}
        )
    endif()

    if(${THIS_LIB_NAME}_ARGS_INTERFACE)
    else(${THIS_LIB_NAME}_ARGS_INTERFACE)
    add_custom_command(TARGET ${THIS_LIB_NAME}
      PRE_LINK
      COMMAND mkdir -p ${LIBRARY_OUTPUT_PATH}
      )
    endif(${THIS_LIB_NAME}_ARGS_INTERFACE)

    # Disabled all special linker flags for APPLE:
    #  - single_module: this is the default anyway
    #  - undefined dynamic_lookup: it seems not to hurt letting the
    #      linker throw an error for undefined symbols.
    #  - flat_namespace: not using the two-level namespace (library+symbol name)
    #      seems to introduce bugs in exception handling with boost::python.
    #
    #if(APPLE)
    #  set_target_properties(${THIS_LIB_NAME}
    #  PROPERTIES
    #    LINK_FLAGS "-single_module -undefined dynamic_lookup -flat_namespace"
    #    )
    #endif(APPLE)

    if(${THIS_LIB_NAME}_ARGS_INTERFACE AND ${CMAKE_VERSION} VERSION_LESS "3.19.0")
    else(${THIS_LIB_NAME}_ARGS_INTERFACE AND ${CMAKE_VERSION} VERSION_LESS "3.19.0")
    if(NOT ${THIS_LIB_NAME}_ARGS_WITHOUT_CCM_HEADERS)
      set_target_properties(${THIS_LIB_NAME}
    PROPERTIES
    COMPILE_FLAGS ""
    )
    endif(NOT ${THIS_LIB_NAME}_ARGS_WITHOUT_CCM_HEADERS)
    if(${THIS_LIB_NAME}_ARGS_COMPILE_FLAGS)
      set_target_properties(${THIS_LIB_NAME}
    PROPERTIES
    COMPILE_FLAGS ${${THIS_LIB_NAME}_ARGS_COMPILE_FLAGS}
    )
    endif(${THIS_LIB_NAME}_ARGS_COMPILE_FLAGS)
    endif(${THIS_LIB_NAME}_ARGS_INTERFACE AND ${CMAKE_VERSION} VERSION_LESS "3.19.0")

    if (${THIS_LIB_NAME}_ARGS_LINK_LIBRARIES)
      target_link_libraries(${THIS_LIB_NAME} ${${THIS_LIB_NAME}_ARGS_LINK_LIBRARIES})
    endif (${THIS_LIB_NAME}_ARGS_LINK_LIBRARIES)

    use_projects(${THIS_LIB_NAME}
      PROJECTS "${${THIS_LIB_NAME}_ARGS_USE_PROJECTS}"
      )
    set(${THIS_LIB_NAME}_PROJECT_DEPENDS ${${THIS_LIB_NAME}_ARGS_USE_PROJECTS} CACHE INTERNAL "Projects needed by library ${THIS_LIB_NAME}")

    use_tools(${THIS_LIB_NAME}
      TOOLS "${${THIS_LIB_NAME}_ARGS_USE_TOOLS}"
      )

    if(NOT ${THIS_LIB_NAME}_ARGS_EXCLUDE_FROM_ALL)
      if(${THIS_LIB_NAME}_ARGS_INSTALL_DESTINATION)
        install(TARGETS ${THIS_LIB_NAME} DESTINATION ${${THIS_LIB_NAME}_ARGS_INSTALL_DESTINATION})
      else(${THIS_LIB_NAME}_ARGS_INSTALL_DESTINATION)
        install(TARGETS ${THIS_LIB_NAME} DESTINATION lib)
      endif(${THIS_LIB_NAME}_ARGS_INSTALL_DESTINATION)
    endif(NOT ${THIS_LIB_NAME}_ARGS_EXCLUDE_FROM_ALL)

    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/doxyfile
      @ONLY
      )

    if(NOT ${THIS_LIB_NAME}_ARGS_NO_DOXYGEN AND DOXYGEN_FOUND)

      add_custom_target(${PROJECT_NAME}-${THIS_LIB_NAME}-doxygen
    COMMAND mkdir -p ${DOXYGEN_OUTPUT_PATH}/${PROJECT_NAME}
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/doxyfile
    )
      add_dependencies(${PROJECT_NAME}-doxygen
    ${PROJECT_NAME}-${THIS_LIB_NAME}-doxygen
    )
      add_dependencies(doxygen ${PROJECT_NAME}-doxygen)
      foreach(l ${${THIS_LIB_NAME}_PROJECT_DEPENDS})
    file(APPEND ${DOXYGEN_OUTPUT_PATH}/.tagfiles/${PROJECT_NAME}.include "${DOXYGEN_OUTPUT_PATH}/.tagfiles/${l}.tag=../${l}\n")
    add_dependencies(${PROJECT_NAME}-${THIS_LIB_NAME}-doxygen ${l}-doxygen)
      endforeach()

    endif(NOT ${THIS_LIB_NAME}_ARGS_NO_DOXYGEN AND DOXYGEN_FOUND)

    if(DPKG_INSTALL_PREFIX)
      set_target_properties(${THIS_LIB_NAME}
    PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    )

    endif(DPKG_INSTALL_PREFIX)

  endif(BUILD_${CCM_PROJECT})
endmacro(ccm_add_library)


macro(ccm_project PROJECT_NAME)
  project(${PROJECT_NAME})
  parse_arguments(ARG
    "PYTHON_DIR;PYTHON_DEST;DOCS_DIR"
    "USE_SETUPTOOLS"
    ${ARGN}
    )

  string(TOUPPER "${PROJECT_NAME}" CCM_PROJECT)
  string(TOUPPER "BUILD_${CCM_PROJECT}" BUILD_PROJECT_OPTION)
  option(${BUILD_PROJECT_OPTION}
    "Build project CCM.${PROJECT_NAME} (prefer using make targets, not this, for building individual libs)"
    ON)

  if(BUILD_${CCM_PROJECT})

    add_custom_target(${PROJECT_NAME}-doxygen)

    #ccm_add_testing_targets()

    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/resources)
      install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/resources/
	  DESTINATION lib/CCMAnalysis/${PROJECT_NAME}/resources/
    PATTERN ".svn" EXCLUDE
    PATTERN ".git" EXCLUDE
    PATTERN "*.py"
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/resources)

    if(ARG_DOCS_DIR)
      file(MAKE_DIRECTORY "${SPHINX_DIR}/source/projects")
      execute_process(COMMAND ln -fsn
        ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_DOCS_DIR}
        ${SPHINX_DIR}/source/projects/${PROJECT_NAME})
    endif(ARG_DOCS_DIR)

    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/public/CCMAnalysis/${PROJECT_NAME} AND INSTALL_HEADERS)
      install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/public/CCMAnalysis/${PROJECT_NAME}
      DESTINATION include/CCMAnalysis
      PATTERN ".git" EXCLUDE
      PATTERN ".svn" EXCLUDE)
    endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/public/CCMAnalysis/${PROJECT_NAME} AND INSTALL_HEADERS)

    if(ARG_PYTHON_DIR AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR})
      if(ARG_USE_SETUPTOOLS)
    colormsg(GREEN "+-- python [setuptools]")

    #
    # do the 'setup.py develop'
    #
    execute_process(COMMAND
      /usr/bin/env PYTHONPATH=${LIBRARY_OUTPUT_PATH}:$ENV{PYTHONPATH}
      ${PYTHON_EXECUTABLE} setup.py -q develop
      --install-dir ${LIBRARY_OUTPUT_PATH}
      --script-dir  ${EXECUTABLE_OUTPUT_PATH}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )

    #
    # Install targets
    #
    add_custom_target(${PROJECT_NAME}-install-to-tarball
      /usr/bin/env PYTHONPATH=${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/lib:$ENV{PYTHONPATH}
      ${PYTHON_EXECUTABLE} setup.py install
      --install-lib ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/lib
      --install-scripts  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/bin
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
    add_dependencies(${PROJECT_NAME}-install-to-tarball
      tarball-install)
    add_dependencies(tarball-finish
      ${PROJECT_NAME}-install-to-tarball)

      else(ARG_USE_SETUPTOOLS)
    if (COPY_PYTHON_DIR)
      colormsg(GREEN "+-- python [directory copy]")
    else (COPY_PYTHON_DIR)
      colormsg(GREEN "+-- python [symlinks]")
    endif (COPY_PYTHON_DIR)

    if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR}/__init__.py)
      message(FATAL_ERROR
        "Project ${PROJECT_NAME} has PYTHON_DIR specified, but the directory contains no file '__init__.py' and will be useless")
    endif (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR}/__init__.py)

    if(NOT ARG_PYTHON_DEST)
      set(ARG_PYTHON_DEST ccm/${PROJECT_NAME})
      string(REPLACE "-" "_" ARG_PYTHON_DEST "ccm/${PROJECT_NAME}")
    endif(NOT ARG_PYTHON_DEST)

    #
    #  Just bare python, no setuptools
    #
    if (COPY_PYTHON_DIR)
      file(GLOB_RECURSE python_components RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${ARG_PYTHON_DIR}/*.py)
      foreach(file ${python_components})
            string(REPLACE ${ARG_PYTHON_DIR}/ "" file ${file})
            configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR}/${file} ${CMAKE_BINARY_DIR}/lib/${ARG_PYTHON_DEST}/${file} COPYONLY)
          endforeach()
        else (COPY_PYTHON_DIR)
      execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR} ${CMAKE_BINARY_DIR}/lib/${ARG_PYTHON_DEST})
      install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR}/
        DESTINATION lib/${ARG_PYTHON_DEST}
        PATTERN ".git" EXCLUDE
        PATTERN ".svn" EXCLUDE)
      execute_process(COMMAND python -m compileall -fq ${CMAKE_BINARY_DIR}/lib/${ARG_PYTHON_DEST} OUTPUT_QUIET)
    endif (COPY_PYTHON_DIR)
      endif(ARG_USE_SETUPTOOLS)

    endif(ARG_PYTHON_DIR AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PYTHON_DIR})

  endif(BUILD_${CCM_PROJECT})
endmacro(ccm_project PROJECT_NAME)


macro(ccm_executable_script THIS_EXECUTABLE_NAME THIS_SCRIPT_NAME)
    if(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_NO_PREFIX)
      set(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME
    ${THIS_EXECUTABLE_NAME})
    else(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_NO_PREFIX)
      set(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME
    ${PROJECT_NAME}-${THIS_EXECUTABLE_NAME})
    endif(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_NO_PREFIX)

    # copy it for local use
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${THIS_SCRIPT_NAME}
        ${EXECUTABLE_OUTPUT_PATH}/${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
        COPYONLY)
    # and install it in the tarball when the time comes
    install(PROGRAMS ${THIS_SCRIPT_NAME} DESTINATION bin
        RENAME ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME})
endmacro(ccm_executable_script THIS_EXECUTABLE_NAME THIS_SCRIPT_NAME)


macro(ccm_executable THIS_EXECUTABLE_NAME)
  if(BUILD_${CCM_PROJECT})
    parse_arguments(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}
      "USE_TOOLS;USE_PROJECTS;LINK_LIBRARIES"
      "NO_PREFIX;WITHOUT_CCM_HEADERS;IWYU"
      ${ARGN}
      )
    no_dotfile_glob(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_SOURCES
      ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_DEFAULT_ARGS})

    if(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_NO_PREFIX)
      set(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME
    ${THIS_EXECUTABLE_NAME})
    else(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_NO_PREFIX)
      set(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME
    ${PROJECT_NAME}-${THIS_EXECUTABLE_NAME})
    endif(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_NO_PREFIX)

    add_executable(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_SOURCES})

    if(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_WITHOUT_CCM_HEADERS)
      message(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_WITHOUT_CCM_HEADERS)
      set(THIS_CCMH_FLAGS "")
    else()
      set(THIS_CCMH_FLAGS "")
    endif()

    add_dependencies(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME} env-check)

    use_projects(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      PROJECTS ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_USE_PROJECTS})
    use_tools(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      TOOLS "${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_USE_TOOLS}"
      )

    if(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_LINK_LIBRARIES)
    target_link_libraries(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_LINK_LIBRARIES})
    endif(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_LINK_LIBRARIES)

    ## FIXME: temporarily force pthread linking
    target_compile_options(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      BEFORE PUBLIC "-pthread")
    set_target_properties(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      PROPERTIES LINK_FLAGS "-pthread")
    ## :FIXME

    set_target_properties(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
      PROPERTIES
      COMPILE_FLAGS "${THIS_CCMH_FLAGS} -DPROJECT=${PROJECT_NAME}"
      )

    if(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_IWYU AND USE_IWYU)
      set_target_properties(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
        PROPERTIES
        CXX_INCLUDE_WHAT_YOU_USE ${IWYU_PROGRAM}
        )
    endif()

    if(APPLE)
      set_target_properties(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
    PROPERTIES
    LINK_FLAGS "-bind_at_load -multiply_defined suppress"
    )
    endif(APPLE)
    install(TARGETS ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME} RUNTIME DESTINATION bin)

    if(DPKG_INSTALL_PREFIX)
      set_target_properties(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_TARGET_NAME}
    PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    )
    endif(DPKG_INSTALL_PREFIX)

  endif(BUILD_${CCM_PROJECT})
endmacro(ccm_executable THIS_EXECUTABLE_NAME)

macro(ccm_test_executable THIS_EXECUTABLE_NAME)
  if (BUILD_${CCM_PROJECT})
    add_test(NAME "${PROJECT_NAME}::${THIS_EXECUTABLE_NAME}" #::${testable_file}/${unittest}"
             WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
             COMMAND ${PROJECT_NAME}-${THIS_EXECUTABLE_NAME} -t 1100 -saf)
    set_tests_properties("${PROJECT_NAME}::${THIS_EXECUTABLE_NAME}"
                     PROPERTIES LABELS BINARY)

    parse_arguments(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}
      "USE_TOOLS;USE_PROJECTS;USE_PYBINDINGS;LINK_LIBRARIES"
      ""
      ${ARGN}
      )
    no_dotfile_glob(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_SOURCES
      ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_DEFAULT_ARGS}
      )

    add_executable(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      EXCLUDE_FROM_ALL
      ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_SOURCES}
      )

    add_dependencies(test-bins ${PROJECT_NAME}-${THIS_EXECUTABLE_NAME})
    use_pybindings(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME})

    ## FIXME: temporarily force pthread linking
    target_compile_options(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      BEFORE PUBLIC "-pthread")
    set_target_properties(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      PROPERTIES LINK_FLAGS "-pthread")
    ## :FIXME

    set_target_properties(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      PROPERTIES
      COMPILE_DEFINITIONS PROJECT=${PROJECT_NAME}
      )
    if(APPLE)
      set_target_properties(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
    PROPERTIES
    LINK_FLAGS "-bind_at_load -multiply_defined suppress"
    )
    endif(APPLE)

    use_projects(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      PROJECTS ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_USE_PROJECTS})

    use_pybindings(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      PROJECTS ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_USE_PYBINDINGS})

    use_tools(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
      TOOLS "${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_USE_TOOLS}")

    set_source_files_properties(${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_SOURCES}
      PROPERTIES
      COMPILE_FLAGS ""
      )

    if(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_LINK_LIBRARIES)
      target_link_libraries(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
    ${${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_LINK_LIBRARIES})
    endif(${PROJECT_NAME}_${THIS_EXECUTABLE_NAME}_LINK_LIBRARIES)

    if(DPKG_INSTALL_PREFIX)
      set_target_properties(${PROJECT_NAME}-${THIS_EXECUTABLE_NAME}
    PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    )
    endif(DPKG_INSTALL_PREFIX)
  endif ()
endmacro(ccm_test_executable THIS_EXECUTABLE_NAME)


macro(ccm_test_scripts)
  no_dotfile_glob(${PROJECT_NAME}_SCRIPTS ${ARGN})
  sort(${PROJECT_NAME}_SCRIPTS_ALPHABETICAL "${${PROJECT_NAME}_SCRIPTS}")
  foreach(script ${${PROJECT_NAME}_SCRIPTS_ALPHABETICAL})
    get_filename_component(script_basename ${script} NAME)
    if(${script} MATCHES "\\.py$")
      add_test(NAME ${PROJECT_NAME}::${script_basename}
               WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
               COMMAND ${PYTHON_EXECUTABLE} ${script})
    else()
      add_test(NAME ${PROJECT_NAME}::${script_basename}
               WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
               COMMAND ${script})
    endif()
  endforeach()
endmacro(ccm_test_scripts)


macro(project_moved_to NEWLOCATION)
  get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME_WE)
  message(STATUS "************************************************************")
  message(STATUS "***")
  message(STATUS "***   Subversion for project ${PROJECT_NAME} has moved:")
  message(STATUS "***")
  message(STATUS "***   The new location is ")
  message(STATUS "***   ${NEWLOCATION}")
  message(STATUS "***")
  message(STATUS "***   Please update your svn:externals accordingly.")
  message(STATUS "***")
  message(STATUS "************************************************************")
  message(FATAL_ERROR "Repository moved.")
endmacro(project_moved_to PROJECT NEWLOCATION)

