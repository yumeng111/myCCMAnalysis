# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(CCMAnalysis)

if(DEFINED ENV{CCMINSTALL})
  #set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "$ENV{CCMINSTALLBuild}")
  #set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "$ENV{CCMINSTALLBuild}")
  #set(CMAKE_INCLUDE_OUTPUT_DIRECTORY "$ENV{CCMINSTALLBuild}")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{CCMINSTALLBuild}")
else()
  message(FATAL_ERROR "You must have CCMINSTALL defined")
endif()

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/root/cmake)
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/share/root/cmake)
#list(APPEND CMAKE_MODULE_PATH /usr/share/root/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#message("${ROOT_USE_FILE}")
include(${ROOT_USE_FILE})
#message("INCLUDE ${ROOT_INCLUDE_DIR} LIBRARY ${ROOT_LIBRARY_DIR}")
#message("LIBRARIES ${ROOT_LIBRARIES}")
#message("XERCESINCLUDE $ENV{XERCESINCLUDE}")
include_directories(${ROOT_INCLUDE_DIR} $ENV{XERCESINCLUDE})
set(LIBRARIES ${LIBRARIES} ${ROOT_LIBRARIES} $ENV{XERCESLIB}/libxerces-c.so)
#message("LIBRARIES ${LIBRARIES}")

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} $ENV{XERCESINCLUD})
add_definitions(${ROOT_CXX_FLAGS} -std=c++14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Werror -std=c++14")

set(CCMAnalysis_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_SOURCE_DIR}/src/ds
  ${CMAKE_SOURCE_DIR}/src/io
  ${CMAKE_SOURCE_DIR}/src/spe
  ${CMAKE_SOURCE_DIR}/src/modules/framework
  ${CMAKE_SOURCE_DIR}/src/modules/reco
  ${CMAKE_SOURCE_DIR}/src/modules/nearline
  ${CMAKE_SOURCE_DIR}/src/modules/plots
  )

#Set the library version
set(CCMAnalysis_MAJOR_VERSION 4)
set(CCMAnalysis_MINOR_VERSION 0)
set(CCMAnalysis_PATCH_VERSION 0)
set(CCMAnalysis_VERSION
  "${CCMAnalysis_MAJOR_VERSION}.${CCMAnalysis_MINOR_VERSION}.${CCMAnalysis_PATCH_VERSION}")
set(CCMAnalysis_LIBRARY_PROPERTIES ${CCMAnalysis_LIBRARY_PROPERTIES}
  VERSION "${CCMAnalysis_VERSION}"
  SOVERSIONI "${CCMAnalysis_MAJOR_VERSION}"
  SUFFIX ".so"
  )


#uninstall target
configure_file(
  "${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY
  )

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# Recurse into the the following subdirectories.
#   src/utils
#   src/ds
#   src/spe
#   modules
# This does not actually cause another cmake executable to run. The same process
# will talk through the project's entire directory structure.
add_subdirectory(src/utils)
add_subdirectory(src/ds)
add_subdirectory(src/io)
add_subdirectory(src/spe)
add_subdirectory(src/modules)

#---Create  a main program using the library
#add_executable(EnergyCalibration energyCalibration.cc)
#target_link_libraries(EnergyCalibration Event)
#add_executable(NearlineDiag nearlineDiag.cc)
#target_link_libraries(NearlineDiag Event)
#add_executable(FindEvents findEvents.cc)
#target_link_libraries(FindEvents Event)
#add_executable(ApplyPreCuts applyPreCuts.cc)
#target_link_libraries(ApplyPreCuts Event)
