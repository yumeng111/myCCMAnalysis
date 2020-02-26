# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(EventTree)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/root/cmake)
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/share/root/cmake)
list(APPEND CMAKE_MODULE_PATH /usr/share/root/cmake)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${ROOT_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS} -std=c++11)

ROOT_GENERATE_DICTIONARY(G__Event PMTInfoMap.h SimplifiedEvent.h SinglePulse.h Pulses.h PMTInformation.h MsgLog.h Utility.h RawData.h SPECalibrationVariables.h MakeSPECalibrationHists.h NearlineSPEDiag.h LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(Event SHARED PMTInfoMap.cxx SimplifiedEvent.cxx SinglePulse.cxx Pulses.cxx PMTInformation.cxx MsgLog.cxx RawData.cxx SPECalibrationVariables.cxx MakeSPECalibrationHists.cxx NearlineSPEDiag.cxx Utility.cxx G__Event.cxx)
target_link_libraries(Event ${ROOT_LIBRARIES})

#---Create  a main program using the library
add_executable(EnergyCalibration energyCalibration.cc)
target_link_libraries(EnergyCalibration Event)
add_executable(NearlineDiag nearlineDiag.cc)
target_link_libraries(NearlineDiag Event)
add_executable(FindEvents findEvents.cc)
target_link_libraries(FindEvents Event)
add_executable(ApplyPreCuts applyPreCuts.cc)
target_link_libraries(ApplyPreCuts Event)
