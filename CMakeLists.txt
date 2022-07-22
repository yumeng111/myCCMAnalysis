# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(CCMAnalysis)

message(STATUS "source dir: ${CMAKE_SOURCE_DIR}")
message(STATUS "binary dir: ${CMAKE_BINARY_DIR}")

set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include/${CMAKE_PROJECT_NAME})

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

if(DEFINED ENV{XERCESINSTALL})
    set(XercesC_INCLUDE_DIRS $ENV{XERCESINSTALL}/include)
    set(XercesC_LIBRARIES xerces-c)
    set(XercesC_LIBRARY_DIRS $ENV{XERCESINSTALL}/lib)
else()
    message("-- XercesC not found in env variables")
    message("-- Checking for XercesC with find_package")
    find_package(XercesC)
    if(XercesC_FOUND)
        message("-- XercesC found by find_package")
    else()
        message(FATAL_ERROR "XercesC not found by find_package")
    endif()
endif()

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR} ${XercesC_INCLUDE_DIRS})
set(LIBRARIES ${LIBRARIES} ${ROOT_LIBRARIES} ${XercesC_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${XercesC_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS} -std=c++14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Werror -std=c++14 -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-unused-result -std=c++14 -Og")

set(CCMAnalysis_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/src/
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

if(NOT DEFINED BUILD_SIMULATION)
    message(STATUS "Setting BUILD_SIMULATION=TRUE by default.")
    set(BUILD_SIMULATION true)
endif()
if(${BUILD_SIMULATION})
    add_subdirectory(simulationCCM)
endif()
#---Create  a main program using the library
#add_executable(EnergyCalibration energyCalibration.cc)
#target_link_libraries(EnergyCalibration Event)
#add_executable(NearlineDiag nearlineDiag.cc)
#target_link_libraries(NearlineDiag Event)
#add_executable(FindEvents findEvents.cc)
#target_link_libraries(FindEvents Event)
#add_executable(ApplyPreCuts applyPreCuts.cc)
#target_link_libraries(ApplyPreCuts Event)

