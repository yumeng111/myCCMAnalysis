# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(simulationCCM)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# R. T. Thornton Added these
#
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${ROOT_LIBRARIES})
set(CCMAnalysis_INCLUDE_DIRS
  $ENV{CCMINSTALLBuild}/src/utils
  $ENV{CCMINSTALLBuild}/src/ds
  $ENV{CCMINSTALLBuild}/src/io
  $ENV{CCMINSTALLBuild}/include
  )


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${CCMAnalysis_INCLUDE_DIRS})

set(LINK_DIRECTORIES /usr/lib64/ ${ROOT_LIBRARY_DIR} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

link_directories(${LINK_DIRECTORIES})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
set(CMAKE_CXX_LINKER_FLAGS "${CMAKE_CXX_LINKER_FLAGS} ${ROOT_LIBRARIES}")
add_executable(simulationCCM fullsim_CCM.cc ${sources} ${headers})
target_link_libraries(simulationCCM ${Geant4_LIBRARIES} CCMIO CCMDS CCMUtils)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(FS_SCRIPTS
  init_vis.mac
  runlaser.mac
  runcalib.mac
  )

foreach(_script ${FS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS simulationCCM DESTINATION bin)
