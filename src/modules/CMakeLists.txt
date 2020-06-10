# The name of our project is "CCMAnalysis". CMakeLists files in this project can
# refer to the root source directory of the project as ${CCMAnalysis_SOURCE_Dir} and
# to the root binary directory of the project as ${CCMAnalysis_BINARY_DIR}.

# Recurse into the the following subdirectories.
#   fracework
#   detsim
#   reco
#   evtDisplay
# This does not actually cause another cmake executable to run. The same process
# will talk through the project's entire directory structure.
add_subdirectory(framework)
add_subdirectory(reco)
add_subdirectory(nearline)
add_subdirectory(plots)


set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR} $ENV{XERCESINCLUDE} ${CCMAnalysis_INCLUDE_DIRS})

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES /usr/lib64/ ${ROOT_LIBRARY_DIR} $ENV{XERCESLIB} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

link_directories(${LINK_DIRECTORIES})

set(CMAKE_CXX_LINKER_FLAGS "${CMAKE_CXX_LINKER_FLAGS} ${ROOT_LIBRARIES} $ENV{XERCESLIB}")

add_executable(CCMAnalysis CCMAnalysis.cxx)
target_link_libraries(CCMAnalysis ${ROOT_LIBRARIES} CCMUtils CCMDS CCMSPE CCMXML CCMReco CCMIO CCMNearline CCMPlots)

