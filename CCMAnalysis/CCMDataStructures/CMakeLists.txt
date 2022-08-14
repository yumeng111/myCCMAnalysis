set(PROJECT_NAME CCMDataStructures)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/Pulses.cxx
    private/RawData.cxx
    private/SimplifiedEvent.cxx
    private/SinglePulse.cxx
    private/Events.cxx
    private/AccumWaveform.cxx
    private/MCTruth.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils
    ROOTCINT
    public/CCMAnalysis/CCMDataStructures/Pulses.h
    public/CCMAnalysis/CCMDataStructures/RawData.h
    public/CCMAnalysis/CCMDataStructures/SimplifiedEvent.h
    public/CCMAnalysis/CCMDataStructures/SinglePulse.h
    public/CCMAnalysis/CCMDataStructures/Events.h
    public/CCMAnalysis/CCMDataStructures/AccumWaveform.h
    public/CCMAnalysis/CCMDataStructures/MCTruth.h
)


# set(LOCAL_HEADERS
#   data_structures.hh
#   Pulses.h
#   RawData.h
#   SimplifiedEvent.h
#   SinglePulse.h
#   Events.h
#   AccumWaveform.h
#   MCTruth.h
#   )
# 
# set(LOCAL_SRCS
#   Pulses.cxx
#   RawData.cxx
#   SimplifiedEvent.cxx
#   SinglePulse.cxx
#   Events.cxx
#   AccumWaveform.cxx
#   MCTruth.cxx
#   )
# 
# ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

