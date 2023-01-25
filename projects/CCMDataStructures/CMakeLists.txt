set(PROJECT_NAME CCMDataStructures)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/Pulses.cxx
    private/RawData.cxx
    private/SimplifiedEvent.cxx
    private/SinglePulse.cxx
    private/Events.cxx
    private/AccumWaveform.cxx
    private/MCTruth.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils serialization icetray
    ROOTCINT
    ../icetray/public/icetray/I3FrameObject.h
    public/CCMAnalysis/CCMDataStructures/Pulses.h
    private/Pulses_v1.h
    public/CCMAnalysis/CCMDataStructures/RawData.h
    private/RawData_v2.h
    public/CCMAnalysis/CCMDataStructures/SimplifiedEvent.h
    public/CCMAnalysis/CCMDataStructures/SinglePulse.h
    public/CCMAnalysis/CCMDataStructures/Events.h
    public/CCMAnalysis/CCMDataStructures/AccumWaveform.h
    private/AccumWaveform_v3.h
    public/CCMAnalysis/CCMDataStructures/MCTruth.h
)
