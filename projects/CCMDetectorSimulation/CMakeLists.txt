set(PROJECT_NAME CCMDetectorSimulation)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/actionInitialization.cc
    private/detectorConstruction.cc
    private/detectorMessenger.cc
    private/eventAction.cc
    private/primaryGenerator.cc
    private/runAction.cc
    private/steppingAction.cc
    USE_TOOLS root geant4
    USE_PROJECTS CCMBinary CCMUtils CCMIO CCMDataStructures
)

ccm_executable(CCMG4Simulation
    apps/CCMDetectorSimulation.cxx
    USE_TOOLS root geant4
    USE_PROJECTS CCMBinary CCMUtils CCMIO CCMDataStructures CCMDetectorSimulation
    NO_PREFIX
)
