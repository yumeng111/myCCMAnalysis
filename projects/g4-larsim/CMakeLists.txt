set(PROJECT_NAME g4-larsim)

if (GEANT4_CONFIG_ERROR OR NOT DEFINED GEANT4_CONFIG_ERROR)
    include(FindGeant4.cmake)
    include(FindCLHEP.cmake)
else (GEANT4_CONFIG_ERROR OR NOT DEFINED GEANT4_CONFIG_ERROR)
    set(GEANT4_FOUND TRUE)
    set(CLHEP_FOUND TRUE)
endif (GEANT4_CONFIG_ERROR OR NOT DEFINED GEANT4_CONFIG_ERROR)

if (GEANT4_FOUND AND CLHEP_FOUND)
    i3_project(${PROJECT_NAME}
        PYTHON_DIR python
        DOCS_DIR resources/docs
    )

    i3_add_library(${PROJECT_NAME}
        private/g4-larsim/g4classes/G4CCMDetectorConstruction.h
        private/g4-larsim/g4classes/G4CCMDetectorMessenger.h
        private/g4-larsim/g4classes/G4CCMPhysicsList.h
        private/g4-larsim/g4classes/G4CCMRunManager.h
        private/g4-larsim/g4classes/G4CCMUserSteppingAction.h
        private/g4-larsim/g4classes/G4Interface.h
        private/g4-larsim/g4classes/G4CCMDetectorConstruction.cxx
        private/g4-larsim/g4classes/G4CCMDetectorMessenger.cxx
        private/g4-larsim/g4classes/G4CCMPhysicsList.cxx
        private/g4-larsim/g4classes/G4CCMRunManager.cxx
        private/g4-larsim/g4classes/G4CCMUserSteppingAction.cxx
        private/g4-larsim/g4classes/G4Interface.cxx
        USE_TOOLS boost python clhep geant4 gsl
        USE_PROJECTS serialization icetray dataclasses phys-services
    )

    include(FindG4VisDrivers.cmake)

else (GEANT4_FOUND AND CLHEP_FOUND)
    message(STATUS "  - Geant4 and/or CLHEP not found - Skipping build of g4-larsim")
endif (GEANT4_FOUND AND CLHEP_FOUND)

