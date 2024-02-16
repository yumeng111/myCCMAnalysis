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
        private/g4-larsim/g4classes/*.hh
        private/g4-larsim/g4classes/*.cxx
        private/g4-larsim/g4classes/*.cc
        USE_TOOLS boost python clhep geant4 gsl
        USE_PROJECTS serialization icetray dataclasses phys-services
    )

    include(FindG4VisDrivers.cmake)

else (GEANT4_FOUND AND CLHEP_FOUND)
    message(STATUS "  - Geant4 and/or CLHEP not found - Skipping build of g4-larsim")
endif (GEANT4_FOUND AND CLHEP_FOUND)

