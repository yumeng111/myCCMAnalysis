set(PROJECT_NAME g4-larsim)

if(GEANT4_FOUND)
    i3_project(${PROJECT_NAME}
        PYTHON_DIR python
        DOCS_DIR resources/docs
    )

    i3_add_library(${PROJECT_NAME}
        private/g4-larsim/CCMSimpleInjector.cxx
        private/g4-larsim/CCM200Response.cxx
        private/g4-larsim/CCMSimulator.cxx
        private/g4-larsim/g4classes/*.h
        private/g4-larsim/g4classes/*.cxx
        USE_TOOLS boost python clhep geant4 gsl
        USE_PROJECTS serialization icetray dataclasses phys-services simclasses
    )

    i3_executable(CCM
        apps/*.cxx
        USE_TOOLS clhep geant4
        USE_PROJECTS serialization icetray dataclasses phys-services simclasses g4-larsim
        NO_PREFIX
    )
else (GEANT4_FOUND)
    message(STATUS "  - Geant4 and/or CLHEP not found - Skipping build of g4-larsim")
endif (GEANT4_FOUND)

