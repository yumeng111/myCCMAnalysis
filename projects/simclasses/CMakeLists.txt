set(PROJECT_NAME simclasses)

i3_project(${PROJECT_NAME}
    PYTHON_DIR python
    DOCS_DIR resources/docs)

SET(LIB_${PROJECT_NAME}_SOURCEFILES
    private/simclasses/CCMMCPE.cxx
    private/simclasses/PhotonSummary.cxx
)


SET(LIB_${PROJECT_NAME}_PROJECTS serialization icetray dataclasses phys-services)
SET(LIB_${PROJECT_NAME}_TOOLS boost python)

if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)
    LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS tableio)
    LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES private/simclasses/converter/*.cxx)
    ADD_DEFINITIONS(-DUSE_TABLEIO)
    MESSAGE(STATUS "     + tableio support")
endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)

i3_add_library(${PROJECT_NAME} ${LIB_${PROJECT_NAME}_SOURCEFILES}
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

i3_add_pybindings(simclasses
    private/pybindings/CCMMCPE.cxx
    private/pybindings/PhotonSummary.cxx
    private/pybindings/module.cxx
    USE_TOOLS boost python
    USE_PROJECTS simclasses icetray tableio
)

#i3_test_scripts(resources/test/*.py)
