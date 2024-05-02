set(PROJECT_NAME photon-propagation)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/${PROJECT_NAME}
    )

i3_add_library(${PROJECT_NAME}
    private/photon-propagation/*.cxx
    IWYU
    USE_TOOLS boost python nlopt
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services
    )

i3_add_pybindings(${PROJECT_NAME}
    private/pybindings/module.cxx
    private/pybindings/*.cxx
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses phys-services photon-propagation
    )
