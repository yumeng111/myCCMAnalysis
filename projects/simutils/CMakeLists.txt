set(PROJECT_NAME simutils)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/simutils
    )

i3_add_library(${PROJECT_NAME}
    private/simutils/PMTResponse.cxx
    private/simutils/OverlayNoise.cxx
    IWYU
    USE_TOOLS boost python nlopt
    USE_PROJECTS icetray serialization dataio dataclasses simclasses CCMBinary phys-services
    )


