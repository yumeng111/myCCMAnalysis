set(PROJECT_NAME marley)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/marley
)

i3_add_library(${PROJECT_NAME}
    private/marley/MarleySimulator.cxx
    IWYU
    USE_TOOLS boost python marley
    USE_PROJECTS icetray serialization dataio dataclasses phys-services CCMBinary
)
