set(PROJECT_NAME CCMBinary)

i3_project(${PROJECT_NAME}
    PYTHON_DIR python
    )

i3_add_library(${PROJECT_NAME}
    private/BinaryFormat.cxx
    private/CCMBinaryReader.cxx
    IWYU
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses
)

i3_add_pybindings(CCMBinary
    private/pybindings/module.cxx
    USE_TOOLS boost python
    USE_PROJECTS icetray dataio CCMBinary
)
