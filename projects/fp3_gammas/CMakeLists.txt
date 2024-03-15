set(PROJECT_NAME fp3_gammas)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
)

i3_add_library(${PROJECT_NAME}
    private/fp3_gammas/FindGammas.cxx
    IWYU
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary daqtools
)
