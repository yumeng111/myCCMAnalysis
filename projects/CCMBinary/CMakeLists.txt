set(PROJECT_NAME CCMBinary)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/BinaryFormat.cxx
    USE_TOOLS boost
    USE_PROJECTS icetray serialization
)
