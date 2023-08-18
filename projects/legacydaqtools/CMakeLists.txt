set(PROJECT_NAME legacydaqtools)

i3_project(${PROJECT_NAME}
)

i3_add_library(${PROJECT_NAME}
    private/legacydaqtools/CCMDerivativeFinderRootOutput.cxx
    IWYU
    USE_TOOLS boost python root
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services
)
