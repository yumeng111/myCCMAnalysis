set(PROJECT_NAME analytic-light-yields)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/analytic_light_yields 
    )

i3_add_library(${PROJECT_NAME}
    private/analytic-light-yields/SodiumVertexDistribution.cxx
    private/analytic-light-yields/YieldsPerPMT.cxx
    private/analytic-light-yields/GenerateExpectation.cxx
    IWYU
    USE_TOOLS boost python nlopt
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services
    )

i3_add_pybindings(analytic_light_yields
    private/pybindings/module.cxx
    private/pybindings/*.cxx
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses phys-services analytic-light-yields
    )
