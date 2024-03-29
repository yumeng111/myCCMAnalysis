set(PROJECT_NAME nearline-reco)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/nearline_reco
)

i3_add_library(${PROJECT_NAME}
    private/nearline-reco/IntervalChargeSum.cxx
    private/nearline-reco/EventFinder.cxx
    private/nearline-reco/CountPulses.cxx
    private/nearline-reco/BaselineCheck.cxx
    IWYU
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses
)
