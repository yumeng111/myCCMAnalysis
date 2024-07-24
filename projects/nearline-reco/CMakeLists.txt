set(PROJECT_NAME nearline-reco)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/nearline_reco
)

i3_add_library(${PROJECT_NAME}
    private/nearline-reco/IntervalChargeSum.cxx
    private/nearline-reco/EventFinder.cxx
    private/nearline-reco/PulseChargeHistogram.cxx
    private/nearline-reco/BaselineCheck.cxx
    private/nearline-reco/CCMEventSplitter.cxx
    private/nearline-reco/PulseChargeTimeHistogram.cxx
    private/nearline-reco/PulseCorrelations.cxx
    private/nearline-reco/NHit.cxx
    private/nearline-reco/LargestPMTFraction.cxx
    private/nearline-reco/BasicVertexReco.cxx
    private/nearline-reco/ChargeUnfolding.cxx
    private/nearline-reco/ComputeChargeContamination.cxx
    IWYU
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses phys-services CCMBinary
)
