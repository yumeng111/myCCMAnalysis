set(MODULENAME CCMPlots)

set(LOCAL_HEADERS
    CCMSumWaveforms.h
    CCMProcessingPlots.h
    CCMNuDMPlots.h
    )

set(LOCAL_SRCS
    CCMSumWaveforms.cxx
    CCMProcessingPlots.cxx
    CCMNuDMPlots.cxx
    )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

