set(MODULENAME CCMIO)

set(LOCAL_HEADERS
    CCMRootIO.h
    CCMRawIO.h
    CCMEventTreeHandle.h
    )

set(LOCAL_SRCS
    CCMRootIO.cxx
    CCMRawIO.cxx
    CCMEventTreeHandle.cxx
    )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

