set(MODULENAME CCMXML)

set(LOCAL_HEADERS
  CCMTaskManager.h
  CCMTaskConfig.h
  CCMConfig.h
  CCMConfigParam.h
  CCMConfigTable.h
  CCMModule.h
  CCMModuleTable.h
  )

set(LOCAL_SRCS
  CCMTaskManager.cxx
  CCMTaskConfig.cxx
  CCMConfig.cxx
  CCMConfigParam.cxx
  CCMConfigTable.cxx
  CCMModule.cxx
  CCMModuleTable.cxx
  )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

