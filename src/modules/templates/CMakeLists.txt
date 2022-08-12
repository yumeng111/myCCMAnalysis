set(MODULENAME CCMTemplates)

set(LOCAL_HEADERS 
  CCMTemplateAna.h
  )

set(LOCAL_SRCS
  CCMTemplateAna.cxx
  )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

