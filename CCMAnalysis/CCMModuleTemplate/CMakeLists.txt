set(PROJECT_NAME CCMModuleTemplate)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/CCMTemplateAna.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMFramework
    ROOTCINT
    public/CCMAnalysis/CCMModuleTemplate/CCMTemplateAna.h
)
