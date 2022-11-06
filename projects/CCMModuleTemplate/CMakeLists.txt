set(PROJECT_NAME CCMModuleTemplate)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMTemplateAna.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMFramework
    ROOTCINT
    public/CCMAnalysis/CCMModuleTemplate/CCMTemplateAna.h
)
