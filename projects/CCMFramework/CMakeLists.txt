set(PROJECT_NAME CCMFramework)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMTaskManager.cxx
    private/CCMTaskConfig.cxx
    private/CCMConfig.cxx
    private/CCMConfigParam.cxx
    private/CCMConfigTable.cxx
    private/CCMModule.cxx
    private/CCMModuleTable.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMIO CCMDataStructures
    ROOTCINT
    public/CCMAnalysis/CCMFramework/CCMTaskManager.h
    public/CCMAnalysis/CCMFramework/CCMTaskConfig.h
    public/CCMAnalysis/CCMFramework/CCMConfig.h
    public/CCMAnalysis/CCMFramework/CCMConfigParam.h
    public/CCMAnalysis/CCMFramework/CCMConfigTable.h
    public/CCMAnalysis/CCMFramework/CCMModule.h
    public/CCMAnalysis/CCMFramework/CCMModuleTable.h
)
