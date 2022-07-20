/*------------------------------------------------------

  CCMModuleTable

  Based on the TPCModuleTable from the NIFFTE Experiment

  Table to hold registered modules for automatic
  registration in a reconstruction job

  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMMODULETABLE_H
#define CCMMODULETABLE_H

#include <map>
#include <string>

#include "CCMAnalysis/modules/framework/CCMModule.h"
#include "CCMAnalysis/utils/MsgLog.h"

class CCMConfig;
class CCMModule;

class CCMModuleTable 
{

public:

  static CCMModuleTable& Instance();
  void   Insert(const char* name, ModuleMaker_t mm);
  ModuleMaker_t Lookup(const char* name);

private:
  typedef std::map<std::string,ModuleMaker_t> ModMakerMap_t;

  CCMModuleTable() {};              //constructor
  ModMakerMap_t          fMakerMap; ///< Map module makers by name
  static CCMModuleTable* fInstance; ///< The sole instance of this class

};

////////////////////////////////////////////////////////////////////////
// MODULE_DECL is a macro to declare the existence of a module to a
// program
//
// Modules *MUST* provide a constructor of the form:
//
// MyModule::MyModule(const char* version) : CCMModule("MyModule", version) 
// {
// ...
// }
//
// where "version" is the version tag of the set of configuration
// parameters used to setup the module.
//
#ifndef __CINT__
#define MODULE_DECL(NAME) \
 static CCMModule* gs_##NAME##_Maker(const char* v) { return new NAME(v); } \
 class _##NAME##_MakerHelper { \
 public: \
   _##NAME##_MakerHelper() { \
     CCMModuleTable::Instance().Insert(#NAME, gs_##NAME##_Maker); \
   } \
 }; \
 static _##NAME##_MakerHelper gs_##NAME##_MakerHelper;
#endif //__CINT__
//
////////////////////////////////////////////////////////////////////////

#endif //CCMMODULE_H
