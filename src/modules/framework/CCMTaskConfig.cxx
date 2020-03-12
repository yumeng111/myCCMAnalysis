/*------------------------------------------------------

  CCMTaskConfig

  Based on the TPCTaskConfig from the NIFFTE Experiment


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <fstream>

#include "RConfigure.h"

//SciBath resources
#include "CCMTaskConfig.h"
#include "MsgLog.h"
#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMConfigTable.h"
#include "Utility.h"
//XML Parsing
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include "DOMCount.hpp"

XERCES_CPP_NAMESPACE_USE

//------------------------------------------------------
CCMTaskConfig::CCMTaskConfig()
{
  //default constructor
  

}

//------------------------------------------------------
CCMTaskConfig::CCMTaskConfig(const CCMTaskConfig& reco)
{
 //copy constructor

}

//------------------------------------------------------
CCMTaskConfig::~CCMTaskConfig() 
{ 
  //destructor
  
}

//------------------------------------------------------
CCMTaskConfig::CCMTaskConfig(std::string configfile,std::vector<std::string> infileList, std::vector<std::string> outfileList)
  : fConfigFile(configfile),fInputFileList(infileList),fOutputFileList(outfileList)
{
  //constructor
  int status = ReadConfigFile();
  if(status == kCCMSuccess) {
    MsgInfo("Configuring with following inputs:");
    Print();
  } else {
    if(status == kCCMWarning) MsgWarning("Problem with configuration.  Proceeding anyway.");
    if(status == kCCMError) { MsgError("Problem with configuration.  Exiting."); exit(0); }
    if(status == kCCMFailure) MsgFatal("Problem with configuration.  Aborting.");
  }

}

//------------------------------------------------------
CCMTaskConfig::CCMTaskConfig(std::string configfile)
  : fConfigFile(configfile),fInputFileList(0),fOutputFileList(0)
{
  //constructor
  int status = ReadConfigFile();
  if(status == kCCMSuccess) {
    MsgInfo("Configuring with following inputs:");
    Print();
  } else {
    if(status == kCCMWarning) MsgWarning("Problem with configuration.  Proceeding anyway.");
    if(status == kCCMError) { MsgError("Problem with configuration.  Exiting."); exit(0); }
    if(status == kCCMFailure) MsgFatal("Problem with configuration.  Aborting.");
  }

}

//------------------------------------------------------
int CCMTaskConfig::ReadConfigFile()
{

  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    MsgError(MsgLog::Form("Error during XML initialization! : %s \n",message));
    XMLString::release(&message);
    return kCCMError;
  }

  XercesDOMParser* parser = new XercesDOMParser();
  parser->setValidationScheme(XercesDOMParser::Val_Always);    
  parser->setDoNamespaces(true);    // optional

  ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
  parser->setErrorHandler(errHandler);

  try {
    parser->parse(fConfigFile.c_str());
  }
  catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    MsgError(MsgLog::Form("XML Error during parsing: %s \n %s \n",fConfigFile.c_str(),message));
    XMLString::release(&message);
    return kCCMError;
  }
  catch (const DOMException& toCatch) {
    char* message = XMLString::transcode(toCatch.msg);
    MsgError(MsgLog::Form("DOM Error during parsing: %s \n %s \n",fConfigFile.c_str(),message));
    XMLString::release(&message);
    return kCCMError;
  }
  catch (...) {
    MsgError(MsgLog::Form("Unexpected Exception during parsing %s\n",fConfigFile.c_str()));
    return kCCMError;
  }

  // Need Analysis Manager, because the xml file might contain information about cuts
  //CCMAnalysisManager * anaMan = CCMAnalysisManager::GetInstance();

  //We got the parser, now let's get the info out of it
  DOMDocument *doc = parser->getDocument();

  if(doc) {
    XMLCh xa[1000];
    char* pEnd;

    //get the CCMConfigTable instance so we can add configurations to it
    CCMConfigTable& cfgTable = CCMConfigTable::Instance();

    //Modules has attributes and a list of items
    XMLString::transcode("Modules",xa,99);
    DOMNodeList* xmlModList = doc->getElementsByTagName(xa);
    for(unsigned int i = 0; i < xmlModList->getLength(); i++) {
      DOMElement* modEle = (DOMElement*)(xmlModList->item(i));
      //      std::string mods = StrX(modEle->getChildNodes()->item(0)->getNodeValue()).localForm();
      char* modchar = XMLString::transcode(modEle->getChildNodes()->item(0)->getNodeValue());
      /*int nmods =*/ Split(modchar," ,\n\t\r",fModuleList);
      XMLString::release(&modchar);
      XMLCh *type = XMLString::transcode("type");
      fProcessType = StrX(modEle->getAttribute(type)).localForm();
      XMLString::release(&type);
      if(!(fProcessType == "DetSim" || fProcessType == "EventDisplay" || fProcessType == "Reco" ||
            fProcessType == "Calibration"))
        MsgFatal(MsgLog::Form("Unrecognized process type %s.  Check your XML config file.",fProcessType.c_str()));
    }

    //Config has attributes and an unknown number of parameters
    XMLString::transcode("Config",xa,99);
    DOMNodeList* xmlCfgList = doc->getElementsByTagName(xa);
    for(unsigned int icfg = 0; icfg < xmlCfgList->getLength(); icfg++) {
      DOMElement* cfgEle = (DOMElement*)(xmlCfgList->item(icfg));
      XMLCh *tag = XMLString::transcode("name");
      std::string cfgname = StrX(cfgEle->getAttribute(tag)).localForm();
      XMLString::release(&tag);
      tag = XMLString::transcode("version");
      std::string cfgver = StrX(cfgEle->getAttribute(tag)).localForm();
      XMLString::release(&tag);

      //make a CCMConfig object for this configuration
      CCMConfig* cfgObj = new CCMConfig(cfgname.c_str(),cfgver.c_str(),fConfigFile.c_str());

      //Param has attributes, a child node of a certain type with a
      //cut value and a comment
      XMLString::transcode("Param",xa,99);
      DOMNodeList* paramlist = ((DOMElement*)cfgEle)->getElementsByTagName(xa);
      for (unsigned int iparam = 0; iparam < paramlist->getLength(); iparam++) {
        DOMElement* parEle = (DOMElement*)(paramlist->item(iparam));
        XMLCh* tag = XMLString::transcode("name");
        std::string parname = StrX(parEle->getAttribute(tag)).localForm();
        XMLString::release(&tag);
        DOMElement* cutEle = (DOMElement*)(parEle->getChildNodes()->item(0));
        std::string cuttype = StrX(cutEle->getTagName()).localForm();
        std::string cutval = StrX(cutEle->getChildNodes()->item(0)->getNodeValue()).localForm();
        std::string comment = "empty";
        if(parEle->getChildNodes()->getLength() > 1) {
          comment = StrX(parEle->getChildNodes()->item(1)->getNodeValue()).localForm();
        }
        //Check that this is one of our allowed types, and the value
        //matches.  If not, exit
        //	CheckTypeValue(cuttype,cutval);

        //Now that we have the cut, add it to the CCMConfig object
        if(cuttype != "bool") {

          if(cuttype == "int") {
            cfgObj->AdoptParam(new CCMConfigParam(parname.c_str(),atoi(cutval.c_str()),comment.c_str()));	
            //if(cfgname.find("Setup") != std::string::npos)
            //  anaMan->SetParameter(parname,atoi(cutval.c_str()));
          }
          if(cuttype == "double" || cuttype == "float") {
            cfgObj->AdoptParam(new CCMConfigParam(parname.c_str(),strtod(cutval.c_str(),&pEnd),comment.c_str()));	
            //if(cfgname.find("Setup") != std::string::npos)
            //  anaMan->SetParameter(parname,strtod(cutval.c_str(), &pEnd));
          }
          if(cuttype == "string") {	    
            //Convert the string to uppercase
            if(parname.find("File") == std::string::npos)
              std::transform(cutval.begin(), cutval.end(), cutval.begin(), toupper);
            cfgObj->AdoptParam(new CCMConfigParam(parname.c_str(),cutval,comment.c_str()));
            //if(cfgname.find("Setup") != std::string::npos)
            //  anaMan->SetParameter(parname,cutval);
          }
        }
      }
      //add the config object to the table
      cfgTable.AdoptConfig(cfgObj);
    }
    if(MsgLog::GetGlobalDebugLevel() > 3) {
      cfgTable.Print();	
    }
  } else {
    MsgError(MsgLog::Form("Error getting XML info out of %s",fConfigFile.c_str()));
    return kCCMError;
  }

  delete parser;
  delete errHandler;

  XMLPlatformUtils::Terminate();

  return kCCMSuccess;

}

//------------------------------------------------------
void CCMTaskConfig::Print()
{
  //Use MsgLog to print out info about the configuration

  MsgInfo(MsgLog::Form("\tConfig File: %s",fConfigFile.c_str()));
  for(unsigned int i = 0; i < fInputFileList.size(); i++) {
    MsgInfo(MsgLog::Form("\tInput file[%d]: %s",i,fInputFileList[i].c_str()));
  }
  MsgInfo(MsgLog::Form("\tJob type: %s",fProcessType.c_str()));
  MsgInfo(MsgLog::Form("\tModules to run:"));
  for(unsigned int i = 0; i< fModuleList.size(); i++) {
    MsgInfo(MsgLog::Form("\t\t %s",fModuleList.at(i).c_str()));
  }
  for(unsigned int i = 0; i < fOutputFileList.size(); i++) {
    MsgInfo(MsgLog::Form("\tOutput file[%d]: %s",i,fOutputFileList[i].c_str()));
  }

}

//------------------------------------------------------
int CCMTaskConfig::Split(const char* line, const char* tok, 
    std::vector<std::string>& fields)
{
  //Split the text line "line" into sub-strings delimited by the
  //character tokens listed in the string "tok".  Returns the number
  //of fields added to the vector "fields"
  // Thanks MippXML/MXMLString class!

  bool        inDQuotedStr = false;
  bool        inSQuotedStr = false;
  bool        specialChar  = false;
  const char* c            = line;
  int         nfields      = 0;
  std::string f;

  for (; *c!='\0'; ++c) {
    int         istoken = 0;
    const char* t;

    // If last character was '\' interpret this character literally
    if (specialChar)  { f+= *c; specialChar = false; continue; }

    // If this character is a '\' raise the special character flag
    if (*c == '\\') { specialChar = true; continue; }

    // Determine if we're entering or leaving a quoted string
    if (*c == '\"') { inDQuotedStr = !inDQuotedStr; continue; }
    if (*c == '\'') { inSQuotedStr = !inSQuotedStr; continue; }

    // In quoted string just copy everything until the close quotes
    if (inDQuotedStr)  { f+= *c; continue; }
    if (inSQuotedStr)  { f+= *c; continue; }

    // Determine if this character is one of the deliminators
    for (t=tok; *t!='\0'; ++t) if (*t==*c) ++istoken;

    // Handle non-token charaters
    if (istoken == 0) { f+= *c; continue; }

    // Handle token charaters
    if (istoken>0) {
      if (f!="") { // Don't allow empty strings
        fields.push_back(f);
        ++nfields;
        f = "";
      }
    }
  }
  // Take care of the last one
  if (f!="") {
    fields.push_back(f);
    ++nfields;
    f = "";
  }
  if (inSQuotedStr || inDQuotedStr) {
    MsgFatal(MsgLog::Form("Unterminated string: '%s'",line));
  }
  return nfields;
}
