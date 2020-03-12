/*------------------------------------------------------

  CCMConfigParam

  Based on TPCConfigParam from the NIFFTE Experiment which is
  based on CfgParam from the E907/MIPP Experiment
  
  Stores a single CCM reconstruction parameter

  A parameter consists of a name, explanation ("comment")
  and a collection of data values, all of which must be
  of a single type

  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMCONFIGPARAM_H
#define CCMCONFIGPARAM_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>

#include "MsgLog.h"

class CCMConfigParam
{

public:
  CCMConfigParam();  //default constructor
  CCMConfigParam(const CCMConfigParam& p); //copy constructor
  ~CCMConfigParam();

  CCMConfigParam(const char* name, const char* comment);
  template <class T>
  CCMConfigParam(const char* name, const T& d, const char* comment);

  //N.B. Do we need this?  Can we use MsgLog?
  friend std::ostream& operator<<(std::ostream& os, const CCMConfigParam& p);

  const char* Name() const { return fName.c_str(); }
  const char* Comment() const { return fComment.c_str(); }
  const char* XMLTag() const;
  void Print(std::ostream& os) const { fDataPrint(os, fData); }

  template <class T> void Get(T& t) const;
  template <class T> void Set(const T& t);
  CCMConfigParam& operator=(const CCMConfigParam& rhs);

private:
  std::string fName;    //!< Name of the parameter
  std::string fComment; //!< Explanation of what parameter is used for
  void* fData;          //!< Data values held by parameter

  // The following methods allow one to correctly treat the void*
  // fData as a class of correct type
  const std::type_info& (*fDataType)(void);
  void                  (*fDataDelete)(void*);
  void                  (*fDataCopyCons)(void** dest, const void* src);
  void                  (*fDataCopy)(void* dest, const void* src);
  void                  (*fDataPrint)(std::ostream& os, const void* d);

};

////////////////////////////////////////////////////////////////////////
// Inlined templated functions
////////////////////////////////////////////////////////////////////////

//------------------------------------------------------
template <class T>
static const std::type_info& gsDataType(void)
{
  //Return the type info for data of type T
  return typeid(T);
}

//------------------------------------------------------
template <class T>
static void gsDataDelete(void* p)
{
  //Delete the data pointed to by p assuming its of type T

  T* pc = (T*)p;
  delete pc;

}

//------------------------------------------------------
template <class T>
static void gsDataCopyCons(void** dest, const void* src)
{
  //Construct a new object of type T at location (*dest) given object
  //of same type at location src

  T* srcc = (T*)src;
  T* newt = new T(*srcc);
  *dest   = (void*)newt;

}

//------------------------------------------------------
template <class T>
static void gsDataCopy(void* dest, const void* src)
{
  //Copy data from location src to location dest assuming that dest
  //and src point to objects of type T

  T* destc = (T*)dest;
  T* srcc  = (T*)src;
  (*destc) = (*srcc);

}

//------------------------------------------------------
template <class T>
static void gsDataPrint(std::ostream& os, const void* p)
{
  //Print the data located at p assuming its of type T

  T* pc = (T*)p;
  os << (*pc);

}

//------------------------------------------------------
template <class T>
static std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  //Print a vector of basic types

  int sz   = v.size();
  int szm1 = sz-1;
  if ( (typeid(T) == typeid(std::string)) || 
       (typeid(T) == typeid(const char*)) ) {
    // Add quotes if data is of string type
    for (int i=0; i<sz; ++i) {
      if (i==szm1) os << "\"" << v[i] << "\"";
      else         os << "\"" << v[i] << "\" ";
    }
  }
  else {
    // Non-string data gets streamed normally
    for (int i=0; i<sz; ++i) {
      if (i==szm1) os << v[i];
      else         os << v[i] << " ";
    }
  }
  return os;

}

//------------------------------------------------------
template <class T>
CCMConfigParam::CCMConfigParam(const char* name, const T& d, const char* comment)
  : fName        ( name    ),
    fComment     ( comment ),
    fData        ( 0 ),
    fDataType    ( 0 ),
    fDataDelete  ( 0 ),
    fDataCopyCons( 0 ),
    fDataCopy    ( 0 ),
    fDataPrint   ( 0 ) 
{
  //Construct a parameter from a name, comment, and data

  this->Set(d);

}

//------------------------------------------------------
template <class T>
void CCMConfigParam::Get(T& t) const
{
  //Retrieve the data held by the parameter into t

  //Check that types match
  const std::type_info& ti = typeid(T);
  if (ti != this->fDataType()) {
    MsgFatal(MsgLog::Form("Attempt to get data of type %s from parameter %s which is of type %s",
			  ti.name(),this->Name(),this->fDataType().name()));
  }

  // Normal case: types match so just make a copy to t
  t = (*(T*)fData);

}

//------------------------------------------------------
template <class T>
void CCMConfigParam::Set(const T& t)
{
  //Set the value of the data held by the parameter to t

  // Delete the data currently held
  if (fData) {
    this->fDataDelete(fData);
    fData = 0;
  }
  
  // Setup the methods for handling this datum
  fDataType     = gsDataType<T>;
  fDataDelete   = gsDataDelete<T>;
  fDataCopyCons = gsDataCopyCons<T>;
  fDataCopy     = gsDataCopy<T>;
  fDataPrint    = gsDataPrint<T>;

  // Create the data from the data passed in
  fDataCopyCons(&fData, &t);
}

#endif // CCMCONFIGPARAM_H
