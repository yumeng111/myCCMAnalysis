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
#include "CCMAnalysis/modules/framework/CCMConfigParam.h"

//------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const CCMConfigParam& p)
{
  p.Print(os);
  return os;

}

//------------------------------------------------------
const char* CCMConfigParam::XMLTag() const
{
  //Return a recommended xml tag for this data type

  static const std::type_info& bool_id   = typeid(bool);
  static const std::type_info& char_id   = typeid(char);
  static const std::type_info& short_id  = typeid(short);
  static const std::type_info& int_id     =typeid(int);
  static const std::type_info& uchar_id   =typeid(unsigned char);
  static const std::type_info& ushort_id  =typeid(unsigned short);
  static const std::type_info& uint_id    =typeid(unsigned int);
  static const std::type_info& float_id   =typeid(float);
  static const std::type_info& double_id  =typeid(double);
  static const std::type_info& string_id  =typeid(std::string);
  static const std::type_info& boolv_id   =typeid(std::vector<bool>);
  static const std::type_info& charv_id   =typeid(std::vector<char>);
  static const std::type_info& shortv_id  =typeid(std::vector<short>);
  static const std::type_info& intv_id    =typeid(std::vector<int>);
  static const std::type_info& ucharv_id  =typeid(std::vector<unsigned char>);
  static const std::type_info& ushortv_id =typeid(std::vector<unsigned short>);
  static const std::type_info& uintv_id   =typeid(std::vector<unsigned int>);
  static const std::type_info& floatv_id  =typeid(std::vector<float>);
  static const std::type_info& doublev_id =typeid(std::vector<double>);
  static const std::type_info& stringv_id =typeid(std::vector<std::string>);

  if (bool_id   == this->fDataType()) return "bool";
  if (char_id   == this->fDataType()) return "char";
  if (short_id  == this->fDataType()) return "short";
  if (int_id    == this->fDataType()) return "int";
  if (uint_id   == this->fDataType()) return "uint";
  if (uchar_id  == this->fDataType()) return "uchar";
  if (ushort_id == this->fDataType()) return "ushort";
  if (uint_id   == this->fDataType()) return "uint";
  if (float_id  == this->fDataType()) return "float";
  if (double_id == this->fDataType()) return "double";
  if (string_id == this->fDataType()) return "string";

  if (boolv_id   == this->fDataType()) return "bool";
  if (charv_id   == this->fDataType()) return "char";
  if (shortv_id  == this->fDataType()) return "short";
  if (intv_id    == this->fDataType()) return "int";
  if (uintv_id   == this->fDataType()) return "uint";
  if (ucharv_id  == this->fDataType()) return "uchar";
  if (ushortv_id == this->fDataType()) return "ushort";
  if (uintv_id   == this->fDataType()) return "uint";
  if (floatv_id  == this->fDataType()) return "float";
  if (doublev_id == this->fDataType()) return "double";
  if (stringv_id == this->fDataType()) return "string";
  
  return "?";

}

//------------------------------------------------------
CCMConfigParam::CCMConfigParam()
  : fName         ( "<null>" ),
    fComment      ( "" ),
    fData         ( 0 ),
    fDataType     ( 0 ),
    fDataDelete   ( 0 ),
    fDataCopyCons ( 0 ),
    fDataCopy     ( 0 ),
    fDataPrint    ( 0 )
{ 
  //Default constructor

}

//------------------------------------------------------
CCMConfigParam::CCMConfigParam(const CCMConfigParam& p)
  : fName         ( p.fName    ),
    fComment      ( p.fComment ),
    fData         ( 0 ),
    fDataType     ( p.fDataType     ),
    fDataDelete   ( p.fDataDelete   ),
    fDataCopyCons ( p.fDataCopyCons ),
    fDataCopy     ( p.fDataCopy     ),
    fDataPrint    ( p.fDataPrint    )
{
  //copy constructor

  fDataCopyCons(&fData, p.fData);

}

//------------------------------------------------------
CCMConfigParam::CCMConfigParam(const char* name, const char* comment)
  : fName         ( name ),
    fComment      ( comment ),
    fData         ( 0 ),
    fDataType     ( 0 ),
    fDataDelete   ( 0 ),
    fDataCopyCons ( 0 ),
    fDataCopy     ( 0 ),
    fDataPrint    ( 0 )
{
  //constructor

}

//------------------------------------------------------
CCMConfigParam::~CCMConfigParam()
{
  //destructor

  if (fData) {
    fDataDelete(fData); 
    fData = 0;
  }

}

//------------------------------------------------------
CCMConfigParam& CCMConfigParam::operator=(const CCMConfigParam& rhs)
{
  //copy operation

  // Avoid self-copy
  if (&rhs == this) return *this;
  
  // Copy the name and comments over
  fName    = rhs.fName;
  fComment = rhs.fComment;
  
  // Copy the data handlers
  fDataType     = rhs.fDataType;
  fDataDelete   = rhs.fDataDelete;
  fDataCopyCons = rhs.fDataCopyCons;
  fDataCopy     = rhs.fDataCopy;
  fDataPrint    = rhs.fDataPrint;
  
  // Delete the data from the left-hand side and reallocate it using
  // the copy constructor
  if (fData) { fDataDelete(fData); fData = 0; }
  fDataCopyCons(&fData, rhs.fData);

  return *this;

}
