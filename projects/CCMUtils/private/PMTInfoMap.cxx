/*!************************************************************************************************
 * \file PMTInfoMap.cxx
 * \brief Functions of the #PMTInfoMap class
 * \author R. T. Thornton (LANL)
 * \date February 24, 2020
 **************************************************************************************************/

#include <sstream>
#include <utility>
#include <algorithm>

#include "CCMAnalysis/CCMUtils/PMTInfoMap.h"
#include "CCMAnalysis/CCMUtils/PMTInformation.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/CSVrow.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

const std::string PMTInfoMap::fgkTreeName = "pmtMap"; ///< Tree name for map
const std::string PMTInfoMap::fgkBranchName = "pmtInfo"; ///< Branch name for map
std::map<int,PMTInformation*> PMTInfoMap::fgPMTInfo; ///< map to all the individual #PMTInformation
std::vector<int> PMTInfoMap::fgHVOffList;

bool PMTInfoMap::fgMapLoaded = false;
bool PMTInfoMap::fgBadListLoaded = false;

int PMTInfoMap::fgEJStart = 249;
int PMTInfoMap::fgEJEnd = 254;
size_t PMTInfoMap::fgMinKey = 0;
size_t PMTInfoMap::fgMaxKey = 240;

/*
detsim/libCCMDetSim.so: undefined reference to `PMTInfoMap::ConvertRowColToKey(int const&, int const&)'
*/

/*!************************************************************************************************
 * \fn void PMTInfoMap::DefaultFillPMTMap()
 * \brief Function that fills the PMTInformation map with default csv files
 **************************************************************************************************/
void PMTInfoMap::DefaultFillPMTMap()
{
  if (fgMapLoaded) {
    return;
  }

  std::ifstream infile;

  
  
  std::string env = std::getenv("CCMINSTALL");
  
  
  std::string eightInchFile = env + "/mappings/2021/mapping_master_8inch.csv";
  std::string oneInchFile = env + "/mappings/2021/mapping_master_1inch.csv";
  
  std::cout << "eightInchFile: " << eightInchFile << "\n";


  MsgInfo(MsgLog::Form("Loading 8in map from %s",eightInchFile.c_str()));
  infile.open(eightInchFile.c_str());
  FillPMTMap(infile);
  infile.close();

  fgMapLoaded = false;

  MsgInfo(MsgLog::Form("Loading 1in map from %s",oneInchFile.c_str()));
  infile.open(oneInchFile.c_str());
  FillPMTMap(infile);
  infile.close();

  fgMapLoaded = true;
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::FillPMTMap(TTree * tree)
 * \brief Function that fills the PMTInformation map
 * \param[in] tree The tree that contains the PMTInformation
 **************************************************************************************************/
void PMTInfoMap::FillPMTMap(TTree * tree)
{
  //std::cout << "PMT INFO MAP" << "\n"; 
  if (fgMapLoaded) {
    //std::cout << "Map Already Loaded" << "\n"; 
    return;
  }

  PMTInformation * pmtInfo = 0;
  tree->SetBranchAddress(fgkBranchName.c_str(),&pmtInfo);

  const long kNEntries = tree->GetEntries();
  for (long entry = 0; entry < kNEntries; ++entry) {
    tree->GetEntry(entry);
    int key = CreateKey(pmtInfo->GetBoard(),pmtInfo->GetBoardChan());
    if (fgPMTInfo.find(key) == fgPMTInfo.end()) {
      if (MsgLog::GetGlobalDebugLevel() >= 1) {
        MsgDebug(1,MsgLog::Form("Added key %3d board %2d chan %2d col %3d row %3d coated %1d name %s",
              key,pmtInfo->GetBoard(),pmtInfo->GetBoardChan(),pmtInfo->GetColumn(),pmtInfo->GetRow(),
              pmtInfo->IsUncoated(),pmtInfo->GetLocName().c_str()));
      }
      fgPMTInfo.insert(std::make_pair(key,new PMTInformation(*pmtInfo)));
    }
  }

  fgMapLoaded = true;

  return;
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::FillPMTMap(std::istream& file)
 * \brief Function that fills the PMTInformation map
 * \param[in] file The csv file to loop through
 **************************************************************************************************/
void PMTInfoMap::FillPMTMap(std::istream& file)
{
  //std::cout << "FILL  MY PMT MAP" << std::endl;  
  if (fgMapLoaded) {
    return;
  }

  if (MsgLog::GetGlobalDebugLevel() >= 2) {
    MsgDebug(1,MsgLog::Form("Filling mapping from file..."));
  }

  fgEJStart = CreateKey(Utility::fgkNumDigitizers-1, 9);
  fgEJEnd = CreateKey(Utility::fgkNumDigitizers-1, 14);

  int hvBoard = 0;
  int hvChannel = 0;
  int flange = 0;
  int ring = 0;
  int ringLoc = 0;
  int adcBoardOrder = 0;
  //int adcBoard = 0; // not used
  int adcCH = 0;
  int row = 0;
  int col = 0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double adcToPE = 0.0;
  double adcToPEDer = 0.0;
  double adcThreshold = 0.0;
  bool coating = false;

  std::size_t counter = 0;
  for (CSVIterator iter(file); iter != CSVIterator(); ++iter, ++counter) {
    if (counter == 0) {
      continue;
    }
    //std::cout << "Counter: " << counter << "\n"; 
  

    std::string currentString = "";
    std::string tempString = "";
    char tempChar = ' ';
    std::size_t findLoc = 0;


    const CSVRow * currRow  = &(*iter);
    std::size_t currRowSize = currRow->size();
    for (size_t index = 0; index < currRowSize; ++index) {
        
      currentString = currRow->operator[](index);
      
      //std::cout << "Index: " << index << "\n";
      //std::cout << "STRING: " << std::stoi(currentString) << "\n";
     // std::cout << "Index: " << index << " String: " << currentString << "\n";
     
      switch(index) {
        case 0: ///Decopler Box 
          break; /// do not use this information for this peace of code
        case 1: ///HV Board
          hvBoard = std::stoi(currentString);
          break;
        case 2: ///HV CH
          hvChannel= std::stoi(currentString);
          break; /// do not use this information for this peace of code
        case 3: ///Flange
          currentString.erase(currentString.begin());
          flange = std::stoi(currentString);
          break;
        case 4: ///Flange Hole
          tempChar = currentString.front();
          switch(tempChar) {
            case 'R': ring = 1; break;
            case 'Y': ring = 2; break;
            case 'B': ring = 3; break;
            case 'G': ring = 4; break;
            case 'I': ring = 5; break;
            default: ring = -1; break;
          }
          currentString.erase(currentString.begin());
          ringLoc = std::stoi(currentString);
          break;
        case 5: ///Flange ID
          break; /// information is gotten from the previous two fields
        case 6: ///ADC Board Order
          adcBoardOrder = std::stoi(currentString);
          break;
        case 7: ///ADC Board
          // not used
          //adcBoard = std::stoi(currentString);
          break;
        case 8: ///ADC CH
          adcCH = std::stoi(currentString);
          break;
        case 9: ///PMT
          break; /// do not use this information for this peace of code
        case 10: ///HV
          break; /// do not use this information for this peace of code
        case 11: ///Position
          tempChar = currentString.front();
          currentString.erase(currentString.begin());
          if (tempChar == 'C') { /// PMT is a detector PMT
            findLoc = currentString.find("R");
            tempString.assign(currentString.begin(),currentString.begin()+findLoc);
            col = std::stoi(tempString);
            tempString.clear();
            currentString.erase(currentString.begin(),currentString.begin()+findLoc+1);
            row = std::stoi(currentString);
          } else if (tempChar == 'V') { /// PMT is a veto PMT
            tempChar = currentString.front();
            if (tempChar == 'B') { /// veto PMT on the bottom
              row = 8;
              currentString.erase(currentString.begin());
            } else if (tempChar == 'C') { /// veto PMT on the side
              currentString.erase(currentString.begin());
              tempChar = currentString.front();
              if (tempChar == 'B') { /// veto PMT is looking up
                row = 7;
              } else if (tempChar == 'T') { /// veto PMT is looking down
                row = -1;
              } else {
                row = -22;
              }
              currentString.erase(currentString.begin());
            } else if (tempChar == 'T') { /// veto PMT on the top
              row = -2;
              currentString.erase(currentString.begin());
            } else {
              row = -21;
            }
            col = std::stoi(currentString);
          } /*else if (tempChar = 'E') {
	    row = -3;
	    findLoc = currentString.find('J');
	    tempString.assign(currentString.begin(),currentString.begin()+findLoc);
            tempString.clear();
            currentString.erase(currentString.begin(),currentString.begin()+findLoc+1);
            col = std::stoi(currentString);
	    } */  else {
            row = -20; 
            col = -20;
          }
          break;
        case 12: ///Coating
          tempChar = currentString.front();
          if (tempChar == '0') {
            coating = false;
          } else {
            coating = true;
          }
          break;
        case 13: /// x position
          x = std::stod(currentString);
          break;
        case 14: /// y position
          y = std::stod(currentString);
          break;
        case 15: /// z position
          z = std::stod(currentString);
          break;
        case 16: /// adcToPE Value
          adcToPE = std::stod(currentString);
          break;
        case 17: /// adcThreshold Value
          adcThreshold = std::stod(currentString);
          break;
        case 18: /// adcToPE Value Calculated from Derivative
          adcToPEDer = std::stod(currentString);
          break;
        default: break;
      } // end switch(index)
    } // end for index < currRow->size()
    int key = CreateKey(adcBoardOrder,adcCH);
    std::cout << "key: " << key << "\n";//" threshold " << adcThreshold << " adcToPE " << adcToPE <<  " adcToPEder " << adcToPEDer << "\n";
    if (fgPMTInfo.find(key) == fgPMTInfo.end()) {
      fgPMTInfo.insert(std::make_pair(key,new PMTInformation()));
      std::map<int,PMTInformation*>::iterator itMap = fgPMTInfo.find(key);
      if (itMap == fgPMTInfo.end()) {
        MsgFatal(MsgLog::Form("pmt key %d was not added tothe map",key));
        abort();
      }

      itMap->second->SetHVBoard(hvBoard);
      itMap->second->SetHVBoardChan(hvChannel);
      itMap->second->SetBoard(adcBoardOrder);
      itMap->second->SetBoardChan(adcCH);
      itMap->second->SetFlange(flange);
      itMap->second->SetRing(ring);
      itMap->second->SetRingLoc(ringLoc);
      itMap->second->SetColumn(col);
      itMap->second->SetRow(row);
      itMap->second->SetUncoated(coating);
      itMap->second->SetPosition(TVector3(x,y,z));
      itMap->second->SetADCToPE(adcToPE);
      itMap->second->SetADCToPEDer(adcToPEDer);
      itMap->second->SetADCThreshold(adcThreshold);
      itMap->second->CreateNames();
      if (MsgLog::GetGlobalDebugLevel() >= 1) {
        MsgDebug(1,MsgLog::Form("Added key %3d hvBoard %2d hvChan %2d board %2d chan %2d col %3d row %3d coated %1d name %s",
              key,hvBoard,hvChannel,adcBoardOrder,adcCH,col,row,coating,itMap->second->GetLocName().c_str()));
      }
      if (MsgLog::GetGlobalDebugLevel() >= 2) {
        MsgDebug(1,MsgLog::Form("Added flange  %3d ring  %2d, ring loc, %2d,flange name %s  position(x,y,z) x: %f, y: %f z:  %f, adctoPE %f adcToPEder %f threshold %f",
             flange,ring,ringLoc,itMap->second->GetFlangeName().c_str(),x,y,z,adcToPE,adcToPEDer,adcThreshold));
        //std::cout << "IsVeto" << itMap->second->IsVeto() << "\n";
      }
    } // end if fgPMTInfo.find(key) == fgPMTInfo.end()

  } // end for CSVIterator

  fgMapLoaded = true;
}

/*!************************************************************************************************
 * \fn int PMTInfoMap::CreateKey(int digit, int channel)
 * \brief Static function that creates the map key
 * \param[in] digit Board number
 * \param[in] channel Channel number
 * \return key value
 **************************************************************************************************/
int PMTInfoMap::CreateKey(int digit, int channel) {
  return digit*16+channel;
}


/*!************************************************************************************************
 * \fn int PMTInfoMap::ConvertHVBoardChanToKey(const int & box, const int & channel, bool veto)
 * \brief Static function that gets the key given HV board and channel numbers
 * \param[in] box HV Board number
 * \param[in] channel HV Channel number
 * \param[in] veto Boolean: true if the pmt is a veto pmt
 * \return key value
 **************************************************************************************************/
int PMTInfoMap::ConvertHVBoardChanToKey(const int & box, const int & channel, bool veto) 
{
  for (const auto & pmtInfo : fgPMTInfo) {
    if (veto && !pmtInfo.second->Is1in()) {
      continue;
    } else if (!veto && pmtInfo.second->Is1in()) {
      continue;
    }
    if (pmtInfo.second->GetHVBoard() == box && pmtInfo.second->GetHVBoardChan() == channel) {
      return CreateKey(pmtInfo.second->GetBoard(),pmtInfo.second->GetBoardChan());
    }
  }

  return -1;
}

/*!************************************************************************************************
 * \fn int PMTInfoMap::ConvertRowColToKey(const int & row, const int & col)
 * \brief Static function that gets the key given HV board and col numbers
 * \param[in] row HV Board number
 * \param[in] col HV Channel number
 * \return key value
 **************************************************************************************************/
int PMTInfoMap::ConvertRowColToKey(const int & row, const int & col) {
  for (const auto & pmtInfo : fgPMTInfo) {
    if (pmtInfo.second->GetRow() == row && pmtInfo.second->GetColumn() == col) {
      return CreateKey(pmtInfo.second->GetBoard(),pmtInfo.second->GetBoardChan());
    }
  }
  return -1;
} 

/*!************************************************************************************************
 * \fn void PMTInfoMap::ConvertKeyToHVBoardChan(cosnt int key, int & box, int & channel)
 * \brief Static function that gets the HV board and channel numbers given the PMT key
 * \param[in] key PMT key number
 * \param[out] box HV Board number
 * \param[out] channel HV Channel number
 **************************************************************************************************/
void PMTInfoMap::ConvertKeyToHVBoardChan(const int key, int & box, int & channel) 
{
  std::map<int,PMTInformation*>::iterator itMap = fgPMTInfo.find(key);
  if (itMap == fgPMTInfo.end()) {
    box = -1;
    channel = -1;
    return;
  }

  box = itMap->second->GetHVBoard();
  channel = itMap->second->GetHVBoardChan();
  return;
}

/*!************************************************************************************************
 * \fn bool PMTInfoMap::IsActive(int key)
 * \brief returns true if the key is active
 * \param[in] key Key value
 * \return True if active
 **************************************************************************************************/
bool PMTInfoMap::IsActive(int key)
{
  std::map<int,PMTInformation*>::iterator itMap = fgPMTInfo.find(key);
  if (itMap == fgPMTInfo.end()) {
    return false;
  }

  if (std::find(fgHVOffList.begin(),fgHVOffList.end(),key) == fgHVOffList.end()) {
    return true;
  }

  return false;

  //if (itMap->second->GetADCToPEDer() <= 0) {
  //  return false;
  //}

  //if (itMap->second->GetBoard() != 9) {
  //  return false;
  //}

  //return itMap->second->GetColumn() >= 0;
  //return ((itMap->second->GetRow() >= 1) && (itMap->second->GetRow() <= 5)) || (itMap->second->GetRow() == 7);
  //return true;
}

/*!************************************************************************************************
 * \fn bool PMTInfoMap::IsActive(int digitizer, int channel)
 * \brief returns true if the channel board combo is active
 * \param[in] digitizer Board number
 * \param[in] channel Board channel number
 * \return True if active
 **************************************************************************************************/
bool PMTInfoMap::IsActive(int digitizer, int channel)
{
  return IsActive(CreateKey(digitizer,channel));
}

/*!************************************************************************************************
 * \fn const PMTInformation * PMTInfoMap::GetPMTInfo(int board, int channel)
 * \brief Returns the PMTInformation object for the key passed, nullptr otherwise
 * \param[in] board Board number
 * \param[in] channel Channel number
 * \return #PMTInformation object or nullptr
 *
 * Returns the result of calling #GetPMTInfo(size_t key)
 **************************************************************************************************/
const PMTInformation * PMTInfoMap::GetPMTInfo(int board, int channel)

{

  return GetPMTInfo(CreateKey(board,channel));
}

/*!************************************************************************************************
 * \fn const PMTInformation * PMTInfoMap::GetPMTInfo(size_t key)
 * \brief Returns the PMTInformation object for the key passed, nullptr otherwise
 * \param[in] key PMT key
 * \return #PMTInformation object or nullptr
 **************************************************************************************************/
const PMTInformation * PMTInfoMap::GetPMTInfo(size_t key)
{

//std::cout << "function 2" << "\n";
 // std::cout << "key : " << key << "\n";

  std::map<int,PMTInformation*>::const_iterator itMap = fgPMTInfo.find(key);
  
// std::cout << "itmap : " << itMap->first <<" " << itMap->second <<"\n";
  
  if (itMap == fgPMTInfo.end()) {
    return nullptr;
  }

  

  return itMap->second;
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::WritePMTMap(TTree *& tree)
 * \brief Fills the tree with the PMTInformation objects
 * \param[in,out] tree TTree to fill
 **************************************************************************************************/
void PMTInfoMap::WritePMTMap(TTree *& tree)
{
  if (gROOT->FindObject(fgkTreeName.c_str())) {
    delete gROOT->FindObject(fgkTreeName.c_str());
  }

  PMTInformation * pmtInfo = 0;

  tree = new TTree(fgkTreeName.c_str(),fgkTreeName.c_str());
  tree->Branch(fgkBranchName.c_str(),&pmtInfo);

  for (auto & mapIndex : fgPMTInfo) {
    pmtInfo = mapIndex.second;
    tree->Fill();
  }
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::ClearMap()
 * \brief Clears the #fgPMTInfo map container
 **************************************************************************************************/
void PMTInfoMap::ClearMap()
{
  for (auto & mapIndex : fgPMTInfo) {
    if (mapIndex.second) {
      delete mapIndex.second;
    }
  }

  if (!fgPMTInfo.empty()) {
    fgPMTInfo.clear();
  }

  if (!fgHVOffList.empty()) {
    fgHVOffList.clear();
  }

  fgMapLoaded = false;
  fgBadListLoaded = false;
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::LoadHVOffList(std::string fileName)
 * \brief Loads the list of PMTs to turn off
 * \param[in] fileName The name of the file to parse
 **************************************************************************************************/
void PMTInfoMap::LoadHVOffList(std::string fileName)
{
  if (fgBadListLoaded) {
    return;
  }

  if (fileName == "default" || fileName == "") {
    std::string env = std::getenv("CCMPROJECT");
    fileName = env + "/calibrationFiles/2021/hv_off_2021.csv";
  }

  MsgInfo(MsgLog::Form("hvOffList %s",fileName.c_str()));
  std::ifstream hvOffList(fileName.c_str());
  int board = 0;
  int chan = 0;
  std::string line;
  while (hvOffList.good()) {
    getline(hvOffList,line);
    if (hvOffList.eof()) {
      break;
    }
    std::stringstream ss(line);
    ss >> board >> chan;
    int key = ConvertHVBoardChanToKey(board,chan);

    if (MsgLog::GetGlobalDebugLevel() >= 1) {
      MsgDebug(1,MsgLog::Form("To Remove HV board %d HV channel %d key %d",board,chan,key));
    }

    if (key < 0) {
      continue;
    }

    if (std::find(fgHVOffList.begin(),fgHVOffList.end(),key) == fgHVOffList.end()) {
      fgHVOffList.push_back(key);
    }
  }
  hvOffList.close();

  fgBadListLoaded = true;
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::LoadCalibrationFile(std::string fileName, bool fixedThreshold, double threshValue, double maxValue)
 * \brief Loads the calibration file
 * \param[in] fileName The name of the file to parse
 * \param[in] fixedThreshold Boolean to choose between the treshold to be the same for all PMTs or to allow them to differ (default = true)
 * \param[in] threshValue The threshold value to use if \p fixedThreshold is true (default is 5)
 * \param[in] maxValue The maximum ADCtoSPE value to allow (default = 40)
 **************************************************************************************************/
void PMTInfoMap::LoadCalibrationFile(std::string fileName, bool fixedThreshold, double threshValue, double maxValue)
{
  //TFile * calibrationFileTop = TFile::Open("root_out_2019ledData_run179_legFix_integral_all_round4_.root","READ");

  if (fileName == "") {
    MsgFatal("No file name passed. Exiting");
  }

  if (fileName == "default" || fileName == "") {
    fileName = std::getenv("CCMPROJECT");
   // fileName += "/calibrationFiles/2019/root_out_2019ledData_run179_legFix_integral_all_round4_.root"
      fileName += "/calibrationFiles/2021/speCalc_run1529to1542.root";
  }

  MsgInfo(MsgLog::Form("Loading file from %s",fileName.c_str()));

  TFile * calibrationFile = TFile::Open(fileName.c_str(),"READ");
  if (!calibrationFile) {
    MsgFatal(MsgLog::Form("Could not open calibration file %s. Exiting",fileName.c_str()));
  }

  if (!calibrationFile->IsOpen()) {
    MsgFatal(MsgLog::Form("Could not open calibration file %s (IsOpen() check). Exiting",fileName.c_str()));
  }

  TTreeReader calibrationTree("spe",calibrationFile);
  TTreeReaderValue<float> speValue(calibrationTree,"speValue");
  TTreeReaderValue<float> speThreshold(calibrationTree,"endNoiseWallFitRangeStart");
  TTreeReaderValue<int> pmtID(calibrationTree,"pmtID");
  std::map<int,PMTInformation*>::iterator itMap = fgPMTInfo.begin();
  while (calibrationTree.Next()) {
    itMap = fgPMTInfo.find(*pmtID);
    if (itMap == fgPMTInfo.end()) {
      continue;
    }

    // change spe calibration value
    if (*speValue > maxValue) {
      itMap->second->SetADCToPE(maxValue);
    } else {
      itMap->second->SetADCToPE(*speValue);
    }

    // change threshold value
    if (fixedThreshold) {
      itMap->second->SetADCThreshold(threshValue);
    } else {
      itMap->second->SetADCThreshold(*speThreshold);
    }

    if (MsgLog::GetGlobalDebugLevel() >= 1) {
      MsgDebug(1,MsgLog::Form("PMT %d SPE %.2f Threshold %.2f",*pmtID,itMap->second->GetADCToPE(),itMap->second->GetADCThreshold()));
    }
    
  }
  delete calibrationFile;
}

//__________________________________________________
void PMTInfoMap::SetParameter(std::string /*name*/, const int /*value*/)
{
  // nothing yet
}

//--------------------------------------------------------------------
void PMTInfoMap::SetParameter(std::string /*name*/, const double /*value*/)
{
  // nothing yet
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::LoadThresholdFile(std::istream& file)
 * \brief Function that changes the adc thresholds to hand generated values
 * \param[in] file The csv file to loop through
 **************************************************************************************************/
void PMTInfoMap::LoadThresholdFile(std::string fileName)
{
  if (fileName == "default" || fileName == "") {
    std::string env = std::getenv("CCMPROJECT");
    fileName = env + "/calibrationFiles/2021/threshold_corrections.csv";
  }

  MsgInfo(MsgLog::Form("thresholdFile %s",fileName.c_str()));
  std::ifstream file(fileName.c_str());

  if (MsgLog::GetGlobalDebugLevel() >= 2) {
    MsgDebug(1,MsgLog::Form("Filling thresholds from file..."));
  }

  double adcThreshold = 0.0;

  std::size_t counter = 0;
  for (CSVIterator iter(file); iter != CSVIterator(); ++iter, ++counter) {
    if (counter == 0) {
      continue;
    }
    //std::cout << "Counter: " << counter << "\n";   

    std::string currentString = "";
    int key;

    const CSVRow * currRow  = &(*iter);
    std::size_t currRowSize = currRow->size();
    for (size_t index = 0; index < currRowSize; ++index) {
        
      currentString = currRow->operator[](index);
      
      switch(index) {
        case 0: ///PMT Key integer
	  key = std::stoi(currentString);
          break; 
        case 1: ///User Input Threshold from file
	  adcThreshold = std::stod(currentString);
          break;
        default: break;
      } // end switch(index)
    } // end for index < currRow->size()
    std::cout << "key: " << key << " threshold " << adcThreshold << "\n";

    std::map<int,PMTInformation*>::iterator itMap = fgPMTInfo.begin();
    itMap = fgPMTInfo.find(key);
    if (itMap == fgPMTInfo.end()) {
      continue;
    }

    // change threshold value
    itMap->second->SetADCThreshold(adcThreshold);
    
    if (MsgLog::GetGlobalDebugLevel() >= 1) {
      MsgDebug(1,MsgLog::Form("Changed key %3d threshold %f", key, adcThreshold));
      //std::cout << "IsVeto" << itMap->second->IsVeto() << "\n";
    }
  } // end for CSVIterator
}

//--------------------------------------------------------------------
//Runs PMTInfoMap setup from configFile for ./CCMAnalysis
void PMTInfoMap::SetParameter(std::string name, std::string value)
{
  if (name.compare("HVOffFile") == 0) {
    LoadHVOffList(value);
  } else if (name.compare("CalibrationFile") == 0) {
    LoadCalibrationFile(value,false);
  } else if (name.compare("MappingFile") == 0) {
    ClearMap();
    if (value.compare("default") == 0) {
      DefaultFillPMTMap();
    } else {
      std::ifstream infile(value);
      FillPMTMap(infile);
      infile.close();
    }
  } else if (name.compare("ThresholdFile") == 0) {
    LoadThresholdFile(value);
  }
}

