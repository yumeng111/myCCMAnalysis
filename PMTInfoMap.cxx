/*!************************************************************************************************
 * \file PMTInfoMap.cxx
 * \brief Functions of the #PMTInfoMap class
 * \author R. T. Thornton (LANL)
 * \date February 24, 2020
 **************************************************************************************************/
#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "MsgLog.h"
#include "CSVrow.h"

#include "TTree.h"
#include "TROOT.h"

#include <algorithm>
#include <utility>

const std::string PMTInfoMap::fgkTreeName = "pmtMap"; ///< Tree name for map
const std::string PMTInfoMap::fgkBranchName = "pmtInfo"; ///< Branch name for map
std::map<int,PMTInformation*> PMTInfoMap::fgPMTInfo; ///< map to all the individual #PMTInformation

/*!************************************************************************************************
 * \fn void PMTInfoMap::DefaultFillPMTMap()
 * \brief Function that fills the PMTInformation map with default csv files
 **************************************************************************************************/
void PMTInfoMap::DefaultFillPMTMap()
{
  std::ifstream infile;

  std::string env = std::getenv("CCMINSTALL");
  std::string eightInchFile = env + "/mapping_master_8inch.csv";
  std::string oneInchFile = env + "/mapping_master_1inch.csv";

  infile.open(eightInchFile.c_str());
  FillPMTMap(infile);
  infile.close();

  infile.open(oneInchFile.c_str());
  FillPMTMap(infile);
  infile.close();
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::FillPMTMap(TTree * tree)
 * \brief Function that fills the PMTInformation map
 * \param[in] tree The tree that contains the PMTInformation
 **************************************************************************************************/
void PMTInfoMap::FillPMTMap(TTree * tree)
{
  PMTInformation * pmtInfo = 0;
  tree->SetBranchAddress(fgkBranchName.c_str(),&pmtInfo);

  const long kNEntries = tree->GetEntries();
  for (long entry = 0; entry < kNEntries; ++entry) {
    tree->GetEntry(entry);
    int key = CreateKey(pmtInfo->GetBoard(),pmtInfo->GetBoardChan());
    if (fgPMTInfo.find(key) == fgPMTInfo.end()) {
      MsgInfo(MsgLog::Form("Added key %3d board %2d chan %2d col %3d row %3d coated %1d name %s",
            key,pmtInfo->GetBoard(),pmtInfo->GetBoardChan(),pmtInfo->GetColumn(),pmtInfo->GetRow(),
            pmtInfo->IsUncoated(),pmtInfo->GetLocName().c_str()));
      fgPMTInfo.insert(std::make_pair(key,new PMTInformation(*pmtInfo)));
    }
  }

  return;
}

/*!************************************************************************************************
 * \fn void PMTInfoMap::FillPMTMap(std::istream& file)
 * \brief Function that fills the PMTInformation map
 * \param[in] file The csv file to loop through
 **************************************************************************************************/
void PMTInfoMap::FillPMTMap(std::istream& file)
{
  int hvBoard = 0;
  int hvChannel = 0;
  int flange = 0;
  int ring = 0;
  int ringLoc = 0;
  int adcBoardOrder = 0;
  int adcBoard = 0;
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

    std::string currentString = "";
    std::string tempString = "";
    char tempChar = ' ';
    std::size_t findLoc = 0;


    const CSVRow * currRow  = &(*iter);
    std::size_t currRowSize = currRow->size();
    for (size_t index = 0; index < currRowSize; ++index) {
      currentString = currRow->operator[](index);
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
          adcBoard = std::stoi(currentString);
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
              row = 7;
              currentString.erase(currentString.begin());
            } else if (tempChar == 'C') { /// veto PMT on the side
              currentString.erase(currentString.begin());
              tempChar = currentString.front();
              if (tempChar == 'B') { /// veto PMT is looking up
                row = 6;
              } else if (tempChar == 'T') { /// veto PMT is looking down
                row = 0;
              } else {
                row = -22;
              }
              currentString.erase(currentString.begin());
            } else if (tempChar == 'T') { /// veto PMT on the top
              row = -1;
              currentString.erase(currentString.begin());
            } else {
              row = -21;
            }
            col = std::stoi(currentString);
          } else {
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
        case 14: /// x position
          y = std::stod(currentString);
          break;
        case 15: /// x position
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
      MsgInfo(MsgLog::Form("Added key %3d hvBoard %2d hvChan %2d board %2d chan %2d col %3d row %3d coated %1d name %s",
            key,hvBoard,hvChannel,adcBoardOrder,adcCH,col,row,coating,itMap->second->GetLocName().c_str()));
    } // end if fgPMTInfo.find(key) == fgPMTInfo.end()
  } // end for CSVIterator
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
 * \fn void PMTInfoMap::DecodeKey(int key, int & digit, int & channel)
 * \brief Static function that decodes the map key
 * \param[in] key The key value
 * \param[out] digit Board number
 * \param[out] channel Channel number
 **************************************************************************************************/
void PMTInfoMap::DecodeKey(int key, int & digit, int & channel) {
  digit = key/16;
  channel = key - digit*10;

  return;
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
 * \fn bool PMTInfoMap::IsActive(int digitizer, int channel)
 * \brief returns true if the channel board combo is active
 * \param[in] digitizer Board number
 * \param[in] channel Board channel number
 * \return True if active
 **************************************************************************************************/
bool PMTInfoMap::IsActive(int digitizer, int channel)
{
  int key = CreateKey(digitizer,channel);
  std::map<int,PMTInformation*>::iterator itMap = fgPMTInfo.find(key);
  if (itMap == fgPMTInfo.end()) {
    return false;
  }

  return true;

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
  std::map<int,PMTInformation*>::const_iterator itMap = fgPMTInfo.find(key);
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
}

