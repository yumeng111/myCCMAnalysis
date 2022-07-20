/*!************************************************************************************************
 * \file PMTInfoMap.h
 * \brief Header file of the #PMTInfoMap class
 * \author R. T. Thornton (LANL)
 * \date February 24, 2020
 **************************************************************************************************/
#ifndef PMTInfoMap_h
#define PMTInfoMap_h

#include <map>
#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>

class TTree;
class PMTInformation;


/*!************************************************************************************************
 * \class PMTInfoMap
 * \brief Contains the map of #PMTInformation
 *
 * An indeividual PMT's information is stored in the #PMTInformation class.
 * #PMTInfoMap is a class containing all the #PMTInformation's that are created.
 * The class also has the ability to read in a csv file containing all the necessary information,
 * write the information to a TFile, and create and decode the unique id of the PMT.
 **************************************************************************************************/
class PMTInfoMap
{
  public:
    static bool IsActive(int key);
    static bool IsActive(int digitizer, int channel);
    static int CreateKey(int digit, int channel);
    template<typename T, typename U>
    static void DecodeKey(T key, U & digit, U & channel);
    static int ConvertHVBoardChanToKey(const int & box, const int & channel, bool veto = false);
    static void ConvertKeyToHVBoardChan(const int key, int & box, int & channel);
    static int ConvertRowColToKey(const int & row, const int & col);
    static void ConvertKeyToRowCol(const int key, int & row, int & col);

    static const PMTInformation * GetPMTInfo(size_t key);
    static const PMTInformation * GetPMTInfo(int board, int channel);

    static void DefaultFillPMTMap();
    static void FillPMTMap(std::istream& file);
    static void FillPMTMap(TTree * tree);
    static void LoadHVOffList(std::string fileName);
    static void LoadCalibrationFile(std::string fileName, bool fixedThreshold = true, double threshValue = 5.0, double maxValue = 40.0);

    static void WritePMTMap(TTree *& tree);

    static size_t GetMaxKey() { return fgMaxKey; }
    static size_t GetMinKey() { return fgMinKey; }

    static std::string TreeName() { return fgkTreeName; }
    static std::string BranchName() { return fgkBranchName; }

    static void ClearMap();

    static void SetParameter(std::string name, const int value);
    static void SetParameter(std::string name, const double value);
    static void SetParameter(std::string name, std::string value);

    static int GetEJStart() {return fgEJStart;};
    static int GetEJEnd() {return fgEJEnd;};

  private:

    static const std::string fgkTreeName; ///< Tree name for map
    static const std::string fgkBranchName; ///< Branch name for map
    static std::map<int,PMTInformation*> fgPMTInfo; ///< map to all the individual #PMTInformation
    static std::vector<int> fgHVOffList;

    static int fgEJStart;
    static int fgEJEnd;

    static size_t fgMaxKey;
    static size_t fgMinKey;

    static bool fgMapLoaded;
    static bool fgBadListLoaded;

};

#include "CCMAnalysis/utils/PMTInfoMap.tcc"

#endif // PMTInfoMap_h
