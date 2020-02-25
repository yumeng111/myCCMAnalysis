/*!**********************************************
 * \file NearlineSPEDiag.h
 * \brief Code that holds the functions for determining the SPE values
 * \author T.J. Schuab and R.T. Thornton
 * \date February 24, 2020
 ***********************************************/
#ifndef NearlineSPEDiag_h
#define NearlineSPEDiag_h

class TChain;
class Pulses;
class TH1D;
class SPECalibrationVariables;

#include <vector>
#include <memory>

/*!**********************************************
 * \class NearlineSPEDiag
 * \brief Contains functions and containers for the SPE values
 *
 * Code that reads in the data and calculates the SPE values
 ***********************************************/
class NearlineSPEDiag{
	private:
    /// Vector of #SPECalibrationVariables in the region of interest
		std::vector<std::shared_ptr<SPECalibrationVariables>> fCal;
    /// Vector of #SPECalibrationVariables prior to the region of interest 
    /// to determine background in the region of interest
		std::vector<std::shared_ptr<SPECalibrationVariables>> fCalBefore;
    /// The size of #fCal
		size_t fSizeOfCal;
    /// The number of triggers used in the calibration
		int fNEntries;
    /// The region of interest start time in the DAQ window
		float fWindowStart;
    /// The region of interest end time in the DAQ window
		float fWindowEnd;
    /// The region of interest total time in the DAQ window
		float fWindow;

	public: 
    /// \fn NearlineSPEDiag()
    /// \brief Default constructor (does nothing)
		NearlineSPEDiag() {};
    /// \fn ~NearlineSPEDiag()
    /// \brief Deconstructor (does nothing)
		~NearlineSPEDiag() {};
		
		void CreatePEHists();
    void SetClassVec(int value);
    void MakeChainFillPulses(const std::vector<std::string> & fileList, bool ledRunFlag, bool numLEDTriggers = -1);
    void GetHistsToAdjust(std::string fileToAdjust, bool ledRunFlag);
    void CalculateRates();
    void FitPEHists();
    void FillHist(std::string pathPrefix, std::string pdfVar,bool fillAllHistFlag,bool rate2DHistFlag,bool spe2DHistFlag,bool peHistsFlag,
        bool adc2DHistFlag,bool flangeHistsFlag,bool spe1DHistFlag,bool rate1DHistFlag,bool noisePeakHistFlag,bool noiseEndHistFlag,
        bool fitEndHistFlag,bool fitRMSHistFlag,bool noiseIntegralHistFlag,bool tailIntegralHistFlag,bool rootOnlyFlag,bool pdfOnlyFlag,
        bool rootAndPDFFlag,bool rootTreeFlag);

    double GetThreshold(int pmt);
    double GetSPE(int pmt);
    double GetRate(int pmt);

};
#endif

