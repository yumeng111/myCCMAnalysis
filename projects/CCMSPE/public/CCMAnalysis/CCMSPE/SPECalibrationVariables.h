/*!**********************************************
 * \file SPECalibrationVariables.h
 * \brief Header file for the #SPECalibrationVariables class
 * \author T. J. Schuab and R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef SPECalibrationVariables_h
#define SPECalibrationVariables_h

#include<string>
#include <memory>

class TChain;
class Pulses;
class TH1D;
class TCanvas;


/*!**********************************************
 * \class SPECalibrationVariables
 * \brief Container for the information obtained while determining the SPE values
 ***********************************************/
class SPECalibrationVariables{
	public:
		SPECalibrationVariables();
	  ~SPECalibrationVariables() {};

		std::shared_ptr<TH1D>   		GetPEHistPtr() { return fPEHistPtr; }
		TH1D* 			GetPEHist() { return fPEHist; }
		TCanvas* 		GetPECanvas() { return fPECanvas; }
		int    			GetPMTID() { return fPMTID; }
		std::string 			GetPMTName() { return fPMTName; }
		int 	 			GetPMTRow() { return fPMTRow; }
		int 	 			GetPMTColumn() { return fPMTColumn; }
		int 	 			GetADCBoard() { return fADCBoard; }
		int 	 			GetADCChannel() { return fADCChannel; }
		int    			GetPMTFlange() { return fPMTFlange; }
		float  			GetSPE() { return fSPE; }
		float    			GetSPEPeakPos() { return fSPEPeakPos; }
		float  			GetIntegralForRate() { return fIntegralForRate; }
		float  			GetRate() { return fRate; }
		float  			GetPEHistMean() { return fPEHistMean; }
		float  			GetPEHistRMS() { return fPEHistRMS; }
		int    			GetPEHistEntries() { return fPEHistEntries; }
		float    			GetNoisePeakPos() { return fNoiseWallPeakPos; }
		float    			GetNoiseEndPos() { return fNoiseWallEndPos; }
		float  			GetNoiseIntegral() { return fNoiseWallIntegral; }
		float  			GetNoiseEndContent() { return fNoiseEndBinContent; } 
		float    			GetFitStartPos() { return fFitStartPos; }
		float 	 			GetFitEndPos() { return fFitEndPos; }
		float  			GetFitEndContent() { return fFitEndContent; }
		float  			GetFitHeight() { return fFitHeight; }
		float  			GetFitRMS() { return fFitRMS; }
		float  			GetFitMean() { return fFitMean; }
		float  			GetFitIntegral() { return fFitIntegral; }
		float  			GetTailIntegral() { return fTailIntegral; }
		bool   			GetIsVeto() { return fIsVeto; }
		int    			GetFlangeRing() { return fRing; }
		int    			GetFlangeRingLoc() { return fRingLoc; }
//------------------------------------------------------------
    
    void   			SetPEHist(TH1D* value) { fPEHist = value; }
    void				SetPEHistPtr(std::shared_ptr<TH1D> value) {fPEHistPtr = value;}
    void   			SetPECanvas(TCanvas* value) { fPECanvas = value; }
    void 	 			SetPMTID(int value) { fPMTID = value;}
    void   			SetPMTName(std::string value) { fPMTName = value;}
    void   			SetPMTRow(int value) { fPMTRow = value; }
    void   			SetPMTColumn(int value) { fPMTColumn = value; }
    void 	 			SetADCBoard(int value) { fADCBoard = value; }
		void 	 			SetADCChannel(int value) { fADCChannel = value; }
		void   			SetPMTFlange(int value) { fPMTFlange = value; }
    void   			SetSPE(float value) { fSPE = value; }
    void   			SetSPEPeakPos(float value) { fSPEPeakPos = value; }
    void 	 			SetIntegralForRate(float value) { fIntegralForRate = value; }
    void   			SetRate(float value) { fRate = value; }
    void   			SetPEHistMean(float value) { fPEHistMean = value; }
		void   			SetPEHistRMS(float value) { fPEHistRMS = value; }
		void   			SetPEHistEntries(int value) { fPEHistEntries = value; }
		void   			SetNoisePeakPos(float value) { fNoiseWallPeakPos = value; }
		void   			SetNoiseEndPos(float value) { fNoiseWallEndPos = value; }
		void   			SetNoiseIntegral(float value) { fNoiseWallIntegral = value; }
		void   			SetNoiseEndContent(float value) { fNoiseEndBinContent = value;} 
		void   			SetFitStartPos(float value) { fFitStartPos = value; }
		void   			SetFitEndPos(float value) { fFitEndPos = value; }
		void   			SetFitEndContent(float value) { fFitEndContent = value; }
		void   			SetFitHeight(float value) { fFitHeight = value; }
		void   			SetFitRMS(float value) { fFitRMS = value; }
		void   			SetFitMean(float value) { fFitMean = value; }
		void   			SetFitIntegral(float value) { fFitIntegral = value; }
		void   			SetTailIntegral(float value) { fTailIntegral = value; }
		void   			SetIsVeto(bool value) { fIsVeto = value; }
		void   			SetFlangeRing(int value) { fRing = value; }
		void   			SetFlangeRingLoc(int value) { fRingLoc = value; }

	private:
    /// A canvas (do not know what it is here)
		TCanvas* 		 fPECanvas; //class does not own
    /// Pointer to the ADC distribution
		std::shared_ptr<TH1D>  fPEHistPtr;//class does not own
    /// Pointer to the ADC distribution (again???)
		TH1D* 			 fPEHist;
    /// The PMT ID
		int     		 fPMTID;
    /// The PMT name
		std::string  fPMTName;
    /// The row the PMT is on
		int 				 fPMTRow;
    /// The column the PMT is on
		int 				 fPMTColumn;
    /// The digitizer board
		int     		 fADCBoard;
    /// The channel on the digitizer board
		int     		 fADCChannel;
    /// The flange
		int     		 fPMTFlange;
    /// The SPE value
		float   		 fSPE;
    /// The SPE value based on peak position
		float     		 fSPEPeakPos;
    /// The integral used for the rate calculation
		float				 fIntegralForRate; 
    /// The rate
		float   		 fRate;
    /// The mean of the ADC distribution
		float   		 fPEHistMean;
    /// The RMS of the ADC distribution
		float   		 fPEHistRMS;
    /// The number of entries in the hist
		int     		 fPEHistEntries;
    /// Where the noise wall peaks
		float     		 fNoiseWallPeakPos;//(pos means position)
    /// Where the noise wall ends
		float     		 fNoiseWallEndPos;
    /// The integral of the noise wall
		float   		 fNoiseWallIntegral;
    /// The bin content at the end of the noise wall
		float   		 fNoiseEndBinContent;
    /// The start position of the fit
		float     		 fFitStartPos;
    /// The end position of the fit
		float     		 fFitEndPos;
    /// The contents in the end position of the fit
		float   		 fFitEndContent; 
    /// The hight of the distribution from the fit
		float   		 fFitHeight;
    /// The RMS of the distribution from the fit
		float   		 fFitRMS;
    /// The Mean of the distribution from the fit
		float   		 fFitMean;
    /// The integral of the distribution from the fit (???)
		float   		 fFitIntegral;
    /// The integral of the tail of the distribution
		float   		 fTailIntegral;
    /// Is the PMT a veto PMT?
		bool    		 fIsVeto;
    /// Which flange ring it is located
		int     		 fRing;
    /// Location in the flange ring
		int     		 fRingLoc;
};
#endif
