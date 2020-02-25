/*!**********************************************
 * \file SPECalibrationVariables.cxx
 * \brief Source code for the #SPECalibrationVariables class
 * \author T. J. Schuab and R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "SPECalibrationVariables.h"

SPECalibrationVariables::SPECalibrationVariables()

{

	//fPECanvas = 0;
		fPEHist = 0;
		fPMTID = 0;
		fPMTName = "";
		fPMTRow = 0;
		fPMTColumnumn = 0;
		fADCBoard = 0;
		fADCChannel = 0;
		fPMTFlange = 0;
		fSPE = 0;
		fSPEPeakPos = 0;
		fIntegralForRate = 0; 
		fRate = 0;
		fPEHistMean = 0;
		fPEHistRMS = 0;
		fPEHistEntries = 0;
		fNoiseWallPeakPos = 0;//(pos means position)
		fNoiseWallEndPos = 0;
		fNoiseWallIntegral = 0;
		fNoiseEndBinContent = 0;
		fFitStartPos = 0;
		fFitEndPos = 0;
		fFitEndContent = 0; 
		fFitHeight = 0;
		fFitRMS = 0;
		fFitMean = 0;
		fFitIntegral = 0;
		fTailIntegral = 0;
		fIsVeto = 0;
		fRing = 0;
		fRingLoc = 0;
}


