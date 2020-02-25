/*!**********************************************
 * \file MakeSPECalibrationHists.h
 * \brief Holds the class MakeSPECalibrationHists
 * \author T.J. Schuab and R. T. Thornton (LANL)
 *
 * Holds the class MakeSPECalibrationHists
 ***********************************************/
#ifndef MakeSPECalibrationHists_h
#define MakeSPECalibrationHists_h

#include "Utility.h"
#include <memory>
#include <map>

class TLine;
class TTree;
class TLatex;
class TGaxis;
class TH2D;
class TH2Poly;
class TH1D;
class TH1;
class TCanvas;
class SPECalibrationVariables;
class TLegend;




/*!**********************************************
 * \class MakeSPECalibrationHists
 * \brief Class that fills the SPE results into hists and saves the information
 *
 * This class holds the functions that saves the SPE results into histograms
 * and then saves the histograms to a ROOT file or a pdf.
 *
 * There is a lot of duplication between the many functions, hence a separate
 * function for each fill and write. We should streamline the code to make
 * it more user friendly.
 ***********************************************/
class MakeSPECalibrationHists{
	private:
		static std::shared_ptr<TH2D> fDetector2DHist;
		static std::shared_ptr<TH2D> fSPE2DHist;
		static std::shared_ptr<TH2D> fADC2DHist;
		static std::shared_ptr<TH2Poly> fFlange1PMTs;
		static std::shared_ptr<TH2Poly> fFlange2PMTs;
		
		static std::shared_ptr<TH1D> fSPEHist;
		static std::shared_ptr<TH1D> fRateHist1D;
		static std::shared_ptr<TH1D> fNoiseWallPeak;
		static std::shared_ptr<TH1D> fEndNoiseWallHist;
		static std::shared_ptr<TH1D> fFitRangeEndHist;
  	static std::shared_ptr<TH1D> fFitRMSHist;
  	static std::shared_ptr<TH1D> fIntegrateNoiseWallHist;
  	static std::shared_ptr<TH1D> fIntegrateTailHist;
  
		static TGaxis *fAxis;
		static int fSizeOfCal;
		
		static std::shared_ptr<TTree> fSPETree; 
		static std::shared_ptr<TTree> fTrigTree;

		static std::map<HistInfo_t,std::shared_ptr<TH1>> fHistMap;
		static std::vector<HistInfo_t> fHistInfoVec;
		static std::map<HistInfo_t,std::shared_ptr<TH1>>::iterator fHistMapIter;

		static std::vector<double> fFlange1XVals;
  	static std::vector<double> fFlange1YVals;
  	static std::vector<std::string>fFlange1Names;
  	static std::vector<int> fFlange1RunNums;
  	static std::vector<int> fFlange1RunLocs;
  	static std::vector<double> fFlange2XVals;
  	static std::vector<double> fFlange2YVals;
  	static std::vector<std::string>fFlange2Names;
  	static std::vector<int> fFlange2RunNums;
  	static std::vector<int> fFlange2RunLocs;
		
	

	public:
    static void Clear();

		static void MakeAllHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, std::string pathPrefix, std::string pdfVar, int nEntries,
        bool FillAllHists,bool peHistsFlag, bool rate2DHistFlag,bool spe2DHistFlag, bool adc2DHistFlag, bool flangeHistsFlag, bool spe1DHistFlag,
        bool rate1DHistFlag,bool noisePeakHistFlag,bool noiseEndHistFlag,bool fitEndHistFlag,bool fitRMSHistFlag, bool noiseIntegralHistFlag,
        bool tailIntegralHistFlag, bool rootOnlyFlag, bool pdfOnlyFlag, bool rootAndPDFFlag, bool rootTreeFlag);
		
		static void FillDetRate2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag,bool rootAndPDFFlag);
		static void WriteDetRate2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillSPE2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag,bool rootAndPDFFlag);
		static void WriteSPE2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillADC2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag,bool rootAndPDFFlag);
		static void WriteADC2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void Arc(int n, double x, double y, double r, double *px, double *py);
		static void FillFlangeHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag,bool rootAndPDFFlag);
		static void WriteFlangeHists(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void MakeLegend1D(std::shared_ptr<TH1> hist);

		static void Fill1DSPEHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void Write1DSPEHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void Fill1DRateHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void Write1DRateHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillNoisePeakHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void WriteNoisePeakHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillNoiseEndHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void WriteNoiseEndHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillFitEndHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void WriteFitEndHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillFitRMSHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void WriteFitRMSHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillNoiseIntegralHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void WriteNoiseIntegralHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void FillTailIntegralHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries,bool rootOnlyFlag, bool rootAndPDFFlag);
		static void WriteTailIntegralHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd);

		static void DrawPEHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, std::string pathPrefix, std::string pdfVar, std::string pdfEnd);
		static void WritePEHistsToRoot(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal);
		
		static void WriteRootTree(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries);


    /// \fn static void SetSizeOfCal(int value)
    /// \brief Store the size of the #SPECalibrationVariables vector
    /// \param[in] value Size of the #SPECalibrationVariables vector
		static void SetSizeOfCal(int value) { fSizeOfCal = value; }

    /// \fn static int GetSizeOfCal()
    /// \brief Return the value #fSizeOfCal is equal to
    /// \return #fSizeOfCal
		static int  GetSizeOfCal() { return fSizeOfCal; }

};
#endif

