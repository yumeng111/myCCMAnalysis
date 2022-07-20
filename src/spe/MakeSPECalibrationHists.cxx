/*!**********************************************
 * \file MakeSPECalibrationHists.cxx
 * \brief Makes all the various histograms to debug the SPE calibraions
 * \author T.J. Schuab and R.T. Thornton
 *
 * This code was originally written by T.J. Schuab and was adapted 
 * by R. T. Thornton to fit the CCM analysis framework
 ***********************************************/

#include <sstream>

#include "CCMAnalysis/spe/MakeSPECalibrationHists.h"

#include "CCMAnalysis/spe/SPECalibrationVariables.h"
#include "CCMAnalysis/utils/MsgLog.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2Poly.h"
#include "TLegend.h"

/////////////////////////////////////////////////
// Define static variables
/////////////////////////////////////////////////

/// 2D Histogram to save the rates for individual PMTs with row, column mapping
std::shared_ptr<TH2D> MakeSPECalibrationHists::fDetector2DHist = nullptr;
/// 2D Histogram to save the SPE for individual PMTs
std::shared_ptr<TH2D> MakeSPECalibrationHists::fSPE2DHist = nullptr;
/// 2D Histogram to save the rates for individual PMTs with digitizer, channel mapping
std::shared_ptr<TH2D> MakeSPECalibrationHists::fADC2DHist = nullptr;
/// 2D Histogram to save the rates for individual PMTs with the flange mapping (flange 1)
std::shared_ptr<TH2Poly> MakeSPECalibrationHists::fFlange1PMTs = nullptr;
/// 2D Histogram to save the rates for individual PMTs with the flange mapping (flange 2)
std::shared_ptr<TH2Poly> MakeSPECalibrationHists::fFlange2PMTs = nullptr;

/// 1D Histogram to show the SPE distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fSPEHist = nullptr;
/// 1D Histogram to show the rate distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fRateHist1D = nullptr;
/// 1D Histogram to show the noise wall peak distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fNoiseWallPeak = nullptr;
/// 1D Histogram to show the noise wall end distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fEndNoiseWallHist = nullptr;
/// 1D Histogram to show the fit range end distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fFitRangeEndHist = nullptr;
/// 1D Histogram to show the fit RMS distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fFitRMSHist = nullptr;
/// 1D Histogram to show the noise wall integral distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fIntegrateNoiseWallHist = nullptr;
/// 1D Histogram to show the tail integral distribution across all PMTs
std::shared_ptr<TH1D> MakeSPECalibrationHists::fIntegrateTailHist = nullptr;

/// Pointer to make a custom axis for a plot
TGaxis *MakeSPECalibrationHists::fAxis = 0;
/// Variable to state the size of the #SPECalibrationVariables vector
int MakeSPECalibrationHists::fSizeOfCal = 0;

/// Pointer to the SPE tree that might be saved
std::shared_ptr<TTree> MakeSPECalibrationHists::fSPETree = nullptr; 
/// Pointer to the trigger tree that might be saved
std::shared_ptr<TTree> MakeSPECalibrationHists::fTrigTree = nullptr;

/// Map to all the various plots that were created.
/// They are stored in the map for easier access
std::map<HistInfo_t,std::shared_ptr<TH1>> MakeSPECalibrationHists::fHistMap;
/// Vector of all the types of histograms that were created
std::vector<HistInfo_t> MakeSPECalibrationHists::fHistInfoVec;
/// Iterator for the #fHistMap
std::map<HistInfo_t,std::shared_ptr<TH1>>::iterator MakeSPECalibrationHists::fHistMapIter;

/// Values for the X points of Flange 1
std::vector<double> MakeSPECalibrationHists::fFlange1XVals;
/// Values for the Y points of Flange 1
std::vector<double> MakeSPECalibrationHists::fFlange1YVals;
/// Values for the names of each feed through on Flange 1
std::vector<std::string> MakeSPECalibrationHists::fFlange1Names;
/// Do not know what this variable is for
std::vector<int> MakeSPECalibrationHists::fFlange1RunNums;
/// Do not know what this variable is for
std::vector<int> MakeSPECalibrationHists::fFlange1RunLocs;
/// Values for the X points of Flange 2
std::vector<double> MakeSPECalibrationHists::fFlange2XVals;
/// Values for the Y points of Flange 2
std::vector<double> MakeSPECalibrationHists::fFlange2YVals;
/// Values for the names of each feed through on Flange 2
std::vector<std::string> MakeSPECalibrationHists::fFlange2Names;
/// Do not know what this variable is for
std::vector<int> MakeSPECalibrationHists::fFlange2RunNums;
/// Do not know what this variable is for
std::vector<int> MakeSPECalibrationHists::fFlange2RunLocs;


/*!**********************************************
 * \fn void MakeSPECalibrationHists::Clear()
 * \brief Delete #fAxis
 ***********************************************/
void MakeSPECalibrationHists::Clear()
{
  if (fAxis) {
    delete fAxis;
  }
}


/*!**********************************************
 * \fn void MakeSPECalibrationHists::MakeAllHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, std::string pdfVar, int nEntries, std::string flags)
 * \brief Draw all the hists that have the corresponding flag true 
 * \param[in] cal The vector of calibrations results to plot 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] nEntries Number of entries used to calculate the SPE values
 * \param[in] flags A list of the flags to turn on
 *
 * This function has a lot of inputs and maybe should be rewritten so that the inputs are
 * set independently ahead of time. I can see a function for \p fillAllHists that will set
 * all the other flags to true or false.
 *
 * Anyway, this function calls the corresponding fill and draw function for the various flags
 * that are passed if that flag is true
 ***********************************************/
void MakeSPECalibrationHists::MakeAllHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, std::string pdfVar,
    int nEntries, std::string flags)
{

  std::string pathPrefix = "";

  std::stringstream ss(flags);
  std::vector<std::string> flagVec;
  std::string tempFlag = "";
  while (ss >> tempFlag) {
    flagVec.push_back(tempFlag);
  }

  std::vector<std::string>::iterator itFlagBegin = flagVec.begin();
  std::vector<std::string>::iterator itFlagEnd = flagVec.end();

  bool pdfOnlyFlag = std::find(itFlagBegin,itFlagEnd,"pdfOnly") != itFlagEnd;
  bool rootOnlyFlag = std::find(itFlagBegin,itFlagEnd,"rootOnly") != itFlagEnd;
  bool rootAndPDFFlag = std::find(itFlagBegin,itFlagEnd,"rootAndPDF") != itFlagEnd;
  bool rootTreeFlag = std::find(itFlagBegin,itFlagEnd,"rootTree") != itFlagEnd;

  SetSizeOfCal(cal.size());
  HistInfo_t myEnum;
  if(rootOnlyFlag == true || rootAndPDFFlag == true || rootTreeFlag == true){
    TFile * outfile = TFile::Open(std::string(pdfVar + ".root").c_str(),"RECREATE");
    outfile->cd();
  }

  for (const auto & flag : flagVec) {

    if(flag.compare("rate2DHist") == 0){
      myEnum = kRate2DHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("spe2DHist") == 0){
      myEnum = kSPE2DHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("adc2DHist") == 0){
      myEnum = kADC2DHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("flangeHists") == 0){
      myEnum = kFlangeHists1;
      fHistInfoVec.push_back(myEnum);
      myEnum = kFlangeHists2;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("spe1DHist") == 0){
      myEnum = kSPE1DHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("rate1DHist") == 0){
      myEnum = kRate1DHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("noisePeakHist") == 0){
      myEnum = kNoisePeakHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("noiseEndHist") == 0){
      myEnum = kNoiseEndHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("fitEndHist") == 0){
      myEnum = kFitEndHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("fitRMSHist") == 0){
      myEnum = kFitRMSHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("noiseIntegralHist") == 0){
      myEnum = kNoiseIntegralHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("tailIntegralHist") == 0){
      myEnum = kTailIntegralHist;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("peHists") == 0){
      myEnum = kPEHists;
      fHistInfoVec.push_back(myEnum);
    }
    if(flag.compare("rootTree") == 0){
      WriteRootTree(cal,nEntries);
    }
  } // end for flagVec

  std::string pdfEnd = "";
  size_t histInfoVecSize = fHistInfoVec.size();
  for(size_t i=0; i<histInfoVecSize; ++i){
    if(i==0){
      pdfEnd = ".pdf(";
    }else if(i==histInfoVecSize-1){
      pdfEnd = ".pdf)";
    }
    else{
      pdfEnd = ".pdf";
    }
    if(histInfoVecSize == 1){
      pdfEnd = ".pdf()";
    }

    if(fHistInfoVec[i]==kRate2DHist){
      FillDetRate2DHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteDetRate2DHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kSPE2DHist){
      FillSPE2DHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteSPE2DHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kADC2DHist){
      FillADC2DHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteADC2DHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kFlangeHists2){
      FillFlangeHists(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteFlangeHists(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kSPE1DHist){
      Fill1DSPEHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        Write1DSPEHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kRate1DHist){
      Fill1DRateHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        Write1DRateHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kNoisePeakHist){
      FillNoisePeakHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteNoisePeakHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kNoiseEndHist){
      FillNoiseEndHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteNoiseEndHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i]==kFitEndHist){
      FillFitEndHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteFitEndHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i] == kFitRMSHist){
      FillFitRMSHist(cal,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteFitRMSHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i] == kNoiseIntegralHist){
      FillNoiseIntegralHist(cal, nEntries,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteNoiseIntegralHist(pathPrefix,pdfVar,pdfEnd);
      }
    }
    else if(fHistInfoVec[i] == kTailIntegralHist){
      FillTailIntegralHist(cal, nEntries,rootOnlyFlag,rootAndPDFFlag);
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        WriteTailIntegralHist(pathPrefix,pdfVar,pdfEnd);
      }//end if
    }//end else if
    else if(fHistInfoVec[i]==kPEHists){
      if(rootOnlyFlag == true || rootAndPDFFlag == true){
        WritePEHistsToRoot(cal);
      }
      if(pdfOnlyFlag == true || rootAndPDFFlag == true){
        DrawPEHists(cal,pathPrefix,pdfVar,pdfEnd);
      }//end if
    }//end else if
  }//end for i		
}//end MakeAllHists function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WritePEHistsToRoot(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal)
 * \brief Writes the ADC distributions to the output ROOT file
 * \param[in] cal The vector containing the SPE calibration results
 *
 * Writes the ADC distributions to the output ROOT file.
 ***********************************************/
void MakeSPECalibrationHists::WritePEHistsToRoot(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal)
{
  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    //cal.at(i)->GetpeHist()->Write();
    cal.at(i)->GetPEHistPtr()->Write();
  }//end for i
}//end WritePEHistsToRoot

/*!**********************************************
 * \fn void MakeSPECalibrationHists::DrawPEHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Draws the ADC distributions to canvas
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::DrawPEHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal,std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{

  TCanvas *c = new TCanvas("c","c",1000,1000);
  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    std::shared_ptr<TH1D> peHist = cal.at(i)->GetPEHistPtr();
    peHist->GetXaxis()->UnZoom();

    std::shared_ptr<TH1D> derPEHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(peHist->Clone("derPEHist")));
    derPEHist->Reset("ICESM");
    std::shared_ptr<TH1D> secDerPEHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(peHist->Clone("secDerPEHist")));
    secDerPEHist->Reset("ICESM");
    if (peHist->Integral() != 0) {
      for (int i=0; i < derPEHist->GetNbinsX(); ++i) {
        if (i-3 < 0 || i+5 >= derPEHist->GetNbinsX()) {
          continue;
        } 
        derPEHist->SetBinContent(i+1, 1/60.0 * (peHist->GetBinContent(i+2) - peHist->GetBinContent(i) + 2.0*peHist->GetBinContent(i+3) - 2.0*peHist->GetBinContent(i-1) +
              3.0*peHist->GetBinContent(i+4) - 3.0*peHist->GetBinContent(i-2) + 4.0*peHist->GetBinContent(i+5) - 4.0*peHist->GetBinContent(i-3)));
      }
      for (int i=0; i < secDerPEHist->GetNbinsX(); ++i) {
        if (i-3 < 0 || i+5 >= secDerPEHist->GetNbinsX()) {
          continue;
        } 
        secDerPEHist->SetBinContent(i+1, 1/60.0 * (derPEHist->GetBinContent(i+2) - derPEHist->GetBinContent(i) + 2.0*derPEHist->GetBinContent(i+3) - 2.0*derPEHist->GetBinContent(i-1) +
              3.0*derPEHist->GetBinContent(i+4) - 3.0*derPEHist->GetBinContent(i-2) + 4.0*derPEHist->GetBinContent(i+5) - 4.0*derPEHist->GetBinContent(i-3)));
      }
    }

    float noiseWallEndPos = cal.at(i)->GetNoiseEndPos();
    float fitRangeEndPos = cal.at(i)->GetFitEndPos();
    float fitEndContent = cal.at(i)->GetFitEndContent();
    float noiseEndContent = cal.at(i)->GetNoiseEndContent();
    float noisePeakPos = cal.at(i)->GetNoisePeakPos();
    float SPEPeakPos = cal.at(i)->GetSPEPeakPos();
    int histEntries = cal.at(i)->GetPEHistEntries();
    float fitHeight= cal.at(i)->GetFitHeight();    
    float fitMean= cal.at(i)->GetFitMean();
    float fitRMS= cal.at(i)->GetFitRMS();



    gStyle->SetOptStat(0);
    TLegend leg3 = TLegend(.70,.50,.94,.94);
    leg3.SetTextFont(62);
    std::string opttitle = cal.at(i)->GetPEHistPtr()->GetName();
    leg3.AddEntry((TObject*)0, Form("PMT ID: %s",opttitle.c_str()),"");
    leg3.AddEntry((TObject*)0, Form("# of Entries = %d", histEntries),"");

    if(fitRangeEndPos - noiseWallEndPos<10) {
      MsgWarning(MsgLog::Form("POTENTIAL BAD FIT | PMT ID = %d",i));
    }

    leg3.AddEntry((TObject*)0, Form("Fit Mean = %.2f", fitMean),"");
    leg3.AddEntry((TObject*)0, Form("Fit Start = %.2f", noiseWallEndPos),"");
    leg3.AddEntry((TObject*)0, Form("Fit End = %.2f", fitRangeEndPos), "");
    leg3.AddEntry((TObject*)0, Form("Fit Height = %.2f", fitHeight),"");
    leg3.AddEntry((TObject*)0, Form("Fit RMS = %.2f", fitRMS),"");
    //leg3->AddEntry((TObject*)0, Form("Overflow = %2f", histOverflow),"");
    float ymax= peHist->GetMaximum();
    TLine line(noiseWallEndPos,0,noiseWallEndPos,ymax);
    TLine line2(fitMean,0,fitMean,ymax);
    TLine line3(noiseWallEndPos,noiseEndContent,fitRangeEndPos,fitEndContent);
    TLine line4(SPEPeakPos,0,SPEPeakPos,ymax);//
    TLine line5(noisePeakPos,0,noisePeakPos,ymax);
    TLine line0(peHist->GetXaxis()->GetXmin(),0,peHist->GetXaxis()->GetXmax(),0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(9);
    line4.SetLineColor(kViolet-2);
    line3.SetLineColor(kAzure+1);
    line.SetLineColor(kGreen);
    line2.SetLineColor(kOrange);
    line5.SetLineColor(kBlack);
    //leg3.AddEntry(&line, "Integral for Rate Start","L");
    //leg3.AddEntry(&line2, "SPE Value (Fit Mean)","L");
    //leg3.AddEntry(&line4, "SPE Peak (@ Max Value)", "L");
    //leg3.AddEntry(&line3, "Bottom of Fit (X1,Y1,X2,Y2)", "L");
    //leg3.AddEntry(&line5, "NoiseWall Peak", "L");
    //leg3->SetBorderSize(0);
    leg3.SetFillStyle(0);

    //peHist->GetXaxis()->SetRangeUser(1,50);
    //derPEHist->GetXaxis()->SetRangeUser(5,50);
    //secDerPEHist->GetXaxis()->SetRangeUser(5,50);

    c->cd();
    c->Clear();
    //c->Divide(1,3);
    //c->cd(1);
    gPad->SetLogx(true);
    gPad->SetLogy(false);
    gPad->SetMargin(0.12,0.05,0.1,0.06);
    TGaxis::SetMaxDigits(3);
    peHist->Draw("");
    leg3.Draw();
    //line.Draw();
    //line2.Draw();
    //line3.Draw();
    //line4.Draw();
    //line5.Draw();
    //c->cd(2);
    //gPad->SetLogx(true);
    //derPEHist->Draw();
    //leg3.Draw();
    //line.Draw();
    //line2.Draw();
    //line3.Draw();
    //line4.Draw();
    //line5.Draw();
    //line0.Draw();
    //c->cd(3);
    //gPad->SetLogx(true);
    //secDerPEHist->Draw();
    //leg3.Draw();
    //line.Draw();
    //line2.Draw();
    //line3.Draw();
    //line4.Draw();
    //line5.Draw();
    //line0.Draw();
    std::string pePDF = "";
    if(i==0){
      if(pdfEnd == ".pdf(" || pdfEnd == ".pdf)"){
        pePDF = ".pdf(";
      }else if (pdfEnd == ".pdf)"){
        pePDF = ".pdf";
      }
      c->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pePDF).c_str(),"pdf");
      cal.at(i)->SetPECanvas(c);
    } 
    else if(i>0 && i<158){  
      c->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + ".pdf").c_str(),"pdf");
      cal.at(i)->SetPECanvas(c);
    } 
    else if(i>157 && i < 160){
      if(pdfEnd == ".pdf)" || pdfEnd == ".pdf)"){
        pePDF = ".pdf)";
      }else if(pdfEnd == ".pdf" || pdfEnd == ".pdf("){
        pePDF = ".pdf";
      }
      c->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pePDF).c_str(),"pdf");
      cal.at(i)->SetPECanvas(c);
    }//end else if
  }// end for i 
  delete c;
}//end FitDrawPEHists function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillDetRate2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Draws the 2D Rate hists
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillDetRate2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fDetector2DHist = std::make_shared<TH2D>("fDetector2DHist","Detector PMT Rate Plot", 24,1,25,9,-1,8);
  fDetector2DHist->GetZaxis()->SetRangeUser(0,6);
  fDetector2DHist->GetYaxis()->SetLabelOffset(999);
  fDetector2DHist->GetYaxis()->SetTickLength(0);
  fDetector2DHist->GetXaxis()->SetTitle("Column");
  fDetector2DHist->GetZaxis()->SetTitle("Pulses per #mus");
  //fDetector2DHist->GetZaxis()->SetTitle("Occupancy (# Pulses/# Triggers)");
  //fDetector2DHist->GetZaxis()->SetTitle("KHz");
  fDetector2DHist->GetZaxis()->SetTitleOffset(.8);
  fDetector2DHist->SetMarkerColor(kWhite);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    int column = cal.at(i)->GetPMTColumn();
    int row = cal.at(i)->GetPMTRow();
    float rate = cal.at(i)->GetRate();
    int rowInverted = 6-row; 
    if(row==-1 || row==0 || row==6){
      rate=rate*64.;
    }

    fDetector2DHist->Fill(column,rowInverted,rate); 

  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fDetector2DHist->Write();
  }//end if
  fHistMap.insert(std::make_pair(kRate2DHist, fDetector2DHist));
}//end FillDetRate2DHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteDetRate2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Writes the 2D Rate hists
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteDetRate2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *rate2DDetectorCanvas = new TCanvas("rate2DDetectorCanvas","rate2DDetectorCanvas", 2500,1500);
  rate2DDetectorCanvas->cd();

  fHistMapIter=fHistMap.find(kRate2DHist);
  if(fHistMapIter==fHistMap.end()){
    MsgWarning("kRate not in map");
  }
  fHistMapIter->second->Draw("COLZTEXT");

  //SetDetRate2D(fDetector2DHist);
  gStyle->SetPalette(55);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetOptStat(0);
  TLine line(0,0,0,0);
  line.SetLineColor(kRed);
  line.SetLineWidth(2);
  line.DrawLine(10,-1,10,8); 
  line.DrawLine(17,-1,17,8);
  line.DrawLine(23,-1,23,8);
  line.DrawLine(4,-1,4,8);
  line.SetLineColor(kBlack);
  //line->DrawLine(7,7,8,7);
  //line->DrawLine(8,7,8,8);
  //line->DrawLine(7,7,7,8); 
  line.SetLineWidth(1);
  line.DrawLine(9,-1,9,0);
  line.DrawLine(9,0,10,0);
  line.DrawLine(10,0,10,-1);

  fAxis = new TGaxis(gPad->GetUxmin()-0,8,gPad->GetUxmin()-.001,-1,-1,8,9,"-"); //x values cant be same for inverted axis
  fAxis->SetLineColor(46);
  fAxis->SetLabelColor(46);
  fAxis->SetTitle("Ring");
  fAxis->CenterTitle();
  fAxis->SetTitleOffset(-1.25);
  fAxis->SetTitleColor(46);
  fAxis->SetLabelOffset(-.02);
  fAxis->CenterLabels();
  TLatex latex;
  latex.SetTextFont(63);
  latex.SetTextSize(30);
  latex.SetTextColor(kGray+1);
  latex.DrawLatex(1,-.8,"Beam Exit");
  latex.DrawLatex(5.5,-.8,"Beam Right");
  latex.DrawLatex(12,-.8,"Beam Enter");
  latex.DrawLatex(18.5,-.8,"Beam Left");
  latex.DrawLatex(23,-.8,"Beam Exit");
  latex.DrawLatex(-1.4,7.75,"VT (x64)");
  latex.DrawLatex(-1.4,6.75,"VCT (x64)");
  latex.DrawLatex(-1.4,.75,"VCB (x64)");
  //latex.DrawLatex(-1.4,7.75,"VT");
  //latex.DrawLatex(-1.4,6.75,"VCT");
  //latex.DrawLatex(-1.4,.75,"VCB");
  latex.DrawLatex(-1.4,-.25,"VB");
  latex.SetTextColor(kRed);
  latex.DrawLatex(13.25,-.4,"#odot");
  latex.DrawLatex(1.25,-.4,"#otimes");
  latex.SetTextColor(kBlack);
  latex.SetTextSize(15);
  latex.SetTextAngle(45);
  //latex.DrawLatex(1.24,5.2,"No SPE");
  //latex.DrawLatex(11.24,5.2,"No SPE");
  //latex.DrawLatex(2.25,3.2,"No SPE");
  //latex.DrawLatex(5.25,4.2,"No SPE");
  latex.DrawLatex(14.25,3.2,"No SPE");
  //latex.DrawLatex(22.25,2.2,"No SPE");
  //latex.DrawLatex(7.25,7.2,"No SPE");
  latex.DrawLatex(9.4,-.65,"OFF");
  //latex.DrawLatex(1.24,5.2,"No SPE");

  //latex.DrawLatex(6.3,2.3,"No HV");//C6R4
  //latex.DrawLatex(7.3,3.3,"No HV");//C7R3
  //latex.DrawLatex(9.3,3.3,"No HV");//C9R3
  //latex.DrawLatex(10.3,1.3,"No HV");//C10R5
  //latex.DrawLatex(10.3,4.3,"No HV");//C10R2
  //latex.DrawLatex(15.3,1.3,"No HV");//C15R5
  //latex.DrawLatex(18.3,4.3,"No HV");//C18R2
  //latex.DrawLatex(19.3,2.3,"No HV");//C19R4
  //latex.DrawLatex(20.3,4.3,"No HV");//C20R2
  //latex.DrawLatex(23.3,5.3,"No HV");//C23R1
  //latex.DrawLatex(23.3,4.3,"No HV");//C23R2
  //latex.DrawLatex(24.3,2.3,"No HV");//C24R4

  fAxis->Draw();
  rate2DDetectorCanvas->Update();
  rate2DDetectorCanvas->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");

  delete rate2DDetectorCanvas;
  delete fAxis;
}//end WriteDetRate2DHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillSPE2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fills the SPE2D Hists
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillSPE2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fSPE2DHist = std::make_shared<TH2D>("fSPE2DHist","Detector SPE Plot", 24,1,25,9,-1,8);
  fSPE2DHist->GetZaxis()->SetRangeUser(0,100); 
  fSPE2DHist->GetZaxis()->SetRangeUser(0,100);
  fSPE2DHist->GetXaxis()->SetTitle("Column");
  fSPE2DHist->GetZaxis()->SetTitle("SPE Value");
  fSPE2DHist->GetZaxis()->SetTitleOffset(.8);
  fSPE2DHist->SetMarkerColor(kWhite);
  fSPE2DHist->GetYaxis()->SetLabelOffset(999);
  fSPE2DHist->GetYaxis()->SetTickLength(0);


  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    int column = cal.at(i)->GetPMTColumn();
    int row = cal.at(i)->GetPMTRow();
    float spe = cal.at(i)->GetFitMean();
    int rowInverted = 6-row; 

    fSPE2DHist->Fill(column,rowInverted,spe); 
  }
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fSPE2DHist->Write();
  }
  fHistMap.insert(std::make_pair(kSPE2DHist, fSPE2DHist));
}//end FillSPE2DHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteSPE2DHist(std::string pathPrefix, std::string pdfVar,std::string pdfEnd)
 * \brief Writes the SPE2D Hists
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteSPE2DHist(std::string pathPrefix, std::string pdfVar,std::string pdfEnd)
{
  TCanvas *spe2DCanvas = new TCanvas("spe2DCanvas", "spe2DCanvas",2500,1500);
  spe2DCanvas->cd();

  fHistMapIter=fHistMap.find(kSPE2DHist);
  if(fHistMapIter==fHistMap.end()){
    MsgWarning("kRate not in map");
  }
  fHistMapIter->second->Draw("COLZTEXT");

  gStyle->SetPalette(55);
  //gStyle->SetPaintTextFormat(".2f");
  spe2DCanvas->Update();
  gStyle->SetOptStat(0);
  TLine *lineSPE = new TLine();
  lineSPE->SetLineColor(kRed);
  lineSPE->SetLineWidth(2);
  lineSPE->DrawLine(10,-1,10,8);
  lineSPE->DrawLine(17,-1,17,8);
  lineSPE->DrawLine(23,-1,23,8);
  lineSPE->DrawLine(4,-1,4,8);
  lineSPE->SetLineColor(kBlack);
  //lineSPE->DrawLine(7,7,8,7);
  //lineSPE->DrawLine(8,7,8,8);
  //lineSPE->DrawLine(7,7,7,8);
  lineSPE->SetLineWidth(1);
  lineSPE->DrawLine(9,-1,9,0);
  lineSPE->DrawLine(9,0,10,0);
  lineSPE->DrawLine(10,0,10,-1);
  //fDetector2DHist->GetXaxis()->CenterLabels();
  //fDetector2DHist->GetZaxis()->SetTitle("KHz");

  gStyle->SetPaintTextFormat(".2f");

  TGaxis *axisSPE = new TGaxis(gPad->GetUxmin()+0,8,gPad->GetUxmin()+.001,-1,-1,8,9,"-"); //x values cant be same for inverted axis
  axisSPE->SetLineColor(46);
  axisSPE->SetLabelColor(46);
  axisSPE->SetTitle("Ring");
  axisSPE->CenterTitle();
  axisSPE->SetTitleOffset(-1.25);
  axisSPE->SetTitleColor(46);
  axisSPE->SetLabelOffset(-.01);
  axisSPE->CenterLabels();


  TLatex *latexSPE = new TLatex();
  latexSPE->SetTextFont(63);
  latexSPE->SetTextSize(30);
  latexSPE->SetTextColor(kGray+1);
  latexSPE->DrawLatex(1,-.8,"Beam Exit");
  latexSPE->DrawLatex(5.5,-.8,"Beam Right");
  latexSPE->DrawLatex(12,-.8,"Beam Enter");
  latexSPE->DrawLatex(18.5,-.8,"Beam Left");
  latexSPE->DrawLatex(23,-.8,"Beam Exit");
  //latexSPE->DrawLatex(-1.4,7.75,"VT (x64)");
  //latexSPE->DrawLatex(-1.4,6.75,"VCT (x64)");
  //latexSPE->DrawLatex(-1.4,.75,"VCB (x64)");
  latexSPE->DrawLatex(-1.4,7.75,"VT");
  latexSPE->DrawLatex(-1.4,6.75,"VCT");
  latexSPE->DrawLatex(-1.4,.75,"VCB");
  latexSPE->DrawLatex(-1.4,-.25,"VB");
  latexSPE->SetTextColor(kRed);
  latexSPE->DrawLatex(13.25,-.4,"#odot");
  latexSPE->DrawLatex(1.25,-.4,"#otimes");
  latexSPE->SetTextColor(kBlack);
  latexSPE->SetTextSize(15);
  latexSPE->SetTextAngle(45);
  //T1->DrawLatex(1.24,5.2,"No SPE");
  //T1->DrawLatex(11.24,5.2,"No SPE");
  //T1->DrawLatex(2.25,3.2,"No SPE");
  //T1->DrawLatex(5.25,4.2,"No SPE");
  latexSPE->DrawLatex(14.25,3.2,"No SPE");
  //T1->DrawLatex(22.25,2.2,"No SPE");
  //T1->DrawLatex(7.25,7.2,"No SPE");
  latexSPE->DrawLatex(9.4,-.65,"OFF");
  //T1->DrawLatex(1.24,5.2,"No SPE");

  //T1->DrawLatex(6.3,2.3,"No HV");//C6R4
  //T1->DrawLatex(7.3,3.3,"No HV");//C7R3
  //T1->DrawLatex(9.3,3.3,"No HV");//C9R3
  //T1->DrawLatex(10.3,1.3,"No HV");//C10R5
  //T1->DrawLatex(10.3,4.3,"No HV");//C10R2
  //T1->DrawLatex(15.3,1.3,"No HV");//C15R5
  //T1->DrawLatex(18.3,4.3,"No HV");//C18R2
  //T1->DrawLatex(19.3,2.3,"No HV");//C19R4
  //T1->DrawLatex(20.3,4.3,"No HV");//C20R2
  //T1->DrawLatex(23.3,5.3,"No HV");//C23R1
  //T1->DrawLatex(23.3,4.3,"No HV");//C23R2
  //T1->DrawLatex(24.3,2.3,"No HV");//C24R4  
  axisSPE->Draw();
  spe2DCanvas->Update();
  spe2DCanvas->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete spe2DCanvas;
  delete lineSPE;
  delete latexSPE;
  delete axisSPE;
}//end WriteSPE2DHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillADC2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill 2D ADC Hist
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillADC2DHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fADC2DHist = std::make_shared<TH2D>("fADC2DHist","ADC Board:Channel Rate Plot", 10,0,10,15,0,15);

  int channelInverted=0;
  int board=0;
  int channel=0;
  float pmtRate = 0.;
  for(int i=0; i<fSizeOfCal; ++i){
    board = cal.at(i)->GetADCBoard();
    channel = cal.at(i)->GetADCChannel();
    pmtRate = cal.at(i)->GetRate();


    if (board==0){
      pmtRate=pmtRate*64;
    }
    if (board==1 && channel<8){
      pmtRate=pmtRate*64;
    }
    channelInverted = 14-channel;
    fADC2DHist->Fill(board,channelInverted,pmtRate); 
  }

  fADC2DHist->GetZaxis()->SetRangeUser(0,6);
  gStyle->SetPalette(55);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetOptStat(0);



  fADC2DHist->GetXaxis()->CenterLabels();
  fADC2DHist->GetYaxis()->SetLabelOffset(999);
  fADC2DHist->GetYaxis()->SetTickLength(0);
  fADC2DHist->GetXaxis()->SetNdivisions(10);
  gPad-> Update();

  fADC2DHist->GetXaxis()->SetTitle("ADC Board");
  fADC2DHist->GetZaxis()->SetTitle("Rate (Pulses per #mus)");
  //fDetector2DHist->GetZaxis()->SetTitle("KHz");
  fADC2DHist->GetZaxis()->SetTitleOffset(.8);
  fADC2DHist->SetMarkerColor(kWhite);
  //fADC2DHist->Draw("colzTEXT");
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fADC2DHist->Write();
  }
  fHistMap.insert(std::make_pair(kADC2DHist, fADC2DHist));
}//end FillADC2DHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteADC2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Writes 2D ADC Hist
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteADC2DHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *rate2DADCCanvas = new TCanvas("rate2DADCCanvas","rate2DADCCanvas",2000,1500);
  rate2DADCCanvas->cd();

  fHistMapIter = fHistMap.find(kADC2DHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("ADC hist not in map");
  }
  fHistMapIter->second->Draw("COLZTEXT");

  TGaxis *axis2 = new TGaxis(gPad->GetUxmin(),15,gPad->GetUxmin()+.001,0,0,15,16,"-"); //x values cant be same for inverted axis
  axis2->SetLineColor(46);
  axis2->SetLabelColor(46);
  axis2->SetTitle("ADC Channel");
  axis2->CenterTitle();
  axis2->SetTitleOffset(-1.25);
  axis2->SetTitleColor(46);
  axis2->SetLabelOffset(-.01);
  axis2->CenterLabels();

  TLatex *latex2 = new TLatex();
  latex2->SetTextFont(63);
  latex2->SetTextSize(30);
  latex2->SetTextColor(kMagenta-9);
  latex2->SetTextAngle(90);
  latex2->DrawLatex(1.05,9.2,"Veto Channels (x64)");
  latex2->SetTextAngle(0);
  latex2->SetTextColor(kBlack);
  latex2->SetTextSize(15);
  //latex2->DrawLatex(.33,1.33,"No SPE");
  //latex2->DrawLatex(2.33,3.33,"No SPE");
  //latex2->DrawLatex(3.33,14.33,"No SPE");
  latex2->DrawLatex(5.45,1.33,"EJ");
  latex2->DrawLatex(5.45,.33,"EJ");
  //latex2->DrawLatex(6.33,6.33,"No SPE");
  latex2->DrawLatex(6.33,4.33,"No SPE");
  //latex2->DrawLatex(7.33,3.33,"No SPE");
  latex2->DrawLatex(9.45,10.33,"EJ");

  TLine *lineADC = new TLine();
  lineADC->SetLineColor(kMagenta-9);
  lineADC->SetLineWidth(4);
  lineADC->DrawLine(0,15,2,15);
  lineADC->DrawLine(0,0,0,15);
  lineADC->DrawLine(2,7,2,15);
  lineADC->DrawLine(1,0,1,7);
  lineADC->DrawLine(1,7,2,7);
  lineADC->DrawLine(0,0,1,0);
  rate2DADCCanvas->Update();

  gPad->Update();
  axis2->Draw();

  gStyle->SetPaintTextFormat(".2f");
  rate2DADCCanvas->Update();
  rate2DADCCanvas->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete rate2DADCCanvas;
  delete lineADC;
  delete axis2;
  delete latex2;
}//end WriteADC2DHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::Arc(int n, double x, double y, double r, double *px, double *py)
 * \brief To create individual circles.
 * \param[in] n Number of points to represent the arc
 * \param[in] x Center of the arc in the X axis
 * \param[in] y Center of the arc in the Y axis
 * \param[in] r The radius of the arc
 * \param[out] px The x corrdinate for each \p n points in the arc
 * \param[out] py The y corrdinate for each \p n points in the arc
 *
 * Arc is a function to create individual circles. Will be used when Creating Flange  TH2Poly Histograms.
 ***********************************************/
void MakeSPECalibrationHists::Arc(int n, double x, double y, double r, double *px, double *py)
{
  double fgPi = TMath::Pi();
  double fg2Pi = 2.0*fgPi;
  // Add points on a arc of circle from point 2 to n-2
  double step = fg2Pi/static_cast<double>(n);
  double angle = 0.0;
  int count = 0;
  while (angle < fg2Pi) {
    px[count] = r*TMath::Cos(angle) +x;
    py[count] = r*TMath::Sin(angle) +y;
    angle += step;
    ++count;
  }//end while
}//end Arc

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillFlangeHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the flange histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillFlangeHists(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fFlange1PMTs = std::make_shared<TH2Poly>(); // TH2Poly containing flange 1
  fFlange1PMTs->SetName("fFlange1PMTs");
  fFlange1PMTs->SetTitle("Flange 1 PMTs");
  fFlange1PMTs->SetZTitle("Pulses per #mus");
  fFlange1PMTs->SetLineColor(kWhite);
  fFlange1PMTs->SetLineWidth(0);
  fFlange1PMTs->SetStats(0);
  fFlange1PMTs->GetZaxis()->SetTitleOffset(0.65);
  fFlange1PMTs->GetZaxis()->SetRangeUser(0,6);
  fFlange1PMTs->AddBin(-2.5,0,-2.5,0);//Adding Bins to Set Range of plot -- TH2Poly Limits are determined by Bins Added:X1,Y1,X2,Y2
  fFlange1PMTs->AddBin(0,2.5,0,2.5);
  fFlange1PMTs->AddBin(2.5,0,2.5,0);
  fFlange1PMTs->AddBin(0,-2.5,0,-2.5);
  fFlange1PMTs->SetMarkerColor(kWhite);
  fFlange1PMTs->SetMarkerSize(.5);

  fFlange2PMTs = std::make_shared<TH2Poly>(); // TH2Poly containing flange 2
  fFlange2PMTs->SetName("fFlange2PMTs");
  fFlange2PMTs->SetTitle("Flange 2 PMTs");
  fFlange2PMTs->SetZTitle("Pulses per #mus");
  fFlange2PMTs->SetLineColor(kWhite);
  fFlange2PMTs->SetLineWidth(0);
  fFlange2PMTs->SetStats(0);
  fFlange2PMTs->GetZaxis()->SetTitleOffset(0.65);
  fFlange2PMTs->GetZaxis()->SetRangeUser(0,6);
  fFlange2PMTs->AddBin(-2.5,0,-2.5,0);//Adding Bins to Set Range of plot -- TH2Poly Limits are determined by Bins Added:X1,Y1,X2,Y2
  fFlange2PMTs->AddBin(0,2.5,0,2.5);
  fFlange2PMTs->AddBin(2.5,0,2.5,0);
  fFlange2PMTs->AddBin(0,-2.5,0,-2.5);
  fFlange2PMTs->SetMarkerColor(kWhite);
  fFlange2PMTs->SetMarkerSize(.5);

  double fgPi = TMath::Pi();
  double fg2Pi = 2.0*fgPi;
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(kBird);

  // defind how many points to have per PMT circle
  const int nPoints = 50;
  double px[nPoints];
  double py[nPoints];
  for (int i=0; i <nPoints; ++i) {
    px[i] = 0;
    py[i] = 0;
  }

  // loop over the map of PMTs and place each one in its corresponding position
  // in both detector and flange space

  for (int p = 0; p < fSizeOfCal; ++p) {
    if((p+1)%16 == 0){
      continue;
    }

    if (cal.at(p)->GetPMTFlange() > 2){
      continue;
    }
    double x = 0;
    double y = 0;

    double radius = 0.0;
    double totalNum = 0.0;
    switch (cal.at(p)->GetFlangeRing()) {
      case 1: radius = 0.5; totalNum = 12.0; break;
      case 2: radius = 1.0; totalNum = 18.0; break;
      case 3: radius = 1.5; totalNum = 24.0; break;
      case 4: radius = 2.0; totalNum = 30.0; break;
      default: break;
    }
    std::string name=cal.at(p)->GetPMTName();
    double step = fg2Pi/totalNum;
    double angle = static_cast<double>(cal.at(p)->GetFlangeRingLoc())*step;
    double realAngle = angle - fgPi;
    x = radius*TMath::Cos(realAngle);
    y = radius*TMath::Sin(realAngle);

    Arc(nPoints,x,y,0.11,px,py);
    if (cal.at(p)->GetPMTFlange() == 1) {
      fFlange1XVals.push_back(x-.065);
      fFlange1YVals.push_back(y-.17);
      fFlange1Names.push_back(name);
      fFlange1RunNums.push_back(cal.at(p)->GetFlangeRing());
      fFlange1RunLocs.push_back(cal.at(p)->GetFlangeRingLoc());  
      fFlange1PMTs->AddBin(nPoints,px,py);
      fFlange1PMTs->Fill(x,y,cal.at(p)->GetRate());

    } else {
      fFlange2XVals.push_back(x-.065);
      fFlange2YVals.push_back(y-.17);
      fFlange2Names.push_back(name);
      fFlange2RunNums.push_back(cal.at(p)->GetFlangeRing());
      fFlange2RunLocs.push_back(cal.at(p)->GetFlangeRingLoc());
      fFlange2PMTs->AddBin(nPoints,px,py);
      fFlange2PMTs->Fill(x,y,cal.at(p)->GetRate());
    }
    //pmt->SetFlangeX(x);
    //pmt->SetFlangeY(y);
    if (std::isnan(x) || std::isnan(y)) {
      MsgInfo(MsgLog::Form("x %.2f y %.2f Ring %d RingLoc %d angle %.2f realAngle %.2f fgPi %.2f Flange %d",
            x,y,cal.at(p)->GetFlangeRing(),cal.at(p)->GetFlangeRingLoc(),angle,realAngle,fgPi,cal.at(p)->GetPMTFlange()))
    }//end if
    // done looping over all pmts

  }// end for p
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fFlange1PMTs->Write();
    fFlange2PMTs->Write();
  }//end if
  fHistMap.insert(std::make_pair(kFlangeHists1, fFlange1PMTs));
  fHistMap.insert(std::make_pair(kFlangeHists2, fFlange2PMTs));
}//end FillFlangeHists

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteFlangeHists(std::string pathPrefix, std::string pdfVar,std::string pdfEnd)
 * \brief Write the flange histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteFlangeHists(std::string pathPrefix, std::string pdfVar,std::string pdfEnd)
{
  TCanvas *flangeCanvas = new TCanvas("flangeCanvas", "flangeCanvas", 2000,1000);
  flangeCanvas->Divide(2,1);
  flangeCanvas->cd(1);

  fHistMapIter = fHistMap.find(kFlangeHists1);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("Flange1 hist not in map");
  }
  fHistMapIter->second->Draw("COLZTEXT");

  gStyle->SetPaintTextFormat(".2f");

  TLatex T3;
  T3.SetTextFont(63);
  T3.SetTextSize(9);
  T3.SetTextColor(kGreen);
  flangeCanvas->cd(2);

  fHistMapIter = fHistMap.find(kFlangeHists2);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("Flange2 hist not in map");
  }
  fHistMapIter->second->Draw("COLZTEXT");

  TLatex T4;
  T4.SetTextFont(63);
  T4.SetTextSize(9);
  //T2.SetTextColor(kGreen);
  const char *nameChar;
  const char *nameChar2;
  std::string namez;
  std::string namez2;
  flangeCanvas->cd(1);
  for(size_t i=0; i<fFlange1XVals.size(); ++i){
    //namez = fFlange1Names[i];
    namez = fFlange1Names[i];
    nameChar = namez.c_str();

    if(fFlange1RunNums[i]==1 && fFlange1RunLocs[i]==12){fFlange1XVals[i]=fFlange1XVals[i]-.1;}
    if(fFlange1RunNums[i]==1 && fFlange1RunLocs[i]==8){fFlange1XVals[i]=fFlange1XVals[i]-.03;}
    if(fFlange1RunNums[i]==1 && fFlange1RunLocs[i]==7){fFlange1XVals[i]=fFlange1XVals[i]-.1;}
    if(fFlange1RunNums[i]==1 && fFlange1RunLocs[i]==6){fFlange1XVals[i]=fFlange1XVals[i]-.2; fFlange1YVals[i]=fFlange1YVals[i]+.05;}
    if(fFlange1RunNums[i]==1 && fFlange1RunLocs[i]==11){fFlange1XVals[i]=fFlange1XVals[i]+.07;}
    if(fFlange1RunNums[i]==1 && fFlange1RunLocs[i]==1){
      double tempX=fFlange1XVals[i]+.2;
      double tempY=fFlange1YVals[i]-.05;
      T3.SetTextColor(kBlack);
      T3.DrawLatex(tempX,tempY,"VB2");

    }
    if(fFlange1RunNums[i]==1){T3.SetTextColor(kRed);}
    else if(fFlange1RunNums[i]==2){T3.SetTextColor(kYellow+1);}
    else if(fFlange1RunNums[i]==3){T3.SetTextColor(kBlue);}
    else if(fFlange1RunNums[i]==4){T3.SetTextColor(kGreen);}
    T3.DrawLatex(fFlange1XVals[i],fFlange1YVals[i],nameChar);
    /////////////////////////////////////////////////////////////////////
  }
  flangeCanvas->cd(2);
  flangeCanvas->Update();
  for(size_t i=0; i<fFlange2XVals.size(); ++i){
    //namez = fFlange1Names[i];
    namez2 = fFlange2Names[i];
    nameChar2 = namez2.c_str();

    if(fFlange2RunNums[i]==1 && fFlange2RunLocs[i]==12){fFlange2XVals[i]=fFlange2XVals[i]-.12;}
    if(fFlange2RunNums[i]==1 && fFlange2RunLocs[i]==1){fFlange2XVals[i]=fFlange2XVals[i]-.1;}
    if(fFlange2RunNums[i]==1 && fFlange2RunLocs[i]==8){fFlange2XVals[i]=fFlange2XVals[i]-.03;}
    if(fFlange2RunNums[i]==1 && fFlange2RunLocs[i]==7){fFlange2XVals[i]=fFlange2XVals[i]-.1;}
    if(fFlange2RunNums[i]==1 && fFlange2RunLocs[i]==6){fFlange2XVals[i]=fFlange2XVals[i]-.2; fFlange2YVals[i]=fFlange2YVals[i]+.05;}
    if(fFlange2RunNums[i]==1 && fFlange2RunLocs[i]==11){fFlange2XVals[i]=fFlange2XVals[i]+.07;}

    if(fFlange2RunNums[i]==1){T4.SetTextColor(kRed);}
    if(fFlange2RunNums[i]==2){T4.SetTextColor(kYellow+1);}
    if(fFlange2RunNums[i]==3){T4.SetTextColor(kBlue);}
    if(fFlange2RunNums[i]==4){T4.SetTextColor(kGreen);}
    flangeCanvas->Update();
    T4.DrawLatex(fFlange2XVals[i],fFlange2YVals[i],nameChar2);

  }//end for fFlange2XVals
  flangeCanvas->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  //delete fVetoPMTs;
  //delete fTankPMTs;
  delete flangeCanvas;

}//end MakeFlangeHists function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::Fill1DSPEHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the 1D SPE distributions
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::Fill1DSPEHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fSPEHist = std::make_shared<TH1D>("fSPEHist","SPE Hist 1D; SPE Value (ADC); Count", 70,-.5,140);
  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fSPEHist->Fill(cal.at(i)->GetFitMean());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fSPEHist->Write();
  }//end if
  fHistMap.insert(std::make_pair(kSPE1DHist, fSPEHist));
}//end Fill1DSPEHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::MakeLegend1D(std::shared_ptr<TH1> hist)
 * \brief Make the legend for the 1D histograms
 * \param[in] hist The pointer to the hitogram that is being plotted
 *
 * Makes the legend for the 1D histograms assuming the histogram is the SPE histogram
 * or at least the legend should appear in the top right corner of the plot
 ***********************************************/
void MakeSPECalibrationHists::MakeLegend1D(std::shared_ptr<TH1> hist)
{	
  TLegend leg = TLegend(.75,.75,.9,.9);
  leg.SetTextFont(62);
  leg.Clear();
  int entries = hist->GetEntries();
  leg.AddEntry((TObject*)0, Form("# of Entries = %d", entries),"");
  float histMean = hist->GetMean();
  leg.AddEntry((TObject*)0, Form("Hist Mean %.2f", histMean),"");
  float histRMS = hist->GetRMS();
  leg.AddEntry((TObject*)0, Form("Hist RMS %.2f", histRMS),"");
  leg.Draw();
}//end MakeLegend1D

/*!**********************************************
 * \fn void MakeSPECalibrationHists::Write1DSPEHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the 1D SPE histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::Write1DSPEHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kSPE1DHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("SPE1D Hist not in map");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}//end Make1DSPEHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::Fill1DRateHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the 1D Rate histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::Fill1DRateHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fRateHist1D = std::make_shared<TH1D>("fRateHist1D", "Rate Hist 1D; Pulses per #mus; Count",60,-.5,6);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fRateHist1D->Fill(cal.at(i)->GetRate());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fRateHist1D->Write();
  }
  fHistMap.insert(std::make_pair(kRate1DHist,fRateHist1D));
}//end Fill1DRateHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::Write1DRateHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the 1D Rate histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::Write1DRateHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kRate1DHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("Rate1DHist not in map");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}//end Make1DSPEHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillNoisePeakHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the noise peak histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillNoisePeakHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fNoiseWallPeak = std::make_shared<TH1D>("fNoiseWallPeak", "Noise Wall Peak Hist 1D; Peak (ADC); Count)",30,0,30);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fNoiseWallPeak->Fill(cal.at(i)->GetNoisePeakPos());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fNoiseWallPeak->Write();
  }
  fHistMap.insert(std::make_pair(kNoisePeakHist, fNoiseWallPeak));
}//end FillNoisePeakHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteNoisePeakHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the noise peak histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteNoisePeakHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();

  fHistMapIter = fHistMap.find(kNoisePeakHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("noise peak 1d hist not in map");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillNoiseEndHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the noise end histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillNoiseEndHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fEndNoiseWallHist = std::make_shared<TH1D>("fEndNoiseWallHist","End NoiseWall/Start of FitRange Hist 1D; End NoiseWall (ADC); Count",60,0,60);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fEndNoiseWallHist->Fill(cal.at(i)->GetNoiseEndPos());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fEndNoiseWallHist->Write();
  }
  fHistMap.insert(std::make_pair(kNoiseEndHist, fEndNoiseWallHist));
}//end FillNoiseEndHist function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteNoiseEndHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the noise end histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteNoiseEndHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kNoiseEndHist);
  if(fHistMapIter==fHistMap.end()){
    MsgWarning("noise end hist not in map");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}//end WriteNoiseEndHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillFitEndHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the fit end histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillFitEndHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fFitRangeEndHist = std::make_shared<TH1D>("fFitRangeEndHist", "Fitrange End Hist 1D; Fit End (ADC); Count)",500,0,500);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fFitRangeEndHist->Fill(cal.at(i)->GetFitEndPos());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fFitRangeEndHist->Write();
  }
  fHistMap.insert(std::make_pair(kFitEndHist, fFitRangeEndHist));
}//end FillFitEndHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteFitEndHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the fit end histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteFitEndHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kFitEndHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning(" fit end hist not in map ");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}//end WriteFitEndHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillFitRMSHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the fit RMS histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillFitRMSHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fFitRMSHist = std::make_shared<TH1D>("fFitRMSHist", "Fit RMS Hist 1D; Fit RMS; Count",30,0,60);
  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fFitRMSHist->Fill(cal.at(i)->GetFitRMS());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fFitRMSHist->Write();
  }
  fHistMap.insert(std::make_pair(kFitRMSHist, fFitRMSHist));
}//end for FillFitRMS function

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteFitRMSHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the fit RMS histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteFitRMSHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{	
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kFitRMSHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning("Fit RMS hist not in map ");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);

  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}//end WriteFitRMSHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillNoiseIntegralHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the noise intergral histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] nEntries Number of entries used to calculate the SPE values
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillNoiseIntegralHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fIntegrateNoiseWallHist = std::make_shared<TH1D>("fIntegrateNoiseWallHist","NoiseWall Integral Hist 1D; Integral (ADC); Count)",nEntries/40,0,nEntries*40);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fIntegrateNoiseWallHist->Fill(cal.at(i)->GetNoiseIntegral());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fIntegrateNoiseWallHist->Write();
  }
  fHistMap.insert(std::make_pair(kNoiseIntegralHist, fIntegrateNoiseWallHist));
}//end for FillNoiseIntegralHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteNoiseIntegralHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the noise integral histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteNoiseIntegralHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kNoiseIntegralHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning(" noise integral hist not in map ");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;
}//end WriteNoiseIntegralHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::FillTailIntegralHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries, bool rootOnlyFlag, bool rootAndPDFFlag)
 * \brief Fill the tail integral histograms
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] nEntries Number of entries used to calculate the SPE values
 * \param[in] rootOnlyFlag Flag to save the hists only in the ROOT file
 * \param[in] rootAndPDFFlag Flag to save the hists in a ROOT file and in a pdf
 ***********************************************/
void MakeSPECalibrationHists::FillTailIntegralHist(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries, bool rootOnlyFlag, bool rootAndPDFFlag)
{
  fIntegrateTailHist = std::make_shared<TH1D>("fIntegrateTailHist","Tail Integral Hist 1D; Tail Integral(ADC DerIntegral); Count)",nEntries/10,0,nEntries*10);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16 == 0){
      continue;
    }
    fIntegrateTailHist->Fill(cal.at(i)->GetTailIntegral());
  }//end for i
  if(rootOnlyFlag == true || rootAndPDFFlag == true){
    fIntegrateTailHist->Write();
  }
  fHistMap.insert(std::make_pair(kTailIntegralHist, fIntegrateTailHist));
}//end for FillTailIntegralHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteTailIntegralHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
 * \brief Write the tail integral histograms
 * \param[in] pathPrefix The prefix of the file name 
 * \param[in] pdfVar The postfix of the file name
 * \param[in] pdfEnd The ending to the pdf (does it contain a parenthesis or not?)
 ***********************************************/
void MakeSPECalibrationHists::WriteTailIntegralHist(std::string pathPrefix, std::string pdfVar, std::string pdfEnd)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
  fHistMapIter = fHistMap.find(kTailIntegralHist);
  if(fHistMapIter == fHistMap.end()){
    MsgWarning(" tail integral hist not in map ");
  }
  fHistMapIter->second->Draw();
  MakeLegend1D(fHistMapIter->second);
  c1->Print(std::string(pathPrefix + "HIsTZOUT_" + pdfVar + pdfEnd).c_str(),"pdf");
  delete c1;	
}//end WriteTailIntegralHist

/*!**********************************************
 * \fn void MakeSPECalibrationHists::WriteRootTree(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries)
 * \brief Write the root tree to the current file
 * \param[in] cal The vector containing the SPE calibration results
 * \param[in] nEntries Number of entries used to calculate the SPE values
 ***********************************************/
void MakeSPECalibrationHists::WriteRootTree(const std::vector<std::shared_ptr<SPECalibrationVariables>> & cal, int nEntries)
{
  int row = 0;
  int column = 0;
  int adcBoard = 0;
  int adcChannel = 0;
  int flange = 0;
  float noiseEnd = 0;
  float fitRangeEnd = 0;
  float fitRMS = 0.;
  float fitMean = 0;
  float rate = 0.;
  int entries = 0;
  std::string name = "";
  int pmtID = 0;

  fSPETree = std::make_shared<TTree>("spe", "Info for Rate Hists");
  fSPETree->Branch("row",&row);
  fSPETree->Branch("column",&column);
  fSPETree->Branch("adcBoard",&adcBoard);
  fSPETree->Branch("adcChannel",&adcChannel);
  fSPETree->Branch("flange", &flange);
  fSPETree->Branch("endNoiseWallFitRangeStart", &noiseEnd);
  fSPETree->Branch("endFitRange", &fitRangeEnd);
  fSPETree->Branch("fitRMS", &fitRMS);
  fSPETree->Branch("speValue", &fitMean);
  fSPETree->Branch("rate", &rate);
  fSPETree->Branch("histEntries", &entries);
  fSPETree->Branch("name",&name);
  fSPETree->Branch("pmtID", &pmtID);

  for(int i=0; i<fSizeOfCal; ++i){
    if((i+1)%16==0){
      continue;
    }
    row = cal.at(i)->GetPMTRow();
    column = cal.at(i)->GetPMTColumn();
    adcBoard = cal.at(i)->GetADCBoard();
    adcChannel = cal.at(i)->GetADCChannel();
    flange = cal.at(i)->GetPMTFlange();
    noiseEnd = cal.at(i)->GetNoiseEndPos();
    fitRangeEnd = cal.at(i)->GetFitEndPos();
    fitRMS = cal.at(i)->GetFitRMS();
    fitMean = cal.at(i)->GetFitMean();
    rate = cal.at(i)->GetRate();
    entries = cal.at(i)->GetPEHistEntries();
    //entries = cal.at(i)->GetEntries();
    name = cal.at(i)->GetPMTName();
    pmtID = cal.at(i)->GetPMTID();
    fSPETree->Fill();	
  }//end for i

  fTrigTree = std::make_shared<TTree>("triggers", "Number of Triggers");
  fTrigTree->Branch("triggers", &nEntries);
  fTrigTree->Fill();
  fSPETree->Write();
  fTrigTree->Write();

}//end WriteRootTree function

