/*!**********************************************
 * \file NearlineSPEDiag.cxx
 * \brief Code that holds the functions for determining the SPE values
 * \author T.J. Schuab and R.T. Thornton
 * \date February 24, 2020
 ***********************************************/
#include "NearlineSPEDiag.h"
#include "MakeSPECalibrationHists.h"
#include "SPECalibrationVariables.h"
#include "Pulses.h"
#include "RawData.h"
#include "MsgLog.h"
#include "PMTInfoMap.h"
#include "Utility.h"
#include "PMTInformation.h"

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include <fstream>
#include <iostream>
#include <cstring>

/*!**********************************************
 * \fn double SPEFitFunction(double * x, double * par)
 * \brief A fit function to determine the SPE calibration
 * \param[in] x The independent value (where on the ADC axis) to calculate the function
 * \param[in] par The assumed parameters of the fit
 * \return The guessed value at \p x
 *
 * Thit fit function is a double Gaussian of the form
 * \f[
 * [0]*\exp\left\{\frac{(x[0] - [1])^2}{2*[2]^2}\right\}+
 * [3]*\exp\left\{\frac{(x[0] - 2*[1])^2}{2*[4]^2}\right\}
 * \f]
 * where the mean of the second Guassian is assumed to be twice the mean 
 * of the first Gaussian
 ***********************************************/
double SPEFitFunction(double * x, double * par)
{
  //if (par[1] > par[4]) {// || par[1] > par[7] || par[4] > par[7]) {
  //  return 999999.0;
  //}

  return par[0]+TMath::Gaus(x[0],par[1],par[2]) + 
    par[3]*TMath::Gaus(x[0],par[1]*2.0,par[4]);
}

/*!**********************************************
 * \fn void NearlineSPEDiag::SetClassVec(int value)
 * \brief Set the size of the #SPECalibrationVariables vectors to \p value
 * \param[in] value The new size of the #SPECalibrationVariables
 *
 * Must be called before #CreatePEHists
 ***********************************************/
void NearlineSPEDiag::SetClassVec(int value)
{
  fCal.resize(value);
  fCalBefore.resize(value);
  for(int i=0; i<value; ++i){
  	fCal.at(i) = std::make_shared<SPECalibrationVariables>();
  	fCalBefore.at(i) = std::make_shared<SPECalibrationVariables>();
  }
  fSizeOfCal = value;
}

/*!**********************************************
 * \fn void NearlineSPEDiag::CreatePEHists()
 * \brief Create the ADC distribution histograms needed to determine the SPE calibraion values
 ***********************************************/
void NearlineSPEDiag::CreatePEHists()
{
  //std::vector<double> bins;
  //for (double exp = -2; exp < 5; ++exp) {
  //  for (double digit = 1; digit < 10; digit+=0.25) {
  //    bins.push_back(digit*std::pow(10.0,exp));
  //  }
  //}
	for(size_t i=0; i<fSizeOfCal; ++i){
		//fCal.at(i)->SetpeHist(new TH1D(Form("peHist%d",i),";PE;Count",500,0,1000));
		//fCal.at(i)->SetPEHistPtr(std::make_shared<TH1D>(Form("peHist%d",i),";PE;Count",bins.size()-1,&bins.front()));
		fCal.at(i)->SetPEHistPtr(std::make_shared<TH1D>(Form("peHist%zu",i),";ADC;Count Above Background",2000,0,1000));
		fCal.at(i)->SetPMTID(i);
		fCalBefore.at(i)->SetPEHistPtr(std::make_shared<TH1D>(Form("peHistBefore%zu",i),";PE;Count",2000,0,1000));
		fCalBefore.at(i)->SetPMTID(i);
	}
}

/*!**********************************************
 * \brief Loop through the files and create the necessary distributions
 * \param[in] fileList The list of files to process
 * \param[in] ledRunFlag A flag to state if to look at LED triggers or all triggers (true = LED triggers only)
 * \param[in] numLEDTriggers The number of LED triggers to look at (default is -1 which means all)
 *
 * This function adds all the files in \p fileList to a chain and then loops through the chain to create
 * the ADC distributions needed to determine the SPE calibration.
 *
 * The \p ledRunFlag allows the user to state if the candidate triggers have to come from LED triggers or 
 * could come from any trigger. The flag also determines the DAQ window to look for events.
 * These numbers are hard coded which should be fixed. If \p ledRunFlag is true, it assumes only the 
 * 11th board (board number 10) has its waveforms saved.
 * 
 * For looking at data with beam triggers, you may not know how many files to loop over, you just want to
 * stop looking at the data once you have so many LED triggers. Set the \p numLEDTriggers to the number of
 * LED triggers you want to look at and the loop will quite once that many triggers has been processed. If
 * you did not pass enough files, then the total number of LED triggers will be less than \p numLEDTriggers.
 ***********************************************/
void NearlineSPEDiag::MakeChainFillPulses(const std::vector<std::string> & fileList, bool ledRunFlag, double windowStart, double windowEnd, int numLEDTriggers)
{
	TChain* chain=0;
	TChain* chainRaw=0;
	Pulses* pulses=0;
  RawData * rawData = 0;
  TTree * tree=0;
  size_t count=0;
  chain = new TChain("pulses","chain"); 
  chainRaw = new TChain("rawData","chainRaw"); 
  for (const auto & fileName : fileList) {
    MsgInfo(MsgLog::Form("Looking at: %s",fileName.c_str()));
    chain->AddFile(fileName.c_str());
    chainRaw->AddFile(fileName.c_str());
    if (count == 0) {
      TFile * file = TFile::Open(fileName.c_str(),"READ");
      file->GetObject(PMTInfoMap::TreeName().c_str(),tree);
      if (tree) {
        PMTInfoMap::FillPMTMap(tree);
      }else {
        PMTInfoMap::DefaultFillPMTMap();
      }
      delete file;
      ++count;
    } //end if count==
  }//end for fileName

  chain->SetBranchAddress("pulses",&pulses);
  chainRaw->SetBranchAddress("rawData",&rawData);
  fNEntries = chain->GetEntries();
  int ledEntries = 0;
  int totalSkipped = 0;
	
	
  fWindowStart = windowStart;
  fWindowEnd = windowEnd;
  fWindow = std::fabs(fWindowStart-fWindowEnd);
  for (long e = 0; e < fNEntries-1; ++e) {
    if (e % 500 == 0) {
      MsgInfo(MsgLog::Form("At event %ld out of %ld (window start %f window end %f)",e,fNEntries,fWindowStart,fWindowEnd));  
    }//end if e%
    int nbytes = chainRaw->GetEntry(e);
    if (nbytes == 0) {
      ++totalSkipped;
      continue;
    }// end if nbytes
      

    if (ledRunFlag) {
      if (ledEntries  >= numLEDTriggers && numLEDTriggers != -1) {
        break;
      }
      if (!rawData->IsTriggerPresent("LED")) {
        ++totalSkipped;
        continue;
      }
      ++ledEntries;
      fNEntries = ledEntries;
    }// end if ledRunFlag 


    nbytes = chain->GetEntry(e);
    if (nbytes == 0) {
      ++totalSkipped;
      continue;
    }// end if nbytes

    if (pulses->GetTriggerTime()*2e-3 - 9.92 <= -0.5) {
      ++totalSkipped;
      continue;
    }

    size_t numPulses = pulses->GetNumPulses();
    //pulses->Sort();  
    for (size_t pulse = 0; pulse < numPulses; ++pulse) {
      //const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(pulses->GetBoard(pulse),pulses->GetChannel(pulse));
      //if (!pmtInfo) {
      //  continue;
      //}
      float startTime = pulses->GetPulseTime(pulse)*2e-3 - 9.92;
      float length = pulses->GetPulseLength(pulse);
      if (length*2.0 < 20) {
        continue;
      }

      if (startTime < fWindowStart) {
        int key = pulses->GetKey(pulse);
        float integral = pulses->GetPulseIntegral(pulse);
        fCalBefore.at(key)->GetPEHistPtr()->Fill(integral);     
        continue;
      }
      if (startTime > fWindowEnd) {
        continue;
      }
      //int endTime = pulses->GetPulseTime(pulse) + length;
      //bool isVeto = pmtInfo->IsVeto();
      //if (!PMTInfoMap::IsActive(pulses->GetBoard(pulse),pulses->GetChannel(pulse))) {
      //  continue;
      //}
      //if (startTime < fWindowStart) {
        //continue;
      //}
      //if (startTime > fWindowEnd) {
        //break;
      //}
      int key = pulses->GetKey(pulse);
      float integral = pulses->GetPulseIntegral(pulse);
      fCal.at(key)->GetPEHistPtr()->Fill(integral);     
    }//end for pulse   
  }// for e < fNEntries
  fNEntries -= totalSkipped; // remove from count the ones that were skipped
}//end MakeChainFillPulses Function

/*!**********************************************
 * \brief Add the pulse that is passed to the correct CalibrationVariables object
 * \param[in] pulse The pulse to add
 * \param[in] windowStart The start time of the window
 * \param[in] windowEnd The end time of the window
 *
 * This function is mainly used by #CCMSPECalc::ProcessEvent. Please use #MakeChainFillPulses to use
 * this class in a standalone function.
 ***********************************************/
void NearlineSPEDiag::FillPulses(const Pulses & pulses, double windowStart, double windowEnd)
{
  fWindowStart = windowStart;
  fWindowEnd = windowEnd;
  fWindow = std::fabs(fWindowStart-fWindowEnd);   

  size_t numPulses = pulses.GetNumPulses();
  for (size_t pulse = 0; pulse < numPulses; ++pulse) {
    float startTime = pulses.GetPulseTime(pulse)*2e-3 - 9.92;
    float length = pulses.GetPulseLength(pulse);
    if (length*2.0 < 20) {
      continue;
    }

    if (startTime < fWindowStart) {
      int key = pulses.GetKey(pulse);
      float integral = pulses.GetPulseIntegral(pulse);
      fCalBefore.at(key)->GetPEHistPtr()->Fill(integral);     
      continue;
    }
    if (startTime > fWindowEnd) {
      continue;
    }
    int key = pulses.GetKey(pulse);
    float integral = pulses.GetPulseIntegral(pulse);
    fCal.at(key)->GetPEHistPtr()->Fill(integral);     
  }//end for pulse   
}//end FillPulses Function

/*!**********************************************
 * \fn void NearlineSPEDiag::GetHistsToAdjust(std::string fileToAdjust, bool ledRunFlag)
 * \brief Read in a previous calibration file and read in the histograms to determine the SPE calibration
 * \param[in] fileToAdjust Name of the file to use
 * \param[in] ledRunFlag A flag to state if to look at LED triggers or all triggers (true = LED triggers only)
 *
 * Instead of reading in the data, read in a file that has SPE calibration values and ADC
 * distributions saved in it. This allows new SPE values to be determined, with maybe an
 * improved fitter, without the time consuming part of recreating the ADC distributions
 * from data.
 *
 * See #MakeChainFillPulses for how \p ledRunFlag and used and its current limitations.
 ***********************************************/
void NearlineSPEDiag::GetHistsToAdjust(std::string fileToAdjust, bool ledRunFlag, double windowStart, double windowEnd)
{ 
  fWindowStart = windowStart;
  fWindowEnd = windowEnd;
  fWindow = std::fabs(fWindowStart-fWindowEnd);   
  TFile * file = TFile::Open(fileToAdjust.c_str(),"READ"); // opens a file as READ only
  //std::shared_ptr<TH1D> hist = 0;
  TH1D * hist=0; // where we store a histogram from a file
  int startNum = 0;
  int itNum = 160;
  //int numTrigTreeEntries = 0;
  //TTree *t=0;
  //int fNEntries = 0;
    
  //TH1D * peHist[160];
  //for (int i=startNum; i <itNum; ++i) {
  //  peHist[i] = 0;
  //}
    
  PMTInfoMap::DefaultFillPMTMap();
  	
  for (int i=startNum; i < itNum; ++i){
      if((i+1)%16==0){continue;}
    file->GetObject(Form("peHist%d",i),hist);     
    if (!hist){
      MsgError(MsgLog::Form("could not get hist %d",i));
    }
   
    /*
    file->GetObject("TrigTree",t);
    if (!t){
      MsgError("could not get trigtree");
    }
    t->SetBranchAddress("Triggers",&fNEntries);
    numTrigTreeEntries = t->GetEntries();
    if(numTrigTreeEntries == 0){
      MsgWarning("trig tree was not filled");
    }
    if(numTrigTreeEntries > 1){
      MsgWarning(MsgLog::Form("more than one trigtree entry, # = %d",numTrigTreeEntries));
    }
    t->GetEntry(0);  
    */
    //fCal.at(i)->SetpeHist(hist);
    std::shared_ptr<TH1D> hist2(hist); //= std::make_shared<TH1D>(hist);
    hist2->SetXTitle("ADC");
    hist2->SetYTitle("Count Above Background");
    //delete hist;
    fCal.at(i)->SetPEHistPtr(hist2);
  }//end GetObjects for loop
  MsgInfo(MsgLog::Form("Triggers in file to adjust = %d",fNEntries));

}//end GetHistsToAdjust function

/*!**********************************************
 * \fn void NearlineSPEDiag::CalculateRates()
 * \brief Find the SPE fit ranges and calculate various rates
 *
 * The following describes the method for finding the SPE peak:
 * -# Calculate the background subtracted ADC distribution.
 *  The background is computed from the time region before the LED trigger.
 *  So when \p ledRunFlag is set to false in #MakeChainFillPulses the background
 *  is empty and therefore the background subtracted distribution is equal to the
 *  unsubtracted distribution.
 * -# Compute the 9-point derivative (formula from Introduction to Instrumentation and
 *  Measurements Second Edition Robert B. Northrop page 677).
 *  \f[
 *  y = (1/60)(4x+3x_{-1}+2x_{-2}+x_{-3}-x_{-5}-2x_{-6}-3x_{-7}-4x_{-8})
 *  \f]
 *  where \f$x\f$ has been shifted to be in the middle of the derivative instead at the end.
 *  The second derivative is also calculated but not used.
 * -# Find where the noise wall ends and store some values.
 *  The terminology is left over when the SPE calculated was done with the unsubtracted distribution.
 *  With the subtracted distribution we are finding when the SPE peak rises and then goes back to 0.
 *  This is calculated by using the first derivative since the first derivative is most sensitive to the
 *  turn on and off than a fixed threshold.
 * -# Based on when the SPE distribution peak is found, determine the fit range and calculate
 *  all the various rates.
 ***********************************************/
void NearlineSPEDiag::CalculateRates()
{
  const double totalTime = static_cast<double>(fNEntries) * fWindow; 
  for (size_t i=0; i< fSizeOfCal; ++i) {
    if((i+1)%16 == 0){
      continue;
    }
    
    fCal.at(i)->SetPMTID(i);
   
    // Calculate background subtracted distribution
    std::shared_ptr<TH1D> peHist = fCal.at(i)->GetPEHistPtr();
    std::shared_ptr<TH1D> peHistBefore = fCalBefore.at(i)->GetPEHistPtr();
    peHistBefore->Scale(std::fabs(fWindowStart - fWindowEnd)/std::fabs(-9.92 - fWindowStart));
    peHist->Add(peHistBefore.get(),-1.0);


    // Calculate the first and second derivative with a 9 point derivative method
    std::shared_ptr<TH1D> derPEHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(peHist->Clone("derPEHist")));
    derPEHist->Reset("ICESM");
    std::shared_ptr<TH1D> secDerPEHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(peHist->Clone("secDerPEHist")));
    secDerPEHist->Reset("ICESM");
    if (peHist->Integral() != 0) {
      for (int i=0; i < derPEHist->GetNbinsX(); ++i) {
        if (i-3 < 0 || i+5 >= derPEHist->GetNbinsX()) {
          continue;
        } 
        derPEHist->SetBinContent(i+1, 1/60.0 * ( 
              peHist->GetBinContent(i+2) - peHist->GetBinContent(i) + 2.0*peHist->GetBinContent(i+3) - 
              2.0*peHist->GetBinContent(i-1) + 3.0*peHist->GetBinContent(i+4) - 3.0*peHist->GetBinContent(i-2) + 
              4.0*peHist->GetBinContent(i+5) - 4.0*peHist->GetBinContent(i-3)));
      }
      for (int i=0; i < secDerPEHist->GetNbinsX(); ++i) {
        if (i-3 < 0 || i+5 >= secDerPEHist->GetNbinsX()) {
          continue;
        } 
        secDerPEHist->SetBinContent(i+1, 1/60.0 * (
              derPEHist->GetBinContent(i+2) - derPEHist->GetBinContent(i) + 2.0*derPEHist->GetBinContent(i+3) - 
              2.0*derPEHist->GetBinContent(i-1) + 3.0*derPEHist->GetBinContent(i+4) - 3.0*derPEHist->GetBinContent(i-2) + 
              4.0*derPEHist->GetBinContent(i+5) - 4.0*derPEHist->GetBinContent(i-3)));
      }
    }

    bool above0 = false;
    double noiseWallPeak = 0;
    double noiseWallEnd = 0;
    int noiseWallPeakBin = 0;
    int noiseWallEndBin = 0;
    for (int i=0; i < derPEHist->GetNbinsX(); ++i) {
      if (derPEHist->GetBinContent(i+1) < 0 && !above0 && noiseWallPeakBin == 0) {
        above0 = true;
      } else if (derPEHist->GetBinContent(i+1)  >= 0 && above0 && noiseWallPeakBin == 0) {
        above0 = false;
        noiseWallPeak = derPEHist->GetXaxis()->GetBinCenter(i+1);
        noiseWallPeakBin = i+1;
        break;
      } else if (secDerPEHist->GetBinContent(i+1) <= 0 && !above0 && noiseWallPeakBin != 0) {
        above0 = true;
        noiseWallEnd = derPEHist->GetXaxis()->GetBinCenter(i+1);
        noiseWallEndBin = i+1;
        break;
      }
    }

    secDerPEHist->GetXaxis()->SetRange(noiseWallPeakBin,secDerPEHist->GetNbinsX());
    noiseWallEndBin = secDerPEHist->GetMaximumBin();
    noiseWallEnd = secDerPEHist->GetXaxis()->GetBinCenter(noiseWallEndBin);
    secDerPEHist->GetXaxis()->UnZoom();
    int highBin = peHist->GetNbinsX()+1;

    peHist->GetXaxis()->SetRange(noiseWallEndBin,highBin);
    int spePeak = peHist->GetMaximumBin();
    double spePeakPos = peHist->GetXaxis()->GetBinCenter(spePeak);
    double mean = peHist->GetMean();
    MsgInfo(MsgLog::Form("Noise Peak %f (bin %d content %f derivative %f), noise end %f (bin %d content %f derivative %f) spePeak %f (bin %d mean %f)",
          noiseWallPeak,noiseWallPeakBin,peHist->GetBinContent(noiseWallPeakBin),derPEHist->GetBinContent(noiseWallPeakBin),
          noiseWallEnd,noiseWallEndBin,peHist->GetBinContent(noiseWallEndBin),derPEHist->GetBinContent(noiseWallEndBin),
          spePeakPos, spePeak, mean
          ));

    peHist->GetXaxis()->UnZoom();
    spePeakPos = std::max(mean,spePeakPos);
    fCal.at(i)->SetNoisePeakPos(noiseWallPeak);
    fCal.at(i)->SetNoiseEndPos(noiseWallEnd);
    fCal.at(i)->SetSPEPeakPos(spePeakPos);

    double noiseEndContent = peHist->GetBinContent(noiseWallEndBin);
    fCal.at(i)->SetNoiseEndContent(noiseEndContent);

    int fitRangeEnd =0;
    for(int x=spePeak; x<highBin; ++x){
      float fitRangeEndContent = peHist->GetBinContent(x);
      float fitRangeEndContent2 = peHist->GetBinContent(x+1);
      float fitRangeEndContent3 = peHist->GetBinContent(x+2);
      if (noiseEndContent==0){
        noiseEndContent=peHist->GetBinContent(noiseWallEndBin+2);
      }
      if(fitRangeEndContent<noiseEndContent && fitRangeEndContent2<noiseEndContent && fitRangeEndContent3<noiseEndContent){
        fitRangeEnd = x+10;
        break;
      }
    }

    if (peHist->GetXaxis()->GetBinLowEdge(fitRangeEnd) < 100) {
      fitRangeEnd = peHist->FindBin(100);
    }
    
    float fitEndContent = peHist->GetBinContent(fitRangeEnd);
    float fitRangeEndPos = peHist->GetBinCenter(fitRangeEnd);
    float tailIntegral = peHist->Integral(fitRangeEnd, highBin);
    float fitIntegral = peHist->Integral(noiseWallEndBin,fitRangeEnd);
    float integralForRate = peHist->Integral(noiseWallEndBin,highBin);
    float noiseWallIntegral = peHist->Integral(0, noiseWallEnd);
    float pmtRate = integralForRate/totalTime;

    fCal.at(i)->SetFitEndContent(fitEndContent);
    fCal.at(i)->SetFitEndPos(fitRangeEndPos);
    fCal.at(i)->SetTailIntegral(tailIntegral);
    fCal.at(i)->SetFitIntegral(fitIntegral);
    fCal.at(i)->SetIntegralForRate(integralForRate);  
    fCal.at(i)->SetNoiseIntegral(noiseWallIntegral);
    fCal.at(i)->SetRate(pmtRate);

    fCal.at(i)->GetPEHistPtr()->GetXaxis()->UnZoom();
    
  }// end for i
         
}// end for CalculateRates function 

/*!**********************************************
 * \fn void NearlineSPEDiag::FitPEHists()
 * \brief Fit the ADC distributions created by #CalculateRates
 *
 * Fit the ADC distributions created by #CalculateRates to determine
 * the SPE values. Currently a single Gaussian is used to determine the
 * SPE value.
 ***********************************************/
void NearlineSPEDiag::FitPEHists()
{
	for(size_t i=0; i<fSizeOfCal; ++i){
		if((i+1)%16 == 0){
      continue;
    }
    std::shared_ptr<TH1D> peHist = fCal.at(i)->GetPEHistPtr();

    // Get the ranges and positions determined in #CalculateRates
		//float noiseWallEndPos = fCal.at(i)->GetNoiseEndPos();
		//float fitRangeEndPos = fCal.at(i)->GetFitEndPos();
		//float fitEndContent = fCal.at(i)->GetFitEndContent();
		//float noiseEndContent = fCal.at(i)->GetNoiseEndContent();
		//float noisePeakPos = fCal.at(i)->GetNoisePeakPos();
		//float spePeakPos = fCal.at(i)->GetSPEPeakPos();

    peHist->GetXaxis()->UnZoom();
    int totalMaxBin = peHist->GetMaximumBin();
    double totalValue = peHist->GetBinContent(totalMaxBin);
    double totalMean = peHist->GetXaxis()->GetBinCenter(totalMaxBin);

    peHist->GetXaxis()->SetRangeUser(6,100);
    int maxBin = peHist->GetMaximumBin();
    double value = peHist->GetBinContent(maxBin);
    double mean = peHist->GetXaxis()->GetBinCenter(maxBin);
    double valueLow = peHist->GetXaxis()->GetBinCenter(maxBin-6);
    double valueHigh = peHist->GetXaxis()->GetBinCenter(maxBin+6);
    peHist->GetXaxis()->UnZoom();

    double start = valueLow;
    double end = valueHigh;

    // Set the fit function string
    // Tried various functions, currently using a single Gaussian
    std::string strFunc = "gaus(0)";//+ gaus(3)";// + gaus(6)";

    MsgInfo("======================================================================================");
    MsgInfo(MsgLog::Form("Looking at PMT info %d",i));
    MsgInfo(MsgLog::Form("TotalMaxBin = %d, MaxBin = %d",totalMaxBin,maxBin));
    MsgInfo(MsgLog::Form("TotalValue = %f, Value = %f",totalValue,value));
    MsgInfo(MsgLog::Form("TotalMean = %f, Mean = %f",totalMean,mean));
    MsgInfo("======================================================================================");
   
    // Fit the distribution
    std::shared_ptr<TF1> f3 = std::make_shared<TF1>("f3", strFunc.c_str(),start,end);
    f3->SetParameters(totalValue,totalMean,std::sqrt(totalMean),value,std::sqrt(mean));
    peHist->Fit("f3","R");
    f3->SetLineColor(kGreen+2);
    
       
    // Get the fit results
    maxBin = peHist->GetMaximumBin();
    float histMean = peHist->GetMean();
    float histRMS = peHist->GetRMS();
    int histEntries = peHist->GetEntries();

    fCal.at(i)->SetPEHistMean(histMean);
    fCal.at(i)->SetPEHistRMS(histRMS);
    fCal.at(i)->SetPEHistEntries(histEntries);  


    int meanLoc = 1;
    float fitHeight= f3->GetParameter(meanLoc-1);
    float fit_mean= f3->GetParameter(meanLoc);
    float fitRMS= f3->GetParameter(meanLoc+1);

    fCal.at(i)->SetFitHeight(fitHeight);
    fCal.at(i)->SetFitMean(fit_mean);
		fCal.at(i)->SetFitRMS(fitRMS);  

    //delete f2;
  }
}
  
/*!**********************************************
 * \fn void NearlineSPEDiag::FillHist(std::string pdfVar, std::string saveParameters)
 * \brief Call the #MakeSPECalibrationHists::MakeAllHists function to create any of the histograms and save the data
 * \param[in] pdfVar The name of the output file
 * \param[in] saveParameters A list of words stating which flags to turn on
 ***********************************************/
void NearlineSPEDiag::FillHist(std::string pdfVar, std::string saveParameters)
{
  MakeSPECalibrationHists::MakeAllHists(fCal,pdfVar,fNEntries,saveParameters);
}

/*!**********************************************
 * \fn double NearlineSPEDiag::GetThreshold(int pmt)
 * \brief Return the threshold found for the \p pmt number passed
 * \param[in] pmt The PMT key
 * \return The threshold found, if the PMT key corresponds to a PMT in the detector
 ***********************************************/
double NearlineSPEDiag::GetThreshold(int pmt)
{
  if (pmt > static_cast<int>(fSizeOfCal)) {
    return 9999;
  }

  if (fCal[pmt] == nullptr) {
    return 9999;
  }

  return fCal[pmt]->GetNoiseEndPos();
}

/*!**********************************************
 * \fn double NearlineSPEDiag::GetSPE(int pmt)
 * \brief Return the SPE value found for the \p pmt number passed
 * \param[in] pmt The PMT key
 * \return The SPE value found, if the PMT key corresponds to a PMT in the detector
 ***********************************************/
double NearlineSPEDiag::GetSPE(int pmt)
{
  if (pmt > static_cast<int>(fSizeOfCal)) {
    return 9999;
  }

  if (fCal[pmt] == nullptr) {
    return 9999;
  }

  if (fCal[pmt]->GetFitMean() > 0) {
    return fCal[pmt]->GetFitMean();
  }

  return fCal[pmt]->GetSPEPeakPos();
}

/*!**********************************************
 * \fn double NearlineSPEDiag::GetRate(int pmt)
 * \brief Return the SPE rate found for the \p pmt number passed
 * \param[in] pmt The PMT key
 * \return The SPE rate found, if the PMT key corresponds to a PMT in the detector
 ***********************************************/
double NearlineSPEDiag::GetRate(int pmt)
{
  if (pmt > static_cast<int>(fSizeOfCal)) {
    return 9999;
  }

  if (fCal[pmt] == nullptr) {
    return 9999;
  }

  return fCal[pmt]->GetRate();
}

