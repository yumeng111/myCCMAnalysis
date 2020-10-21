#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include "TText.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include <iostream>
#include <math.h>
#include "TLine.h"

#include "TVector.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <array>
#include <map>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include <memory>

#include "myfunc.C"
#include "MsgLog.h"

#include <sstream>
#include <numeric>

#include "TROOT.h"
#include "TH1.h"

/*
// read in all ntuple files, make one single ratio plot as f(days since oct 15)
// for each day of data make an average ratio.  plot that ratio as the day's ratio
// then take an average of daily plot (so average of a bunch of averages)
                                                                                                                                                                                   
********* there are 2 lines you need to change, if you want to include events that probably                                                                                        
    have no triplet, or ignore them.                                                                                                                                               
 

How to run through ROOT:

.x singletTripletPlots.cc++("listofinputs_plot.txt")


*/


//-------------------------------------------------------------------------------------------------
int singletTripletPlots(char *infile)
{

  // *******************************************************************************
  // open all ntuple root files, get ntuple
  // *******************************************************************************
  //  TFile *inFile = new TFile("stNtuple_10-21.root");
  
  TChain *chain = new TChain("tree"); // name of ntuple in root file is "tree"
  //TTree *chain = (TTree *) inFile -> Get("tree");

  
  std::ifstream infileStream(infile);
  std::string line= "";
  while (infileStream >> line) {
    std::cout << line << std::endl;
    chain -> Add(line.data());
  }
  infileStream.close();
  
  // print out list of files included in chain to check everything was done right
  std::cout << "list of input files: " << std::endl << std::endl;
  chain -> GetListOfFiles() -> ls();
  
  if (!chain) {
    MsgError("s/t ntuple tree was not loaded");
  }
  
  // *******************************************************************************
  // declare vars for ntuple and ntuple here
  // *******************************************************************************
  long evt;
  double singInt, tripInt, ratio, thistime;
  int singletLowBin, singletHighBin, tripletLowBin, tripletHighBin;
  int notrip;

  chain -> SetBranchAddress ("evt", &evt);
  
  chain -> SetBranchAddress ("singInt", &singInt);
  chain -> SetBranchAddress ("tripInt", &tripInt);
  chain -> SetBranchAddress ("ratio", &ratio);
  chain -> SetBranchAddress ("thistime", &thistime);

  chain -> SetBranchAddress ("singletLowBin", &singletLowBin);
  chain -> SetBranchAddress ("singletHighBin", &singletHighBin);
  chain -> SetBranchAddress ("tripletLowBin", &tripletLowBin);
  chain -> SetBranchAddress ("tripletHighBin", &tripletHighBin);
  chain -> SetBranchAddress ("notrip", &notrip);

  
  // *******************************************************************************
  // global vars
  // *******************************************************************************
  // for time conversions and daily average
  double octtime = 1571097600;
  double secToDays = 1.0 / (60 * 60 * 24); // event time in sec, 60 sec/min, 60 min/hr, 24 hr/day

  int firstevent = 0;
  int thisday = 0;
  int firstnewday;
  
  // for all events where we calculate s/t ratio
  double this_ratio = 0.0;
  double num_ratios = 0.0;

  double daily_avg = 0.0;
  double cum_avg = 0.0;
  double total_days = 0.0;
  
  int mybin;
  
  int nentries = (int) chain -> GetEntries();
  std::cout << nentries << std::endl;

  int evtsused = 0;

  int days_count[150];
  for (int i = 0; i < 150; i++) days_count[i] = 0;
  
  // *******************************************************************************
  // declare histogram for ratio, days since oct 15
  // *******************************************************************************
 
  TH1F *h_STratio = new TH1F("h_STratio", "Days since Oct 15 (x) vs Daily Avg Ratio (y)", 150, 0.5, 150.5); 

  
  // ##############################################################################
  // loop over all entries in ntuple
  // ##############################################################################
 
  for (int i = 0; i < nentries; i++) {

     chain -> GetEntry(i);

     std::cout << i << std::endl;

    // ------------------------------------------------------------------------------
    // want S/T ratio as a function of global time in DAYS since 10/15, 12:00 am
    // EPOCH time for 10/15/2019, near start of good data runs

     // some oddity with thistime being set to 0.  investigate why this happens
     // however only want events starting at Oct 15th
     if (thistime < octtime) continue;
     
    double time_sinceOct = thistime - octtime;
    double days_sinceOct = time_sinceOct * secToDays;

    // ------------------------------------------------------------------------------
    // do daily average
    int dayssince = (int) days_sinceOct;

    // std::cout << "octtime " << octtime << " thistime " << thistime << " time_since " << time_sinceOct
    //	      << " days since = " << dayssince << std::endl;
      

    // set thisday to the days since oct 15 of our first event
    // need to do this for the first event we process
    if (firstevent == 0) {
      // std::cout << "before: thisday " << thisday << " dayssince " << dayssince << std::endl;
      thisday = dayssince;
      // std::cout << "after: thisday " << thisday << " dayssince " << dayssince << std::endl;
      firstevent = 1;
    }

    
    // ------------------------------------------------------------------------------
    // for every event check if the days since oct 15 is same as first event (same day)
    // when it changes we need to clear ratio & start fresh, change value of thisday

    // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    // !!!! these are for ALL events making it into the ntuple, including ones with no apparent triplet
    //      probably easiest to just re-run code uncommenting this line, changing output pdf name

    if (notrip == 1) {
      // std::cout << "entry: " << i << " no good triplet.  days since oct 15: " << dayssince << std::endl;
     continue;
    }

    // std::cout << "entry " << i << " good triplet" << std::endl;
    // shit what if the last event in the file is notrip!  we don't print out the daily total
    // how do I get around that??
    // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    evtsused += 1;
    days_count[dayssince] += 1;
    
    if (dayssince == thisday) {
      // fill counters to calculate daily average
      this_ratio += ratio;
      num_ratios += 1;

      
      // for very last event need to also fill histo, add to global avg
      if (i == nentries - 2) {
	daily_avg = this_ratio / num_ratios;
	
	mybin = h_STratio -> GetBin(thisday); 
	h_STratio -> SetBinContent(mybin, daily_avg);
	
	// add previous day average to global counter, increment data pts by 1 so we can do a final s/t avg ratio across days
	cum_avg += daily_avg;
	std::cout << total_days << std::endl;
	total_days += 1;
	std::cout << "last entry: day " << thisday << " strobe events " << num_ratios << " daily avg ratio " << daily_avg << std::endl;
      } // if on last event in input files
      
    } // when still on same delta days since Oct 15

    
    // ------------------------------------------------------------------------------
    // called first event where we move to a new day of data
    else if (dayssince != thisday) {

      // calculate average of previous days ratio, fill plot with value
      daily_avg = this_ratio / num_ratios;

      mybin = h_STratio -> GetBin(thisday); 
      h_STratio -> SetBinContent(mybin, daily_avg);

      // print values for this day
      std::cout << "first entry: this day " << thisday << " days since " << dayssince << " strobe events " << num_ratios << " daily avg ratio " << daily_avg << std::endl;

      // add previous day average to global counter, increment data pts by 1 so we can do a final s/t avg ratio across days
      cum_avg += daily_avg;
      total_days += 1;
            
      // need to reset ratio and counters so we can do sums for current day
      daily_avg = 0;
      this_ratio = 0;
      num_ratios = 0;

      // now fill counters to calculate daily average first time we hit this new day
      this_ratio += ratio;
      num_ratios += 1;

      // change thisday to be the new number of days since oct 15
      thisday = dayssince; // changes thisday to current delta days for following events in file
     
    } // new number of days since Oct 15



    /*
    // debug purposes
    std::cout.precision(15);
    // in file I see values for ent
    // why are values printed different from what are in the ntuple?
    std::cout << "*** Ntuple Entry: " << i << " File Entry: " << evt
	      << " oct time " << octtime 
	      << " time " << thistime << " conversion " << secToDays
	      << " days since oct 15, 2019 = " << days_sinceOct 
	      << " integral of singlet = " << singInt << " integral of triplet = " << tripInt 
	      << " ratio = " << ratio << std::endl;
    */
    
    
  } // loop over all ntuple entries

    
 
  std::cout << std::endl;
  std::cout << "### Total Strobe events found in file: " << nentries << std::endl;
  std::cout << "### Strobe events used in analysis: " << evtsused << std::endl;

  for (int i = 0; i < 150; i++) {
    std::cout << "days since Oct 15: " << i << " evts used in ratio " << days_count[i] << std::endl;
  }
  
  double finalavg = cum_avg / (total_days);
  std::cout << " S/T Avg Ratio for this range of dates = " << finalavg << " (ratio / num days) = " << cum_avg << " / " << (total_days) << std::endl;

 
  // *******************************************************************************
  // format all histograms before drawing
  // *******************************************************************************
  
  prettyHist(h_STratio, "Days Since Oct 15, 2019", "Daily Avg Singlet/Triplet Ratio", 0);

  
  // *******************************************************************************
  // print all histograms to pdf file
  // *******************************************************************************
 
  gROOT -> SetStyle("Plain");
  
  gStyle -> SetTitleBorderSize(0);
  gStyle -> SetTitleX(0.08);
  gStyle -> SetTitleY(0.985);

  // if we print the stat box always use the same, with over and underflow info
  gStyle -> SetOptStat("eMRuo");
  gStyle -> SetStatX(0.975);
  gStyle -> SetStatY(0.98);
  gStyle -> SetStatH(0.08);

  // *******************************************************************************
  // analysis plot of S/T 

  // make horizontal line at my y axis average value
  TLine *stavg = new TLine(0.5, finalavg, 150.5, finalavg);
  stavg->SetLineColor(kRed);
  stavg->SetLineWidth(5);
   
  // print average on ratio plot
  char thisrat[50];
  sprintf(thisrat, "%f", finalavg);
  TText *t = new TText;
  t->SetNDC();
  t->SetTextFont(2);
  t->SetTextColor(2);
  t->SetTextSize(0.05);
  t->SetTextAlign(23); // centered in width, at top
  
  TCanvas *c1 = new TCanvas("c1", " ", 800, 800);

  c1 -> SetGrid();
  c1 -> SetGrayscale(kFALSE);

  h_STratio -> GetZaxis() -> SetLabelSize(0.01);
  h_STratio -> SetMarkerStyle(20);
  h_STratio -> SetMarkerSize(0.85);
  h_STratio -> Draw("hist p");
 
  t->DrawText(0.5, 0.5, thisrat);

  stavg->Draw();


  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // !!!! these are for ALL events making it into the ntuple, including ones with no apparent triplet
  //      probably easiest to just re-run code uncommenting this line, changing output pdf name

  // call them "all" and "tripFlag" to differentiate
  c1 -> Print("ratioPlot_tripFlag.pdf", "pdf");

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  // --------------------------------------------------------------------
  // if we want to draw p.e. threshold lines on integral plots
  // make horizontal line at 10 pe
  TLine *line0 = new TLine(-9920, 10, 6280, 10);
  line0->SetLineColor(kRed);
  
  // make horizontal line at 40 pe
  TLine *line1 = new TLine(-9920, 40, 6280, 40);
  line1->SetLineColor(kRed);
  
  // want to cut at 100, how many S/T events will this remove?
  TLine *line2 = new TLine(-9920, 90, 6280, 90);
  line2->SetLineColor(kRed);

  // *******************************************************************************
  // *******************************************************************************
 
  return EXIT_SUCCESS;

}// end of main

