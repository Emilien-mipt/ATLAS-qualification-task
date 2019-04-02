#include <cstdarg>
//#include "Riostream.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <vector>
#include <string>

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPaveText.h"

#define TEST(cond) if (cond) { std::cout << "Exiting due to errors" << std::endl; return 1; } else
    
using namespace std;

//extern double RMS_errors[48][13];

//Fuction for reading integrals
int readint(const char* filename, const char* partition, float integrals[65][48][13], int problems[65][48][13]){

  memset(integrals,0,sizeof(float)*65*48*13);
  memset(problems,0,sizeof(int)*65*48*13);
  
  TFile f(filename,"READ");
  TDirectory *d = f.GetDirectory("");
  TTree *t = (TTree*)d->Get("integrals");
  if (!t) {
    cout << "No such tree: integrals in " << filename << endl; 
    return 1;
  }

  stringstream ssbranch;
  ssbranch << partition;
  TBranch *br1 = t->FindBranch(ssbranch.str().c_str());
  if (!br1) {
    cout << "No such branch: " << ssbranch.str() << " in " << filename << endl; 
    return 1;
  }
  br1->SetAddress(integrals);

  ssbranch << "prob";
  TBranch *br2 = t->FindBranch(ssbranch.str().c_str());
  if (!br2) {
    cout << "No such branch: " << ssbranch.str() << " in " << filename << endl; 
    return 1;
  }
  br2->SetAddress(problems);

  if (t->GetEntries()==1) {                                                    
    t->GetEntry(0);
  } else {
    cout << "Summary tree should contain 1 entry!!!" << endl; 
    return 1;
  }

  return 0;
}

//High Voltage
int readhv(const char* filename, const char* partition, int module, double hv[48], double temperature[7])
{
  memset(hv,0,sizeof(double)*48);
  memset(temperature,0,sizeof(double)*48);

  TFile f(filename,"READ");
  f.cd("RawData");

  stringstream ssname;
  ssname << "DCS.HV." << partition << setfill('0') << setw(2) << module;
  TTree *t = (TTree*) gDirectory->Get(ssname.str().c_str());
  if (!t) {
    cout << "No tree " << ssname.str() << " in " << filename << endl; 
    return 1;
  }

  TBranch *b1 = t->FindBranch("pmt");
  b1->SetAddress(hv);
  TBranch *b2 = t->FindBranch("temp");
  b2->SetAddress(temperature);

  if (t->GetEntries()==1) {                                                    
    t->GetEntry(0);
    //cout << "Read HV for " << ssname.str().substr(7) << " from file " << filename << endl; 
  } else {
      //cout << "ERROR - HV TTree for " << ssname.str().substr(7) << " from file " << filename << " has " << t->GetEntries() << " entries" << endl; 
    return 1;
  }

  return 0;
}

int run2year(int run) {
  int year=2000;
  if      (run<=  196) year=2008;
  else if (run<= 2258) year=2009;
  else if (run<= 3603) year=2010;
  else if (run<= 5082) year=2011;
  else if (run<= 6302) year=2012;
  else if (run<= 8059) year=2013;
  else if (run<= 9292) year=2014;
  else if (run<= 9846) year=2015;
  else if (run<=10396) year=2016;
  else if (run<=10526) year=2017;
  else if (run<=12000) year=2018;
  return year;
}

string data_file(int run) 
{
  stringstream ss;
  ss << "/data/cs/" << run2year(run) << "/data/cs" << setfill('0') << setw(5) << run << ".root";
  return ss.str();
}

string integral_file(int run) 
{
  stringstream ss;
  ss << "/data/cs/" << run2year(run) << "/int/summary/integrals_" << setfill('0') << setw(5) << run << ".root";
  return ss.str();
}

string modname0(const char *partition, int mod0) 
{
  stringstream ss;
  ss << partition << setfill('0') << setw(2) << mod0+1;
  return ss.str();
}







//
///////////////////////////////////////////////////SHOWRATIO///////////////////////////////////////////////////////////
int showratio(const char *partition, int run_1, int run_2, int row, vector<double> & rms)
{
  rms.clear();

  bool EB=(partition[0]=='E'); //Whether extended barrel or not

  //Struction for runs 
  struct Run
  {
    float integrals[65][48][13];
    int problems[65][48][13];
    double hv[65][48];
    double temp[65][7];
    bool good[65][48][13];
  };

  Run R1, R2; //Run_1 and Run_2
  string file;

  //Reading Run_1
  file = integral_file(run_1);
  TEST(readint(file.c_str(), partition, R1.integrals, R1.problems)) {
      cout << "Integrals have been read from file " << file << endl;
  }

  //Reading Run_2
  file = integral_file(run_2);
  TEST(readint(file.c_str(), partition, R2.integrals, R2.problems)) {
      cout << "Integrals have been read from file " << file << endl;
  }

  //Array which contains the ratio of values 
  float int_ratio[65][48][13];
  bool goodR1[65];
  bool goodR2[65];
  bool goodR[65];
  memset(int_ratio,0,sizeof(float)*65*48*13);
  memset(goodR1,0,sizeof(bool)*65);
  memset(goodR2,0,sizeof(bool)*65);
  memset(goodR,0,sizeof(bool)*65);
  for (int imod = 0; imod < 65; ++imod){
    for (int ip = 0; ip < 48; ++ip){
      for (int row = 0; row < 13; ++row){
        R1.good[imod][ip][row] = (R1.problems[imod][ip][row]<16 && R1.integrals[imod][ip][row]>300);
        R2.good[imod][ip][row] = (R2.problems[imod][ip][row]<16 && R2.integrals[imod][ip][row]>300);
        if (R1.good[imod][ip][row]) goodR1[imod] = true;
        if (R2.good[imod][ip][row]) goodR2[imod] = true;
        if (R1.good[imod][ip][row] && R2.good[imod][ip][row]) {
          goodR[imod] = true;
          int_ratio[imod][ip][row] = R1.integrals[imod][ip][row]/R2.integrals[imod][ip][row];
        }
      }
    }
    if (goodR[imod]) {
        //cout << modname0(partition,imod);
      for (int ip = 0; ip < 48; ++ip){
          //cout << " " << setw(5) << int_ratio[imod][ip][row];
      }
      //cout << endl;
    } else {
        //cout << modname0(partition,imod);
        //if (!goodR1[imod]) cout << " no data in run "<< run_1;
        //if (!goodR2[imod]) cout << " no data in run "<< run_2;
        //cout << endl;
    }
  }

  string file1 = data_file(run_1);
  for (int imod = 0; imod < 65; ++imod){
    if (goodR1[imod]) readhv(file1.c_str(),partition,imod+1,R1.hv[imod],R1.temp[imod]);
  }
  string file2 = data_file(run_2);
  for (int imod = 0; imod < 65; ++imod){
    if (goodR2[imod]) readhv(file2.c_str(),partition,imod+1,R2.hv[imod],R2.temp[imod]);
  }
  
  // book and fill all histograms
  //gDirectory->GetList()->Delete(); // delete old histograms first

  TH1F *hist[48];
  TString stitle, sname;

  for (int ip=0; ip<48; ip++) {
    sname.Form("h%2.2d",ip+1);
    stitle.Form("%s PMT %2.2d integrals ratio",partition,ip+1);
    TObject * h1 = gROOT->FindObject(sname);
    if (h1) h1->Delete();
    hist[ip] = new TH1F(sname, stitle, 200, -10., 10.);
    for(int imod=0; imod<65; imod++){
      if (goodR[imod]) {
        hist[ip]->Fill((int_ratio[imod][ip][row]-1.)*100.);
      }
    }
    if(hist[ip]->GetEntries()>0){
       rms.push_back(hist[ip]->GetRMS());
    }
    else{
       rms.push_back(-1.);
    }
  }    

  
  return 0;
}
///////////////////////////////////////////////////SHOWRATIO//////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////SHOWRATIO_ROW//////////////////////////////////////////////////////////////////
int showratio_row(const char *partition, vector<int> & runs, int row1, int row2, vector<double> & mean, vector<double> & rms)
{
  mean.clear();
  rms.clear();

  bool EB=(partition[0]=='E'); //Whether extended barrel or not

  //Struction for runs 
  struct Run
  {
    float integrals[65][48][13];
    int problems[65][48][13];
    double hv[65][48];
    double temp[65][7];
    bool good[65][48][13];
  };

  float int_ratio[65][48];
  bool goodR[65][48];
  memset(int_ratio,0,sizeof(float)*65*48);
  memset(goodR,0,sizeof(bool)*65*48);

  // for (size_t ir=0; ir<std::min(3,runs.size()); ++ir) {
  // for (size_t ir=max(3,runs.size())-3; ir<runs.size(); ++ir) {
  for (size_t ir=0; ir<runs.size(); ++ir) {
    int run=runs[ir];
      
    //Reading Run
    Run R1;
    string file = integral_file(run);
    TEST(readint(file.c_str(), partition, R1.integrals, R1.problems)) {
        cout << "Integrals have been read from file " << file << endl;
    }
    
    //Array which contains the ratio of values 
    for (int imod = 0; imod < 65; ++imod){
      for (int ip = 0; ip < 48; ++ip){
        for (int row = 0; row < 13; ++row){
          R1.good[imod][ip][row] = (R1.problems[imod][ip][row]<16 && R1.integrals[imod][ip][row]>300);
        }
        if (R1.good[imod][ip][row1] && R1.good[imod][ip][row2]) {
       //   float r1=R1.integrals[imod][ip][row1]/R1.integrals[imod][ip][row2];
         // if (goodR[imod][ip]) {
         //   float v1 = (r1-1.);
         //   float v2 = (int_ratio[imod][ip]-1.);
         //   if (fabs(v1)<fabs(v2)) {
         //     int_ratio[imod][ip] = r1;
         //   } 
         // } else {
         //   int_ratio[imod][ip] = r1;
         //}
          int_ratio[imod][ip] = R1.integrals[imod][ip][row1]/R1.integrals[imod][ip][row2];
          goodR[imod][ip] = true;
        }
      }//PMT loop
    }//Modules loop
  }//Runs loop
  
  // book and fill all histograms
  TH1F *hist[48];
  TString stitle, sname;

  for (int ip=0; ip<48; ip++) {
    sname.Form("h%2.2d",ip+1);
    stitle.Form("%s PMT %2.2d integrals ratio",partition,ip+1);
    TObject * h1 = gROOT->FindObject(sname);
    if (h1) h1->Delete();
    hist[ip] = new TH1F(sname, stitle, 60, -30., 30.);
    for(int imod=0; imod<65; imod++){
      if (goodR[imod][ip]) {
        hist[ip]->Fill((int_ratio[imod][ip]-1.)*100.);
      }
    }
    //cout << "FILLING the hists" << endl;
    mean.push_back(hist[ip]->GetMean());
    if (hist[ip]->GetEntries() > 0) {
      rms.push_back(hist[ip]->GetRMS());
      //cout << "PMT: " << ip+1 << " RMS: " << hist[ip]->GetRMS() << endl;
    } else {
      rms.push_back(-0.1);
    }
   }
  
  //Drawing
  //create all canvas
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(111110);
    gStyle->SetOptFit(1111);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.15);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetPalette(1);

    sname.Form("All_pmt_ratio_mb%d", 1);
    stitle.Form("All PMTs Ratio motherboard %d", 1);

    TCanvas* C = new TCanvas(sname,stitle,10,10,1000,800);
    C->Divide(4,3);
  
    vector<int> EB_filled;
 
    for(int ip=0; ip<48; ip++){
	if(hist[ip]->GetEntries() > 0) EB_filled.push_back(ip);
    }
     
    cout << EB_filled.size() << endl;
 
    for (int cnt=0; cnt<10; ++cnt) {
        C->cd(cnt+1);
 	hist[EB_filled[cnt]]->Draw();
        hist[EB_filled[cnt]]->Fit("gaus");
    }
    cout << "Saving Canvas with hists..." << endl;
    sname.Form("c_%s_run_%5.5d_r%d_r%d.png",partition, runs[runs.size()-1], row1+1, row2+1);
    C->SaveAs(sname);
    delete C;
  return 0;
}
///////////////////////////////////////////////////SHOWRATIO_ROW//////////////////////////////////////////////////////////////////
/*

http://zenis.dnp.fmph.uniba.sk/tile.html

Barrel:

0 / 2 - cell A  (10 cells, 20 PMTs)
3 / 5 - cell B  (9 cells, 18 PMTs)
6 / 8 - cell C  (8 cells, 16 PMTs)

Ext.barrel

0 / 2 - cell A  (5 cells, 10 PMTs)
3 / 5 - cell B  (5 cells, 10 PMTs)
4 / 6 - cell B  (5 cells, 10 PMTs)

*/
