#include <cstdarg>

#include "Riostream.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLatex.h"

#include "showratio.cpp"
#include "csdate.cpp"
//#include "CountRMS.cpp"
#include "ratio_rms.cpp"

#include "AtlasUtils.h"
#include "AtlasLabels.h"
//#include "AtlasStyle.h"
#include "AtlasStyle.C"
#include "AtlasLabels.C"
#include "AtlasUtils.C"


#define iLBA 0
#define iLBC 1
#define iEBA 2
#define iEBC 3
#define nPART 4

//double RMS_errors[48][13];

const int nsetsize = 80;

//using namespace std;

int pad_pos(int n, int nx, int ny) 
{
    int n1 = (n-1)%ny;
    int n2 = (n-1)/ny;
    return n1*nx+n2+1;
}

void CountROW(const char *partition, int row_1, int row_2, int iopt=1000)
{
  //CountRMS("EBA", "EBA.list");
  const char *partname[nPART] = {"LBA","LBC","EBA","EBC"};
  int ipart=0;
  for ( ; ipart<nPART; ++ipart) {
    if (strcmp(partition,partname[ipart]) == 0) break;
  }
  if (ipart==nPART) {
    std::cout << "Wrong partitition name " << partition << std::endl;
    return;
  }

  static int runday[nPART][nsetsize], runmonth[nPART][nsetsize], runyear[nPART][nsetsize], nrun[nPART][nsetsize][6],
      solstat[nPART][nsetsize], torstat[nPART][nsetsize], badstat[nPART][nsetsize], empty[nPART][nsetsize];
  memset(runday,0,sizeof(runday));
  memset(runmonth,0,sizeof(runmonth));
  memset(runyear,0,sizeof(runyear));
  memset(nrun,0,sizeof(nrun));
  memset(solstat,0,sizeof(solstat));
  memset(torstat,0,sizeof(torstat));
  memset(badstat,0,sizeof(badstat));
  memset(empty,0,sizeof(empty));

  std::ifstream runlist;
  const char * lbrunlist = "lb_run.list"; // list of runs in LB  (the same runs for LBC and LBA)
  const char * ebrunlist = "eb_run.list"; // list of runs in EBs (different runs for EBC and EBA)

  int nmbLB = 0; 
  runlist.open(lbrunlist);
  if (!runlist.is_open()) {
    std::cout << "Can't open " << lbrunlist << std::endl;
  }
  else {
    int first;
    while (runlist >> first) {
      int m = (nmbLB)/6;
      int n = (nmbLB)%6;
      if (n == 0) {
        runday[iLBA][m] = first;
        runlist >> runmonth[iLBA][m] >> runyear[iLBA][m] >> nrun[iLBA][m][n] >> solstat[iLBA][m] >> torstat[iLBA][m] >> badstat[iLBA][m];

        nrun[iLBC][m][n] = nrun[iLBA][m][n];
        runday[iLBC][m]  = runday[iLBA][m];
        runmonth[iLBC][m]= runmonth[iLBA][m];
        runyear[iLBC][m] = runyear[iLBA][m];
        solstat[iLBC][m] = solstat[iLBA][m];
        torstat[iLBC][m] = torstat[iLBA][m];
        badstat[iLBC][m] = badstat[iLBA][m];
      }
      else {
        nrun[iLBC][m][n] = nrun[iLBA][m][n] = first;
      }
      ++nmbLB;
    }
  }
  runlist.close();

  int nmbEB = 0; 
  runlist.open(ebrunlist);
  if (!runlist.is_open()) {
    std::cout << "Can't open " << ebrunlist << std::endl;
  }
  else {
    int first, second;
    while (runlist >> first >> second) {
      int m = (nmbEB)/6;
      int n = (nmbEB)%6;
      if (n == 0) {
        runday[iEBA][m] = first;
        runmonth[iEBA][m] = second;
        runlist >> runyear[iEBA][m] >> nrun[iEBC][m][n] >> nrun[iEBA][m][n] >> solstat[iEBA][m] >> torstat[iEBA][m] >> badstat[iEBA][m];

        runday[iEBC][m] = runday[iEBA][m];
        runmonth[iEBC][m] = runmonth[iEBA][m];
        runyear[iEBC][m] = runyear[iEBA][m];
        solstat[iEBC][m] = solstat[iEBA][m];
        torstat[iEBC][m] = torstat[iEBA][m];
        badstat[iEBC][m] = badstat[iEBA][m];
      }
      else {
        nrun[iEBC][m][n] = first;   //  first column - run number for EBC
        nrun[iEBA][m][n] = second;  // second column - run number for EBA
      }
      ++nmbEB;
    }
  }
  runlist.close();

  gDirectory->GetList()->Delete(); // delete old histograms first

  //SetAtlasStyle();
  static TStyle* atlasStyle = 0;
  std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  if ( atlasStyle==0 ) atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  gStyle->SetTimeOffset(day0().Convert());
  const char * timeformat = "#splitline{%b}{%Y}";

  static const char *cellNameLB[48] =
    {"D0","A1","B1","B1","A1","A2","B2","B2","A2","A3","A3","B3",   
     "B3","D1","D1","A4","B4","B4","A4","A5","A5","B5","B5","A6",  
     "A6","D2","D2","A7","B6","B6","A7","xx","xx","A8","B7","B7",  
     "A8","A9","A9","D3","B8","B8","D3","xx","B9","B9","A10","A10"};

  static const char *cellNameEB[48] =
         {"E3", "E4","D4", "D4", "C10","C10","A12","A12","B11","B11","A13","A13",
          "E1", "E2","B12","B12","D5", "D5", "xx", "xx", "A14","A14","B13","B13",
          "xx", "xx","xx", "xx","A15","A15", "xx", "xx", "B14","B14","xx", "xx",
          "D6","D6", "xx", "xx","A16","A16","B15", "B15","xx", "xx","xx", "xx"};

  std::vector<const char **> cellnames; 
  cellnames.push_back(cellNameLB); // LBA
  cellnames.push_back(cellNameLB); // LBC
  cellnames.push_back(cellNameEB); // EBA
  cellnames.push_back(cellNameEB); // EBC
  
  //Defining histograms
  TH1F* hist_rms[48];
  TH1F* hist_mean[48];
  TProfile* hprof_rms[48];
  TGraphErrors* gr_diff[48];

  TString stitle, sname;

  for (int ip=0; ip<48; ++ip){
    //Average RMS vs date
    sname.Form( "%s_pmt%2.2d_rms_vs_date",partition,ip+1);
    stitle.Form("%s PMT %2.2d rms vs date",partition,ip+1);
    hprof_rms[ip] = new TProfile(sname, stitle, 1000, 0, cstime(1,1,20), 0., 10.);
    hprof_rms[ip]->GetXaxis()->SetTitle("Date [month and year]");
    hprof_rms[ip]->GetXaxis()->SetTimeDisplay(1);
    hprof_rms[ip]->GetXaxis()->SetTimeFormat(timeformat); 
    hprof_rms[ip]->GetXaxis()->SetLabelOffset(0.02);
    hprof_rms[ip]->GetYaxis()->SetTitle("Ratio RMS [%]");
    //Average RMS
    sname.Form( "%s_pmt%2.2d_rms",partition,ip+1);
    stitle.Form("%s PMT %2.2d RMS",partition,ip+1);
    hist_rms[ip] = new TH1F(sname, stitle, 240, -1., 11.); 
    hist_rms[ip]->GetXaxis()->SetTitle("Ratio RMS [%]");
    //Average mean
    sname.Form( "%s_pmt%2.2d_mean",partition,ip+1);
    stitle.Form("%s PMT %2.2d Mean",partition,ip+1);
    hist_mean[ip] = new TH1F(sname, stitle, 200, -10., 10.);
    hist_mean[ip]->GetXaxis()->SetTitle("Mean diff [%]");
  }

  std::cout << "After booking the hists" << std::endl;

  float mean_first[48];
  memset(mean_first,0,sizeof(float)*48);
  bool m_first[48];
  memset(m_first,-1,sizeof(bool)*48);
  std::vector<float> mean_all[48];
  std::vector<float> rms_all[48];
  std::vector<float> date_all[48];

  //Filling histograms
  for(int i=0; i<nsetsize; ++i){
    double dtim = cstime(runday[ipart][i], runmonth[ipart][i], runyear[ipart][i]);
    std::vector<double> mean;
    std::vector<double> rms;
    std::vector<int> runs;
    for (int j=0; j<6; ++j) {
      if (nrun[ipart][i][j]>0) runs.push_back(nrun[ipart][i][j]);
    }
    if (runs.size()>0) {
      showratio_row(partition, runs, row_1-1, row_2-1, mean, rms);
      for (int ip=0; ip<48; ++ip){
        if (rms[ip]>=0) {
          if (m_first[ip]) {
            m_first[ip] = false;
            mean_first[ip] = mean[ip];
          }
          hist_mean[ip]->Fill(mean[ip]-mean_first[ip]);
          hist_rms[ip]->Fill(rms[ip]);
          hprof_rms[ip]->Fill(dtim,rms[ip]);
          mean_all[ip].push_back(mean[ip]-mean_first[ip]);
          rms_all[ip].push_back(rms[ip]);
          date_all[ip].push_back(dtim);
        }
      }
    }
  }

  std::vector<float> err_x;
  std::vector<float> err_y;
  err_x.resize(nsetsize,0.0);

  for (int ip=0; ip<48; ++ip){
    sname.Form ("%s_pmt%2.2d_r%d_r%d",partition,ip+1,row_1,row_2);
    stitle.Form("%s PMT %2.2d row %d to row %d difference vs date",partition,ip+1,row_1,row_2);
    if (date_all[ip].size()) {
      std::cout << sname << " Mean ";
      for(size_t i=0; i<mean_all[ip].size();++i) std::cout << " " << mean_all[ip][i];
      std::cout << std::endl;
      std::cout << sname << " RMS  ";
      for(size_t i=0; i<rms_all[ip].size();++i) std::cout << " " << rms_all[ip][i];
      std::cout << std::endl << std::endl;
    }
    err_y.clear();
    err_y.resize(nsetsize, ratio_rms(iEBA,ip,row_1-1,row_2-1));
    TGraphErrors *gr = new TGraphErrors(date_all[ip].size(),&(date_all[ip][0]),&(mean_all[ip][0]),&(err_x[0]),&(err_y[0]));
    gr->SetName(sname);
    gr->SetTitle(stitle);
    gr->GetXaxis()->SetTitle("Date [month and year]");
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat(timeformat); 
    gr->GetXaxis()->SetLabelOffset(0.02); 
    gr->GetYaxis()->SetTitle("Mean diff [%]");
    gr_diff[ip] = gr;
  }

  std::cout << "After filling the histograms" << std::endl;

  //Output

  //create all canvas
  std::vector<TCanvas *> c1;
  std::vector<TCanvas *> c2;
  std::vector<TCanvas *> c3;
  std::vector<TCanvas *> c4;
  TCanvas *c = 0;
  
  int nx=3;
  int ny=2;
  int nxy=nx*ny;
  int ic;
  int id;

  if ((iopt>>0)%10) {
    ic=0;
    id=999;
    for (int ip=0; ip<48; ++ip){
      if (date_all[ip].size()) {
        if (id>=nxy) {
          ic+=1;
          id=0;
          sname.Form("All_pmt_ratio_p%d",ic);
          stitle.Form("All PMTs Ratio page %d",ic);
          c = new TCanvas(sname,stitle,10,10,1000,800);
          c->Divide(nx,ny);
          c1.push_back(c);
        }
        c->cd(pad_pos(++id,nx,ny));
        hist_mean[ip]->Draw();
      }
    }
  }

  if ((iopt>>1)%10) {
    ic=0;
    id=999;
    for (int ip=0; ip<48; ++ip){
      if (date_all[ip].size()) {
        if (id>=nxy) {
          ic+=1;
          id=0;
          sname.Form("All_pmt_RMS_p%d",ic);
          stitle.Form("All PMTs RMS page %d",ic);
          c = new TCanvas(sname,stitle,10,10,1000,800);
          c->Divide(nx,ny);
          c2.push_back(c);
        }
        c->cd(pad_pos(++id,nx,ny));
        hist_rms[ip]->Draw();
      }
    }
  }

  if ((iopt>>2)%10) {
    ic=0;
    id=999;
    for (int ip=0; ip<48; ++ip){
      if (date_all[ip].size()) {
        if (id>=nxy) {
          ic+=1;
          id=0;
          sname.Form("All_pmt_RMS_vs_date_p%d",ic);
          stitle.Form("All PMTs RMS vs date page %d",ic);
          c = new TCanvas(sname,stitle,10,10,1000,800);
          c->Divide(nx,ny);
          c3.push_back(c);
        }
        c->cd(pad_pos(++id,nx,ny));
        TProfile * h=hprof_rms[ip];
        h->SetMarkerStyle(20);
        h->SetMarkerColor(2);
        h->SetMarkerSize(1.0);
        h->SetMaximum(10.0);
        h->Draw("P");
        
      }
    }
  }
  
  if ((iopt>>3)%10) {
    ic=0;
    id=999;
    for (int ip=0; ip<48; ++ip){
      if (date_all[ip].size()) {
        if (id>=nxy) {
          ic+=1;
          id=0;
          sname.Form("All_pmt_ratio_vs_date_p%d",ic);
          stitle.Form("All PMTs ratio vs date page %d",ic);
          c = new TCanvas(sname,stitle,10,10,1000,800);
          c->Divide(nx,ny);
          c4.push_back(c);
        }
        c->cd(pad_pos(++id,nx,ny));
        TGraphErrors * h=gr_diff[ip];
        if (date_all[ip][date_all[ip].size()-1]-date_all[ip][0] < 9*365*24*3600) 
          h->GetXaxis()->SetNdivisions(505); 
        h->SetMarkerStyle(20);
        h->SetMarkerColor(2);
        h->SetMarkerSize(1.0);
        h->SetMinimum(-6.0);
        h->SetMaximum(4.0);
        h->Draw("AP");

        //ATLASLabel(0.40, 0.86, "Preliminary", 1);
        ATLASLabel(0.40, 0.86, "   Internal", 1);
        myText(0.40, 0.81, 1, "Tile Calorimeter");
        sname.Form("%s %s pmt %d",partition,cellnames[ipart][ip],ip+1);
        myText(0.21, 0.26, 1, sname);
        sname.Form("Row %d - Row %d",row_1,row_2);
        myText(0.21, 0.21, 1, sname);
      }
    }
  }

  //save all canvas
  std::cout << "Saving Canvas..." << std::endl;
  for (size_t i=0; i<c1.size(); ++i) {
    sname.Form("c1_%s_r%d_r%d_p%lu.png",partition,row_1,row_2,i+1);
    c1[i]->SaveAs(sname);
  }
  for (size_t i=0; i<c2.size(); ++i) {
    sname.Form("c2_%s_r%d_r%d_p%lu.png",partition,row_1,row_2,i+1);
    c2[i]->SaveAs(sname);
  }
  for (size_t i=0; i<c3.size(); ++i) {
    sname.Form("c3_%s_r%d_r%d_p%lu.png",partition,row_1,row_2,i+1);
    c3[i]->SaveAs(sname);
  }
  for (size_t i=0; i<c4.size(); ++i) {
    sname.Form("c4_%s_r%d_r%d_p%lu.png",partition,row_1,row_2,i+1);
    c4[i]->SaveAs(sname);
  }

  return;
}
