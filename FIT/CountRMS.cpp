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

#include "showratio.cpp"

using namespace std;

double RMS_errors[48][13];

void CountRMS( const char * partition, const char * filename){

    //Reading text file with runs
    vector<int> run_1, run_2;
    int buf_1, buf_2;

    fstream runs;
    runs.open(filename, ios::in);
    runs >> buf_1 >> buf_2;

    while(!runs.eof()){
        run_1.push_back(buf_1);
        run_2.push_back(buf_2);
        runs >> buf_1 >> buf_2;
    }

    runs.close();

    cout << "run_1";
    for(int s:run_1) cout << " " << s;
    cout << endl;

    cout << "run_2";
    for(int p:run_2) cout << " " << p;
    cout << endl;

    int runmin = min( *min_element(run_1.begin(), run_1.end()), *min_element(run_2.begin(), run_2.end()) ) - 10;
    int runmax = max( *max_element(run_1.begin(), run_1.end()), *max_element(run_2.begin(), run_2.end()) ) + 10;
    int nruns = run_1.size();

    cout << "N runs: " << nruns << "\tMIN Run: " << runmin << "\tMAX Run: " << runmax << endl;
    if (run_1.size() != run_2.size()) {
        cout << "Inconsistent length of run_1 and run_2 vectors " << run_1.size() << " vs " << run_2.size() << endl;
        return;
    }
    if (nruns == 0) {
        cout << "Nothing to do" << endl;
        return;
    }

    gDirectory->GetList()->Delete(); // delete old histograms first

    //Defining histograms
    TH1F* hist_rms[48][13];
    
    TString stitle, sname;
   
   for(int row=0; row<13; ++row){
       for (int ip=0; ip<48; ++ip){
           //Average RMS
           sname.Form( "%s_pmt%2.2d_rms%2.2d_row",partition,ip+1,row+1);
           stitle.Form("%s PMT %2.2d RMS %2.2d ROW",partition,ip+1,row+1);
           hist_rms[ip][row] = new TH1F(sname, stitle, 120, -1., 5.); 
           hist_rms[ip][row]->GetXaxis()->SetTitle("Ratio RMS [%]");
       }
    }
    cout << "TEST after booking the hists" << endl;

    //Filling histograms
    for(int row=0; row<13; row++){
       for(int i=0; i<nruns; ++i){
	  //vector<double> mean;
          vector<double> rms;
          showratio(partition, run_1[i], run_2[i], row, rms);
          for (int ip=0; ip<48; ++ip){
	     if(rms[ip]>0.){
	        hist_rms[ip][row]->Fill(rms[ip]);
	     }
          }

       }
       for(int ip=0; ip<48; ip++){
          double mean = hist_rms[ip][row]->GetMean();
          if (mean > 0.) {
             TFitResultPtr r = hist_rms[ip][row]->Fit("gaus","SQ","",0.,3.*mean);
             double par1 = r->Parameter(1);
             if (par1>mean/3. && par1<mean) mean=par1;
          }
          RMS_errors[ip][row] = mean;   //Getting Mean of RMS distributions
          //hist_rms[ip]->Reset("ICESM");
       }
    }

    cout << "TEST after filling the histograms" << endl;
    hist_rms[40][0]->Draw();   
    
    //Print
    cout << "Print values for RMS errors for partition " << partition << ": " << endl; 
    for(int ip=0; ip<48; ip++){
       for(int row=0; row<13; row++){
          cout << RMS_errors[ip][row] << ", ";
       }
    cout << endl;
    }

    return;
}
