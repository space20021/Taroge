#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>// for usleep(usec);
#include <iomanip> //setw
#include <cmath>
#include <algorithm>
#include <time.h>
#include "datalib.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "Math/Interpolator.h"
#include "TVirtualFFT.h"
using namespace std;

union short_bytes {
    short shorty[1000];
    char bytes[2000];
};

int Read()
{
    const int length=1000;
    double ch1_fft[length], ch1_fft_ave[300000], max_ele;
    int event_rate_bin;
    int num_tot;
    char root_file[80] = "";
    long long time, time_start;
    int time_ind=1;
    double time_sto[3000000];
    TDatime time_raw;
    int time_mday, time_hour, time_min;
    bool chk;
    
    TChain* chain = new TChain("t1");
    for (int i=1840;i<=2264;i++) //1840~2264
    {
        sprintf(root_file,"/Users/judi/Desktop/Taroge_Root/r24_%d.root",i);
        chain->Add(root_file);
    }
    for (int i=1;i<=135;i++) //1~135
    {
        sprintf(root_file,"/Users/judi/Desktop/Taroge_Root/r49_%d.root",i);
        chain->Add(root_file);
    }
    num_tot=chain->GetEntries(); //
    chain->SetBranchAddress("time",&time);
    chain->SetBranchAddress("ch1_fft[1]",&ch1_fft);
    
    TFile *f1 = new TFile("Forced_Trigger_spectra.root","recreate");
    TCanvas* c1 = new TCanvas("c1","canvas",1600,1200);
    //c1->Divide(1,2);
    //c1->cd(1);
    
    //////////////////////////
    //Start editing from here
    ////////////
    time_start=1426608381;
    event_rate_bin=0;
    for (int k=0;k<500;k++)
        ch1_fft_ave[k]=0;
    //TH2F* h1 = new TH2F("h1","Spectra -- Time",6314,1426608381,1426987221,500,0,500); //378840
    TH2F* h1 = new TH2F("h1","Time variation of spectra (dBm)",9350,1426608381,1427169436,500,0,500); //561055
    for (int i=0;i<num_tot;i++)
    {
        chain->GetEntry(i);
        time=time%100000000;
        time_mday=time/1000000;
        time=time%1000000;
        time_hour=time/10000;
        time=time%10000;
        time_min=time/100;
        time=time%100;
        time_raw.Set(115,3,time_mday,time_hour,time_min,time);
        time_sto[i]=time_raw.Convert();
        //printf( "%10f\n", time_sto[i] );
        if (time_sto[i]>time_start+60*time_ind)
        {
            if (event_rate_bin==0)
            {
                for (int k=0;k<500;k++)
                    h1->SetBinContent(time_ind,k,0);
            }
            else
            {
                for (int k=0;k<500;k++)
                {
                    ch1_fft_ave[k]/=event_rate_bin;
                    h1->SetBinContent(time_ind,k,ch1_fft_ave[k]+30.); //+30 to convert to dBm
                }
            }
            time_ind++;
            event_rate_bin=0;
        }
        for (int k=0;k<500;k++)
            ch1_fft_ave[k]+=ch1_fft[k];
        event_rate_bin++;
        
        if (i%100000==0)
            cout << "Reading the " << i << " th event...\n"; //Just to track progress
    }
    cout << "num_tot: " << num_tot << "\n";
    h1->GetXaxis()->SetTimeDisplay(1);
    h1->GetXaxis()->SetTimeFormat("%m/%d %H:%M%F1995-01-02 08:00:00");
    h1->GetXaxis()->SetTitle("Date / Time");
    h1->GetXaxis()->SetLabelSize(0.022);
    h1->GetYaxis()->SetTitle("Frequency (MHz)");
    h1->GetYaxis()->SetRange(100,320);
    h1->GetZaxis()->SetRangeUser(-75, -25);
    h1->SetStats(false);
    h1->Draw("colz");
    c1->Update();
    //h1->Write();
    //f1->Close();
}