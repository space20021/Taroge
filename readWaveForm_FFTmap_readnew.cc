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

long TimeAdd(long time_now, int inc)
{
    int inc_day=inc/86400;
    inc=inc%86400;
    int inc_hr=inc/3600;
    inc=inc%3600;
    int inc_min=inc/60;
    inc=inc%60;
    int inc_sec=inc;
    long output=time_now+inc_day*1000000+inc_hr*10000+inc_min*100+inc_sec;
    if (output%100>=60)
        output+=40;
    if (output%10000>=6000)
        output+=4000;
    if (output%1000000>=240000)
        output+=760000;
    return output;
}

int Read()
{
    const int length=1000;
    double ch1_fft[length], ch1_fft_143[3000000], max_ele, ch1_fft_143_red[1000000];
    int red_ind=0;
    int event_rate[3000000];
    int num_tot;
    char root_file[80] = "";
    long long time;
    double time_sto[3000000], time_sto_red[1000000];
    TDatime time_raw;
    int time_mday, time_hour, time_min;
    bool chk;
    int event_am_143=0, event_am_tot=0;
    int event_pm_143=0, event_pm_tot=0;
    
    TChain* chain = new TChain("t1");
    for (int i=10;i<=59;i++) //10~59
    {
        sprintf(root_file,"/Users/judi/Desktop/Taroge_Root/r57_%d.root",i);
        chain->Add(root_file);
    }
    num_tot=chain->GetEntries(); //250050
    chain->SetBranchAddress("time",&time);
    chain->SetBranchAddress("ch1_fft[1]",&ch1_fft);
    
    TCanvas* c1 = new TCanvas("c1","canvas",1600,1200);
    c1->Divide(1,2);
    c1->cd(1);
    TH1F* h1 = new TH1F("h1","Event Rate -- Time",193678,1427462911,1427656589);
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
        h1->Fill(time_sto[i]);
        max_ele=ch1_fft[140];
        for (int k=1;k<=5;k++)
        {
            if (ch1_fft[140+k]>max_ele)
                max_ele=ch1_fft[140+k];
        }
        ch1_fft_143[i]=max_ele;
        chk=1;
        for (int k=0;k<length/2-1;k++)
        {
            if (ch1_fft_143[i]<ch1_fft[k])
            {
                chk=0;
                break;
            }
        }
        if (time_hour>=23||time_hour<11) // 7~19 -8 hours = -1~11
        {
            event_am_tot++;
            if (chk==1)
            {
                event_am_143++;
                time_sto_red[red_ind]=time_sto[i];
                ch1_fft_143_red[red_ind]=max_ele;
                red_ind++;
            }
        }
        else
        {
            event_pm_tot++;
            if (chk==1)
            {
                event_pm_143++;
                time_sto_red[red_ind]=time_sto[i];
                ch1_fft_143_red[red_ind]=max_ele;
                red_ind++;
            }
        }
        if (i%100000==0)
            cout << "Reading the " << i << " th event...\n"; //Just to track progress
    }
    cout << "event_am_143: " << event_am_143 << "\n";
    cout << "event_am_tot: " << event_am_tot << "\n";
    cout << "event_pm_143: " << event_pm_143 << "\n";
    cout << "event_pm_tot: " << event_pm_tot << "\n";
    cout << "num_tot: " << num_tot << "\n";
    h1->SetMarkerStyle(24);
    h1->SetMarkerColorAlpha(1,0.004);
    h1->GetXaxis()->SetTimeDisplay(1);
    h1->GetXaxis()->SetTimeFormat("%m/%d %H:%M%F1995-01-02 08:00:00");
    h1->GetXaxis()->SetTitle("Date / Time");
    h1->GetXaxis()->SetLabelSize(0.022);
    h1->GetYaxis()->SetTitle("(Hz)");
    h1->Draw("P");
    c1->cd(2);
    TGraph* gr1 = new TGraph(num_tot,time_sto,ch1_fft_143);
    gr1->SetTitle("141~145 MHz Peak Magnitude -- Time");
    gr1->SetMarkerStyle(1);
//    gr1->SetMarkerSize(0.08);
    gr1->GetXaxis()->SetTimeDisplay(1);
    gr1->GetXaxis()->SetTimeFormat("%m/%d %H:%M%F1995-01-02 08:00:00");
    gr1->GetXaxis()->SetLimits(1427462911,1427656589);
    gr1->GetXaxis()->SetTitle("Date / Time");
    gr1->GetXaxis()->SetLabelSize(0.022);
    gr1->GetYaxis()->SetTitle("Power (dBm)");
    gr1->Draw("AP");
    TGraph* gr2 = new TGraph(red_ind,time_sto_red,ch1_fft_143_red);
    gr2->SetMarkerStyle(1);
    gr2->SetMarkerColor(2);
    gr2->GetXaxis()->SetTimeDisplay(1);
    gr2->GetXaxis()->SetTimeFormat("%m/%d %H:%M%F1995-01-02 08:00:00");
    gr2->GetXaxis()->SetLimits(1427462911,1427656589);
    gr2->GetXaxis()->SetTitle("Date / Time");
    gr2->GetXaxis()->SetLabelSize(0.022);
    gr2->Draw("P");
    c1->Update();
}