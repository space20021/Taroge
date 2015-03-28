#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>// for usleep(usec);
#include <iomanip> //setw
#include <cmath>
#include <algorithm>
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
    double ch1_fft[length];
    
    TFile *f1 = new TFile("/Users/judi/Desktop/Taroge_Root/r49_45.root","read");
    TTree *t1 = (TTree*)f1->Get("t1");
    t1->SetBranchAddress("ch1_fft",&ch1_fft);
    t1->GetEntry(100);
    TH1F *h1 = new TH1F("h1","ch1_fft",length/2,0.0e6,500.0e6);
    for (int k=0;k<500;k++)
        h1->SetBinContent(k,ch1_fft[k]);
    TCanvas *c1;
    h1->Draw();
    //TH1F *h1 = (Th1F*)f1->Get()
}