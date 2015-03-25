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

//int abc(int argc, char* argv[])
//int Corr(int target, char* path, char* filename)
int FFTmap(char* path)
{
    TFile *f = new TFile("test2.root","UPDATE");
    
    cl_data dataBlock;
    ifstream infile;
    bool readStart=true;
    char char_date[80] = "";
    char PdfnameBin[80] = "";
    char Pdfname[80] = "";
    char name[80] = "";
    char PdfnameEnd[80] = "";
    char fftname[80]="";
    const int length = 1000;
    int infile_ind = 0;
    double ch1[12][length], t[12][length], grid[10000], ch1_grid[12][10000];
    double average=0;
    double ch1_fftRe[12][length], ch1_fftIm[12][length], ch1_fftMAG[12][length];
    double ch1_fftPow[12][length], ch1_fftdBm[12][length], ch1_fftBack[12][length], ch1_fftBuf1[length], thres;
    int sort_ind[length];
    short_bytes buffer;
    double delay_time[12]={0,0.000000000,-0.000000035,0,0,0,0,0,-0.000000025,0.00000000,0.00000000,0.00000000};
    double corr=0, corr_inc=0, corr_self[12];
    int Nb=0;
    int chk=0;
    double map_time[5000][length]; //Bin every 10 seconds
    int num_time[5000], map_ind=0;
    long start_time=20150323000032, time;
    int FileNames[500]={45 , 46 , 47 , 48 , 49 , 50 , 51 , 52 , 53 , 54 , 55 , 56 , 57 , 58 , 59 , 60 , 61 , 62 , 63 , 64 , 65 , 66 , 67 , 68 , 69 , 70 , 71 , 72 , 73 , 74 , 75 , 76 , 77 , 78 , 79 , 80 , 81 , 82 , 83 , 84 , 85 , 86 , 87 , 88 , 89 , 90 , 91 , 92 , 93 , 94 , 95 , 96 , 97 , 98 , 99 , 100 , 101 , 102 , 103 , 104 , 105 , 106 , 107 , 108 , 109 , 110 , 111 , 112 , 113 , 114 , 115 , 116 , 117 , 118 , 119 , 120 , 121 , 122 , 123 , 124 , 125 , 126 , 127 , 128 , 129 , 130 , 131 , 132 , 133 , 134 , 135}, FileInd=0;

    TCanvas* c =new TCanvas("WaveForm","data",1600,1200);
    c->Divide(2,2);

    sprintf(char_date,"%sr49_%d.bin", path, FileNames[FileInd]);
    cout <<"opening file: " <<char_date <<endl; // show the name of data
    infile.open (char_date, ios::in | ios::binary);// open biniary file
    if (infile.good()==false)
    {
        cout << "Cannot open file" << endl;
        return 0;
    }
    
    //Init
    for (int k=0;k<length;k++)
        ch1_fftBuf1[k]=0;
    
    while(true)
    {
        //Open next file when the previous one ends
        if (infile.eof())
        {
            infile.close();
            if (FileNames[++FileInd]==0)
                break;
            sprintf(char_date,"%sr49_%d.bin", path, FileNames[FileInd]);
            cout <<"opening file: " <<char_date <<endl; // show the name of data
            infile.open (char_date, ios::in | ios::binary);
            if (infile.good()==false)
            {
                cout << "Cannot open file" << endl;
                return 0;
            }
        }
        
        for(int j=0; j<12; j++) // read channels data for biniary file , ch 0-11 (12 channels);
        {
            infile.read ( (char *)(&dataBlock), sizeof(dataBlock) ); // read jth channel data
            if (j==0||j==2||j==3||j==4||j==5||j==6||j==7||j==8||j==9||j==10||j==11)
                continue; //Pick only the horizontally polarized antennae
            //if (infile_ind==target)
            time=dataBlock.globaltime;
            if (time>TimeAdd(start_time,60*map_ind))
            {
                cout << "map_ind: " << map_ind << ", infile_ind: " << infile_ind << endl;
                //cout << "time: " << time << endl;
                num_time[map_ind]=infile_ind;
                if (infile_ind!=0)
                {
                    //if (map_ind>=100) cout << "ccc" << endl;
                    for (int k=0;k<length;k++)
                    {
                        map_time[map_ind][k]=ch1_fftBuf1[k]/infile_ind;
                        ch1_fftBuf1[k]=0;
                    }
                }
                else
                {
                    //if (map_ind>=100) cout << "ddd" << endl;
                    for (int k=0;k<length;k++)
                        map_time[map_ind][k]=0;
                }
                infile_ind=0;
                map_ind++;
            }
            infile_ind++;
            //cout << start_time+20 << endl;
            //cout << setprecision(14) << dataBlock.globaltime <<endl;
            for (int k=0; k<length; k++)
            {
                buffer.bytes[2*k+1] = dataBlock.WaveForm[2*k];//===========================WAVEFORM
                buffer.bytes[2*k] = dataBlock.WaveForm[2*k+1];//===========================WAVEFORM
                ch1[j][k] = (float)buffer.shorty[k]*0.04*dataBlock.v_Scale;
                t[j][k] = (float)k *dataBlock.h_scale*1E-2 +dataBlock.h_Position+delay_time[j];
            }
            //Subtract average
            average=0;
            for (int k=0; k<length; k++)
                average+=ch1[j][k];
            average/=length;
            for (int k=0; k<length; k++)
                ch1[j][k]-=average;
            
            TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &length, "R2C M K");
            fft_own->SetPoints(ch1[j]);
            fft_own->Transform();
            fft_own->GetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
            for (int k=0;k<length/2;k++)
            {
                ch1_fftMAG[j][k]=sqrt((pow(ch1_fftRe[j][k],2)+pow(ch1_fftIm[j][k],2))/length);
                ch1_fftPow[j][k]=pow(ch1_fftMAG[j][k],2)/50;
                ch1_fftdBm[j][k]=10*TMath::Log10(ch1_fftPow[j][k]);
                ch1_fftBuf1[k]+=ch1_fftdBm[j][k];
            }
            delete fft_own;
        }
    }
    
    c->cd(1);
    TH1F *h1 = new TH1F("h1","1st bin",length/2,1.0e6,500e6);
    for (int k=0;k<length/2;k++)
        h1->SetBinContent(k,map_time[1][k]);
    h1->Draw();
    
    c->cd(2);
    TH1F *h2 = new TH1F("h2","50th bin",length/2,1.0e6,500e6);
    for (int k=0;k<length/2;k++)
        h2->SetBinContent(k,map_time[50][k]);
    h2->Draw();
    
    c->cd(3);
    TH2F *h3 = new TH2F("h3","Frequency - Time map",length/2,1.0e6,500e6,5000,0,5000);
    for (int k=0;k<length/2;k++)
    {
        for (int l=0;l<5000;l++)
            h3->SetBinContent(k,l,map_time[l][k]);
    }
    h3->GetZaxis()->SetRangeUser(-70, -20);
    h3->Draw("colz");
    
    c->cd(4);
    TH1F *h4 = new TH1F("h4","Number of events in each bin",5000,0,5000);
    for (int l=0;l<5000;l++)
        h4->SetBinContent(l,num_time[l]);
    h4->Draw();
    
    c->Update();
    
    c->Write();
    
    c->Print(PdfnameEnd,"pdf");
    cout << "+-----------------------------+" << endl;
    cout << "+-----------------------------+" << endl;
}