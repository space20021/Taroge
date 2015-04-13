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
int Tree()
{
    cl_data dataBlock;
    ifstream infile;
    bool readStart=true;
    char path[80] = "/Users/judi/Desktop/Taroge_Data/";
    char char_date[80] = "", char_date_root[80] = "";
    char name[80] = "";
    char fftname[80]="";
    const int length = 1000;
    int infile_ind = 0;
    double ch1[12][length], t[12][length], grid[10000], ch1_grid[12][10000];
    double average=0;
    double ch1_fftRe[12][length], ch1_fftIm[12][length], ch1_fftMAG[12][length];
    double ch1_fftPow[12][length], ch1_fftdBm[12][length], ch1_fftBack[12][length], thres;
    int sort_ind[length];
    short_bytes buffer;
    double delay_time[12]={0,0.000000000,-0.000000035,0,0,0,0,0,-0.000000025,0.00000000,0.00000000,0.00000000};
    double corr=0, corr_inc=0, corr_self[12];
    int Nb=0;
    int chk=0;
    double map_time[5000][length]; //Bin every 10 seconds
    int num_time[5000], map_ind=0;
    long start_time=20150323000032, time;
    Double_t ch1_Buf[length], ch1_fftBuf[length];
    
    char FileNames[200][80]={"r57_10","r57_11","r57_12","r57_13","r57_14","r57_15","r57_16","r57_17","r57_18","r57_19","r57_20","r57_21","r57_22","r57_23","r57_24","r57_25","r57_26","r57_27","r57_28","r57_29","r57_30","r57_31","r57_32","r57_33","r57_34","r57_35","r57_36","r57_37","r57_38","r57_39","r57_40","r57_41","r57_42","r57_43","r57_44","r57_45","r57_46","r57_47","r57_48","r57_49","r57_50","r57_51","r57_52","r57_53","r57_54","r57_55","r57_56","r57_57","r57_58","r57_59"};
    int FileInd=0;

    sprintf(char_date,"%s%s.bin", path, FileNames[FileInd]);
    cout <<"opening file: " <<char_date <<endl; // show the name of data
    infile.open (char_date, ios::in | ios::binary);// open biniary file
    if (infile.good()==false)
    {
        cout << "Cannot open file" << endl;
        return 0;
    }
    
    //Init
    for (int k=0;k<length;k++)
        ch1_fftBuf[k]=0;
    sprintf(char_date_root,"/Users/judi/Desktop/Taroge_Root/%s.root", FileNames[0]);
    TFile *f1 = new TFile(char_date_root,"recreate");
    TTree *t1 = new TTree("t1","TTree for each event");
    t1->Branch("time",&time,"time/L");
    t1->Branch("length",&length,"length/I");
    t1->Branch("ch1[1]",&ch1_Buf,"ch1_Buf[length]/D");
    t1->Branch("ch1_fft[1]",&ch1_fftBuf,"ch1_fftBuf[length]/D");
    
    while(true)
    {
        //Open next file when the previous one ends
        if (infile.eof())
        {
            infile.close();
            
            t1->Write();
            delete t1;
            f1->Close();
            
            delete f1;
            //f1->Open(char_date_root);
            
            FileInd++;
            
            sprintf(char_date,"%s%s.bin", path, FileNames[FileInd]);
            cout <<"opening file: " <<char_date <<endl; // show the name of data
            infile.open (char_date, ios::in | ios::binary);
            if (infile.good()==false)
            {
                cout << "Cannot open file" << endl;
                return 0;
            }
            
            sprintf(char_date_root,"/Users/judi/Desktop/Taroge_Root/%s.root", FileNames[FileInd]);
            TFile *f1 = new TFile(char_date_root,"recreate");
            TTree *t1 = new TTree("t1","TTree for each event");
            t1->Branch("time",&time,"time/L");
            t1->Branch("length",&length,"length/I");
            t1->Branch("ch1[1]",&ch1_Buf,"ch1_Buf[length]/D");
            t1->Branch("ch1_fft[1]",&ch1_fftBuf,"ch1_fftBuf[length]/D");
        }
        
        for(int j=0; j<12; j++) // read channels data for biniary file , ch 0-11 (12 channels);
        {
            infile.read ( (char *)(&dataBlock), sizeof(dataBlock) ); // read jth channel data
            if (j==0||j==2||j==3||j==4||j==5||j==6||j==7||j==8||j==9||j==10||j==11)
                continue; //Pick only the horizontally polarized antennae
            //if (infile_ind==target)
            time=dataBlock.globaltime;
            
            
            //cout << start_time+20 << endl;
            //cout << setprecision(14) << dataBlock.globaltime <<endl;
            for (int k=0; k<length; k++)
            {
                buffer.bytes[2*k+1] = dataBlock.WaveForm[2*k];//===========================WAVEFORM
                buffer.bytes[2*k] = dataBlock.WaveForm[2*k+1];//===========================WAVEFORM
                ch1[j][k] = (float)buffer.shorty[k]*0.04*dataBlock.v_Scale;
                t[j][k] = (float)k *dataBlock.h_scale*1E-2 +dataBlock.h_Position+delay_time[j];
            }
            for (int k=0; k<length; k++)
                ch1_Buf[k]=ch1[j][k];
            
            TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &length, "R2C M K");
            fft_own->SetPoints(ch1[j]);
            fft_own->Transform();
            fft_own->GetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
            for (int k=0;k<length/2;k++)
            {
                ch1_fftMAG[j][k]=sqrt((pow(ch1_fftRe[j][k],2)+pow(ch1_fftIm[j][k],2)))/length;
                ch1_fftPow[j][k]=pow(ch1_fftMAG[j][k],2)/50;
                ch1_fftdBm[j][k]=10*TMath::Log10(ch1_fftPow[j][k]);
                ch1_fftBuf[k]=ch1_fftdBm[j][k];
            }
            t1->Fill();
            delete fft_own;
        }
    }
    
    cout << "+-----------------------------+" << endl;
}