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

//int abc(int argc, char* argv[])
int Corr(int target, char* path, char* filename)
//int Corr(char* path, char* filename)
{
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
    double average=0, peak=0, rms=0;
    double ch1_fftRe[12][length], ch1_fftIm[12][length], ch1_fftMAG[12][length];
    double ch1_fftPow[12][length], ch1_fftdBm[12][length], ch1_fftBack[12][length], ch1_fftBuf1[length], thres;
    int sort_ind[length];
    short_bytes buffer;
    double delay_time[12]={0,0.000000000,-0.000000035,0,0,0,0,0,-0.000000025,0.00000000,0.00000000,0.00000000};
    double corr=0, corr_inc=0, corr_self[12];
    int Nb=0;
    int chk=0;
    double map_time[100][length]; //Bin every 10 seconds
    long time, map_ind=0;
    long start_time=20150228010000;
    float jizz_t[250], jizz_ch1[250];

    TCanvas* c =new TCanvas("WaveForm","data",1600,1200);
    c->Divide(2,2);
    c->cd(1);
    TGraph* gr1[12],gr2[12];

    sprintf(char_date,"%s%s.bin", path, filename );

    sprintf(PdfnameBin,"Corr%s.pdf(", filename );
    sprintf(Pdfname,"Corr%s.pdf", filename );
    sprintf(PdfnameEnd,"Corr%s.pdf)", filename ); //cout<<"PDFnameEnd="<<PdfnameEnd<<endl;

    cout <<"opened file is " <<char_date <<endl; // show the name of data 
    //usleep(1.E6);

    infile.open (char_date, ios::in | ios::binary);// open biniary file
    cout << char_date <<endl;
    if (infile.good()==false){
    cout << "Cannot open file" << endl;
    return 0;
    }
    int eventCounter=1;

    while(!infile.eof())
    {
        infile_ind++;
        for(int j=0; j<12; j++) // read channels data for biniary file , ch 0-11 (12 channels);
        {
            infile.read ( (char *)(&dataBlock), sizeof(dataBlock) ); // read jth channel data
            time=dataBlock.globaltime;
            if (j==0||j==3||j==4||j==5||j==6||j==7||j==9||j==10||j==11)
                continue; //Pick only the horizontally polarized antennae
            if (infile_ind==target)
            //if ((time>start_time+10*map_ind)&&(time<start_time+10*map_ind+10))
            {
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
                //plot WaveForm figure to j-th pannels
                gr1[j] = new TGraph(length, t[j], ch1[j]);
                
                //Calculate Peak-to-Average
                peak=0;
                for (int k=0;k<length;k++)
                {
                    if (ch1[j][k]>peak)
                        peak=ch1[j][k];
                }
                cout << "peak: " << peak << "\n";
                rms=0;
                for (int k=0;k<length;k++)
                    rms+=pow(ch1[j][k],2);
                rms/=length;
                cout << "rms: " << rms << "\n";
                cout << "Peak-to-Average: " << pow(peak,2)/rms << "\n";
                
        
                
                if(j<4)// give labels
                {
                            sprintf(name,"Global Time: %ld", time );//sprintf(name,"GEN141421 ANT%d", j );
                }
                gr1[j]->SetTitle(name);
                gr1[j]->SetMaximum(2);
                gr1[j]->SetMinimum(-2);
                if (j==1)
                {
                  gr1[j]->SetLineColor(4);
                  gr1[j]->Draw("AL");
                  /*
                  for (int k=50;k<300;k++)
                  {
                      jizz_t[k-50]=t[j][k];
                      jizz_ch1[k-50]=ch1[j][k];
                  }
                  
                  gr2[j] = new TGraph(250,jizz_t,jizz_ch1);
                  gr2[j]->SetLineColor(2);
                  gr2[j]->Draw("AL");8/
                }
                /*else if (j==8)
                {
                  gr1[j]->SetLineColor(2);
                  gr1[j]->Draw("L");
                }
                else
                {
                  gr1[j]->SetLineColor(4);
                  gr1[j]->Draw("L");
                }*/
                
                //gr1[j]->GetXaxis()->SetLimits(-1E-6,1E-6 ); //Set the range of X axis
                
                
            }
        }
        
    }

    //FFT
    for (int j=0; j<12; j++)
    {
        if (j==0||j==3||j==4||j==5||j==6||j==7||j==9||j==10||j==11)
            continue; //Pick only the horizontally polarized antennae
        TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &length, "R2C M K");
        fft_own->SetPoints(ch1[j]);
        fft_own->Transform();
        fft_own->GetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
        cout <<"Zzsfzsdf\n"
        for (int i=0;i<length/2;i++)
        {
            ch1_fftMAG[j][i]=sqrt((pow(ch1_fftRe[j][i],2)+pow(ch1_fftIm[j][i],2))/length);
            ch1_fftPow[j][i]=pow(ch1_fftMAG[j][i],2)/50;
            ch1_fftdBm[j][i]=10*TMath::Log10(ch1_fftPow[j][i]);
        }
        if (j==1)
        {
            c->cd(2);
            TH1F *h1 = new TH1F("h1","Forward FFT",length/2,1.0e6,500e6);
            for (int i=0;i<length/2;i++)
                h1->SetBinContent(i,ch1_fftdBm[j][i]);
            h1->Draw();
        }
    }
    c->Update();
    cout <<"adsfasdf\n"
    //Filter spectrum
    //3/9: Use high-pass filter. Also compare the CW frequencies across all recorded events
    for (int j=0; j<12; j++)
    {
        if (j==0||j==3||j==4||j==5||j==6||j==7||j==9||j==10||j==11)
            continue; //Pick only the horizontally polarized antennae
        if (j==1)
        {
            c->cd(2*j+1);
            TH1F *h2 = new TH1F("h2","jadfasdklfjalsdkfj",length,1.0e6,1000e6);
            for (int i=0;i<length;i++)
                h2->SetBinContent(i,ch1_fftMAG[j][i]);
            //h2->Draw();
        }
        //Backward FFT
        TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &length, "C2R M K");
        fft_back->SetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
        fft_back->Transform();
        fft_back->GetPoints(ch1_fftBack[j]);
        for (int i=0;i<length;i++)
            ch1_fftBack[j][i]/=length;
        if (j==1)
        {
            c->cd(2*j+2);
            TH1F *h3 = new TH1F("h3","Backward FFT",length,1.0e-9,1.0e-6);
            for (int i=0;i<length;i++)
                h3->SetBinContent(i,ch1_fftBack[j][i]);
            h3->Draw();
        }
    }
    c->Update();
    
    //Re-grid
    for (int i=0; i<10000; i++)
        grid[i]=double(i)*0.0000000002;
    for (int j=0; j<12; j++)
    {
        if (j==0||j==2||j==3||j==4||j==5||j==6||j==7||j==8||j==9||j==10||j==11)
            continue; //Pick only the horizontally polarized antennae
        ROOT::Math::Interpolator interp(length);
        interp.SetData(length,t[j],ch1[j]);
        for (int i=0; i<10000; i++)
        {
            if (grid[i]<t[j][0]||grid[i]>t[j][length-1]):
                ch1_grid[j][i]=0;
            else
                ch1_grid[j][i]=interp.Eval(grid[i]);
        }
    }
    /*TGraph* gr2 = new TGraph(10000, grid, ch1_grid[8]);
    gr2->SetLineColor(4);
    gr2->Draw("*");*/
    
    /*
    //Calculate self-correlation: Integrate (sqrt(v_i^2))
    for (int j=0; j<12; j++)
    {
        if (j==0||j==3||j==4||j==5||j==6||j==7||j==9||j==10||j==11)
            continue; //Pick only the horizontally polarized antennae
        corr_self[j]=0;
        for (int i=0; i<10000; i++)
            corr_self[j]+=ch1_grid[j][i]*ch1_grid[j][i];
        corr_self[j]=sqrt(corr_self[j]);
    }
    //Calculate cross-correlation M(r) = Sum( C_ij(r) ) / Nb
    //                            C_ij(r) = (v_i * v_j) / (sqrt(v_i^2) * sqrt(v_j^2))
    for (int j=0; j<12; j++)
    {
        if (j==0||j==3||j==4||j==5||j==6||j==7||j==9||j==10||j==11)
            continue; //Pick only the horizontally polarized antennae
        for (int k=j+1; k<12; k++)
        {
            if (k==0||k==3||k==4||k==5||k==6||k==7||k==9||k==10||k==11)
                continue; //Pick only the horizontally polarized antennae
            Nb+=1;
            corr_inc=0;
            for (int i=0; i<10000; i++)
                corr_inc+=ch1_grid[j][i]*ch1_grid[k][i];
            corr_inc/=corr_self[j];
            corr_inc/=corr_self[k];
            corr+=corr_inc;
        }
    }
    corr/=Nb;
    */
    infile.close();
    cout << "+-----------------------------+" << endl;
    //cout << "| Cross-Correlation = " << corr << " |" << endl;
    cout << "+-----------------------------+" << endl;
}