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

long TimeAdd(long time_now, int inc, int steps)
{
    int inc_min=inc*steps/60;
    int inc_sec=inc*steps%60;
    long output=time_now+inc_min*100+inc_sec;
    if (output%100>=60)
        output+=40;
    return output;
}

//int abc(int argc, char* argv[])
//int Corr(int target, char* path, char* filename)
int FFTbackground(char* path, char* filename)
{
    //TFile *f = new TFile("test.root","UPDATE");
    
    cl_data dataBlock;
    ifstream infile;
    ifstream infile2;
    bool readStart=true;
    char char_date[80] = "";
    char PdfnameBin[80] = "";
    char Pdfname[80] = "";
    char name[80] = "";
    char PdfnameEnd[80] = "";
    char fftname[80]="";
    const int length = 1000, length2 = 250;
    int infile_ind = 0;
    double ch1[12][length], t[12][length], grid[10000], ch1_grid[12][10000];
    double ch2[12][length2];
    double average=0, peak, rms;
    double ch1_fftRe[12][length], ch1_fftIm[12][length], ch1_fftMAG[12][length];
    double ch1_fftPow[12][length], ch1_fftdBm[12][length], ch1_fftBack[12][length], ch1_fftBuf1[length], thres;
    int sort_ind[length];
    short_bytes buffer;
    double delay_time[12]={0,0.000000000,-0.000000035,0,0,0,0,0,-0.000000025,0.00000000,0.00000000,0.00000000};
    double corr=0, corr_inc=0, corr_self[12];
    int Nb=0;
    int chk=0;
    double map_time[2000][length]; //Bin every 10 seconds
    int num_time[400], map_ind=0;
    long start_time=20150323171540, time;


    //Read all the response curves
    char* input = new char[1000];
    char* p;
    
    infile2.open("/Users/judi/Desktop/Measurement/13-0079.S2P",ios::in);
    double LNA[2][401]];
    int inp_ind=0, dat_ind=0;
    while (infile2.getline(input,1000,'\n'))
    {
        p=strtok(input,"\t");
        while (p)
        {
            if (inp_ind%9==0)
                LNA[0][dat_ind]=atof(p);
            else if (inp_ind%9==3)
                LNA[1][dat_ind++]=atof(p);
            p = strtok(NULL, "\t");
            inp_ind++;
        }
    }
    infile2.close();
    ROOT::Math::Interpolator LNA_inter(401);
    LNA_inter.SetData(401,LNA[0],LNA[1]);

    
    infile2.open("/Users/judi/Desktop/Measurement/A8.csv",ios::in);
    double FILTER[2][201];
    int inp_ind=0, dat_ind=0;
    while (infile2.getline(input,1000,'\n'))
    {
        p=strtok(input,",");
        while (p)
        {
            if (inp_ind%2==0)
                FILTER[0][dat_ind]=atof(p);
            else
                FILTER[1][dat_ind++]=atof(p);
            p = strtok(NULL, ",");
            inp_ind++;
        }
    }
    infile2.close();
    ROOT::Math::Interpolator FILTER_inter(201);
    FILTER_inter.SetData(201,FILTER[0],FILTER[1]);
    
    
    infile2.open("/Users/judi/Desktop/Measurement/RG-316-CableLoss.csv",ios::in);
    double RG316[2][401];
    int inp_ind=0, dat_ind=0;
    while (infile2.getline(input,1000,'\n'))
    {
        p=strtok(input,",");
        while (p)
        {
            if (inp_ind%2==0)
                RG316[0][dat_ind]=atof(p);
            else
                RG316[1][dat_ind++]=atof(p);
            p = strtok(NULL, ",");
            inp_ind++;
        }
    }
    infile2.close();
    ROOT::Math::Interpolator RG316_inter(401);
    RG316_inter.SetData(401,RG316[0],RG316[1]);
    
    
    infile2.open("/Users/judi/Desktop/Measurement/RG-400-CableLoss.csv",ios::in);
    double RG400[2][401];
    int inp_ind=0, dat_ind=0;
    while (infile2.getline(input,1000,'\n'))
    {
        p=strtok(input,",");
        while (p)
        {
            if (inp_ind%2==0)
                RG400[0][dat_ind]=atof(p);
            else
                RG400[1][dat_ind++]=atof(p);
            p = strtok(NULL, ",");
            inp_ind++;
        }
    }
    infile2.close();
    ROOT::Math::Interpolator RG400_inter(401);
    RG400_inter.SetData(401,RG400[0],RG400[1]);
    
    
    infile2.open("/Users/judi/Desktop/Measurement/Splitter-Loss-ZX10-2-12-S1.csv",ios::in);
    double SPLITTER[2][137];
    int inp_ind=0, dat_ind=0;
    while (infile2.getline(input,1000,'\n'))
    {
        p=strtok(input,",");
        while (p)
        {
            if (inp_ind%2==0)
                SPLITTER[0][dat_ind]=atof(p)*1.0e6;
            else
                SPLITTER[1][dat_ind++]=atof(p)*-1; // Spliiter Loss has to change sign
            p = strtok(NULL, ",");
            inp_ind++;
        }
    }
    infile2.close();
    ROOT::Math::Interpolator SPLITTER_inter(137);
    SPLITTER_inter.SetData(137,SPLITTER[0],SPLITTER[1]);
    
    
    infile2.open("/Users/judi/Desktop/Measurement/SLP-300+.S2P",ios::in);
    double SLP[2][401];
    int inp_ind=0, dat_ind=0;
    while (infile2.getline(input,1000,'\n'))
    {
        p=strtok(input,"\t");
        while (p)
        {
            if (inp_ind%9==0)
                SLP[0][dat_ind]=atof(p);
            else if (inp_ind%9==3)
                SLP[1][dat_ind++]=atof(p);
            p = strtok(NULL, "\t");
            inp_ind++;
        }
    }
    infile2.close();
    ROOT::Math::Interpolator SLP_inter(401);
    SLP_inter.SetData(401,SLP[0],SLP[1]);
    
    
    float total[length];
    for (int k=0;k<50;k++)
        total[k]=0;
    for (int k=50;k<500;k++)
        total[k] = LNA_inter.Eval(k*1.0e6)+FILTER_inter.Eval(k*1.0e6)+RG316_inter.Eval(k*1.0e6)*12.55+RG400_inter.Eval(k*1.0e6)*1.2+SPLITTER_inter.Eval(k*1.0e6)+SLP_inter.Eval(k*1.0e6)-16.5*pow((k/1.0e3),2);

    
    TCanvas* c =new TCanvas("WaveForm","data",1600,1200);
    c->Divide(2,2);
    c->cd(2);

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

    //Init
    for (int k=0;k<length;k++)
        ch1_fftBuf1[k]=0;
    
    while(!infile.eof())
    {
        for(int j=0; j<12; j++) // read channels data for biniary file , ch 0-11 (12 channels);
        {
            infile.read ( (char *)(&dataBlock), sizeof(dataBlock) ); // read jth channel data
            if (j==0||j==2||j==3||j==4||j==5||j==6||j==7||j==8||j==9||j==10||j==11)
                continue; //Pick only the horizontally polarized antennae
            //if (infile_ind==target)
            time=dataBlock.globaltime;
            if (time>TimeAdd(start_time,10,map_ind))
            {
                //if (map_ind>=100) cout << "bbb" << map_ind << endl;
                num_time[map_ind]=infile_ind;
                if (infile_ind!=0)
                {
                    //if (map_ind>=100) cout << "ccc" << endl;
                    for (int k=0;k<length2/2;k++)
                    {
                        map_time[map_ind][k]=ch1_fftBuf1[k]/infile_ind-total[k*4+2];
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
            //if (map_ind>=100) cout << "eee" << infile_ind << endl;
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
            /*gr1[j] = new TGraph(length, t[j], ch1[j]);
            if(j<4)// give labels
            {
                        sprintf(name,"Global Time: %ld", start_time+20 );//sprintf(name,"GEN141421 ANT%d", j );
            }
            gr1[j]->SetTitle(name);
            gr1[j]->SetMaximum(2);
            gr1[j]->SetMinimum(-2);
            if (j==1)
            {
              gr1[j]->SetLineColor(4);
              gr1[j]->Draw("AL");
            }
            else if (j==8)
            {
              gr1[j]->SetLineColor(2);
              gr1[j]->Draw("L");
            }
            else
            {
              gr1[j]->SetLineColor(4);
              gr1[j]->Draw("L");
            }
            
            //gr1[j]->GetXaxis()->SetLimits(-1E-6,1E-6 ); //Set the range of X axis
            
            c->Update();
            
            if(readStart == true)
            { 
                c->Print(PdfnameBin,"pdf");
                readStart == false;
            }
            */
            
            //Calculate Peak-to-Average
            peak=0;
            for (int k=0;k<length;k++)
            {
                if (ch1[j][k]>peak)
                    peak=ch1[j][k];
            }
            rms=0;
            for (int k=0;k<length;k++)
                rms+=pow(ch1[j][k],2);
            rms/=length;
            //cout << "Peak-to-Average: " << pow(peak,2)/rms << "\n";
            
            if (pow(peak,2)/rms>20)
            {
                infile_ind++;
                for (int k=50;k<300;k++)
                    ch2[j][k-50]=ch1[j][k];
                TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &length2, "R2C M K");
                fft_own->SetPoints(ch2[j]);
                fft_own->Transform();
                fft_own->GetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
                for (int k=0;k<length2/2;k++)
                {
                    ch1_fftMAG[j][k]=sqrt((pow(ch1_fftRe[j][k],2)+pow(ch1_fftIm[j][k],2))/length2);
                    ch1_fftPow[j][k]=pow(ch1_fftMAG[j][k],2)/50;
                    ch1_fftdBm[j][k]=10*TMath::Log10(ch1_fftPow[j][k]);
                    ch1_fftBuf1[k]+=ch1_fftdBm[j][k];
                }
                delete fft_own;
            }
        }
    }
    
    c->cd(1);
    TH1F *h1 = new TH1F("h1","560~570 seconds",length2/2,1.0e6,500e6);
    for (int k=0;k<length2/2;k++)
        h1->SetBinContent(k,map_time[56][k]);
    h1->Draw();
    
    c->cd(2);
    TH1F *h2 = new TH1F("h2","total system response",500,1.0e6,500e6);
    for (int k=0;k<500;k++)
        h2->SetBinContent(k,total[k]);
    h2->Draw();
    
    c->cd(3);
    TH2F *h3 = new TH2F("h3","Frequency - Time map",length2/2,1.0e6,500e6,87,0,87);
    for (int k=0;k<length2/2;k++)
    {
        for (int l=0;l<87;l++)
            h3->SetBinContent(k,l,map_time[l][k]);
    }
    h3->GetXaxis()->SetTitle("Frequency (Hz)");
    h3->GetXaxis()->SetTitleSize(0.05);
    h3->GetXaxis()->SetLabelSize(0.05);
    h3->GetYaxis()->SetTitle("Time (bin #)");
    h3->GetZaxis()  ->SetRangeUser(-105, -60);
    h3->Draw("colz");
    
    c->cd(4);
    TH1F *h4 = new TH1F("h4","Number of events in each bin",90,0,90);
    for (int l=0;l<90;l++)
        h4->SetBinContent(l,num_time[l]);
    h4->Draw();
    
    c->Update();
    
    c->Write();
    
    c->Print(PdfnameEnd,"pdf");
    infile.close();
    cout << "+-----------------------------+" << endl;
    cout << "+-----------------------------+" << endl;
}