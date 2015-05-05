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
#include "TSpectrum.h"
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
int TreeNext()
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
    double ch1_fftPow[12][length], thres;
    int sort_ind[length];
    short_bytes buffer;
    double delay_time[12]={0,0.000000000,-0.000000035,0,0,0,0,0,0,0,0,0};
    double corr=0, corr_inc=0, corr_self[12];
    int Nb=0;
    int chk=0;
    double map_time[5000][length]; //Bin every 10 seconds
    int num_time[5000], map_ind=0;
    long start_time=20150323000032, time;
    Double_t ch1_fftdBm[12][length], max_ele, ch1_back[12][length];
    
    Double_t avg=0;
    Float_t background[500];
    
    Bool_t max_143;
    Double_t ch1_1[1000], ch1_fft1[1000], ch1_back1[1000];
    Double_t ch1_3[1000], ch1_fft3[1000], ch1_back3[1000];
    Double_t ch1_4[1000], ch1_fft4[1000], ch1_back4[1000];
    Double_t ch1_8[1000], ch1_fft8[1000], ch1_back8[1000];
    Double_t ch1_9[1000], ch1_fft9[1000], ch1_back9[1000];
    Double_t ch1_11[1000], ch1_fft11[1000], ch1_back11[1000];
    
    char FileNames[1000][80]={"r24_1924"};
    
    //char FileNames[1000][80]={"r24_1840","r24_1841","r24_1842","r24_1843","r24_1844","r24_1845","r24_1846","r24_1847","r24_1848","r24_1849","r24_1850","r24_1851","r24_1852","r24_1853","r24_1854","r24_1855","r24_1856","r24_1857","r24_1858","r24_1859","r24_1860","r24_1861","r24_1862","r24_1863","r24_1864","r24_1865","r24_1866","r24_1867","r24_1868","r24_1869","r24_1870","r24_1871","r24_1872","r24_1873","r24_1874","r24_1875","r24_1876","r24_1877","r24_1878","r24_1879","r24_1880","r24_1881","r24_1882","r24_1883","r24_1884","r24_1885","r24_1886","r24_1887","r24_1888","r24_1889","r24_1890","r24_1891","r24_1892","r24_1893","r24_1894","r24_1895","r24_1896","r24_1897","r24_1898","r24_1899","r24_1900","r24_1901","r24_1902","r24_1903","r24_1904","r24_1905","r24_1906","r24_1907","r24_1908","r24_1909","r24_1910","r24_1911","r24_1912","r24_1913","r24_1914","r24_1915","r24_1916","r24_1917","r24_1918","r24_1919","r24_1920","r24_1921","r24_1922","r24_1923","r24_1924","r24_1925","r24_1926","r24_1927","r24_1928","r24_1929","r24_1930","r24_1931","r24_1932","r24_1933","r24_1934","r24_1935","r24_1936","r24_1937","r24_1938","r24_1939","r24_1940","r24_1941","r24_1942","r24_1943","r24_1944","r24_1945","r24_1946","r24_1947","r24_1948","r24_1949","r24_1950","r24_1951","r24_1952","r24_1953","r24_1954","r24_1955","r24_1956","r24_1957","r24_1958","r24_1959","r24_1960","r24_1961","r24_1962","r24_1963","r24_1964","r24_1965","r24_1966","r24_1967","r24_1968","r24_1969","r24_1970","r24_1971","r24_1972","r24_1973","r24_1974","r24_1975","r24_1976","r24_1977","r24_1978","r24_1979","r24_1980","r24_1981","r24_1982","r24_1983","r24_1984","r24_1985","r24_1986","r24_1987","r24_1988","r24_1989","r24_1990","r24_1991","r24_1992","r24_1993","r24_1994","r24_1995","r24_1996","r24_1997","r24_1998","r24_1999","r24_2000","r24_2001","r24_2002","r24_2003","r24_2004","r24_2005","r24_2006","r24_2007","r24_2008","r24_2009","r24_2010","r24_2011","r24_2012","r24_2013","r24_2014","r24_2015","r24_2016","r24_2017","r24_2018","r24_2019","r24_2020","r24_2021","r24_2022","r24_2023","r24_2024","r24_2025","r24_2026","r24_2027","r24_2028","r24_2029","r24_2030","r24_2031","r24_2032","r24_2033","r24_2034","r24_2035","r24_2036","r24_2037","r24_2038","r24_2039","r24_2040","r24_2041","r24_2042","r24_2043","r24_2044","r24_2045","r24_2046","r24_2047","r24_2048","r24_2049","r24_2050","r24_2051","r24_2052","r24_2053","r24_2054","r24_2055","r24_2056","r24_2057","r24_2058","r24_2059","r24_2060","r24_2061","r24_2062","r24_2063","r24_2064","r24_2065","r24_2066","r24_2067","r24_2068","r24_2069","r24_2070","r24_2071","r24_2072","r24_2073","r24_2074","r24_2075","r24_2076","r24_2077","r24_2078","r24_2079","r24_2080","r24_2081","r24_2082","r24_2083","r24_2084","r24_2085","r24_2086","r24_2087","r24_2088","r24_2089","r24_2090","r24_2091","r24_2092","r24_2093","r24_2094","r24_2095","r24_2096","r24_2097","r24_2098","r24_2099","r24_2100","r24_2101","r24_2102","r24_2103","r24_2104","r24_2105","r24_2106","r24_2107","r24_2108","r24_2109","r24_2110","r24_2111","r24_2112","r24_2113","r24_2114","r24_2115","r24_2116","r24_2117","r24_2118","r24_2119","r24_2120","r24_2121","r24_2122","r24_2123","r24_2124","r24_2125","r24_2126","r24_2127","r24_2128","r24_2129","r24_2130","r24_2131","r24_2132","r24_2133","r24_2134","r24_2135","r24_2136","r24_2137","r24_2138","r24_2139","r24_2140","r24_2141","r24_2142","r24_2143","r24_2144","r24_2145","r24_2146","r24_2147","r24_2148","r24_2149","r24_2150","r24_2151","r24_2152","r24_2153","r24_2154","r24_2155","r24_2156","r24_2157","r24_2158","r24_2159","r24_2160","r24_2161","r24_2162","r24_2163","r24_2164","r24_2165","r24_2166","r24_2167","r24_2168","r24_2169","r24_2170","r24_2171","r24_2172","r24_2173","r24_2174","r24_2175","r24_2176","r24_2177","r24_2178","r24_2179","r24_2180","r24_2181","r24_2182","r24_2183","r24_2184","r24_2185","r24_2186","r24_2187","r24_2188","r24_2189","r24_2190","r24_2191","r24_2192","r24_2193","r24_2194","r24_2195","r24_2196","r24_2197","r24_2198","r24_2199","r24_2200","r24_2201","r24_2202","r24_2203","r24_2204","r24_2205","r24_2206","r24_2207","r24_2208","r24_2209","r24_2210","r24_2211","r24_2212","r24_2213","r24_2214","r24_2215","r24_2216","r24_2217","r24_2218","r24_2219","r24_2220","r24_2221","r24_2222","r24_2223","r24_2224","r24_2225","r24_2226","r24_2227","r24_2228","r24_2229","r24_2230","r24_2231","r24_2232","r24_2233","r24_2234","r24_2235","r24_2236","r24_2237","r24_2238","r24_2239","r24_2240","r24_2241","r24_2242","r24_2243","r24_2244","r24_2245","r24_2246","r24_2247","r24_2248","r24_2249","r24_2250","r24_2251","r24_2252","r24_2253","r24_2254","r24_2255","r24_2256","r24_2257","r24_2258","r24_2259","r24_2260","r24_2261","r24_2262","r24_2263","r24_2264"};
    
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
    sprintf(char_date_root,"/Users/judi/Desktop/Taroge_Root/%s.root", FileNames[0]);
    TFile *f1 = new TFile(char_date_root,"recreate");
    TTree *t1 = new TTree("t1","TTree for each event");
    t1->Branch("time",&time,"time/L");
    t1->Branch("length",&length,"length/I");
    t1->Branch("ch1_1",ch1_1,"ch1_1[length]/D");
    t1->Branch("ch1_3",ch1_3,"ch1_3[length]/D");
    t1->Branch("ch1_4",ch1_4,"ch1_4[length]/D");
    t1->Branch("ch1_8",ch1_8,"ch1_8[length]/D");
    t1->Branch("ch1_9",ch1_9,"ch1_9[length]/D");
    t1->Branch("ch1_11",ch1_11,"ch1_11[length]/D");
    t1->Branch("ch1_fft1",ch1_fft1,"ch1_fft1[length]/D");
    t1->Branch("ch1_fft3",ch1_fft3,"ch1_fft3[length]/D");
    t1->Branch("ch1_fft4",ch1_fft4,"ch1_fft4[length]/D");
    t1->Branch("ch1_fft8",ch1_fft8,"ch1_fft8[length]/D");
    t1->Branch("ch1_fft9",ch1_fft9,"ch1_fft9[length]/D");
    t1->Branch("ch1_fft11",ch1_fft11,"ch1_fft11[length]/D");
    t1->Branch("max_143",&max_143,"max_143/O");
    t1->Branch("ch1_back1",ch1_back1,"ch1_back1[length]/D");
    t1->Branch("ch1_back3",ch1_back3,"ch1_back3[length]/D");
    t1->Branch("ch1_back4",ch1_back4,"ch1_back4[length]/D");
    t1->Branch("ch1_back8",ch1_back8,"ch1_back8[length]/D");
    t1->Branch("ch1_back9",ch1_back9,"ch1_back9[length]/D");
    t1->Branch("ch1_back11",ch1_back11,"ch1_back11[length]/D");
    
    int jizz=1;
    int jiz=1;
    
    while(true)
    {
        jiz+=1;
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
            t1->Branch("ch1_1",ch1_1,"ch1_1[length]/D");
            t1->Branch("ch1_3",ch1_3,"ch1_3[length]/D");
            t1->Branch("ch1_4",ch1_4,"ch1_4[length]/D");
            t1->Branch("ch1_8",ch1_8,"ch1_8[length]/D");
            t1->Branch("ch1_9",ch1_9,"ch1_9[length]/D");
            t1->Branch("ch1_11",ch1_11,"ch1_11[length]/D");
            t1->Branch("ch1_fft1",ch1_fft1,"ch1_fft1[length]/D");
            t1->Branch("ch1_fft3",ch1_fft3,"ch1_fft3[length]/D");
            t1->Branch("ch1_fft4",ch1_fft4,"ch1_fft4[length]/D");
            t1->Branch("ch1_fft8",ch1_fft8,"ch1_fft8[length]/D");
            t1->Branch("ch1_fft9",ch1_fft9,"ch1_fft9[length]/D");
            t1->Branch("ch1_fft11",ch1_fft11,"ch1_fft11[length]/D");
            t1->Branch("max_143",&max_143,"max_143/O");
            t1->Branch("ch1_back1",ch1_back1,"ch1_back1[length]/D");
            t1->Branch("ch1_back3",ch1_back3,"ch1_back3[length]/D");
            t1->Branch("ch1_back4",ch1_back4,"ch1_back4[length]/D");
            t1->Branch("ch1_back8",ch1_back8,"ch1_back8[length]/D");
            t1->Branch("ch1_back9",ch1_back9,"ch1_back9[length]/D");
            t1->Branch("ch1_back11",ch1_back11,"ch1_back11[length]/D");
        }
        
        for(int j=0; j<12; j++) // read channels data for biniary file , ch 0-11 (12 channels);
        {
            infile.read ( (char *)(&dataBlock), sizeof(dataBlock) ); // read jth channel data
            if (j==0||j==2||j==4||j==5||j==6||j==7||j==8||j==9||j==10||j==11)
                continue;   //Pick the top and bottom H-pol antenna: 1 & 3
                            //H-pol: 0, 1, 2, 3, 9, 10
                            //V-pol: 4, 5, 6, 7, 8, 11
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
            TVirtualFFT* fft_own = TVirtualFFT::FFT(1, &length, "R2C M K");
            fft_own->SetPoints(ch1[j]);
            fft_own->Transform();
            fft_own->GetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
            for (int k=0;k<length/2;k++)
            {
                ch1_fftMAG[j][k]=sqrt((pow(ch1_fftRe[j][k],2)+pow(ch1_fftIm[j][k],2)))/length*2;
                ch1_fftPow[j][k]=pow(ch1_fftMAG[j][k],2)/50;
                ch1_fftdBm[j][k]=10*TMath::Log10(ch1_fftPow[j][k])+30;
            }
            
            
            //Cut off 0MHz, 115MHz, and 124MHz
            /*ch1_fftRe[j][0]=0;
            ch1_fftRe[j][115]=0;
            ch1_fftRe[j][123]=0;
            ch1_fftRe[j][124]=0;
            ch1_fftRe[j][125]=0;
            ch1_fftIm[j][0]=0;
            ch1_fftIm[j][115]=0;
            ch1_fftIm[j][123]=0;
            ch1_fftIm[j][124]=0;
            ch1_fftIm[j][125]=0;
            ch1_fftdBm[j][0]=-100;
            ch1_fftdBm[j][115]=-100;
            ch1_fftdBm[j][123]=-100;
            ch1_fftdBm[j][124]=-100;
            ch1_fftdBm[j][125]=-100;*/
            
            //Check whether 143 MHz is the largest component (139~145)
            if (j==1)
            {
                max_ele=ch1_fftdBm[j][139];
                for (int k=1;k<=6;k++)
                {
                    if (ch1_fftdBm[j][139+k]>max_ele)
                        max_ele=ch1_fftdBm[j][139+k];
                }
                max_143=1;
                for (int k=100;k<=320;k++)
                {
                    if (max_ele<ch1_fftdBm[j][k])
                    {
                        max_143=0;
                        break;
                    }
                }
            }
            
            if (max_143==0&&jiz>100)
            {
                //Apply peak-finding algorithm
                //Smooth
                TH1F* h1 = new TH1F("h1","ch1_fftdBm",500,0,500);
                for (int k=0;k<length/2;k++)
                {
                    if (ch1_fftdBm[j][k]<-60)
                        h1->SetBinContent(k,-60);
                    else
                        h1->SetBinContent(k,ch1_fftdBm[j][k]);
                }
                TH1F* h2 = (TH1F*) h1->Clone(); //h2:"smooth"
                h2->Smooth(3);
                //Subtract average (HAS TO BE LARGER THAN 2.5~3)
                TH1F* h3 = (TH1F*) h2->Clone(); //h3:"average"
                for (int k=0;k<500;k++)
                {
                    avg=0;
                    avg+=ch1_fftdBm[j][k];
                    if (k==0)
                    {
                        avg=avg+ch1_fftdBm[j][k+1]+ch1_fftdBm[j][k+2]+ch1_fftdBm[j][k+3];
                        avg/=4;
                    }
                    else if (k==1)
                    {
                        avg=avg+ch1_fftdBm[j][k-1]+ch1_fftdBm[j][k+1]+ch1_fftdBm[j][k+2]+ch1_fftdBm[j][k+3];
                        avg/=5;
                    }
                    else if (k==2)
                    {
                        avg=avg+ch1_fftdBm[j][k-2]+ch1_fftdBm[j][k-1]+ch1_fftdBm[j][k+1]+ch1_fftdBm[j][k+2]+ch1_fftdBm[j][k+3];
                        avg/=6;
                    }
                    else if (k==497)
                    {
                        avg=avg+ch1_fftdBm[j][k-3]+ch1_fftdBm[j][k-2]+ch1_fftdBm[j][k-1]+ch1_fftdBm[j][k+1]+ch1_fftdBm[j][k+2];
                        avg/=6;
                    }
                    else if (k==498)
                    {
                        avg=avg+ch1_fftdBm[j][k-3]+ch1_fftdBm[j][k-2]+ch1_fftdBm[j][k-1]+ch1_fftdBm[j][k+1];
                        avg/=5;
                    }
                    else if (k==499)
                    {
                        avg=avg+ch1_fftdBm[j][k-3]+ch1_fftdBm[j][k-2]+ch1_fftdBm[j][k-1];
                        avg/=4;
                    }
                    else
                    {
                        avg=avg+ch1_fftdBm[j][k-3]+ch1_fftdBm[j][k-2]+ch1_fftdBm[j][k-1]+ch1_fftdBm[j][k+1]+ch1_fftdBm[j][k+2]+ch1_fftdBm[j][k+3];
                        avg/=7;
                    }
                    h3->SetBinContent(k,avg);
                }
                TH1F* h4 = new TH1F("h4","smooth-average",500,0,500);
                for (int k=0;k<500;k++)
                    h4->SetBinContent(k,h2->GetBinContent(k)-h3->GetBinContent(k));
                TSpectrum* s1 = new TSpectrum();
                TCanvas* c1 = new TCanvas();
                TH1F* h5 = (TH1F*) h2->Clone(); //h5:"background"
                for (int k=0;k<500;k++)
                    background[k]=h2->GetBinContent(k);
                s1->Background(background,500,5,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing3,kFALSE);
                for (int k=0;k<500;k++)
                    h5->SetBinContent(k,background[k]);
                delete s1;
                TH1F* h6 = new TH1F("h6","smooth-background",500,0,500);
                for (int k=0;k<500;k++)
                    h6->SetBinContent(k,h2->GetBinContent(k)-background[k]);
                TH1F* h7 = new TH1F("h7","(smooth-average)+(smooth-background)",500,0,500);
                for (int k=0;k<500;k++)
                    h7->SetBinContent(k,h4->GetBinContent(k)+h6->GetBinContent(k));
                TSpectrum* s1 = new TSpectrum();
                Int_t npeaks=s1->Search(h7,2.5,"goff",0.5);
                Float_t* xpeaks=s1->GetPositionX();
                Float_t* ypeaks=s1->GetPositionY();
                delete s1;
                TH1F* h8 = new TH1F("h8","just for painting",500,0,500);
                for (int k=0;k<500;k++)
                    h8->SetBinContent(k,-100);
                for (int k=0;k<npeaks;k++)
                {
                    if (ypeaks[k]>0)
                        h5->SetBinContent(h5->GetBin(xpeaks[k]),h1->GetBinContent(h1->GetBin(xpeaks[k])));
                }
                h1->SetLineColor(4);
                h1->Draw();
                h3->SetLineColor(2);
                h3->Draw("same");
                h4->SetLineColor(2);
                h4->Draw("same");
                h8->SetMarkerStyle(29);
                h8->SetMarkerSize(2);
                h8->Draw("same P");
                for (int k=0;k<npeaks;k++)
                    cout << xpeaks[k] << "\n";
                return 0;
            }
            
            
            //Backward FFT
            TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &length, "C2R M K");
            fft_back->SetPointsComplex(ch1_fftRe[j],ch1_fftIm[j]);
            fft_back->Transform();
            fft_back->GetPoints(ch1_back[j]);
            for (int k=0;k<length;k++)
                ch1_back[j][k]/=length;
            
            delete fft_own;
            delete fft_back;
        }
        
        //Store to temporary array and fill to tree
        for (int k=0;k<length;k++)
        {
            ch1_1[k]=ch1[1][k];
            ch1_3[k]=ch1[3][k];
            ch1_4[k]=ch1[4][k];
            ch1_8[k]=ch1[8][k];
            ch1_9[k]=ch1[9][k];
            ch1_11[k]=ch1[11][k];
            ch1_fft1[k]=ch1_fftdBm[1][k];
            ch1_fft3[k]=ch1_fftdBm[3][k];
            ch1_fft4[k]=ch1_fftdBm[4][k];
            ch1_fft8[k]=ch1_fftdBm[8][k];
            ch1_fft9[k]=ch1_fftdBm[9][k];
            ch1_fft11[k]=ch1_fftdBm[11][k];
            ch1_back1[k]=ch1_back[1][k];
            ch1_back3[k]=ch1_back[3][k];
            ch1_back4[k]=ch1_back[4][k];
            ch1_back8[k]=ch1_back[8][k];
            ch1_back9[k]=ch1_back[9][k];
            ch1_back11[k]=ch1_back[11][k];
        }
        if (jizz==1)
        {
            t1->Fill();
            jizz=0;
        }
    }
}