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
    //char FileNames[200][80]={"r49_1","r49_2","r49_3","r49_4","r49_5","r49_6","r49_7","r49_8","r49_9","r49_10","r49_11","r49_12","r49_13","r49_14","r49_15","r49_16","r49_17","r49_18","r49_19","r49_20","r49_21","r49_22","r49_23","r49_24","r49_25","r49_26","r49_27","r49_28","r49_29","r49_30","r49_31","r49_32","r49_33","r49_34","r49_35","r49_36","r49_37","r49_38","r49_39","r49_40","r49_41","r49_42","r49_43","r49_44","r49_45","r49_46","r49_47","r49_48","r49_49","r49_50","r49_51","r49_52","r49_53","r49_54","r49_55","r49_56","r49_57","r49_58","r49_59","r49_60","r49_61","r49_62","r49_63","r49_64","r49_65","r49_66","r49_67","r49_68","r49_69","r49_70","r49_71","r49_72","r49_73","r49_74","r49_75","r49_76","r49_77","r49_78","r49_79","r49_80","r49_81","r49_82","r49_83","r49_84","r49_85","r49_86","r49_87","r49_88","r49_89","r49_90","r49_91","r49_92","r49_93","r49_94","r49_95","r49_96","r49_97","r49_98","r49_99","r49_100","r49_101","r49_102","r49_103","r49_104","r49_105","r49_106","r49_107","r49_108","r49_109","r49_110","r49_111","r49_112","r49_113","r49_114","r49_115","r49_116","r49_117","r49_118","r49_119","r49_120","r49_121","r49_122","r49_123","r49_124","r49_125","r49_126","r49_127","r49_128","r49_129","r49_130","r49_131","r49_132","r49_133","r49_134","r49_135"};
    char FileNames[500][80]={"r24_1840","r24_1841","r24_1842","r24_1843","r24_1844","r24_1845","r24_1846","r24_1847","r24_1848","r24_1849","r24_1850","r24_1851","r24_1852","r24_1853","r24_1854","r24_1855","r24_1856","r24_1857","r24_1858","r24_1859","r24_1860","r24_1861","r24_1862","r24_1863","r24_1864","r24_1865","r24_1866","r24_1867","r24_1868","r24_1869","r24_1870","r24_1871","r24_1872","r24_1873","r24_1874","r24_1875","r24_1876","r24_1877","r24_1878","r24_1879","r24_1880","r24_1881","r24_1882","r24_1883","r24_1884","r24_1885","r24_1886","r24_1887","r24_1888","r24_1889","r24_1890","r24_1891","r24_1892","r24_1893","r24_1894","r24_1895","r24_1896","r24_1897","r24_1898","r24_1899","r24_1900","r24_1901","r24_1902","r24_1903","r24_1904","r24_1905","r24_1906","r24_1907","r24_1908","r24_1909","r24_1910","r24_1911","r24_1912","r24_1913","r24_1914","r24_1915","r24_1916","r24_1917","r24_1918","r24_1919","r24_1920","r24_1921","r24_1922","r24_1923","r24_1924","r24_1925","r24_1926","r24_1927","r24_1928","r24_1929","r24_1930","r24_1931","r24_1932","r24_1933","r24_1934","r24_1935","r24_1936","r24_1937","r24_1938","r24_1939","r24_1940","r24_1941","r24_1942","r24_1943","r24_1944","r24_1945","r24_1946","r24_1947","r24_1948","r24_1949","r24_1950","r24_1951","r24_1952","r24_1953","r24_1954","r24_1955","r24_1956","r24_1957","r24_1958","r24_1959","r24_1960","r24_1961","r24_1962","r24_1963","r24_1964","r24_1965","r24_1966","r24_1967","r24_1968","r24_1969","r24_1970","r24_1971","r24_1972","r24_1973","r24_1974","r24_1975","r24_1976","r24_1977","r24_1978","r24_1979","r24_1980","r24_1981","r24_1982","r24_1983","r24_1984","r24_1985","r24_1986","r24_1987","r24_1988","r24_1989","r24_1990","r24_1991","r24_1992","r24_1993","r24_1994","r24_1995","r24_1996","r24_1997","r24_1998","r24_1999","r24_2000","r24_2001","r24_2002","r24_2003","r24_2004","r24_2005","r24_2006","r24_2007","r24_2008","r24_2009","r24_2010","r24_2011","r24_2012","r24_2013","r24_2014","r24_2015","r24_2016","r24_2017","r24_2018","r24_2019","r24_2020","r24_2021","r24_2022","r24_2023","r24_2024","r24_2025","r24_2026","r24_2027","r24_2028","r24_2029","r24_2030","r24_2031","r24_2032","r24_2033","r24_2034","r24_2035","r24_2036","r24_2037","r24_2038","r24_2039","r24_2040","r24_2041","r24_2042","r24_2043","r24_2044","r24_2045","r24_2046","r24_2047","r24_2048","r24_2049","r24_2050","r24_2051","r24_2052","r24_2053","r24_2054","r24_2055","r24_2056","r24_2057","r24_2058","r24_2059","r24_2060","r24_2061","r24_2062","r24_2063","r24_2064","r24_2065","r24_2066","r24_2067","r24_2068","r24_2069","r24_2070","r24_2071","r24_2072","r24_2073","r24_2074","r24_2075","r24_2076","r24_2077","r24_2078","r24_2079","r24_2080","r24_2081","r24_2082","r24_2083","r24_2084","r24_2085","r24_2086","r24_2087","r24_2088","r24_2089","r24_2090","r24_2091","r24_2092","r24_2093","r24_2094","r24_2095","r24_2096","r24_2097","r24_2098","r24_2099","r24_2100","r24_2101","r24_2102","r24_2103","r24_2104","r24_2105","r24_2106","r24_2107","r24_2108","r24_2109","r24_2110","r24_2111","r24_2112","r24_2113","r24_2114","r24_2115","r24_2116","r24_2117","r24_2118","r24_2119","r24_2120","r24_2121","r24_2122","r24_2123","r24_2124","r24_2125","r24_2126","r24_2127","r24_2128","r24_2129","r24_2130","r24_2131","r24_2132","r24_2133","r24_2134","r24_2135","r24_2136","r24_2137","r24_2138","r24_2139","r24_2140","r24_2141","r24_2142","r24_2143","r24_2144","r24_2145","r24_2146","r24_2147","r24_2148","r24_2149","r24_2150","r24_2151","r24_2152","r24_2153","r24_2154","r24_2155","r24_2156","r24_2157","r24_2158","r24_2159","r24_2160","r24_2161","r24_2162","r24_2163","r24_2164","r24_2165","r24_2166","r24_2167","r24_2168","r24_2169","r24_2170","r24_2171","r24_2172","r24_2173","r24_2174","r24_2175","r24_2176","r24_2177","r24_2178","r24_2179","r24_2180","r24_2181","r24_2182","r24_2183","r24_2184","r24_2185","r24_2186","r24_2187","r24_2188","r24_2189","r24_2190","r24_2191","r24_2192","r24_2193","r24_2194","r24_2195","r24_2196","r24_2197","r24_2198","r24_2199","r24_2200","r24_2201","r24_2202","r24_2203","r24_2204","r24_2205","r24_2206","r24_2207","r24_2208","r24_2209","r24_2210","r24_2211","r24_2212","r24_2213","r24_2214","r24_2215","r24_2216","r24_2217","r24_2218","r24_2219","r24_2220","r24_2221","r24_2222","r24_2223","r24_2224","r24_2225","r24_2226","r24_2227","r24_2228","r24_2229","r24_2230","r24_2231","r24_2232","r24_2233","r24_2234","r24_2235","r24_2236","r24_2237","r24_2238","r24_2239","r24_2240","r24_2241","r24_2242","r24_2243","r24_2244","r24_2245","r24_2246","r24_2247","r24_2248","r24_2249","r24_2250","r24_2251","r24_2252","r24_2253","r24_2254","r24_2255","r24_2256","r24_2257","r24_2258","r24_2259","r24_2260","r24_2261","r24_2262","r24_2263","r24_2264"};
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