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
    Double_t ch1_Buf[length], ch1_fftBuf[length], max_ele;
    Bool_t max_115, max_143, max_156;
    
    //char FileNames[1000][80]={"r59_14","r59_15","r59_16","r59_17","r59_18","r59_19","r59_20","r59_21","r59_22","r59_23","r59_24","r59_25","r59_26","r59_27","r59_28","r59_29","r59_30","r59_31","r59_32","r59_33","r59_34","r59_35","r59_36","r59_37","r59_38","r59_39","r59_40","r59_41","r59_42","r59_43","r59_44","r59_45","r59_46","r59_47","r59_48","r59_49","r59_50","r59_51","r59_52","r59_53","r59_54","r59_55","r59_56","r59_57","r59_58","r59_59","r59_60","r59_61","r59_62","r59_63","r59_64","r59_65","r59_66","r59_67","r59_68","r59_69","r59_70","r59_71","r59_72","r59_73","r59_74","r59_75","r59_76","r59_77","r59_78","r59_79","r59_80","r59_81","r59_82","r59_83","r59_84","r59_85","r59_86","r59_87","r21_801","r21_802","r21_803","r21_804","r21_805","r21_806","r21_807","r21_808","r21_809","r21_810","r21_811","r21_812","r21_813","r21_814","r21_815","r21_816","r21_817","r21_818","r21_819","r21_820","r21_821","r21_822","r21_823","r21_824","r21_825","r21_826","r21_827","r21_828","r21_829","r21_830","r21_831","r21_832","r21_833","r21_834","r21_835","r21_836","r21_837","r21_838","r21_839","r21_840","r21_841","r21_842","r21_843","r21_844","r21_845","r21_846","r21_847","r21_848","r21_849","r21_850","r21_851","r21_852","r21_853","r21_854","r21_855","r21_856","r21_857","r21_858","r21_859","r21_860","r21_861","r21_862","r21_863","r21_864","r21_865","r21_866","r21_867","r21_868","r21_869","r21_870","r21_871","r21_872","r21_873","r21_874","r21_875","r21_876","r21_877","r21_878","r21_879","r21_880","r21_881","r21_882","r21_883","r21_884","r21_885","r21_886","r21_887","r21_888","r21_889","r21_890","r21_891","r21_892","r21_893","r21_894","r21_895","r21_896","r21_897","r21_898","r21_899","r21_900","r21_901","r21_902","r21_903","r21_904","r21_905","r21_906","r21_907","r21_908","r21_909","r21_910","r21_911","r21_912","r21_913","r21_914","r21_915","r21_916","r21_917","r21_918","r21_919","r21_920","r21_921","r21_922","r21_923","r21_924","r21_925","r21_926","r21_927","r21_928","r21_929","r21_930","r21_931","r21_932","r21_933","r21_934","r21_935","r21_936","r21_937","r21_938","r21_939","r21_940","r21_941","r21_942","r21_943","r21_944","r21_945","r21_946","r21_947","r21_948","r21_949","r21_950","r21_951","r21_952","r21_953","r21_954","r21_955","r21_956","r21_957","r21_958","r21_959","r21_960","r21_961","r21_962","r21_963","r21_964","r21_965","r21_966","r21_967","r21_968","r21_969","r21_970","r21_971","r21_972","r21_973","r21_974","r21_975","r21_976","r21_977","r21_978","r21_979","r21_980","r21_981","r21_982","r21_983","r21_984","r21_985","r21_986","r21_987","r21_988","r21_989","r21_990","r21_991","r21_992","r21_993","r21_994","r24_1840","r24_1841","r24_1842","r24_1843","r24_1844","r24_1845","r24_1846","r24_1847","r24_1848","r24_1849","r24_1850","r24_1851","r24_1852","r24_1853","r24_1854","r24_1855","r24_1856","r24_1857","r24_1858","r24_1859","r24_1860","r24_1861","r24_1862","r24_1863","r24_1864","r24_1865","r24_1866","r24_1867","r24_1868","r24_1869","r24_1870","r24_1871","r24_1872","r24_1873","r24_1874","r24_1875","r24_1876","r24_1877","r24_1878","r24_1879","r24_1880","r24_1881","r24_1882","r24_1883","r24_1884","r24_1885","r24_1886","r24_1887","r24_1888","r24_1889","r24_1890","r24_1891","r24_1892","r24_1893","r24_1894","r24_1895","r24_1896","r24_1897","r24_1898","r24_1899","r24_1900","r24_1901","r24_1902","r24_1903","r24_1904","r24_1905","r24_1906","r24_1907","r24_1908","r24_1909","r24_1910","r24_1911","r24_1912","r24_1913","r24_1914","r24_1915","r24_1916","r24_1917","r24_1918","r24_1919","r24_1920","r24_1921","r24_1922","r24_1923","r24_1924","r24_1925","r24_1926","r24_1927","r24_1928","r24_1929","r24_1930","r24_1931","r24_1932","r24_1933","r24_1934","r24_1935","r24_1936","r24_1937","r24_1938","r24_1939","r24_1940","r24_1941","r24_1942","r24_1943","r24_1944","r24_1945","r24_1946","r24_1947","r24_1948","r24_1949","r24_1950","r24_1951","r24_1952","r24_1953","r24_1954","r24_1955","r24_1956","r24_1957","r24_1958","r24_1959","r24_1960","r24_1961","r24_1962","r24_1963","r24_1964","r24_1965","r24_1966","r24_1967","r24_1968","r24_1969","r24_1970","r24_1971","r24_1972","r24_1973","r24_1974","r24_1975","r24_1976","r24_1977","r24_1978","r24_1979","r24_1980","r24_1981","r24_1982","r24_1983","r24_1984","r24_1985","r24_1986","r24_1987","r24_1988","r24_1989","r24_1990","r24_1991","r24_1992","r24_1993","r24_1994","r24_1995","r24_1996","r24_1997","r24_1998","r24_1999","r24_2000","r24_2001","r24_2002","r24_2003","r24_2004","r24_2005","r24_2006","r24_2007","r24_2008","r24_2009","r24_2010","r24_2011","r24_2012","r24_2013","r24_2014","r24_2015","r24_2016","r24_2017","r24_2018","r24_2019","r24_2020","r24_2021","r24_2022","r24_2023","r24_2024","r24_2025","r24_2026","r24_2027","r24_2028","r24_2029","r24_2030","r24_2031","r24_2032","r24_2033","r24_2034","r24_2035","r24_2036","r24_2037","r24_2038","r24_2039","r24_2040","r24_2041","r24_2042","r24_2043","r24_2044","r24_2045","r24_2046","r24_2047","r24_2048","r24_2049","r24_2050","r24_2051","r24_2052","r24_2053","r24_2054","r24_2055","r24_2056","r24_2057","r24_2058","r24_2059","r24_2060","r24_2061","r24_2062","r24_2063","r24_2064","r24_2065","r24_2066","r24_2067","r24_2068","r24_2069","r24_2070","r24_2071","r24_2072","r24_2073","r24_2074","r24_2075","r24_2076","r24_2077","r24_2078","r24_2079","r24_2080","r24_2081","r24_2082","r24_2083","r24_2084","r24_2085","r24_2086","r24_2087","r24_2088","r24_2089","r24_2090","r24_2091","r24_2092","r24_2093","r24_2094","r24_2095","r24_2096","r24_2097","r24_2098","r24_2099","r24_2100","r24_2101","r24_2102","r24_2103","r24_2104","r24_2105","r24_2106","r24_2107","r24_2108","r24_2109","r24_2110","r24_2111","r24_2112","r24_2113","r24_2114","r24_2115","r24_2116","r24_2117","r24_2118","r24_2119","r24_2120","r24_2121","r24_2122","r24_2123","r24_2124","r24_2125","r24_2126","r24_2127","r24_2128","r24_2129","r24_2130","r24_2131","r24_2132","r24_2133","r24_2134","r24_2135","r24_2136","r24_2137","r24_2138","r24_2139","r24_2140","r24_2141","r24_2142","r24_2143","r24_2144","r24_2145","r24_2146","r24_2147","r24_2148","r24_2149","r24_2150","r24_2151","r24_2152","r24_2153","r24_2154","r24_2155","r24_2156","r24_2157","r24_2158","r24_2159","r24_2160","r24_2161","r24_2162","r24_2163","r24_2164","r24_2165","r24_2166","r24_2167","r24_2168","r24_2169","r24_2170","r24_2171","r24_2172","r24_2173","r24_2174","r24_2175","r24_2176","r24_2177","r24_2178","r24_2179","r24_2180","r24_2181","r24_2182","r24_2183","r24_2184","r24_2185","r24_2186","r24_2187","r24_2188","r24_2189","r24_2190","r24_2191","r24_2192","r24_2193","r24_2194","r24_2195","r24_2196","r24_2197","r24_2198","r24_2199","r24_2200","r24_2201","r24_2202","r24_2203","r24_2204","r24_2205","r24_2206","r24_2207","r24_2208","r24_2209","r24_2210","r24_2211","r24_2212","r24_2213","r24_2214","r24_2215","r24_2216","r24_2217","r24_2218","r24_2219","r24_2220","r24_2221","r24_2222","r24_2223","r24_2224","r24_2225","r24_2226","r24_2227","r24_2228","r24_2229","r24_2230","r24_2231","r24_2232","r24_2233","r24_2234","r24_2235","r24_2236","r24_2237","r24_2238","r24_2239","r24_2240","r24_2241","r24_2242","r24_2243","r24_2244","r24_2245","r24_2246","r24_2247","r24_2248","r24_2249","r24_2250","r24_2251","r24_2252","r24_2253","r24_2254","r24_2255","r24_2256","r24_2257","r24_2258","r24_2259","r24_2260","r24_2261","r24_2262","r24_2263","r24_2264"};
    char FileNames[200][80]={"r24_2264"};
    
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
    t1->Branch("max_115[1]",&max_115,"max_115/O");
    t1->Branch("max_143[1]",&max_143,"max_143/O");
    t1->Branch("max_156[1]",&max_156,"max_156/O");
    
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
                ch1_fftMAG[j][k]=sqrt((pow(ch1_fftRe[j][k],2)+pow(ch1_fftIm[j][k],2)))/length*2;
                ch1_fftPow[j][k]=pow(ch1_fftMAG[j][k],2)/50;
                ch1_fftdBm[j][k]=10*TMath::Log10(ch1_fftPow[j][k])+30;
                ch1_fftBuf[k]=ch1_fftdBm[j][k];
            }
            
            max_143=0;
            max_156=0;
            
            //Check whether 115 MHz is the largest component
            max_ele=ch1_fftBuf[113];
            for (int k=1;k<=4;k++)
            {
                if (ch1_fftBuf[113+k]>max_ele)
                    max_ele=ch1_fftBuf[113+k];
            }
            max_115=1;
            for (int k=0;k<length/2-1;k++)
            {
                if (max_ele<ch1_fftBuf[k])
                {
                    max_115=0;
                    break;
                }
            }
            
            //Check whether 143 MHz is the largest component
            if (max_115==0)
            {
                max_ele=ch1_fftBuf[141];
                for (int k=1;k<=4;k++)
                {
                    if (ch1_fftBuf[141+k]>max_ele)
                        max_ele=ch1_fftBuf[141+k];
                }
                max_143=1;
                for (int k=0;k<length/2-1;k++)
                {
                    if (max_ele<ch1_fftBuf[k])
                    {
                        max_143=0;
                        break;
                    }
                }
            }
            
            //Check whether 156 MHz is the largest component
            if (max_115==0&&max_143==0)
            {
                max_ele=ch1_fft[154];
                for (int k=1;k<=4;k++)
                {
                    if (ch1_fft[154+k]>max_ele)
                        max_ele=ch1_fft[154+k];
                }
                max_156=1;
                for (int k=0;k<length/2-1;k++)
                {
                    if (max_ele<ch1_fft[k])
                    {
                        max_156=0;
                        break;
                    }
                }
            }
            
            t1->Fill();
            delete fft_own;
        }
    }
    
    cout << "+-----------------------------+" << endl;
}