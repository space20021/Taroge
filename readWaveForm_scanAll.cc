#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>// for usleep(usec);
#include <iomanip> //setw
#include "datalib.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMultiGraph.h"
using namespace std;

union short_bytes {
    short shorty[1000];
    char bytes[2000];
};

//int abc(int argc, char* argv[])
int abc(char* path, char* filename)
{
  cl_data dataBlock;
  ifstream infile;
  bool readStart=true;
  char char_date[80] = "";
  char PdfnameBin[80] = "";
  char Pdfname[80] = "";
  char name[80] = "";
  char PdfnameEnd[80] = "";
  const int length = 1000;
  double ch1[length], t[length];
  short_bytes buffer;

  TCanvas* c =new TCanvas("WaveForm","data",1600,1200);
  c->Divide(4,3);
  TGraph* gr1[12];

  sprintf(char_date,"%s%s.bin", path, filename );

  sprintf(PdfnameBin,"Run%s.pdf(", filename );
  sprintf(Pdfname,"Run%s.pdf", filename );
  sprintf(PdfnameEnd,"Run%s.pdf)", filename ); //cout<<"PDFnameEnd="<<PdfnameEnd<<endl;

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
  //while(eventCounter<300)
  {
	cout <<" writeing " <<eventCounter <<"th events"<<endl; eventCounter++;
	for(int j=0; j<12; j++) // read channels data for biniary file , ch 0-11 (12 channels);
	{
		c->cd(12-j); //(go to j TCanvas)

		infile.read ( (char *)(&dataBlock), sizeof(dataBlock) ); // read jth channel data

//		cout <<"ANt #"<< dataBlock.ant_N <<endl;
		cout << "globaltime =  "<<  setprecision(14) << dataBlock.globaltime <<endl;
		cout <<"Format"<< dataBlock.format <<endl;
//		cout <<"mem Length"<< dataBlock.memLength<<endl;
//		cout << dataBlock.IntpDist <<endl;
//		cout << dataBlock.trig_Add <<endl;
//		cout << dataBlock.trig_level <<endl;
//		cout <<" dataBlock.v_Unit, "  <<dataBlock.v_Unit <<endl;
//		cout <<" dataBlock.v_Unit_div, " << dataBlock.v_Unit_div <<endl;
//		cout <<" dataBlock.v_Unit_ext"<< dataBlock.v_Unit_ext <<endl;
//		cout <<" dataBlock.h_Position"<< dataBlock.h_Position <<endl;
//		cout << "v_scale" <<dataBlock.v_Scale <<endl;
//		cout<<j<<endl;
		//cout <<"time stamp = "<< dataBlock.timeStamp[0]<<endl;
		//cout <<"time stamp = "<< dataBlock.timeStamp[1]<<endl;
		//cout <<"time stamp = "<< dataBlock.timeStamp[2]<<endl<<endl;
		for(int i = 0; i<19; i++) {//cout << dataBlock.date[i];
	} //cout<<endl;

		for(int k=0; k<length; k++)
		{
			
			buffer.bytes[2*k+1] = dataBlock.WaveForm[2*k];//===========================WAVEFORM
			buffer.bytes[2*k] = dataBlock.WaveForm[2*k+1];//===========================WAVEFORM
			ch1[k] = (float)buffer.shorty[k]*0.04*dataBlock.v_Scale;
			t[k] = (float)k *dataBlock.h_scale*1E-2 +dataBlock.h_Position; 
		}

		gr1[j] = new TGraph(length, t, ch1); //plot WaveForm figure to j-th pannels

 		if(j<4)// give labels
		{
                	sprintf(name,"GEN222043, Tier 3, ANT%d", j+1 );//sprintf(name,"GEN141421 ANT%d", j );
		}
		else if(j<8)
		{
                	sprintf(name,"GEN170272, Tier 2, ANT%d", j-3 );
		}
		else
		{
                	sprintf(name,"GEN141421, Tier 1, ANT%d", j-7 );//sprintf(name,"GEN222043 ANT%d", j );
		}

		gr1[j]->SetTitle(name);		
		gr1[j]->SetMaximum(3);
		gr1[j]->SetMinimum(-3);
		gr1[j]->Draw("AL");
		//gr1[j]->GetXaxis()->SetLimits(-1E-6,1E-6 ); //Set the range of X axis
	}

	c->Update();

	if(readStart == true)
	{ 
		c->Print(PdfnameBin,"pdf");
		readStart == false;
	}
	else 
	{c->Print(Pdfname,"pdf");}

 }
  c->Print(PdfnameEnd,"pdf");
  infile.close();
}
