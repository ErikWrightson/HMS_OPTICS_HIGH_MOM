#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCutG.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;


void comp_yfp_xfp(Int_t FileID=-2) {
  
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.14);
  const Int_t Number=3;
  Double_t Red[Number] = { 1.0,0.0,0.0};
 Double_t Blue[Number] = { 1.0,0.0,1.0};
  Double_t Green[Number] = { 0.0,1.0,0.0};
 Double_t Len[Number] = { 0.0,.5,1.0};
 Int_t nb=50;
 TColor::CreateGradientColorTable(Number,Len,Red,Green,Blue,nb);
  //  Get info for that optics run
 TString OpticsFile = "list_of_optics_run.dat";
   ifstream file_optics(OpticsFile.Data());
 TString opticsline;
 vector<TString > OpticsID;
 vector<Double_t > Angle;
 vector<Double_t > ymis;
 vector<Double_t > xbeam;
 vector<vector<Double_t > > ztar;
 vector<vector<Double_t > > ytarCalc;
  TString temp;
 //
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (!file_optics.eof()  ) {
      temp.ReadToDelim(file_optics,',');
      temp.ReadToDelim(file_optics,',');
      //     cout << temp << endl;
      OpticsID.push_back(temp);
      temp.ReadToDelim(file_optics,',');
      //    cout << temp << endl;
      Angle.push_back(temp.Atof());
      temp.ReadToDelim(file_optics,',');
      //  cout << temp << endl;
      Int_t nftot = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      temp.ReadToDelim(file_optics,',');
      temp.ReadToDelim(file_optics);
      ymis.push_back(temp.Atof());
      vector<Double_t > ftemp;
     for (Int_t nf=0;nf<nftot;nf++) {
        temp.ReadToDelim(file_optics,',');
        cout << temp << endl;
 	ftemp.push_back(temp.Atof());
      }
      ztar.push_back(ftemp);
      file_optics >> temp;
      //      cout << temp << endl;
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }
  //
  //
  Int_t nftot=OpticsID.size();
  cout << " n files = " << nftot << endl;
  ytarCalc.resize(nftot);
  xbeam.resize(nftot);
  TString inputroot;
   TFile *fhistroot;
   TH2F* hYFpXFp[nftot];
   TH1F* hxbpm_tar[nftot];
   TH1F* hYtar[nftot];
   TH1F* hZtar[nftot];
   TH2F* hYtarDelta[nftot];
   TH2F* hZtarDelta[nftot];
   TString outputpdf;
   for (Int_t nf=0;nf<nftot;nf++) {
   inputroot=Form("hist/Optics_%s_%d_hist.root",OpticsID[nf].Data(),FileID);
   outputpdf = Form("plots/comp_shms_xfp_yfp_%d",FileID);
     cout << " infile root = " << inputroot << endl;
     fhistroot =  new TFile(inputroot);
     TString title = Form("%s Angle = %4.1f Fit = %d",OpticsID[nf].Data(),Angle[nf],FileID);
     hYFpXFp[nf] = (TH2F*)fhistroot->Get("hYFpXFp_all");
     hYFpXFp[nf]->SetTitle(title);
     hYtar[nf] = (TH1F*)fhistroot->Get("hytar");
     hYtar[nf]->SetTitle(title);
     hxbpm_tar[nf] = (TH1F*)fhistroot->Get("hxbpm_tar");
     hxbpm_tar[nf]->SetTitle(title);
     if (hxbpm_tar[nf]) xbeam[nf] =-hxbpm_tar[nf]->GetMean(); 
     hZtar[nf] = (TH1F*)fhistroot->Get("hztar");
     hZtar[nf]->SetTitle(title);
     hYtarDelta[nf] = (TH2F*)fhistroot->Get("hYtarDelta");
     hYtarDelta[nf]->SetTitle(title);
     hZtarDelta[nf] = (TH2F*)fhistroot->Get("hZtarDelta");
     hZtarDelta[nf]->SetTitle(title);
   }
 //
  for (Int_t nf=0;nf<nftot;nf++) {
     Int_t nfoils = ztar[nf].size();
     ytarCalc[nf].resize(nfoils);
   }  
   for (Int_t nf=0;nf<nftot;nf++) {
   for (Int_t nc=0;nc<ztar[nf].size();nc++) {
     //       ytarCalc[nf][nc]= +ztar[nf][nc]*sin(Angle[nf]/57.3)-0.23*cos(Angle[nf]/57.3);   
     ytarCalc[nf][nc]= ztar[nf][nc]*sin(Angle[nf]/57.3)+xbeam[nf]*cos(Angle[nf]/57.3)- ymis[nf];   
   }}
 //
   TCanvas* cFp = new TCanvas("cFp","cFp",700,700);
   cFp->Divide(2,2);
   for (Int_t nf=0;nf<nftot;nf++) {
     cFp->cd(nf+1);
 	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz();   
	if (hYFpXFp[nf]) {
    hYFpXFp[nf]->Draw("colz");
     hYFpXFp[nf]->SetMinimum(10);
	}
   }  
   cFp->Print(outputpdf+".pdf("); 
 //
     //
   TText* ytext;
   TLine* yline;
   TCanvas* cZtar = new TCanvas("cZtar","cZtar",700,700);
   cZtar->Divide(2,3);
   Int_t cind=1;
   for (Int_t nf=0;nf<nftot;nf++) {
     cZtar->cd(cind++);
     hZtarDelta[nf]->Draw("colz");
     hZtarDelta[nf]->SetMinimum(10);
   for (Int_t nc=0;nc<ztar[nf].size();nc++) {
     yline = new TLine(ztar[nf][nc],-10,ztar[nf][nc],10);
     yline->Draw("same");
   }
     cZtar->cd(cind++);    
     hZtar[nf]->Draw();
   for (Int_t nc=0;nc<ztar[nf].size();nc++) {
     yline = new TLine(ztar[nf][nc],0,ztar[nf][nc],hZtar[nf]->GetMaximum());
     yline->Draw("same");
   }
   }   
   cZtar->Print(outputpdf+".pdf"); 
     //
   TCanvas* cYtar = new TCanvas("cYtar","cYtar",700,700);
   cYtar->Divide(2,3);
   cind=1;
   for (Int_t nf=0;nf<nftot;nf++) {
     cYtar->cd(cind++);
     hYtarDelta[nf]->Draw("colz");
     hYtarDelta[nf]->SetMinimum(10);
       Double_t doff=1;
   for (Int_t nc=0;nc<ztar[nf].size();nc++) {
     yline = new TLine(ytarCalc[nf][nc],-10,ytarCalc[nf][nc],10);
       if(doff >0) {
	 ytext = new TText(ytarCalc[nf][nc],-12.5,Form("%3.1f",ztar[nf][nc]));
       } else {
	 ytext = new TText(ytarCalc[nf][nc],-12.5,Form("%3.1f",ztar[nf][nc]));
       }
       ytext->SetTextAlign(22);
     ytext->Draw("same");
     yline->Draw("same");
       doff=pow(-1,(nc+1));
   }
     cYtar->cd(cind++);    
     hYtar[nf]->Draw();
   for (Int_t nc=0;nc<ztar[nf].size();nc++) {
     yline = new TLine(ytarCalc[nf][nc],0,ytarCalc[nf][nc],hYtar[nf]->GetMaximum());
     yline->Draw("same");
     ytext = new TText(ytarCalc[nf][nc],hYtar[nf]->GetMaximum(),Form("%3.1f",ztar[nf][nc]));
     ytext->Draw("same");
   }
   }   
   cYtar->Print(outputpdf+".pdf)"); 
   //
}
