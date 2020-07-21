#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <TGraphErrors.h>
using namespace std;
void comp_hms_pointing(Int_t fid = -5) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_hms_pointing.pdf";
    //
   ifstream HMSSurveyFile("survey_hms_data.dat");
  vector <Double_t> surveyhms_angles;
  vector <Double_t> surveyhms_ymis;
  vector <Double_t> surveyhms_xmis;
    if (HMSSurveyFile.fail()) exit(1);
    while( !HMSSurveyFile.eof()) {
      Double_t tymis;
      Double_t txmis;
      Double_t tang;
      HMSSurveyFile >> tang >> tymis >> txmis ;
      surveyhms_angles.push_back(tang);
      surveyhms_ymis.push_back(tymis);
      surveyhms_xmis.push_back(txmis);
    }
    	TGraph *gSurveyYmis_Angle= new TGraph(surveyhms_angles.size(), &surveyhms_angles[0], &surveyhms_ymis[0]);
	gSurveyYmis_Angle->SetMarkerStyle(24);
	gSurveyYmis_Angle->SetMarkerColor(6);
    //
   ifstream OldHMSPointingFile("hms_eric_mispointing.dat");
  vector <Double_t> oldhms_angles;
  vector <Double_t> oldhms_ymis;
    if (OldHMSPointingFile.fail()) exit(1);
    while( !OldHMSPointingFile.eof()) {
      Double_t tmis;
      Double_t tang;
      OldHMSPointingFile >> tang >> tmis ;
      oldhms_angles.push_back(tang);
      oldhms_ymis.push_back(tmis);
    }
    	TGraph *gOldYmis_Angle= new TGraph(oldhms_angles.size(), &oldhms_angles[0], &oldhms_ymis[0]);
	gOldYmis_Angle->SetMarkerStyle(23);
	gOldYmis_Angle->SetLineColor(2);
  //  Get info pointing scan
 TString PointingFileName = "list_of_pointing_run.dat";
   ifstream file_Pointing(PointingFileName.Data());
  vector <Int_t> nrun;
  vector <Double_t> angles;
  vector <TString > OpticsID;
   cout << " Open file = " << PointingFileName << endl;
    if (file_Pointing.fail()) exit(1);
    while( !file_Pointing.eof()) {
      Int_t tn;
      Double_t tang;
      TString tOpt;
      file_Pointing >> tn >> tang >>tOpt ;
      nrun.push_back(tn);
      angles.push_back(tang);
      OpticsID.push_back(tOpt);
    }
    //
    vector <Double_t> x_mis(nrun.size());    
    vector <Double_t> y_mis(nrun.size());    
    for (Int_t n=0;n<nrun.size();n++) {
      y_mis[n] = 0.1*(0.52-0.012*angles[n]+0.002*angles[n]*angles[n]); // cm
      x_mis[n] = 0.1*(2.37-0.086*angles[n]+0.0012*angles[n]*angles[n]); //cm
      //      cout << nrun[n] << " " << angles[n] << " "<< y_mis[n] << " "<< x_mis[n] << endl;
    }
    //
    //
    TF1*   Ymis_fit = new TF1("Ymis_fit","(0.52-0.012*x+0.002*x*x)",12.,20.);

    vector<TH1F* > Hxbeam(nrun.size());
    vector<TH1F* > Hybeam(nrun.size());
    vector<TH1F* > Hytar(nrun.size());
    vector <Double_t> xbeam(nrun.size());    
    vector <Double_t> ybeam(nrun.size());    
     //
    for (Int_t n=0;n<nrun.size();n++) {
      TString HistFileName = Form("hist/hms_replay_matrixopt_%s_%d_pointing_hist.root",OpticsID[n].Data(),fid);
      TFile * HistRoot =  new TFile(HistFileName);
      Hxbeam[n] = (TH1F*)HistRoot->Get("hxbeam");
      Hybeam[n] = (TH1F*)HistRoot->Get("hybeam");
      Hytar[n] = (TH1F*)HistRoot->Get("hytar");
      xbeam[n] = Hxbeam[n]->GetMean();
       ybeam[n] = Hybeam[n]->GetMean();
   }    
    //
	TCanvas* canytar;
	TF1* fytar[nrun.size()];
	Double_t angles_err[nrun.size()];
	Double_t ytar_cent[nrun.size()];
	Double_t ytarymis_cent[nrun.size()];
	Double_t ytar_centerr[nrun.size()];
	Double_t ytar_sigma[nrun.size()];
	Double_t ymis_calc[nrun.size()];
	Double_t ymis_calc2[nrun.size()];
	Double_t ymis_calcerr[nrun.size()];
	Double_t ymis_err[nrun.size()];
	canytar = new TCanvas(Form("Canytar"),"ytar", 700,700);
	canytar->Divide(2,3);
	for  (Int_t nc=0;nc<nrun.size();nc++) {
	  Hytar[nc]->GetXaxis()->SetRangeUser(-1,1);
	  fytar[nc] = new TF1(Form("fytar_%d",nc),"gaus");
	      fytar[nc]->SetParameter(1,Hytar[nc]->GetBinCenter(Hytar[nc]->GetMaximumBin()));
	      fytar[nc]->SetParameter(0,Hytar[nc]->GetMaximum());
	      fytar[nc]->SetParameter(2,Hytar[nc]->GetRMS());
	  canytar->cd(nc+1);
	  Hytar[nc]->GetXaxis()->SetRangeUser(-3,3);
	  Hytar[nc]->Draw();
	  Hytar[nc]->Fit(Form("fytar_%d",nc),"Q","",-1,1);
	  Double_t min=fytar[nc]->GetParameter(1)-3*fytar[nc]->GetParameter(2);
	  Double_t max=fytar[nc]->GetParameter(1)+3*fytar[nc]->GetParameter(2);
	  Hytar[nc]->Fit(Form("fytar_%d",nc),"Q","",min,max);
	  angles_err[nc] = 0.0001;
	  ytarymis_cent[nc] = 10*fytar[nc]->GetParameter(1)+10*y_mis[nc];
	  ytar_cent[nc] = 10*fytar[nc]->GetParameter(1);
	  ytar_centerr[nc] = 10*fytar[nc]->GetParError(1);
	  ymis_calc[nc] =10*(-xbeam[nc]*cos(angles[nc]/180.*3.14159) -fytar[nc]->GetParameter(1));
	  ymis_calc2[nc] =10*(-xbeam[nc]*cos(angles[nc]/180.*3.14159) -fytar[nc]->GetParameter(1));
	  ymis_calcerr[nc] = 10*fytar[nc]->GetParError(1);
	  y_mis[nc]=10.*y_mis[nc];
	  ymis_err[nc] = 0.0001;
	  ytar_sigma[nc] = 10*fytar[nc]->GetParameter(2);
	  Double_t YPred = 10*(xbeam[nc]*cos(angles[nc]/180.*3.14159) -y_mis[nc]); 
	  cout << nrun[nc] << " " << angles[nc] << " "<< y_mis[nc] << " "<< ymis_calc[nc] << " "<< xbeam[nc] << " " << ytar_cent[nc]<< " " << ytar_centerr[nc] << " " << ytar_sigma[nc]<< endl;
	}
	//
	TF1 *fYtarYmis = new TF1("fYtarYmis","-2.6*cos(x/180.*3.14159) +  [0]*sin(x/180.*3.14159)");
	TF1 *fYtar = new TF1("fYtar","[0]*sin(x/180.*3.14159) -2.6*cos(x/180.*3.14159)+0.9");
    	TGraphErrors *gYmis_Angle= new TGraphErrors(nrun.size(), &angles[0], &y_mis[0], &angles_err[0], &ymis_err[0]);
	gYmis_Angle->SetMarkerStyle(21);
	gYmis_Angle->SetLineColor(1);
	gYmis_Angle->SetTitle(" HMS Pointing scan" );
	gYmis_Angle->GetXaxis()->SetTitle("HMS Angle (deg)");
	gYmis_Angle->GetYaxis()->SetTitle("YMis  (mm)");
	gYmis_Angle->GetYaxis()->SetRangeUser(-3.,3.);
    	TGraphErrors *gYmisCalc2_Angle= new TGraphErrors(nrun.size(), &angles[0], &ymis_calc2[0], &angles_err[0], &ymis_calcerr[0]);
	gYmisCalc2_Angle->SetMarkerStyle(23);
	gYmisCalc2_Angle->SetMarkerColor(3);
	gYmisCalc2_Angle->SetTitle(" HMS Pointing scan" );
	gYmisCalc2_Angle->GetXaxis()->SetTitle("HMS Angle (deg)");
	gYmisCalc2_Angle->GetYaxis()->SetTitle("YMis Calc (mm)");
	gYmisCalc2_Angle->GetYaxis()->SetRangeUser(-5.,5.);
    	TGraphErrors *gYmisCalc_Angle= new TGraphErrors(nrun.size(), &angles[0], &ymis_calc[0], &angles_err[0], &ymis_calcerr[0]);
	gYmisCalc_Angle->SetMarkerStyle(22);
	gYmisCalc_Angle->SetMarkerColor(2);
	gYmisCalc_Angle->SetTitle(" HMS Pointing scan" );
	gYmisCalc_Angle->GetXaxis()->SetTitle("HMS Angle (deg)");
	gYmisCalc_Angle->GetYaxis()->SetTitle("YMis Calc (mm)");
    	TGraphErrors *gYtar_Angle= new TGraphErrors(nrun.size(), &angles[0], &ytar_cent[0], &angles_err[0], &ytar_centerr[0]);
	gYtar_Angle->SetMarkerStyle(22);
	gYtar_Angle->SetTitle(" HMS Pointing scan" );
	gYtar_Angle->GetXaxis()->SetTitle("HMS Angle (deg)");
	gYtar_Angle->GetYaxis()->SetTitle("Ytar (mm)");
    	TGraphErrors *gYtarYmis_Angle= new TGraphErrors(nrun.size(), &angles[0], &ytarymis_cent[0], &angles_err[0], &ytar_centerr[0]);
	gYtarYmis_Angle->SetMarkerStyle(20);
	gYtarYmis_Angle->SetTitle(" HMS Pointing scan" );
	gYtarYmis_Angle->GetXaxis()->SetTitle("HMS Angle (deg)");
	gYtarYmis_Angle->GetYaxis()->SetTitle("Ytar + Y_mis (mm)");
	TCanvas* canytarfit;
	canytarfit = new TCanvas(Form("Canytarfit"),"ytarfit", 700,700);
	canytarfit->Divide(2,2);
	canytarfit->cd(1);
	gYtarYmis_Angle->Draw("AP");
	gYtarYmis_Angle->Fit("fYtarYmis");
	canytarfit->cd(2);
	gYtar_Angle->Draw("AP");
	gYtar_Angle->Fit("fYtar");
	canytarfit->cd(3);
	//gYmis_Angle->Draw("AP");
	//gYmisCalc_Angle->Draw("AP");
	gYmisCalc2_Angle->Draw("AP");
	gOldYmis_Angle->Draw("L same");
	gSurveyYmis_Angle->Draw("P same");
	Ymis_fit->Draw("L same");
	Ymis_fit->SetLineColor(3);
	//	gYtar_Angle->Fit("fYtar2","","",5.,20.);
	//gYtar_Angle->Fit("fYtar3","","",20.,40.);
    //
}
