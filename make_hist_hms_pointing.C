#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;
void make_hist_hms_pointing(TString tnrun,Int_t nrun,Int_t ID=-2) {
  TString basename=Form("hms_replay_matrixopt_%s_%d",tnrun.Data(),ID);
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist="hist/"+basename+"_pointing_hist.root";
 TObjArray HList(0);
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 Double_t  beamx;
  tsimc->SetBranchAddress("H.rb.raster.fr_xbpm_tar",&beamx);
 Double_t  beamy;
   tsimc->SetBranchAddress("H.rb.raster.fr_ybpm_tar",&beamy);
 Double_t  rastx;
   tsimc->SetBranchAddress("H.rb.raster.fr_xa",&rastx);
 Double_t  rasty;
   tsimc->SetBranchAddress("H.rb.raster.fr_ya",&rasty);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("H.gtr.y",&e_ytar);
   Double_t  npeSum;
   tsimc->SetBranchAddress("H.ngcer.npeSum",&npeSum);
 Double_t  ysieve;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xsieve);
 Double_t  e_xtar;
   tsimc->SetBranchAddress("H.gtr.x",&e_xtar);
 Double_t  gindex;
   tsimc->SetBranchAddress("H.gtr.index",&gindex);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etracknorm);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("H.react.z",&e_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("H.gtr.dp",&e_delta);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("H.gtr.th",&e_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&e_xpfp);
 Double_t  W;
   tsimc->SetBranchAddress("H.kin.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("H.kin.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("H.kin.scat_ang_rad",&ThScat);
   //
     TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etot ;Counts",nrun), 300, -.5,2.0);
    HList.Add(hetot);
     TH1F *hxrast = new TH1F("hxrast",Form("Run %d ; Xrast (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hxrast);
     TH1F *hysieve = new TH1F("hysieve",Form("Run %d ; Ysieve (cm);Counts",nrun), 300, -15,15.);
    HList.Add(hysieve);
     TH1F *hxsieve = new TH1F("hxsieve",Form("Run %d ; Xsieve (cm);Counts",nrun), 300, -15,15.);
    HList.Add(hxsieve);
    TH1F *hysieve_cent = new TH1F("hysieve_cent",Form("Ztar = 0 Run %d ; Ysieve (cm);Counts",nrun), 300, -15,15.);
    HList.Add(hysieve_cent);
    TH2F *hxsieve_ytar_yscent = new TH2F("hxsieve_ytar_yscent",Form("Ztar = 0 Run %d ; Xsieve (cm); Ytar",nrun), 100, -15,15.,100,-2,2);
    HList.Add(hxsieve_ytar_yscent);
     TH1F *hxsieve_cent = new TH1F("hxsieve_cent",Form("Ztar = 0 Run %d ; Xsieve (cm);Counts",nrun), 300, -15,15.);
    HList.Add(hxsieve);
   TH1F *hyrast = new TH1F("hyrast",Form("Run %d ; Yrast (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hyrast);
    TH1F *hxbeam = new TH1F("hxbeam",Form("Run %d ; Xbeam (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hxbeam);
    TH1F *hybeam = new TH1F("hybeam",Form("Run %d ; Ybeam (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hybeam);
    TH1F *hxtar = new TH1F("hxtar",Form("Run %d ; Xtar (cm);Counts",nrun), 300, -.3,.3);
    HList.Add(hxtar);
    TH1F *hytar = new TH1F("hytar",Form("Run %d ; Ytar (cm);Counts",nrun), 150, -3.,3.);
    HList.Add(hytar);
    TH1F *hztar = new TH1F("hztar",Form("Run %d ; Ztar (cm);Counts",nrun), 150, -15.,15.);
    HList.Add(hztar);
    //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
                hetot->Fill(etracknorm);
		Bool_t delta_test = abs(e_delta)<4.;
		//if (nrun > 8416) delta_test = abs(e_delta)>-3.;
		if (etracknorm >0.5 && delta_test) {
                hxrast->Fill(rastx);
                hyrast->Fill(rasty);
                hxbeam->Fill(beamx);
                hybeam->Fill(beamy);
		if (abs(ysieve) < .7) hytar->Fill(e_ytar);
		if (abs(ysieve) < .7) hxsieve_ytar_yscent->Fill(xsieve,e_ytar);
		//hytar->Fill(e_ytar);
              hysieve->Fill(ysieve);
              hxsieve->Fill(xsieve);
	      if (abs(e_reactz)< 5) {
              hysieve_cent->Fill(ysieve);
              hxsieve_cent->Fill(xsieve);
	      }
                hxtar->Fill(e_xtar);
                if (abs(ysieve) < .7) hztar->Fill(e_reactz);
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
   //
}
