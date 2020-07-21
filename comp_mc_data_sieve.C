#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TLine.h>
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


void comp_mc_data_sieve(Int_t nrun=1814,Int_t FileID=-2) {
  
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
  TString OpticsID="";
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Double_t ymis =0.0;
  Int_t NumFoil=0;
  TString temp;
//
  vector <Double_t> ztar_foil;
  Int_t ndelcut=-1;
  vector<Double_t > delcut;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      SieveFlag = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      ndelcut = temp.Atoi();
      temp.ReadToDelim(file_optics);
      ymis = temp.Atof();
      for (Int_t nf=0;nf<NumFoil-1;nf++) {
        temp.ReadToDelim(file_optics,',');
	ztar_foil.push_back(temp.Atof());
      }
        temp.ReadToDelim(file_optics);
	ztar_foil.push_back(temp.Atof());
      for (Int_t nd=0;nd<ndelcut;nd++) {
        temp.ReadToDelim(file_optics,',');
	delcut.push_back(temp.Atof());
        cout << " nd = " << nd << " " << delcut[nd] << endl;
      }
        temp.ReadToDelim(file_optics);
	delcut.push_back(temp.Atof());
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << endl;
  if (NumFoil==0) return;
    CentAngle=CentAngle*3.14159/180.;
 //
  Double_t ys_offset= 0.0;
  vector <Double_t> ys_cent;
  vector <Double_t> xs_cent;
  for (Int_t nys=0;nys<9;nys++) {
    Double_t pos=(nys-4)*0.6*2.54;
    ys_cent.push_back(pos);
  }
 //
   vector<vector<Double_t> > AxisRange;
   AxisRange.resize(ndelcut);
   TString AxisRangeFileName = Form("AxisRange_ypfp_yfp_%d.dat",nrun);
   ifstream AxisRangeFile(AxisRangeFileName.Data());
   for  (Int_t nd=0;nd<ndelcut;nd++) { AxisRange[nd].resize(4) ;}
  TString temp1;
  if (file_optics.is_open()) {
    cout << " Axis range file = " << AxisRangeFileName << endl;
   for  (Int_t nd=0;nd<ndelcut;nd++) { 
     temp1.ReadToDelim(AxisRangeFile,',');
     AxisRange[nd][0] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeFile,',');
     AxisRange[nd][1] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeFile,',');
     AxisRange[nd][2] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeFile);
     AxisRange[nd][3] = temp1.Atof();
      }
  } else {
    cout << " no file = " << AxisRangeFileName << endl;
    return;
  }
 //
 //
   vector<vector<Double_t> > AxisRangeXfp;
   AxisRangeXfp.resize(ndelcut);
   TString AxisRangeXfpFileName = Form("AxisRange_xpfp_xfp_%d.dat",nrun);
   ifstream AxisRangeXfpFile(AxisRangeXfpFileName.Data());
   for  (Int_t nd=0;nd<ndelcut;nd++) { AxisRangeXfp[nd].resize(4) ;}
   if (file_optics.is_open()) {
    cout << " Axis range file = " << AxisRangeXfpFileName << endl;
   for  (Int_t nd=0;nd<ndelcut;nd++) { 
     temp1.ReadToDelim(AxisRangeXfpFile,',');
     AxisRangeXfp[nd][0] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeXfpFile,',');
     AxisRangeXfp[nd][1] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeXfpFile,',');
     AxisRangeXfp[nd][2] = temp1.Atof();
     temp1.ReadToDelim(AxisRangeXfpFile);
     AxisRangeXfp[nd][3] = temp1.Atof();
      }
  } else {
    cout << " no file = " << AxisRangeXfpFileName << endl;
    return;
  }
 //
//
  TString outCutFile;
  TFile *fcut;
  vector<vector<vector<TCutG*> > > ypfp_ypfp_cut;
  vector<vector<vector<Int_t> > > ypfp_ypfp_cut_flag;
  ypfp_ypfp_cut.resize(NumFoil);
  ypfp_ypfp_cut_flag.resize(NumFoil);
	for  (Int_t nf=0;nf<NumFoil;nf++) {
          ypfp_ypfp_cut[nf].resize(ndelcut);
          ypfp_ypfp_cut_flag[nf].resize(ndelcut);
   }
  Bool_t cflag=kTRUE;
  if (cflag) {
   outCutFile=Form("cuts/YpFpYFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcut = new TFile(outCutFile);
    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<9;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	  ypfp_ypfp_cut[nf][nd].push_back(tempg);
      } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  ypfp_ypfp_cut[nf][nd].push_back(tempg);
      }
	}}}
  }
//
//
  TString xpfp_xfp_outCutFile;
  TFile *xpfp_xfp_fcut;
  vector<vector<vector<TCutG*> > > xpfp_xfp_cut;
  vector<vector<vector<Int_t> > > xpfp_xfp_cut_flag;
  xpfp_xfp_cut.resize(NumFoil);
  xpfp_xfp_cut_flag.resize(NumFoil);
	for  (Int_t nf=0;nf<NumFoil;nf++) {
          xpfp_xfp_cut[nf].resize(ndelcut);
          xpfp_xfp_cut_flag[nf].resize(ndelcut);
   }
  if (cflag) {
    xpfp_xfp_outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    xpfp_xfp_fcut = new TFile(xpfp_xfp_outCutFile);
    cout << "xpfp_xfp_ Cut file = " << xpfp_xfp_outCutFile << endl;
    xpfp_xfp_fcut->cd();
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<9;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	  xpfp_xfp_cut[nf][nd].push_back(tempg);
      } else {
	    //cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  xpfp_xfp_cut[nf][nd].push_back(tempg);
      }
	}}}
  }
//

  //
  TString inputroot;
   TFile *fhistroot;
   inputroot=Form("mc_rootfiles/hms_run%d_hist.root",nrun);
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);

	vector<vector<TH2F*> > hYsXs_DelCut;
	vector<vector<TH2F*> > hYpFpYFp_DelCut;
	vector<vector<TH2F*> > hXpFpXFp_DelCut;
	  vector<TH2F*> temp2d;
	TH2F* th;
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  th =  (TH2F*)fhistroot->Get(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd));
	  if (!th) {
	    cout << " no hist : "<< Form("hYsXs_Foil_%d_DelCut_%d",nc,nd) << endl;
	    return;
	  }
	  temp2d.push_back(th);
	}	  
	hYsXs_DelCut.push_back(temp2d);
	  temp2d.clear();
	}	
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  th =  (TH2F*)fhistroot->Get(Form("hYpFpYFp_%d_DelCut_%d",nc,nd));
	  if (!th) {
	    cout << " no hist : "<< Form("hYpFpYFp_%d_DelCut_%d",nc,nd)<< endl;
	    return;
	  }
	  temp2d.push_back(th);
	}	  
	hYpFpYFp_DelCut.push_back(temp2d);
	  temp2d.clear();	
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  th =  (TH2F*)fhistroot->Get(Form("hXpFpXFp_%d_DelCut_%d",nc,nd));
	  if (!th) {
	    cout << " no hist : "<< Form("hXpFpXFp_%d_DelCut_%d",nc,nd)<< endl;
	    return;
	  }
	  temp2d.push_back(th);
	}	  
	hXpFpXFp_DelCut.push_back(temp2d);
	  temp2d.clear();	
	}
	}
	//
   TFile *fhistroot_data;
   inputroot=Form("hist/Optics_%s_%d_hist.root",OpticsID.Data(),FileID);
   TString outputpdf =Form("Optics_%s_%d_mc_comp",OpticsID.Data(),FileID);
     cout << "data  infile root = " << inputroot << endl;
   fhistroot_data =  new TFile(inputroot);
	vector<vector<TH2F*> > hYsXs_data_DelCut;
	vector<vector<TH2F*> > hYpFpYFp_data_DelCut;
	vector<vector<TH2F*> > hXpFpXFp_data_DelCut;
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  th =  (TH2F*)fhistroot_data->Get(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd));
	  if (!th) {
	    cout << " no hist : "<< Form("hYsXs_Foil_%d_DelCut_%d",nc,nd) << endl;
	    return;
	  }
	  temp2d.push_back(th);
	}	  
	hYsXs_data_DelCut.push_back(temp2d);
	  temp2d.clear();
	}	
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  th =  (TH2F*)fhistroot_data->Get(Form("hYpFpYFp_%d_DelCut_%d",nc,nd));
	  if (!th) {
	    cout << " no hist : "<< Form("hYpFpYFp_%d_DelCut_%d",nc,nd) << endl;
	    return;
	  }
	  temp2d.push_back(th);
	}	  
	hYpFpYFp_data_DelCut.push_back(temp2d);
	  temp2d.clear();	
	}
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  th =  (TH2F*)fhistroot_data->Get(Form("hXpFpXFp_%d_DelCut_%d",nc,nd));
	  if (!th) {
	    cout << " no hist : "<< Form("hXpFpXFp_%d_DelCut_%d",nc,nd) << endl;
	    return;
	  }
	  temp2d.push_back(th);
	}	  
	hXpFpXFp_data_DelCut.push_back(temp2d);
	  temp2d.clear();	
	}

 //
	string oldcoeffsfilename="hms_recon_coeff_opt2018.dat";
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  vector<Double_t> xfpcoeffs;
  vector<Double_t> xpfpcoeffs;
  vector<Double_t> yfpcoeffs;
  vector<Double_t> ypfpcoeffs;
  vector<Double_t> lencoeffs;
  vector<Int_t> xtarexpon;
  vector<Int_t> xptarexpon;
  vector<Int_t> ytarexpon;
  vector<Int_t> yptarexpon;
  vector<Int_t> deltaexpon;
  TString currentline;
  int num_recon_terms;
  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){    num_recon_terms++;}
  cout << num_recon_terms;
  oldcoeffsfile.close();
  oldcoeffsfile.open(oldcoeffsfilename.c_str());
  xfpcoeffs.resize(num_recon_terms);
   xpfpcoeffs.resize(num_recon_terms);
  yfpcoeffs.resize(num_recon_terms);
   ypfpcoeffs.resize(num_recon_terms);
   lencoeffs.resize(num_recon_terms);
   xtarexpon.resize(num_recon_terms);
   xptarexpon.resize(num_recon_terms);
   ytarexpon.resize(num_recon_terms);
   yptarexpon.resize(num_recon_terms);
   deltaexpon.resize(num_recon_terms);
   num_recon_terms=0;
 while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){

    TString sc1(currentline(1,14));
    TString sc2(currentline(15,14));
    TString sc3(currentline(29,14));
    TString sc4(currentline(43,14));
    TString sc5(currentline(57,14));
    
    xfpcoeffs[num_recon_terms]=sc1.Atof();
    xpfpcoeffs[num_recon_terms]=sc2.Atof();
    yfpcoeffs[num_recon_terms]=sc3.Atof();
    ypfpcoeffs[num_recon_terms]=sc4.Atof();
    lencoeffs[num_recon_terms]=sc5.Atof();
    
    int expontemp[6];

    for(int expon=0; expon<6; expon++){
      TString stemp(currentline(72+expon,1));
      expontemp[expon] = stemp.Atoi();
    }

  

    xtarexpon[num_recon_terms]=expontemp[0];
    xptarexpon[num_recon_terms]=expontemp[1];
    ytarexpon[num_recon_terms]=expontemp[2];
    yptarexpon[num_recon_terms]=expontemp[3];
    deltaexpon[num_recon_terms]=expontemp[5];
    //cout << xfpcoeffs[num_recon_terms] << " " <<  xpfpcoeffs[num_recon_terms]<< " " <<  yfpcoeffs[num_recon_terms] << " " <<  ypfpcoeffs[num_recon_terms]<< " " << xtarexpon[num_recon_terms] << xptarexpon[num_recon_terms] << ytarexpon[num_recon_terms] << yptarexpon[num_recon_terms]<< expontemp[4]<< deltaexpon[num_recon_terms]  << endl;
   num_recon_terms++;
    
  }
  //
	Double_t xtar=0,xptar=0;
	vector<vector<vector<Double_t> > > vyfp;
	vector<vector<vector<Double_t> > > vypfp;
	vector<vector<vector<TLine*> > > vline;
	vyfp.resize(9);
	vypfp.resize(9);
	vline.resize(9);
	for  (Int_t nd=0;nd<9;nd++) {
	vyfp[nd].resize(9);
	vypfp[nd].resize(9);
	vline[nd].resize(9);
	  for  (Int_t nf=0;nf<ndelcut;nf++) {
	    vyfp[nd][nf].resize(NumFoil);
	    vypfp[nd][nf].resize(NumFoil);
	    vline[nd][nf].resize(NumFoil);
	  }
        }
	Double_t zdis_sieve = 168.0;
        Double_t ytar,yptar;
	Double_t yfpfirst,ypfpfirst;
	for  (Int_t nyscol=0;nyscol<9;nyscol++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	Double_t DelCent = (delcut[nd+1]+delcut[nd])/2;
	yptar = (ys_cent[nyscol]-ztar_foil[nf]*TMath::Sin(CentAngle))/(zdis_sieve-ztar_foil[nf]*TMath::Cos(CentAngle));
	ytar = +ztar_foil[nf]*(TMath::Sin(CentAngle)-yptar*TMath::Cos(CentAngle));
	yfpfirst = -1.89*ytar-0.16*yptar-0.09*DelCent;
	ypfpfirst = (-3.0*ytar-0.77*yptar+0.13*DelCent)/1000.;
	//	cout << " nf = " << nf << " nd = " << nd << " ytar = "  << " yptar = "  << yptar << " " << ys_cent[nyscol] << " yfp = "  << yfpfirst<< " ypfp = "  << ypfpfirst<< endl;
	//
          Double_t yfp = 0.0,ypfp=0.0,xpfp=0.0,xfp=0.0;
          Double_t yfplo = 0.0,ypfplo=0.0;
          Double_t yfphi = 0.0,ypfphi=0.0;
          Double_t etemp;
          for( int icoeffold=0; icoeffold<num_recon_terms; icoeffold++ ){
        	etemp= 
	  pow( xtar, xtarexpon[icoeffold] ) * 
	  pow( xptar, xptarexpon[icoeffold] ) * 
	  pow( ytar, ytarexpon[icoeffold] ) * 
	  pow( yptar*1000, yptarexpon[icoeffold] ) * 
	  pow( DelCent, deltaexpon[icoeffold] );
        	xfp += xfpcoeffs[icoeffold] * etemp;
        	xpfp += xpfpcoeffs[icoeffold] * etemp;
	        yfp += yfpcoeffs[icoeffold] * etemp;
	        ypfp += ypfpcoeffs[icoeffold] * etemp;
        	etemp= 
	  pow( xtar, xtarexpon[icoeffold] ) * 
	  pow( xptar, xptarexpon[icoeffold] ) * 
	  pow( ytar, ytarexpon[icoeffold] ) * 
	  pow( yptar*1000, yptarexpon[icoeffold] ) * 
	  pow(delcut[nd] , deltaexpon[icoeffold] );
	        yfplo += yfpcoeffs[icoeffold] * etemp;
	        ypfplo += ypfpcoeffs[icoeffold] *etemp; 
        	etemp= 
	  pow( xtar, xtarexpon[icoeffold] ) * 
	  pow( xptar, xptarexpon[icoeffold] ) * 
	  pow( ytar, ytarexpon[icoeffold] ) * 
	  pow( yptar*1000, yptarexpon[icoeffold] ) * 
	  pow(delcut[nd+1] , deltaexpon[icoeffold] );
	        yfphi += yfpcoeffs[icoeffold] * etemp;
	        ypfphi += ypfpcoeffs[icoeffold] *etemp; 
	  } // for icoeffold loop
	  vyfp[nyscol][nd][nf]= yfp;
	  vypfp[nyscol][nd][nf] = ypfp;
	  vline[nyscol][nd][nf]= new TLine(ypfplo/1000.,yfplo,ypfphi/1000.,yfphi);
		  vline[nyscol][nd][nf]->SetLineColor(2);
		  vline[nyscol][nd][nf]->SetLineWidth(3);
		  // cout << " full matrix yfp = " << vyfp[nyscol][nd][nf]<< " ypfp = " << vypfp[nyscol][nd][nf] << endl;
	//
	}}} 

	//
	  TLine* ys_line[9];
	  TText* ys_text[9];
	  TLine* xs_line[9];
	  TText* xs_text[9];
  for (Int_t nys=0;nys<9;nys++) {
    Double_t pos=(nys-4)*0.6*2.54;
    ys_line[nys]= new TLine(pos,-12.5 ,pos,12.5);
    ys_text[nys]= new TText(pos,-13.5,Form("%d",nys));
     ys_text[nys]->SetTextColor(2);
     ys_line[nys]->SetLineColor(2);
     ys_line[nys]->SetLineWidth(1);
    pos=(nys-4)*2.54;
    xs_line[nys]= new TLine(-7.,pos,7.,pos);
    xs_text[nys]= new TText(-7.5,pos,Form("%d",nys));
     xs_text[nys]->SetTextColor(2);
     xs_line[nys]->SetLineColor(2);
     xs_line[nys]->SetLineWidth(1);
  }
	//
	TCanvas* can2d[3][10];
	Int_t nctot=0;
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	    for  (Int_t nd=0;nd<ndelcut;nd++) {
	      can2d[nc][nd] = new TCanvas(Form("Can2d_%d_%d",nc,nd),Form("Foil %d Del %d",nc,nd), 700,700);
	  can2d[nc][nd]->Divide(2,3);
	  can2d[nc][nd]->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	//gPad->SetLogz();    
	  hYpFpYFp_data_DelCut[nc][nd]->Draw("colz");
	  hYpFpYFp_data_DelCut[nc][nd]->SetMinimum(1);
	  hYpFpYFp_data_DelCut[nc][nd]->GetYaxis()->SetRangeUser(AxisRange[nd][0],AxisRange[nd][1]);
	  hYpFpYFp_data_DelCut[nc][nd]->GetXaxis()->SetRangeUser(AxisRange[nd][2],AxisRange[nd][3]);
	   for (Int_t nys=0;nys<9;nys++) {
	     if (ypfp_ypfp_cut[nc][nd][nys]) ypfp_ypfp_cut[nc][nd][nys]->Draw("same");
	     //	     if (ypfp_ypfp_cut[nc][nd][nys]) ypfp_ypfp_cut[nc][nd][nys]->SetLineColor(2);
	     	     if (ypfp_ypfp_cut[nc][nd][nys]) ypfp_ypfp_cut[nc][nd][nys]->SetLineWidth(2);
		}
	  can2d[nc][nd]->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	//gPad->SetLogz();    
	  hYpFpYFp_DelCut[nc][nd]->Draw("colz");
	  hYpFpYFp_DelCut[nc][nd]->SetMinimum(1);
	  hYpFpYFp_DelCut[nc][nd]->GetYaxis()->SetRangeUser(AxisRange[nd][0],AxisRange[nd][1]);
	  hYpFpYFp_DelCut[nc][nd]->GetXaxis()->SetRangeUser(AxisRange[nd][2],AxisRange[nd][3]);
	  can2d[nc][nd]->cd(3);
	gPad->SetGridx();
	gPad->SetGridy();
	//gPad->SetLogz();    
	  hXpFpXFp_data_DelCut[nc][nd]->Draw("colz");
	  hXpFpXFp_data_DelCut[nc][nd]->SetMinimum(1);
	  hXpFpXFp_data_DelCut[nc][nd]->GetYaxis()->SetRangeUser(AxisRangeXfp[nd][0],AxisRangeXfp[nd][1]);
	  hXpFpXFp_data_DelCut[nc][nd]->GetXaxis()->SetRangeUser(AxisRangeXfp[nd][2],AxisRangeXfp[nd][3]);
	  	   for (Int_t nys=0;nys<9;nys++) {
	     if (xpfp_xfp_cut[nc][nd][nys]) xpfp_xfp_cut[nc][nd][nys]->Draw("same");
	     	     if (xpfp_xfp_cut[nc][nd][nys]) xpfp_xfp_cut[nc][nd][nys]->SetLineWidth(2);
	  	}
	  can2d[nc][nd]->cd(4);
	gPad->SetGridx();
	gPad->SetGridy();
	//gPad->SetLogz();    
	  hXpFpXFp_DelCut[nc][nd]->Draw("colz");
	  hXpFpXFp_DelCut[nc][nd]->SetMinimum(1);
	  hXpFpXFp_DelCut[nc][nd]->GetYaxis()->SetRangeUser(AxisRangeXfp[nd][0],AxisRangeXfp[nd][1]);
	  hXpFpXFp_DelCut[nc][nd]->GetXaxis()->SetRangeUser(AxisRangeXfp[nd][2],AxisRangeXfp[nd][3]);
	  can2d[nc][nd]->cd(5);
	  hYsXs_data_DelCut[nc][nd]->Draw("colz");
	  hYsXs_data_DelCut[nc][nd]->SetMinimum(1);
          for (Int_t nys=0;nys<9;nys++) { ys_line[nys]->Draw();}
          for (Int_t nys=0;nys<9;nys++) { ys_text[nys]->Draw();}
          for (Int_t nys=0;nys<9;nys++) { xs_line[nys]->Draw();}
          for (Int_t nys=0;nys<9;nys++) { xs_text[nys]->Draw();}
	  can2d[nc][nd]->cd(6);
	  hYsXs_DelCut[nc][nd]->Draw("colz");
	  hYsXs_DelCut[nc][nd]->SetMinimum(1);
          for (Int_t nys=0;nys<9;nys++) { ys_line[nys]->Draw();}
          for (Int_t nys=0;nys<9;nys++) { ys_text[nys]->Draw();}
          for (Int_t nys=0;nys<9;nys++) { xs_line[nys]->Draw();}
          for (Int_t nys=0;nys<9;nys++) { xs_text[nys]->Draw();}
	  TString end = ".pdf";
	  if (nc==0 && nd==0) end=".pdf(";
	  if (nc==NumFoil-1 && nd==ndelcut-1) end=".pdf)";
	  can2d[nc][nd]->Print("plots/"+outputpdf+end);
	}
	}
}
