#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TChain.h"

#include <iostream>

void draw_ZDC_HF_HiForward_210607(int nevt=-1, int f=3)
{

	using namespace std;


	TString filter[9] = {"2Th1","2Th2","2Th3","2Th4","2Th5","2Th1p5","2Th2p5","2Th3p5","2Th4p5"};
	TString fName=filter[f];
	TString s="phfCoincFilter";
	TString title=s+fName;
	cout << "fName : " << fName.Data() << endl;
	TString DATE="210609";
	gSystem->mkdir(Form("figs/%s", DATE.Data()),kTRUE);

	gStyle->SetOptStat(00000);
	gROOT->ForceStyle();
	TString fin = "./roots/HiForward_000*.root";
	TChain *evttree = new TChain("evttree");
	TChain *hlttree = new TChain("hlttree");
	TChain *skimtree = new TChain("skimtree");

	evttree->Add(fin.Data());
	hlttree->Add(fin.Data());
	skimtree->Add(fin.Data());

	//TTree* evttree= (TTree*)fin->Get("evttree");
	//TTree* zdctree= (TTree*)fin->Get("rechitanalyzerpp/zdcrechit");
	//TTree* hlttree= (TTree*)fin->Get("hlttree");
	//TTree* skimtree= (TTree*)fin->Get("skimtree");

	//zdctree->AddFriend(evttree);
	evttree->AddFriend(hlttree);
	evttree->AddFriend(skimtree);

	//### define using variables

	int noareas=20;

	//	int nolines=nolines_region1+nolines_region2+nolines_region3+nolines_region4;
	int nPixBins_P, nZDCBins_P, nTrkBins_P, nHFBins_P;
	int nPixBins_F, nZDCBins_F, nTrkBins_F, nHFBins_F;
	double npix_max_P, zdc_max_P, ntrk_max_P, hf_max_P, npix_min_P, zdc_min_P, ntrk_min_P, hf_min_P;
		nPixBins_P=50;
		nZDCBins_P=50;
		nTrkBins_P=50;
		nHFBins_P =50;
		npix_max_P=110000.0;
		zdc_max_P=450000.0;  
		ntrk_max_P=4500;
		hf_max_P=6500;

		npix_min_P=0000.0;
		zdc_min_P=0000.0;
		ntrk_min_P=0000.0;
		hf_min_P=0000.0;

	double npix_max_F, zdc_max_F, ntrk_max_F, hf_max_F, npix_min_F, zdc_min_F, ntrk_min_F, hf_min_F;
		nPixBins_F=100;
		nZDCBins_F=500;
		nTrkBins_F=50;
		nHFBins_F =50;
		npix_max_F=600; 
		zdc_max_F=4500.0; 
		ntrk_max_F=50;
		hf_max_F=65.0;
		npix_min_F=0000.0;
		zdc_min_F=0000.0;
		ntrk_min_F=0000.0;
		hf_min_F=0000.0;

	// call variables in the input trees
	Int_t evt;
	int hiNpix, hiNpixPlus, hiNpixMinus, hiNtracks, hiBin, n;
	float hiHF;
	float e[18];
	int zside[18];
	int pprimaryVertexFilter, phfCoincFilter2Th1, phfCoincFilter2Th2, phfCoincFilter2Th3, phfCoincFilter2Th4, phfCoincFilter2Th5,phfCoincFilter2Th1p5,phfCoincFilter2Th2p5,phfCoincFilter2Th3p5,phfCoincFilter2Th4p5, pclusterCompatibilityFilter, HLT_HIZeroBias_v1;
	int CoincFilter2Th1, CoincFilter2Th2, CoincFilter2Th3, CoincFilter2Th4, CoincFilter2Th5, CoincFilter2Th1p5, CoincFilter2Th2p5, CoincFilter2Th3p5, CoincFilter2Th4p5;
	int CoincFilter2[9]={CoincFilter2Th1, CoincFilter2Th2, CoincFilter2Th3, CoincFilter2Th4, CoincFilter2Th5, CoincFilter2Th1p5, CoincFilter2Th2p5, CoincFilter2Th3p5, CoincFilter2Th4p5};

	evttree->SetBranchAddress("evt",&evt);
	evttree->SetBranchAddress("hiNpix",&hiNpix);
	evttree->SetBranchAddress("hiNtracks",&hiNtracks);
	evttree->SetBranchAddress("hiBin",&hiBin);
	evttree->SetBranchAddress("hiHF",&hiHF);

	//zdctree->SetBranchAddress("e",&e);
	//zdctree->SetBranchAddress("zside",zside);
	//zdctree->SetBranchAddress("zside",&zside);

	skimtree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	skimtree->SetBranchAddress("phfCoincFilter2Th1"  ,&phfCoincFilter2Th1);
	skimtree->SetBranchAddress("phfCoincFilter2Th2"  ,&phfCoincFilter2Th2);
	skimtree->SetBranchAddress("phfCoincFilter2Th3"  ,&phfCoincFilter2Th3);
	skimtree->SetBranchAddress("phfCoincFilter2Th4"  ,&phfCoincFilter2Th4);
	skimtree->SetBranchAddress("phfCoincFilter2Th5"  ,&phfCoincFilter2Th5);
	skimtree->SetBranchAddress("phfCoincFilter2Th1p5",&phfCoincFilter2Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th2p5",&phfCoincFilter2Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th3p5",&phfCoincFilter2Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th4p5",&phfCoincFilter2Th4p5);
	skimtree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

	hlttree->SetBranchAddress("HLT_HIZeroBias_v1", &HLT_HIZeroBias_v1);

	TH2D* hHF_Pixels_P = new TH2D("hHF_Pixels_P",";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
	TH2D* hHF_Ntracks_P = new TH2D("hHF_Tracks_P",";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
	TH2D* hPix_Trk_P = new TH2D("hPix_Trk_P",";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

	TH2D* hHF_Pixels_F = new TH2D("hHF_Pixels_F",";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
	TH2D* hHF_Ntracks_F = new TH2D("hHF_Tracks_F",";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
	TH2D* hPix_Trk_F = new TH2D("hPix_Trk_F",";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);


	TH1F *hHF_P = new TH1F("hHF_P","",100,0,300);
	TH1F *hTrk_P = new TH1F("hTrk_P","",100,0,300);
	TH1F *hPix_P = new TH1F("hPix_P","",100,0,300);
	TH1F *hHF_F = new TH1F("hHF_F","",100,0,300);
	TH1F *hTrk_F = new TH1F("hTrk_F","",100,0,300);
	TH1F *hPix_F = new TH1F("hPix_F","",100,0,300);

	if (nevt==-1) nevt = evttree->GetEntries();
	for (int iev=0;iev<nevt;iev++) {
		if(iev%10000000==0) cout << ">>>>> EVENT " << iev << " / " << evttree->GetEntries() <<  " ("<<(int)(100.*iev/evttree->GetEntries()) << "%)" << endl;
		evttree->GetEntry(iev);

	  if (f==0)CoincFilter2[f]=phfCoincFilter2Th1; else if(f==1) CoincFilter2[f]=phfCoincFilter2Th2;	else if(f==2) CoincFilter2[f]=phfCoincFilter2Th3; else if(f==3) CoincFilter2[f]=phfCoincFilter2Th4; else if(f==4) CoincFilter2[f]=phfCoincFilter2Th5;
		else if(f==5)CoincFilter2[f]=phfCoincFilter2Th1p5; else if(f==6) CoincFilter2[f]=phfCoincFilter2Th2p5; else if(f==7) CoincFilter2[f]=phfCoincFilter2Th3p5; else if(f==8) CoincFilter2[f]=phfCoincFilter2Th4p5; else if(f==9) CoincFilter2[f]=phfCoincFilter2Th4p5;

		if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[f])==1 ){  
			hHF_Pixels_P->Fill(hiNpix,hiHF);
			hHF_Ntracks_P->Fill(hiNtracks,hiHF);
			hPix_Trk_P->Fill(hiNpix,hiNtracks);
		}
		if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[f])==0 ){  
			hHF_Pixels_F->Fill(hiNpix,hiHF);
			hHF_Ntracks_F->Fill(hiNtracks,hiHF);
			hPix_Trk_F->Fill(hiNpix,hiNtracks);
		}
		if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter2[f]==1) ){
			hHF_P->Fill(hiHF);
			hTrk_P->Fill(hiNtracks);
			hPix_P->Fill(hiNpix);
		}
		if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[f])==0 ){
			hHF_F->Fill(hiHF);
			hTrk_F->Fill(hiNtracks);
			hPix_F->Fill(hiNpix);
		}
	}
	TCanvas* ctest = new TCanvas("ctest","",700,700);
	ctest->cd();
	ctest->SetLogy();
	hTrk_P->Draw();
	TCanvas* cHF1D = new TCanvas("cHF1D","",700,700);
	cHF1D->cd();
	cHF1D->SetLogy();
	hHF_P->SetLineColor(kBlue+2);
	//hHF_P->SetMarkerSize(2);
	hHF_F->SetLineColor(kRed+2);
	//hHF_F->SetMarkerSize(3);
	//hHF_F->SetMarkerStyle(33);
	hHF_P->GetXaxis()->SetTitle("HF");
	hHF_P->GetYaxis()->SetTitle("Entries");
	hHF_F->GetXaxis()->SetTitle("HF");
	hHF_F->GetYaxis()->SetTitle("Entries");
	hHF_P->GetYaxis()->SetRangeUser(0.5,1e+5);
	hHF_F->GetYaxis()->SetRangeUser(0.5,1e+5);
	hHF_P->GetYaxis()->CenterTitle();
	hHF_F->GetYaxis()->CenterTitle();
	hHF_P->GetXaxis()->CenterTitle();
	hHF_F->GetXaxis()->CenterTitle();
	hHF_P->GetYaxis()->SetTitleOffset(1.3);
	hHF_F->GetYaxis()->SetTitleOffset(1.3);
	hHF_P->SetTitle(Form("%s",title.Data()));
	hHF_F->SetTitle(Form("%s",title.Data()));
	hHF_P->Draw("");
	hHF_F->Draw("same");

  TLegend *leg1 = new TLegend(0.58,0.55,0.8,0.8);
	leg1->SetFillColor(0);
  leg1->SetFillStyle(4000);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.2);
  leg1->SetTextSize(0.029);
  leg1->SetTextFont(42);
	leg1->AddEntry(hHF_P,"Passed Condition HF");
	leg1->AddEntry(hHF_F,"Failed Condition HF");
	leg1->Draw("same");
	cHF1D->SaveAs(Form("figs/%s/1D_HF_HiForward_%s.png",DATE.Data(),fName.Data()));

	TCanvas* cTrk1D = new TCanvas("cTrk1D","",700,700);
	cTrk1D->cd();
	cTrk1D->SetLogy();
	hTrk_P->SetLineColor(kBlue+2);
	hTrk_F->SetLineColor(kRed+2);
	hTrk_P->GetXaxis()->SetTitle("Num. of Trackers");
	hTrk_P->GetYaxis()->SetTitle("Entries");
	hTrk_F->GetXaxis()->SetTitle("Num. of Trackers");
	hTrk_F->GetYaxis()->SetTitle("Entries");
	hTrk_P->GetYaxis()->SetRangeUser(0.5,1e+5);
	hTrk_F->GetYaxis()->SetRangeUser(0.5,1e+5);
	hTrk_P->GetYaxis()->CenterTitle();
	hTrk_F->GetYaxis()->CenterTitle();
	hTrk_P->GetXaxis()->CenterTitle();
	hTrk_F->GetXaxis()->CenterTitle();
	hTrk_P->GetYaxis()->SetTitleOffset(1.3);
	hTrk_F->GetYaxis()->SetTitleOffset(1.3);
	hTrk_P->SetTitle(Form("%s",title.Data()));
	hTrk_F->SetTitle(Form("%s",title.Data()));
	hTrk_P->Draw("");
	hTrk_F->Draw("same");

  TLegend *leg2 = new TLegend(0.58,0.55,0.8,0.8);
	leg2->SetFillColor(0);
  leg2->SetFillStyle(4000);
  leg2->SetBorderSize(0);
  leg2->SetMargin(0.2);
  leg2->SetTextSize(0.029);
  leg2->SetTextFont(42);
	leg2->AddEntry(hHF_P,"Passed Condition Tracker");
	leg2->AddEntry(hHF_F,"Failed Condition Tracker");
	leg2->Draw("same");
	cTrk1D->SaveAs(Form("figs/%s/1D_Trk_HiForward_%s.png",DATE.Data(),fName.Data()));

	TCanvas* cPix1D = new TCanvas("cPix1D","",700,700);
	cPix1D->cd();
	cPix1D->SetLogy();
	hPix_P->SetLineColor(kBlue+2);
	hPix_F->SetLineColor(kRed+2);
	hPix_P->GetXaxis()->SetTitle("Num. of Pixels");
	hPix_P->GetYaxis()->SetTitle("Entries");
	hPix_F->GetXaxis()->SetTitle("Num. of Pixels");
	hPix_F->GetYaxis()->SetTitle("Entries");
	hPix_P->GetYaxis()->SetRangeUser(0.5,1e+5);
	hPix_F->GetYaxis()->SetRangeUser(0.5,1e+5);
	hPix_P->GetYaxis()->CenterTitle();
	hPix_F->GetYaxis()->CenterTitle();
	hPix_P->GetXaxis()->CenterTitle();
	hPix_F->GetXaxis()->CenterTitle();
	hPix_P->GetYaxis()->SetTitleOffset(1.3);
	hPix_F->GetYaxis()->SetTitleOffset(1.3);
	hPix_P->SetTitle(Form("%s",title.Data()));
	hPix_F->SetTitle(Form("%s",title.Data()));
	hPix_P->Draw("");
	hPix_F->Draw("same");

  TLegend *leg3 = new TLegend(0.58,0.55,0.8,0.8);
	leg3->SetFillColor(0);
  leg3->SetFillStyle(4000);
  leg3->SetBorderSize(0);
  leg3->SetMargin(0.2);
  leg3->SetTextSize(0.029);
  leg3->SetTextFont(42);
	leg3->AddEntry(hHF_P,"Passed Condition Pixel");
	leg3->AddEntry(hHF_F,"Failed Condition Pixel");
	leg3->Draw("same");
	cPix1D->SaveAs(Form("figs/%s/1D_Pix_HiForward_%s.png",DATE.Data(),fName.Data()));

	TCanvas* c2 = new TCanvas("c2","",1600,800);
	c2->Divide(2,1);
	c2->cd(1);
	hHF_Pixels_P->GetYaxis()->SetMaxDigits(4);
	hHF_Pixels_P->GetYaxis()->CenterTitle();
	hHF_Pixels_P->GetYaxis()->SetLabelSize(0.03);
	hHF_Pixels_P->GetYaxis()->SetTitleOffset(1.5);
	hHF_Pixels_P->GetXaxis()->SetMaxDigits(4);
	hHF_Pixels_P->GetXaxis()->CenterTitle();
	hHF_Pixels_P->GetXaxis()->SetLabelSize(0.03);
	hHF_Pixels_P->SetTitle(Form("HF vs Number of Pixels Passed Cond. (%s)",title.Data()));
	hHF_Pixels_P->Draw("colz");
	c2->cd(2);
	hHF_Ntracks_P->GetYaxis()->SetMaxDigits(4);
	hHF_Ntracks_P->GetYaxis()->CenterTitle();
	hHF_Ntracks_P->GetYaxis()->SetLabelSize(0.03);
	hHF_Ntracks_P->GetYaxis()->SetTitleOffset(1.5);
	hHF_Ntracks_P->GetXaxis()->SetMaxDigits(4);
	hHF_Ntracks_P->GetXaxis()->CenterTitle();
	hHF_Ntracks_P->GetXaxis()->SetLabelSize(0.03);
	hHF_Ntracks_P->GetXaxis()->SetNdivisions(511);
	hHF_Ntracks_P->SetTitle(Form("HF vs Number of Tracks Passed Cond. (%s)",title.Data()));
	hHF_Ntracks_P->Draw("colz");
	c2->SaveAs(Form("figs/%s/HF_pixel_track_HiForward_Passed_%s.png",DATE.Data(),fName.Data()));

	TCanvas* c1 = new TCanvas("c1","",800,800);
	c1->cd();
	hPix_Trk_P->GetYaxis()->SetMaxDigits(3);
	hPix_Trk_P->GetYaxis()->CenterTitle();
	hPix_Trk_P->GetYaxis()->SetLabelSize(0.03);
	hPix_Trk_P->GetYaxis()->SetTitleOffset(1.5);
	hPix_Trk_P->GetXaxis()->SetMaxDigits(3);
	hPix_Trk_P->GetXaxis()->CenterTitle();
	hPix_Trk_P->GetXaxis()->SetLabelSize(0.03);
	hPix_Trk_P->SetTitle(Form("Pixel vs Tracks Passed Cond. (%s)",title.Data()));
	hPix_Trk_P->Draw("colz");
	c1->SaveAs(Form("figs/%s/Npix_Ntrk_HiForward_Passed_%s.png",DATE.Data(),fName.Data()));

	TCanvas* c3 = new TCanvas("c3","",1600,800);
	c3->Divide(2,1);
	c3->cd(1);
	hHF_Pixels_F->GetYaxis()->SetMaxDigits(4);
	hHF_Pixels_F->GetYaxis()->CenterTitle();
	hHF_Pixels_F->GetYaxis()->SetLabelSize(0.03);
	hHF_Pixels_F->GetYaxis()->SetTitleOffset(1.5);
	hHF_Pixels_F->GetXaxis()->SetMaxDigits(4);
	hHF_Pixels_F->GetXaxis()->CenterTitle();
	hHF_Pixels_F->GetXaxis()->SetLabelSize(0.03);
	hHF_Pixels_F->SetTitle(Form("HF vs Number of Pixels Failed Cond.(%s)",title.Data()));
	hHF_Pixels_F->Draw("colz");
	c3->cd(2);
	hHF_Ntracks_F->GetYaxis()->SetMaxDigits(4);
	hHF_Ntracks_F->GetYaxis()->CenterTitle();
	hHF_Ntracks_F->GetYaxis()->SetLabelSize(0.03);
	hHF_Ntracks_F->GetYaxis()->SetTitleOffset(1.5);
	hHF_Ntracks_F->GetXaxis()->SetMaxDigits(4);
	hHF_Ntracks_F->GetXaxis()->CenterTitle();
	hHF_Ntracks_F->GetXaxis()->SetLabelSize(0.03);
	hHF_Ntracks_F->GetXaxis()->SetNdivisions(511);
	hHF_Ntracks_F->SetTitle(Form("HF vs Number of Tracks Failed Cond. (%s)",title.Data()));
	hHF_Ntracks_F->Draw("colz");
	c3->SaveAs(Form("figs/%s/HF_pixel_track_HiForward_Failed_%s.png",DATE.Data(),fName.Data()));

	TCanvas* c4 = new TCanvas("c4","",800,800);
	c4->cd();
	hPix_Trk_F->GetYaxis()->SetMaxDigits(3);
	hPix_Trk_F->GetYaxis()->CenterTitle();
	hPix_Trk_F->GetYaxis()->SetLabelSize(0.03);
	hPix_Trk_F->GetYaxis()->SetTitleOffset(1.5);
	hPix_Trk_F->GetXaxis()->SetMaxDigits(3);
	hPix_Trk_F->GetXaxis()->CenterTitle();
	hPix_Trk_F->GetXaxis()->SetLabelSize(0.03);
	hPix_Trk_F->SetTitle(Form("Pixel vs Tracks Failed Cond. (%s)",title.Data()));
	hPix_Trk_F->Draw("colz");
	c4->SaveAs(Form("figs/%s/Npix_Ntrk_HiForward_Failed_%s.png",DATE.Data(),fName.Data()));
}
