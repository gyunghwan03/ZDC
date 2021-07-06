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
#include "TStopwatch.h"

#include <iostream>

void draw_ZDC_HF_HiForward_v3(int nevt=-1) 
{

	using namespace std;

	TStopwatch *t = new TStopwatch;
	t -> Start();

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
	nPixBins_P=500;
	nZDCBins_P=500;
	nTrkBins_P=500;
	nHFBins_P =500;
	npix_max_P=110000.0;
	zdc_max_P=450000.0;  
	ntrk_max_P=4500;
	hf_max_P=6500;

	npix_min_P=0000.0;
	zdc_min_P=0000.0;
	ntrk_min_P=0000.0;
	hf_min_P=0000.0;

	double npix_max_F, zdc_max_F, ntrk_max_F, hf_max_F, npix_min_F, zdc_min_F, ntrk_min_F, hf_min_F;
	nPixBins_F=500;
	nZDCBins_F=500;
	nTrkBins_F=500;
	nHFBins_F =500;
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
	//int CoincFilter2[9]={phfCoincFilter2Th1, phfCoincFilter2Th2, phfCoincFilter2Th3, phfCoincFilter2Th4, phfCoincFilter2Th5,phfCoincFilter2Th1p5,phfCoincFilter2Th2p5,phfCoincFilter2Th3p5,phfCoincFilter2Th4p5}

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

	int f;
	TString filter[9] = {"2Th1","2Th2","2Th3","2Th4","2Th5","2Th1p5","2Th2p5","2Th3p5","2Th4p5"};
	TString s="phfCoincFilter";
	//cout << "fName : " << fName.Data() << endl;

	TH2D* hHF_Pixels_P[9]; 
	TH2D* hHF_Ntracks_P[9];
	TH2D* hPix_Trk_P[9];

	TH2D* hHF_Pixels_F[9]; 
	TH2D* hHF_Ntracks_F[9];
	TH2D* hPix_Trk_F[9];


	TH1F *hHF_P[9]; 
	TH1F *hTrk_P[9];
	TH1F *hPix_P[9];
	TH1F *hHF_F[9];
	TH1F *hTrk_F[9];
	TH1F *hPix_F[9];


	for(int j=0;j<9;j++){
		hHF_Pixels_P[j]  = new TH2D(Form("hHF_Pixels_P_%d",j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hHF_Ntracks_P[j] = new TH2D(Form("hHF_Tracks_P_%d",j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hPix_Trk_P[j]    = new TH2D(Form("hPix_Trk_P_%d",j),";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

		hHF_Pixels_F[j]  = new TH2D(Form("hHF_Pixels_F_%d",j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hHF_Ntracks_F[j] = new TH2D(Form("hHF_Tracks_F_%d",j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hPix_Trk_F[j]    = new TH2D(Form("hPix_Trk_F_%d",j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);


		hHF_P[j]  = new TH1F(Form("hHF_P_%d",j),"",100,0,300);
		hHF_F[j]  = new TH1F(Form("hHF_F_%d",j),"",100,0,300);
		hTrk_P[j] = new TH1F(Form("hTrk_P_%d",j),"",100,0,300);
		hPix_P[j] = new TH1F(Form("hPix_P_%d",j),"",100,0,300);
		hTrk_F[j] = new TH1F(Form("hTrk_F_%d",j),"",100,0,300);
		hPix_F[j] = new TH1F(Form("hPix_F_%d",j),"",100,0,300);
	}
	if (nevt==-1) nevt = evttree->GetEntries();
	for (int iev=0;iev<nevt;iev++) {
		if(iev%10000000==0) cout << ">>>>> EVENT " << iev << " / " << evttree->GetEntries() <<  " ("<<(int)(100.*iev/evttree->GetEntries()) << "%)" << endl;
		evttree->GetEntry(iev);
		hlttree->GetEntry(iev);
		skimtree->GetEntry(iev);
		for (int j=0;j<9;j++){
			if (j==0)CoincFilter2[j]=phfCoincFilter2Th1; else if(j==1) CoincFilter2[j]=phfCoincFilter2Th2;	else if(j==2) CoincFilter2[j]=phfCoincFilter2Th3; else if(j==3) CoincFilter2[j]=phfCoincFilter2Th4; else if(j==4) CoincFilter2[j]=phfCoincFilter2Th5;
			else if(j==5)CoincFilter2[j]=phfCoincFilter2Th1p5; else if(j==6) CoincFilter2[j]=phfCoincFilter2Th2p5; else if(j==7) CoincFilter2[j]=phfCoincFilter2Th3p5; else if(j==8) CoincFilter2[j]=phfCoincFilter2Th4p5; else if(j==9) CoincFilter2[j]=phfCoincFilter2Th4p5;
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[j])==1 ){  
				hHF_Pixels_P[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_P[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_P[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[j])==0 ){  
				hHF_Pixels_F[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_F[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_F[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter2[j]==1) ){
				hHF_P[j]->Fill(hiHF);
				hTrk_P[j]->Fill(hiNtracks);
				hPix_P[j]->Fill(hiNpix);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[j])==0 ){
				hHF_F[j]->Fill(hiHF);
				hTrk_F[j]->Fill(hiNtracks);
				hPix_F[j]->Fill(hiNpix);
			}
		}
	}
	//TCanvas* ctest = new TCanvas("ctest","",700,700);
	//ctest->cd();
	//ctest->SetLogy();
	//hTrk_P->Draw();
	TCanvas* cHF1D_P = new TCanvas("cHF1D_P","",700,700);
	cHF1D_P->cd();
	cHF1D_P->SetLogy();
	TLegend *leg1_P = new TLegend(0.58,0.55,0.8,0.8);
	leg1_P->SetFillColor(0);
	leg1_P->SetFillStyle(4000);
	leg1_P->SetBorderSize(0);
	leg1_P->SetMargin(0.2);
	leg1_P->SetTextSize(0.029);
	leg1_P->SetTextFont(42);
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		hHF_P[j]->SetLineColor(j+1);
		hHF_P[j]->GetXaxis()->SetTitle("HF");
		hHF_P[j]->GetYaxis()->SetTitle("Entries");
		hHF_P[j]->GetYaxis()->SetRangeUser(1e+2,1e+7);
		hHF_P[j]->GetYaxis()->CenterTitle();
		hHF_P[j]->GetXaxis()->CenterTitle();
		hHF_P[j]->GetYaxis()->SetTitleOffset(1.3);
		hHF_P[j]->SetTitle("Passed Condition HF");
		hHF_P[j]->Draw("same");
		leg1_P->AddEntry(hHF_P[j],Form("%s",title.Data()));
		leg1_P->Draw("same");
	}
	cHF1D_P->SaveAs(Form("figs/%s/1D_HF_HiForward_Passed.png",DATE.Data()));

	TCanvas* cHF1D_F = new TCanvas("cHF1D_F","",700,700);
	cHF1D_F->cd();
	cHF1D_F->SetLogy();
	TLegend *leg1_F = new TLegend(0.58,0.55,0.8,0.8);
	leg1_F->SetFillColor(0);
	leg1_F->SetFillStyle(4000);
	leg1_F->SetBorderSize(0);
	leg1_F->SetMargin(0.2);
	leg1_F->SetTextSize(0.029);
	leg1_F->SetTextFont(42);
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		hHF_F[j]->SetLineColor(j+1);
		hHF_F[j]->GetXaxis()->SetTitle("HF");
		hHF_F[j]->GetYaxis()->SetTitle("Entries");
		hHF_F[j]->GetYaxis()->SetRangeUser(1e+2,1e+7);
		hHF_F[j]->GetYaxis()->CenterTitle();
		hHF_F[j]->GetXaxis()->CenterTitle();
		hHF_F[j]->GetYaxis()->SetTitleOffset(1.3);
		hHF_F[j]->SetTitle("Failed Condition HF");
		hHF_F[j]->Draw("same");
		leg1_F->AddEntry(hHF_F[j],Form("%s",title.Data()));
		leg1_F->Draw("same");
	}
	cHF1D_F->SaveAs(Form("figs/%s/1D_HF_HiForward_Falied.png",DATE.Data()));

	TCanvas* cTrk1D_P = new TCanvas("cTrk1D_P","",700,700);
	cTrk1D_P->cd();
	cTrk1D_P->SetLogy();
	TLegend *leg2_P = new TLegend(0.58,0.55,0.8,0.8);
	leg2_P->SetFillColor(0);
	leg2_P->SetFillStyle(4000);
	leg2_P->SetBorderSize(0);
	leg2_P->SetMargin(0.2);
	leg2_P->SetTextSize(0.029);
	leg2_P->SetTextFont(42);
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		hTrk_P[j]->SetLineColor(j+1);
		hTrk_P[j]->GetXaxis()->SetTitle("Num. of Trackers");
		hTrk_P[j]->GetYaxis()->SetTitle("Entries");
		hTrk_P[j]->GetYaxis()->SetRangeUser(1e+2,1e+7);
		hTrk_P[j]->GetYaxis()->CenterTitle();
		hTrk_P[j]->GetXaxis()->CenterTitle();
		hTrk_P[j]->GetYaxis()->SetTitleOffset(1.3);
		hTrk_P[j]->SetTitle("Passed Condition Trk");
		hTrk_P[j]->Draw("same");
		leg2_P->AddEntry(hHF_P[j],Form("%s",title.Data()));
		leg2_P->Draw("same");
	}
	cTrk1D_P->SaveAs(Form("figs/%s/1D_Trk_HiForward_Passed.png",DATE.Data()));

	TCanvas* cTrk1D_F = new TCanvas("cTrk1D_F","",700,700);
	cTrk1D_F->cd();
	cTrk1D_F->SetLogy();
	TLegend *leg2_F = new TLegend(0.58,0.55,0.8,0.8);
	leg2_F->SetFillColor(0);
	leg2_F->SetFillStyle(4000);
	leg2_F->SetBorderSize(0);
	leg2_F->SetMargin(0.2);
	leg2_F->SetTextSize(0.029);
	leg2_F->SetTextFont(42);
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		hTrk_F[j]->SetLineColor(j+1);
		hTrk_F[j]->GetXaxis()->SetTitle("Num. of Trackers");
		hTrk_F[j]->GetYaxis()->SetTitle("Entries");
		hTrk_F[j]->GetYaxis()->SetRangeUser(0.5,1e+5);
		hTrk_F[j]->GetYaxis()->CenterTitle();
		hTrk_F[j]->GetXaxis()->CenterTitle();
		hTrk_F[j]->GetYaxis()->SetTitleOffset(1.3);
		hTrk_F[j]->SetTitle("Failed Condition Trk");
		hTrk_F[j]->Draw("same");
		leg2_F->AddEntry(hHF_F[j],Form("%s",title.Data()));
		leg2_F->Draw("same");
	}
	cTrk1D_F->SaveAs(Form("figs/%s/1D_Trk_HiForward_Failed.png",DATE.Data()));

	TCanvas* cPix1D_P = new TCanvas("cPix1D_P","",700,700);
	cPix1D_P->cd();
	cPix1D_P->SetLogy();
	TLegend *leg3_P = new TLegend(0.58,0.55,0.8,0.8);
	leg3_P->SetFillColor(0);
	leg3_P->SetFillStyle(4000);
	leg3_P->SetBorderSize(0);
	leg3_P->SetMargin(0.2);
	leg3_P->SetTextSize(0.029);
	leg3_P->SetTextFont(42);
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		hPix_P[j]->SetLineColor(j+1);
		hPix_P[j]->GetXaxis()->SetTitle("Num. of Pixels");
		hPix_P[j]->GetYaxis()->SetTitle("Entries");
		hPix_P[j]->GetYaxis()->SetRangeUser(0.5,1e+5);
		hPix_P[j]->GetYaxis()->CenterTitle();
		hPix_P[j]->GetXaxis()->CenterTitle();
		hPix_P[j]->GetYaxis()->SetTitleOffset(1.3);
		hPix_P[j]->SetTitle("Passed Condition Pixel");
		hPix_P[j]->Draw("same");

		leg3_P->AddEntry(hHF_P[j],Form("%s",title.Data()));
		leg3_P->Draw("same");
	}
	cPix1D_P->SaveAs(Form("figs/%s/1D_Pix_HiForward_Passed.png",DATE.Data()));

	TCanvas* cPix1D_F = new TCanvas("cPix1D_F","",700,700);
	cPix1D_F->cd();
	cPix1D_F->SetLogy();
	TLegend *leg3_F = new TLegend(0.58,0.55,0.8,0.8);
	leg3_F->SetFillColor(0);
	leg3_F->SetFillStyle(4000);
	leg3_F->SetBorderSize(0);
	leg3_F->SetMargin(0.2);
	leg3_F->SetTextSize(0.029);
	leg3_F->SetTextFont(42);
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		hPix_F[j]->SetLineColor(j+1);
		hPix_F[j]->GetXaxis()->SetTitle("Num. of Pixels");
		hPix_F[j]->GetYaxis()->SetTitle("Entries");
		hPix_F[j]->GetYaxis()->SetRangeUser(0.5,1e+5);
		hPix_F[j]->GetYaxis()->CenterTitle();
		hPix_F[j]->GetXaxis()->CenterTitle();
		hPix_F[j]->GetYaxis()->SetTitleOffset(1.3);
		hPix_F[j]->SetTitle("Failed Condition Pixel");
		hPix_F[j]->Draw("same");

		leg3_F->AddEntry(hHF_F[j],Form("%s",title.Data()));
		leg3_F->Draw("same");
	}
	cPix1D_F->SaveAs(Form("figs/%s/1D_Pix_HiForward_Failed.png",DATE.Data()));

	TCanvas* c2[9];
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		c2[j] = new TCanvas(Form("c2_%d",j),"",1600,800);
		c2[j]->Divide(2,1);
		c2[j]->cd(1);
		hHF_Pixels_P[j]->GetYaxis()->SetMaxDigits(4);
		hHF_Pixels_P[j]->GetYaxis()->CenterTitle();
		hHF_Pixels_P[j]->GetYaxis()->SetLabelSize(0.03);
		hHF_Pixels_P[j]->GetYaxis()->SetTitleOffset(1.5);
		hHF_Pixels_P[j]->GetXaxis()->SetMaxDigits(4);
		hHF_Pixels_P[j]->GetXaxis()->CenterTitle();
		hHF_Pixels_P[j]->GetXaxis()->SetLabelSize(0.03);
		hHF_Pixels_P[j]->SetTitle(Form("HF vs Number of Pixels Passed Cond. (%s)",title.Data()));
		hHF_Pixels_P[j]->Draw("colz");
		c2[j]->cd(2);
		hHF_Ntracks_P[j]->GetYaxis()->SetMaxDigits(4);
		hHF_Ntracks_P[j]->GetYaxis()->CenterTitle();
		hHF_Ntracks_P[j]->GetYaxis()->SetLabelSize(0.03);
		hHF_Ntracks_P[j]->GetYaxis()->SetTitleOffset(1.5);
		hHF_Ntracks_P[j]->GetXaxis()->SetMaxDigits(4);
		hHF_Ntracks_P[j]->GetXaxis()->CenterTitle();
		hHF_Ntracks_P[j]->GetXaxis()->SetLabelSize(0.03);
		hHF_Ntracks_P[j]->GetXaxis()->SetNdivisions(511);
		hHF_Ntracks_P[j]->SetTitle(Form("HF vs Number of Tracks Passed Cond. (%s)",title.Data()));
		hHF_Ntracks_P[j]->Draw("colz");
		c2[j]->SaveAs(Form("figs/%s/HF_pixel_track_HiForward_Passed_%s.png",DATE.Data(),fName.Data()));
	}

	TCanvas* c1[9];
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		c1[j] = new TCanvas(Form("c1_%d",j),"",800,800);
		c1[j]->cd();
		hPix_Trk_P[j]->GetYaxis()->SetMaxDigits(3);
		hPix_Trk_P[j]->GetYaxis()->CenterTitle();
		hPix_Trk_P[j]->GetYaxis()->SetLabelSize(0.03);
		hPix_Trk_P[j]->GetYaxis()->SetTitleOffset(1.5);
		hPix_Trk_P[j]->GetXaxis()->SetMaxDigits(3);
		hPix_Trk_P[j]->GetXaxis()->CenterTitle();
		hPix_Trk_P[j]->GetXaxis()->SetLabelSize(0.03);
		hPix_Trk_P[j]->SetTitle(Form("Pixel vs Tracks Passed Cond. (%s)",title.Data()));
		hPix_Trk_P[j]->Draw("colz");
		c1[j]->SaveAs(Form("figs/%s/Npix_Ntrk_HiForward_Passed_%s.png",DATE.Data(),fName.Data()));
	}

	TCanvas* c3[9];
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		c3[j]= new TCanvas(Form("c3_%d",j),"",1600,800);
		c3[j]->Divide(2,1);
		c3[j]->cd(1);
		hHF_Pixels_F[j]->GetYaxis()->SetMaxDigits(4);
		hHF_Pixels_F[j]->GetYaxis()->CenterTitle();
		hHF_Pixels_F[j]->GetYaxis()->SetLabelSize(0.03);
		hHF_Pixels_F[j]->GetYaxis()->SetTitleOffset(1.5);
		hHF_Pixels_F[j]->GetXaxis()->SetMaxDigits(4);
		hHF_Pixels_F[j]->GetXaxis()->CenterTitle();
		hHF_Pixels_F[j]->GetXaxis()->SetLabelSize(0.03);
		hHF_Pixels_F[j]->SetTitle(Form("HF vs Number of Pixels Failed Cond.(%s)",title.Data()));
		hHF_Pixels_F[j]->Draw("colz");
		c3[j]->cd(2);
		hHF_Ntracks_F[j]->GetYaxis()->SetMaxDigits(4);
		hHF_Ntracks_F[j]->GetYaxis()->CenterTitle();
		hHF_Ntracks_F[j]->GetYaxis()->SetLabelSize(0.03);
		hHF_Ntracks_F[j]->GetYaxis()->SetTitleOffset(1.5);
		hHF_Ntracks_F[j]->GetXaxis()->SetMaxDigits(4);
		hHF_Ntracks_F[j]->GetXaxis()->CenterTitle();
		hHF_Ntracks_F[j]->GetXaxis()->SetLabelSize(0.03);
		hHF_Ntracks_F[j]->GetXaxis()->SetNdivisions(511);
		hHF_Ntracks_F[j]->SetTitle(Form("HF vs Number of Tracks Failed Cond. (%s)",title.Data()));
		hHF_Ntracks_F[j]->Draw("colz");
		c3[j]->SaveAs(Form("figs/%s/HF_pixel_track_HiForward_Failed_%s.png",DATE.Data(),fName.Data()));
	}

	TCanvas* c4[9];
	for(int j=0;j<9;j++)
	{
		TString fName=filter[j];
		TString title=s+fName;
		c4[j]= new TCanvas(Form("c4_%d",j),"",800,800);
		c4[j]->cd();
		hPix_Trk_F[j]->GetYaxis()->SetMaxDigits(3);
		hPix_Trk_F[j]->GetYaxis()->CenterTitle();
		hPix_Trk_F[j]->GetYaxis()->SetLabelSize(0.03);
		hPix_Trk_F[j]->GetYaxis()->SetTitleOffset(1.5);
		hPix_Trk_F[j]->GetXaxis()->SetMaxDigits(3);
		hPix_Trk_F[j]->GetXaxis()->CenterTitle();
		hPix_Trk_F[j]->GetXaxis()->SetLabelSize(0.03);
		hPix_Trk_F[j]->SetTitle(Form("Pixel vs Tracks Failed Cond. (%s)",title.Data()));
		hPix_Trk_F[j]->Draw("colz");
		c4[j]->SaveAs(Form("figs/%s/Npix_Ntrk_HiForward_Failed_%s.png",DATE.Data(),fName.Data()));
	}

	t -> Stop();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());
}
