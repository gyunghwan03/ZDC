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

void fill_hist_HiForward_0000(int nevt=-1) 
{

	using namespace std;

	TStopwatch *t = new TStopwatch;
	t -> Start();

	TString DATE="210609";
	gSystem->mkdir(Form("figs/%s", DATE.Data()),kTRUE);

	gStyle->SetOptStat(00000);
	gROOT->ForceStyle();
	TString fin = "./roots/HiForward_0000.root";
	TChain *evttree = new TChain("evttree");
	TChain *hlttree = new TChain("hlttree");
	TChain *skimtree = new TChain("skimtree");

	TFile *outFile = new TFile("./roots/hist_HiForward_0000.root","recreate");

	evttree->Add(fin.Data());
	hlttree->Add(fin.Data());
	skimtree->Add(fin.Data());

	evttree->AddFriend(hlttree);
	evttree->AddFriend(skimtree);

	//### define using variables

	int noareas=20;

	int nPixBins_P, nZDCBins_P, nTrkBins_P, nHFBins_P;
	int nPixBins_F, nZDCBins_F, nTrkBins_F, nHFBins_F;
	double npix_max_P, zdc_max_P, ntrk_max_P, hf_max_P, npix_min_P, zdc_min_P, ntrk_min_P, hf_min_P;
	nPixBins_P=100;
	nZDCBins_P=100;
	nTrkBins_P=100;
	nHFBins_P =100;
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
	nZDCBins_F=100;
	nTrkBins_F=100;
	nHFBins_F =100;
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
	int phfCoincFilterTh1, phfCoincFilterTh2, phfCoincFilterTh3, phfCoincFilterTh4, phfCoincFilterTh5, phfCoincFilterTh1p5, phfCoincFilterTh2p5, phfCoincFilterTh3p5, phfCoincFilterTh4p5;
	int phfCoincFilter3Th1, phfCoincFilter3Th2, phfCoincFilter3Th3, phfCoincFilter3Th4, phfCoincFilter3Th5, phfCoincFilter3Th1p5, phfCoincFilter3Th2p5, phfCoincFilter3Th3p5, phfCoincFilter3Th4p5;
	int phfCoincFilter4Th1, phfCoincFilter4Th2, phfCoincFilter4Th3, phfCoincFilter4Th4, phfCoincFilter4Th5, phfCoincFilter4Th1p5, phfCoincFilter4Th2p5, phfCoincFilter4Th3p5, phfCoincFilter4Th4p5;
	int phfCoincFilter5Th1, phfCoincFilter5Th2, phfCoincFilter5Th3, phfCoincFilter5Th4, phfCoincFilter5Th5, phfCoincFilter5Th1p5, phfCoincFilter5Th2p5, phfCoincFilter5Th3p5, phfCoincFilter5Th4p5;
	int CoincFilterTh1, CoincFilterTh2, CoincFilterTh3, CoincFilterTh4, CoincFilterTh5, CoincFilterTh1p5, CoincFilterTh2p5, CoincFilterTh3p5, CoincFilterTh4p5;
	int CoincFilter2Th1, CoincFilter2Th2, CoincFilter2Th3, CoincFilter2Th4, CoincFilter2Th5, CoincFilter2Th1p5, CoincFilter2Th2p5, CoincFilter2Th3p5, CoincFilter2Th4p5;
	int CoincFilter3Th1, CoincFilter3Th2, CoincFilter3Th3, CoincFilter3Th4, CoincFilter3Th5, CoincFilter3Th1p5, CoincFilter3Th2p5, CoincFilter3Th3p5, CoincFilter3Th4p5;
	int CoincFilter4Th1, CoincFilter4Th2, CoincFilter4Th3, CoincFilter4Th4, CoincFilter4Th5, CoincFilter4Th1p5, CoincFilter4Th2p5, CoincFilter4Th3p5, CoincFilter4Th4p5;
	int CoincFilter5Th1, CoincFilter5Th2, CoincFilter5Th3, CoincFilter5Th4, CoincFilter5Th5, CoincFilter5Th1p5, CoincFilter5Th2p5, CoincFilter5Th3p5, CoincFilter5Th4p5;

	evttree->SetBranchAddress("evt",&evt);
	evttree->SetBranchAddress("hiNpix",&hiNpix);
	evttree->SetBranchAddress("hiNtracks",&hiNtracks);
	evttree->SetBranchAddress("hiBin",&hiBin);
	evttree->SetBranchAddress("hiHF",&hiHF);

	skimtree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	skimtree->SetBranchAddress("phfCoincFilterTh1"  ,&phfCoincFilterTh1);
	skimtree->SetBranchAddress("phfCoincFilterTh2"  ,&phfCoincFilterTh2);
	skimtree->SetBranchAddress("phfCoincFilterTh3"  ,&phfCoincFilterTh3);
	skimtree->SetBranchAddress("phfCoincFilterTh4"  ,&phfCoincFilterTh4);
	skimtree->SetBranchAddress("phfCoincFilterTh5"  ,&phfCoincFilterTh5);
	skimtree->SetBranchAddress("phfCoincFilterTh1p5",&phfCoincFilterTh1p5);
	skimtree->SetBranchAddress("phfCoincFilterTh2p5",&phfCoincFilterTh2p5);
	skimtree->SetBranchAddress("phfCoincFilterTh3p5",&phfCoincFilterTh3p5);
	skimtree->SetBranchAddress("phfCoincFilterTh4p5",&phfCoincFilterTh4p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th1"  ,&phfCoincFilter2Th1);
	skimtree->SetBranchAddress("phfCoincFilter2Th2"  ,&phfCoincFilter2Th2);
	skimtree->SetBranchAddress("phfCoincFilter2Th3"  ,&phfCoincFilter2Th3);
	skimtree->SetBranchAddress("phfCoincFilter2Th4"  ,&phfCoincFilter2Th4);
	skimtree->SetBranchAddress("phfCoincFilter2Th5"  ,&phfCoincFilter2Th5);
	skimtree->SetBranchAddress("phfCoincFilter2Th1p5",&phfCoincFilter2Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th2p5",&phfCoincFilter2Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th3p5",&phfCoincFilter2Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th4p5",&phfCoincFilter2Th4p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th1"  ,&phfCoincFilter3Th1);
	skimtree->SetBranchAddress("phfCoincFilter3Th2"  ,&phfCoincFilter3Th2);
	skimtree->SetBranchAddress("phfCoincFilter3Th3"  ,&phfCoincFilter3Th3);
	skimtree->SetBranchAddress("phfCoincFilter3Th4"  ,&phfCoincFilter3Th4);
	skimtree->SetBranchAddress("phfCoincFilter3Th5"  ,&phfCoincFilter3Th5);
	skimtree->SetBranchAddress("phfCoincFilter3Th1p5",&phfCoincFilter3Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th2p5",&phfCoincFilter3Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th3p5",&phfCoincFilter3Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th4p5",&phfCoincFilter3Th4p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th1"  ,&phfCoincFilter4Th1);
	skimtree->SetBranchAddress("phfCoincFilter4Th2"  ,&phfCoincFilter4Th2);
	skimtree->SetBranchAddress("phfCoincFilter4Th3"  ,&phfCoincFilter4Th3);
	skimtree->SetBranchAddress("phfCoincFilter4Th4"  ,&phfCoincFilter4Th4);
	skimtree->SetBranchAddress("phfCoincFilter4Th5"  ,&phfCoincFilter4Th5);
	skimtree->SetBranchAddress("phfCoincFilter4Th1p5",&phfCoincFilter4Th1p5);
	//skimtree->SetBranchAddress("phfCoincFilter4Th2p5",&phfCoincFilter4Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th3p5",&phfCoincFilter4Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th4p5",&phfCoincFilter4Th4p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th1"  ,&phfCoincFilter5Th1);
	skimtree->SetBranchAddress("phfCoincFilter5Th2"  ,&phfCoincFilter5Th2);
	skimtree->SetBranchAddress("phfCoincFilter5Th3"  ,&phfCoincFilter5Th3);
	skimtree->SetBranchAddress("phfCoincFilter5Th4"  ,&phfCoincFilter5Th4);
	skimtree->SetBranchAddress("phfCoincFilter5Th5"  ,&phfCoincFilter5Th5);
	skimtree->SetBranchAddress("phfCoincFilter5Th1p5",&phfCoincFilter5Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th2p5",&phfCoincFilter5Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th3p5",&phfCoincFilter5Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th4p5",&phfCoincFilter5Th4p5);
	skimtree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);

	hlttree->SetBranchAddress("HLT_HIZeroBias_v1", &HLT_HIZeroBias_v1);

	int f;
	//cout << "fName : " << fName.Data() << endl;

	// ThX
	TH2D* hHF_Pixels_P[5][9]; 
	TH2D* hHF_Ntracks_P[5][9];
	TH2D* hPix_Trk_P[5][9];

	TH2D* hHF_Pixels_F[5][9]; 
	TH2D* hHF_Ntracks_F[5][9];
	TH2D* hPix_Trk_F[5][9];

	TH1F *hHF_P[5][9]; 
	TH1F *hTrk_P[5][9];
	TH1F *hPix_P[5][9];
	TH1F *hHF_F[5][9];
	TH1F *hTrk_F[5][9];
	TH1F *hPix_F[5][9];

	for(int i=1;i<6;i++){
		for(int j=0;j<9;j++){
			// ThX
			hHF_Pixels_P[i][j]  = new TH2D(Form("hHF_Pixels_P_%dTh_%d",i,j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
			hHF_Ntracks_P[i][j] = new TH2D(Form("hHF_Tracks_P_%dTh_%d",i,j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
			hPix_Trk_P[i][j]    = new TH2D(Form("hPix_Trk_P_%dTh_%d",i,j),";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

			hHF_Pixels_F[i][j]  = new TH2D(Form("hHF_Pixels_F_%dTh_%d",i,j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
			hHF_Ntracks_F[i][j] = new TH2D(Form("hHF_Tracks_F_%dTh_%d",i,j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
			hPix_Trk_F[i][j]    = new TH2D(Form("hPix_Trk_F_%dTh_%d",i,j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);

			hHF_P[i][j]  = new TH1F(Form("hHF_P_%dTh_%d",i,j),"",100,0,300);
			hHF_F[i][j]  = new TH1F(Form("hHF_F_%dTh_%d",i,j),"",100,0,300);
			hTrk_P[i][j] = new TH1F(Form("hTrk_P_%dTh_%d",i,j),"",100,0,300);
			hPix_P[i][j] = new TH1F(Form("hPix_P_%dTh_%d",i,j),"",100,0,300);
			hTrk_F[i][j] = new TH1F(Form("hTrk_F_%dTh_%d",i,j),"",100,0,300);
			hPix_F[i][j] = new TH1F(Form("hPix_F_%dTh_%d",i,j),"",100,0,300);
		}
	}


	if (nevt==-1) nevt = evttree->GetEntries();

	///// Start Loop
	for (int iev=0;iev<nevt;iev++) {
		if(iev%10000000==0) cout << ">>>>> EVENT " << iev << " / " << evttree->GetEntries() <<  " ("<<(int)(100.*iev/evttree->GetEntries()) << "%)" << endl;
		evttree->GetEntry(iev);
		hlttree->GetEntry(iev);
		skimtree->GetEntry(iev);
		int CoincFilter[5][9] = {
			{phfCoincFilterTh1, phfCoincFilterTh2, phfCoincFilterTh3, phfCoincFilterTh4, phfCoincFilterTh5, phfCoincFilterTh1p5, phfCoincFilterTh2p5, phfCoincFilterTh3p5, phfCoincFilterTh4p5},
			{phfCoincFilter2Th1, phfCoincFilter2Th2, phfCoincFilter2Th3, phfCoincFilter2Th4, phfCoincFilter2Th5, phfCoincFilter2Th1p5, phfCoincFilter2Th2p5, phfCoincFilter2Th3p5, phfCoincFilter2Th4p5},
			{phfCoincFilter3Th1, phfCoincFilter3Th2, phfCoincFilter3Th3, phfCoincFilter3Th4, phfCoincFilter3Th5, phfCoincFilter3Th1p5, phfCoincFilter3Th2p5, phfCoincFilter3Th3p5, phfCoincFilter3Th4p5},
			{phfCoincFilter4Th1, phfCoincFilter4Th2, phfCoincFilter4Th3, phfCoincFilter4Th4, phfCoincFilter4Th5, phfCoincFilter4Th1p5, phfCoincFilter4Th2p5, phfCoincFilter4Th3p5, phfCoincFilter4Th4p5},
			{phfCoincFilter5Th1, phfCoincFilter5Th2, phfCoincFilter5Th3, phfCoincFilter5Th4, phfCoincFilter5Th5, phfCoincFilter5Th1p5, phfCoincFilter5Th2p5, phfCoincFilter5Th3p5, phfCoincFilter5Th4p5}
		};
		for (int i=1;i<6;i++){
			for (int j=0;j<9;j++){
				/*
				if(i==1) {
					if (j==0)CoincFilter[i][j]=phfCoincFilterTh1; else if(j==1) CoincFilter[i][j]=phfCoincFilterTh2;	else if(j==2) CoincFilter[i][j]=phfCoincFilterTh3; else if(j==3) CoincFilter[i][j]=phfCoincFilterTh4; else if(j==4) CoincFilter[i][j]=phfCoincFilterTh5;
					else if(j==5)CoincFilter[i][j]=phfCoincFilterTh1p5; else if(j==6) CoincFilter[i][j]=phfCoincFilterTh2p5; else if(j==7) CoincFilter[i][j]=phfCoincFilterTh3p5; else if(j==8) CoincFilter[i][j]=phfCoincFilterTh4p5; else if(j==9) CoincFilter[i][j]=phfCoincFilterTh4p5;
				}
				if(i==2) {
					if (j==0)CoincFilter[i][j]=phfCoincFilter2Th1; else if(j==1) CoincFilter[i][j]=phfCoincFilter2Th2;	else if(j==2) CoincFilter[i][j]=phfCoincFilter2Th3; else if(j==3) CoincFilter[i][j]=phfCoincFilter2Th4; else if(j==4) CoincFilter[i][j]=phfCoincFilter2Th5;
					else if(j==5)CoincFilter[i][j]=phfCoincFilter2Th1p5; else if(j==6) CoincFilter[i][j]=phfCoincFilter2Th2p5; else if(j==7) CoincFilter[i][j]=phfCoincFilter2Th3p5; else if(j==8) CoincFilter[i][j]=phfCoincFilter2Th4p5; else if(j==9) CoincFilter[i][j]=phfCoincFilter2Th4p5;
				}
				if(i==3) {
					if (j==0)CoincFilter[i][j]=phfCoincFilter3Th1; else if(j==1) CoincFilter[i][j]=phfCoincFilter3Th2;	else if(j==2) CoincFilter[i][j]=phfCoincFilter3Th3; else if(j==3) CoincFilter[i][j]=phfCoincFilter3Th4; else if(j==4) CoincFilter[i][j]=phfCoincFilter3Th5;
					else if(j==5)CoincFilter[i][j]=phfCoincFilter3Th1p5; else if(j==6) CoincFilter[i][j]=phfCoincFilter3Th2p5; else if(j==7) CoincFilter[i][j]=phfCoincFilter3Th3p5; else if(j==8) CoincFilter[i][j]=phfCoincFilter3Th4p5; else if(j==9) CoincFilter[i][j]=phfCoincFilter3Th4p5;
				}
				if(i==4) {
					if (j==0)CoincFilter[i][j]=phfCoincFilter4Th1; else if(j==1) CoincFilter[i][j]=phfCoincFilter4Th2;	else if(j==2) CoincFilter[i][j]=phfCoincFilter4Th3; else if(j==3) CoincFilter[i][j]=phfCoincFilter4Th4; else if(j==4) CoincFilter[i][j]=phfCoincFilter4Th5;
					else if(j==5)CoincFilter[i][j]=phfCoincFilter4Th1p5; else if(j==6) CoincFilter[i][j]=phfCoincFilter4Th2p5; else if(j==7) CoincFilter[i][j]=phfCoincFilter4Th3p5; else if(j==8) CoincFilter[i][j]=phfCoincFilter4Th4p5; else if(j==9) CoincFilter[i][j]=phfCoincFilter4Th4p5;
				}
				if(i==5) {
					if (j==0)CoincFilter[i][j]=phfCoincFilter5Th1; else if(j==1) CoincFilter[i][j]=phfCoincFilter5Th2;	else if(j==2) CoincFilter[i][j]=phfCoincFilter5Th3; else if(j==3) CoincFilter[i][j]=phfCoincFilter5Th4; else if(j==4) CoincFilter[i][j]=phfCoincFilter5Th5;
					else if(j==5)CoincFilter[i][j]=phfCoincFilter5Th1p5; else if(j==6) CoincFilter[i][j]=phfCoincFilter5Th2p5; else if(j==7) CoincFilter[i][j]=phfCoincFilter5Th3p5; else if(j==8) CoincFilter[i][j]=phfCoincFilter5Th4p5; else if(j==9) CoincFilter[i][j]=phfCoincFilter5Th4p5;
				}*/
				if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter[i][j])==1 ){  
					hHF_Pixels_P[i][j]->Fill(hiNpix,hiHF);
					hHF_Ntracks_P[i][j]->Fill(hiNtracks,hiHF);
					hPix_Trk_P[i][j]->Fill(hiNpix,hiNtracks);
				}
				if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter[i][j])==0 ){  
					hHF_Pixels_F[i][j]->Fill(hiNpix,hiHF);
					hHF_Ntracks_F[i][j]->Fill(hiNtracks,hiHF);
					hPix_Trk_F[i][j]->Fill(hiNpix,hiNtracks);
				}
				if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter[i][j]==1) ){
					hHF_P[i][j]->Fill(hiHF);
					hTrk_P[i][j]->Fill(hiNtracks);
					hPix_P[i][j]->Fill(hiNpix);
				}
				if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter[i][j])==0 ){
					hHF_F[i][j]->Fill(hiHF);
					hTrk_F[i][j]->Fill(hiNtracks);
					hPix_F[i][j]->Fill(hiNpix);
				}
			}
		}
	}//// End Loop

	outFile->cd();
	//// Write 
	for (int i=1;i<6;i++){
		for (int j=0;j<9;j++)
		{
			//ThX
			hHF_Pixels_P[i][j]->Write(); 
			hHF_Ntracks_P[i][j]->Write();
			hPix_Trk_P[i][j]->Write();

			hHF_Pixels_F[i][j]->Write(); 
			hHF_Ntracks_F[i][j]->Write();
			hPix_Trk_F[i][j]->Write();

			hHF_P[i][j]->Write(); 
			hTrk_P[i][j]->Write();
			hPix_P[i][j]->Write();
			hHF_F[i][j]->Write();
			hTrk_F[i][j]->Write();
			hPix_F[i][j]->Write();

		}
	}
	outFile->Close();
	t -> Stop();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());
}
