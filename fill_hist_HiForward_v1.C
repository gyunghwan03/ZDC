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
	int CoincFilter[9]={CoincFilterTh1, CoincFilterTh2, CoincFilterTh3, CoincFilterTh4, CoincFilterTh5, CoincFilterTh1p5, CoincFilterTh2p5, CoincFilterTh3p5, CoincFilterTh4p5};
	int CoincFilter2Th1, CoincFilter2Th2, CoincFilter2Th3, CoincFilter2Th4, CoincFilter2Th5, CoincFilter2Th1p5, CoincFilter2Th2p5, CoincFilter2Th3p5, CoincFilter2Th4p5;
	int CoincFilter2[9]={CoincFilter2Th1, CoincFilter2Th2, CoincFilter2Th3, CoincFilter2Th4, CoincFilter2Th5, CoincFilter2Th1p5, CoincFilter2Th2p5, CoincFilter2Th3p5, CoincFilter2Th4p5};
	int CoincFilter3Th1, CoincFilter3Th2, CoincFilter3Th3, CoincFilter3Th4, CoincFilter3Th5, CoincFilter3Th1p5, CoincFilter3Th2p5, CoincFilter3Th3p5, CoincFilter3Th4p5;
	int CoincFilter3[9]={CoincFilter3Th1, CoincFilter3Th2, CoincFilter3Th3, CoincFilter3Th4, CoincFilter3Th5, CoincFilter3Th1p5, CoincFilter3Th2p5, CoincFilter3Th3p5, CoincFilter3Th4p5};
	int CoincFilter4Th1, CoincFilter4Th2, CoincFilter4Th3, CoincFilter4Th4, CoincFilter4Th5, CoincFilter4Th1p5, CoincFilter4Th2p5, CoincFilter4Th3p5, CoincFilter4Th4p5;
	int CoincFilter4[9]={CoincFilter4Th1, CoincFilter4Th2, CoincFilter4Th3, CoincFilter4Th4, CoincFilter4Th5, CoincFilter4Th1p5, CoincFilter4Th2p5, CoincFilter4Th3p5, CoincFilter4Th4p5};
	int CoincFilter5Th1, CoincFilter5Th2, CoincFilter5Th3, CoincFilter5Th4, CoincFilter5Th5, CoincFilter5Th1p5, CoincFilter5Th2p5, CoincFilter5Th3p5, CoincFilter5Th4p5;
	int CoincFilter5[9]={CoincFilter5Th1, CoincFilter5Th2, CoincFilter5Th3, CoincFilter5Th4, CoincFilter5Th5, CoincFilter5Th1p5, CoincFilter5Th2p5, CoincFilter5Th3p5, CoincFilter5Th4p5};
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

	// 2ThX
	TH2D* hHF_Pixels_P_2Th[9]; 
	TH2D* hHF_Ntracks_P_2Th[9];
	TH2D* hPix_Trk_P_2Th[9];

	TH2D* hHF_Pixels_F_2Th[9]; 
	TH2D* hHF_Ntracks_F_2Th[9];
	TH2D* hPix_Trk_F_2Th[9];

	TH1F *hHF_P_2Th[9]; 
	TH1F *hTrk_P_2Th[9];
	TH1F *hPix_P_2Th[9];
	TH1F *hHF_F_2Th[9];
	TH1F *hTrk_F_2Th[9];
	TH1F *hPix_F_2Th[9];

	// 3ThX
	TH2D* hHF_Pixels_P_3Th[9]; 
	TH2D* hHF_Ntracks_P_3Th[9];
	TH2D* hPix_Trk_P_3Th[9];

	TH2D* hHF_Pixels_F_3Th[9]; 
	TH2D* hHF_Ntracks_F_3Th[9];
	TH2D* hPix_Trk_F_3Th[9];

	TH1F *hHF_P_3Th[9]; 
	TH1F *hTrk_P_3Th[9];
	TH1F *hPix_P_3Th[9];
	TH1F *hHF_F_3Th[9];
	TH1F *hTrk_F_3Th[9];
	TH1F *hPix_F_3Th[9];

	// 4ThX
	TH2D* hHF_Pixels_P_4Th[9]; 
	TH2D* hHF_Ntracks_P_4Th[9];
	TH2D* hPix_Trk_P_4Th[9];

	TH2D* hHF_Pixels_F_4Th[9]; 
	TH2D* hHF_Ntracks_F_4Th[9];
	TH2D* hPix_Trk_F_4Th[9];

	TH1F *hHF_P_4Th[9]; 
	TH1F *hTrk_P_4Th[9];
	TH1F *hPix_P_4Th[9];
	TH1F *hHF_F_4Th[9];
	TH1F *hTrk_F_4Th[9];
	TH1F *hPix_F_4Th[9];

	// 5ThX
	TH2D* hHF_Pixels_P_5Th[9]; 
	TH2D* hHF_Ntracks_P_5Th[9];
	TH2D* hPix_Trk_P_5Th[9];

	TH2D* hHF_Pixels_F_5Th[9]; 
	TH2D* hHF_Ntracks_F_5Th[9];
	TH2D* hPix_Trk_F_5Th[9];

	TH1F *hHF_P_5Th[9]; 
	TH1F *hTrk_P_5Th[9];
	TH1F *hPix_P_5Th[9];
	TH1F *hHF_F_5Th[9];
	TH1F *hTrk_F_5Th[9];
	TH1F *hPix_F_5Th[9];

	for(int j=0;j<9;j++){
		// ThX
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
		// 2ThX
		hHF_Pixels_P_2Th[j]  = new TH2D(Form("hHF_Pixels_P_2Th%d",j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hHF_Ntracks_P_2Th[j] = new TH2D(Form("hHF_Tracks_P_2Th%d",j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hPix_Trk_P_2Th[j]    = new TH2D(Form("hPix_Trk_P_2Th%d",j),";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

		hHF_Pixels_F_2Th[j]  = new TH2D(Form("hHF_Pixels_F_2Th%d",j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hHF_Ntracks_F_2Th[j] = new TH2D(Form("hHF_Tracks_F_2Th%d",j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hPix_Trk_F_2Th[j]    = new TH2D(Form("hPix_Trk_F_2Th%d",j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);

		hHF_P_2Th[j]  = new TH1F(Form("hHF_P_2Th%d",j),"",100,0,300);
		hHF_F_2Th[j]  = new TH1F(Form("hHF_F_2Th%d",j),"",100,0,300);
		hTrk_P_2Th[j] = new TH1F(Form("hTrk_P_2Th%d",j),"",100,0,300);
		hPix_P_2Th[j] = new TH1F(Form("hPix_P_2Th%d",j),"",100,0,300);
		hTrk_F_2Th[j] = new TH1F(Form("hTrk_F_2Th%d",j),"",100,0,300);
		hPix_F_2Th[j] = new TH1F(Form("hPix_F_2Th%d",j),"",100,0,300);
		// 3ThX
		hHF_Pixels_P_3Th[j]  = new TH2D(Form("hHF_Pixels_P_3Th%d",j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hHF_Ntracks_P_3Th[j] = new TH2D(Form("hHF_Tracks_P_3Th%d",j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hPix_Trk_P_3Th[j]    = new TH2D(Form("hPix_Trk_P_3Th%d",j),";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

		hHF_Pixels_F_3Th[j]  = new TH2D(Form("hHF_Pixels_F_3Th%d",j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hHF_Ntracks_F_3Th[j] = new TH2D(Form("hHF_Tracks_F_3Th%d",j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hPix_Trk_F_3Th[j]    = new TH2D(Form("hPix_Trk_F_3Th%d",j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);

		hHF_P_3Th[j]  = new TH1F(Form("hHF_P_3Th%d",j),"",100,0,300);
		hHF_F_3Th[j]  = new TH1F(Form("hHF_F_3Th%d",j),"",100,0,300);
		hTrk_P_3Th[j] = new TH1F(Form("hTrk_P_3Th%d",j),"",100,0,300);
		hPix_P_3Th[j] = new TH1F(Form("hPix_P_3Th%d",j),"",100,0,300);
		hTrk_F_3Th[j] = new TH1F(Form("hTrk_F_3Th%d",j),"",100,0,300);
		hPix_F_3Th[j] = new TH1F(Form("hPix_F_3Th%d",j),"",100,0,300);
		// 4ThX
		hHF_Pixels_P_4Th[j]  = new TH2D(Form("hHF_Pixels_P_4Th%d",j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hHF_Ntracks_P_4Th[j] = new TH2D(Form("hHF_Tracks_P_4Th%d",j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hPix_Trk_P_4Th[j]    = new TH2D(Form("hPix_Trk_P_4Th%d",j),";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

		hHF_Pixels_F_4Th[j]  = new TH2D(Form("hHF_Pixels_F_4Th%d",j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hHF_Ntracks_F_4Th[j] = new TH2D(Form("hHF_Tracks_F_4Th%d",j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hPix_Trk_F_4Th[j]    = new TH2D(Form("hPix_Trk_F_4Th%d",j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);

		hHF_P_4Th[j]  = new TH1F(Form("hHF_P_4Th%d",j),"",100,0,300);
		hHF_F_4Th[j]  = new TH1F(Form("hHF_F_4Th%d",j),"",100,0,300);
		hTrk_P_4Th[j] = new TH1F(Form("hTrk_P_4Th%d",j),"",100,0,300);
		hPix_P_4Th[j] = new TH1F(Form("hPix_P_4Th%d",j),"",100,0,300);
		hTrk_F_4Th[j] = new TH1F(Form("hTrk_F_4Th%d",j),"",100,0,300);
		hPix_F_4Th[j] = new TH1F(Form("hPix_F_4Th%d",j),"",100,0,300);
		// 5ThX
		hHF_Pixels_P_5Th[j]  = new TH2D(Form("hHF_Pixels_P_5Th%d",j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hHF_Ntracks_P_5Th[j] = new TH2D(Form("hHF_Tracks_P_5Th%d",j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,nHFBins_P,hf_min_P,hf_max_P);
		hPix_Trk_P_5Th[j]    = new TH2D(Form("hPix_Trk_P_5Th%d",j),";No. of pixels;No. of tracks",nPixBins_P,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

		hHF_Pixels_F_5Th[j]  = new TH2D(Form("hHF_Pixels_F_5Th%d",j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hHF_Ntracks_F_5Th[j] = new TH2D(Form("hHF_Tracks_F_5Th%d",j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
		hPix_Trk_F_5Th[j]    = new TH2D(Form("hPix_Trk_F_5Th%d",j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);

		hHF_P_5Th[j]  = new TH1F(Form("hHF_P_5Th%d",j),"",100,0,300);
		hHF_F_5Th[j]  = new TH1F(Form("hHF_F_5Th%d",j),"",100,0,300);
		hTrk_P_5Th[j] = new TH1F(Form("hTrk_P_5Th%d",j),"",100,0,300);
		hPix_P_5Th[j] = new TH1F(Form("hPix_P_5Th%d",j),"",100,0,300);
		hTrk_F_5Th[j] = new TH1F(Form("hTrk_F_5Th%d",j),"",100,0,300);
		hPix_F_5Th[j] = new TH1F(Form("hPix_F_5Th%d",j),"",100,0,300);
	}


	if (nevt==-1) nevt = evttree->GetEntries();

	///// Start Loop
	for (int iev=0;iev<nevt;iev++) {
		if(iev%10000000==0) cout << ">>>>> EVENT " << iev << " / " << evttree->GetEntries() <<  " ("<<(int)(100.*iev/evttree->GetEntries()) << "%)" << endl;
		evttree->GetEntry(iev);
		hlttree->GetEntry(iev);
		skimtree->GetEntry(iev);
		for (int j=0;j<9;j++){
			//ThX
			if (j==0)CoincFilter[j]=phfCoincFilterTh1; else if(j==1) CoincFilter[j]=phfCoincFilterTh2;	else if(j==2) CoincFilter[j]=phfCoincFilterTh3; else if(j==3) CoincFilter[j]=phfCoincFilterTh4; else if(j==4) CoincFilter[j]=phfCoincFilterTh5;
			else if(j==5)CoincFilter[j]=phfCoincFilterTh1p5; else if(j==6) CoincFilter[j]=phfCoincFilterTh2p5; else if(j==7) CoincFilter[j]=phfCoincFilterTh3p5; else if(j==8) CoincFilter[j]=phfCoincFilterTh4p5; else if(j==9) CoincFilter[j]=phfCoincFilterTh4p5;
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter[j])==1 ){  
				hHF_Pixels_P[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_P[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_P[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter[j])==0 ){  
				hHF_Pixels_F[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_F[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_F[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter[j]==1) ){
				hHF_P[j]->Fill(hiHF);
				hTrk_P[j]->Fill(hiNtracks);
				hPix_P[j]->Fill(hiNpix);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter[j])==0 ){
				hHF_F[j]->Fill(hiHF);
				hTrk_F[j]->Fill(hiNtracks);
				hPix_F[j]->Fill(hiNpix);
			}
			//2ThX
			if (j==0)CoincFilter2[j]=phfCoincFilter2Th1; else if(j==1) CoincFilter2[j]=phfCoincFilter2Th2;	else if(j==2) CoincFilter2[j]=phfCoincFilter2Th3; else if(j==3) CoincFilter2[j]=phfCoincFilter2Th4; else if(j==4) CoincFilter2[j]=phfCoincFilter2Th5;
			else if(j==5)CoincFilter2[j]=phfCoincFilter2Th1p5; else if(j==6) CoincFilter2[j]=phfCoincFilter2Th2p5; else if(j==7) CoincFilter2[j]=phfCoincFilter2Th3p5; else if(j==8) CoincFilter2[j]=phfCoincFilter2Th4p5; else if(j==9) CoincFilter2[j]=phfCoincFilter2Th4p5;
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[j])==1 ){  
				hHF_Pixels_P_2Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_P_2Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_P_2Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[j])==0 ){  
				hHF_Pixels_F_2Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_F_2Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_F_2Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter2[j]==1) ){
				hHF_P_2Th[j]->Fill(hiHF);
				hTrk_P_2Th[j]->Fill(hiNtracks);
				hPix_P_2Th[j]->Fill(hiNpix);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter2[j])==0 ){
				hHF_F_2Th[j]->Fill(hiHF);
				hTrk_F_2Th[j]->Fill(hiNtracks);
				hPix_F_2Th[j]->Fill(hiNpix);
			}
			//3ThX
			if (j==0)CoincFilter3[j]=phfCoincFilter3Th1; else if(j==1) CoincFilter3[j]=phfCoincFilter3Th2;	else if(j==2) CoincFilter3[j]=phfCoincFilter3Th3; else if(j==3) CoincFilter3[j]=phfCoincFilter3Th4; else if(j==4) CoincFilter3[j]=phfCoincFilter3Th5;
			else if(j==5)CoincFilter3[j]=phfCoincFilter3Th1p5; else if(j==6) CoincFilter3[j]=phfCoincFilter3Th2p5; else if(j==7) CoincFilter3[j]=phfCoincFilter3Th3p5; else if(j==8) CoincFilter3[j]=phfCoincFilter3Th4p5; else if(j==9) CoincFilter3[j]=phfCoincFilter3Th4p5;
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter3[j])==1 ){  
				hHF_Pixels_P_3Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_P_3Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_P_3Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter3[j])==0 ){  
				hHF_Pixels_F_3Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_F_3Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_F_3Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter3[j]==1) ){
				hHF_P_3Th[j]->Fill(hiHF);
				hTrk_P_3Th[j]->Fill(hiNtracks);
				hPix_P_3Th[j]->Fill(hiNpix);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter3[j])==0 ){
				hHF_F_3Th[j]->Fill(hiHF);
				hTrk_F_3Th[j]->Fill(hiNtracks);
				hPix_F_3Th[j]->Fill(hiNpix);
			}
			//4ThX
			if (j==0)CoincFilter4[j]=phfCoincFilter4Th1; else if(j==1) CoincFilter4[j]=phfCoincFilter4Th2;	else if(j==2) CoincFilter4[j]=phfCoincFilter4Th3; else if(j==3) CoincFilter4[j]=phfCoincFilter4Th4; else if(j==4) CoincFilter4[j]=phfCoincFilter4Th5;
			else if(j==5)CoincFilter4[j]=phfCoincFilter4Th1p5; else if(j==6) CoincFilter4[j]=phfCoincFilter4Th2p5; else if(j==7) CoincFilter4[j]=phfCoincFilter4Th3p5; else if(j==8) CoincFilter4[j]=phfCoincFilter4Th4p5; else if(j==9) CoincFilter4[j]=phfCoincFilter4Th4p5;
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter4[j])==1 ){  
				hHF_Pixels_P_4Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_P_4Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_P_4Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter4[j])==0 ){  
				hHF_Pixels_F_4Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_F_4Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_F_4Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter4[j]==1) ){
				hHF_P_4Th[j]->Fill(hiHF);
				hTrk_P_4Th[j]->Fill(hiNtracks);
				hPix_P_4Th[j]->Fill(hiNpix);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter4[j])==0 ){
				hHF_F_4Th[j]->Fill(hiHF);
				hTrk_F_4Th[j]->Fill(hiNtracks);
				hPix_F_4Th[j]->Fill(hiNpix);
			}
			//5ThX
			if (j==0)CoincFilter5[j]=phfCoincFilter5Th1; else if(j==1) CoincFilter5[j]=phfCoincFilter5Th2;	else if(j==2) CoincFilter5[j]=phfCoincFilter5Th3; else if(j==3) CoincFilter5[j]=phfCoincFilter5Th4; else if(j==4) CoincFilter5[j]=phfCoincFilter5Th5;
			else if(j==5)CoincFilter5[j]=phfCoincFilter5Th1p5; else if(j==6) CoincFilter5[j]=phfCoincFilter5Th2p5; else if(j==7) CoincFilter5[j]=phfCoincFilter5Th3p5; else if(j==8) CoincFilter5[j]=phfCoincFilter5Th4p5; else if(j==9) CoincFilter5[j]=phfCoincFilter5Th4p5;
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter5[j])==1 ){  
				hHF_Pixels_P_5Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_P_5Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_P_5Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter5[j])==0 ){  
				hHF_Pixels_F_5Th[j]->Fill(hiNpix,hiHF);
				hHF_Ntracks_F_5Th[j]->Fill(hiNtracks,hiHF);
				hPix_Trk_F_5Th[j]->Fill(hiNpix,hiNtracks);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1 && CoincFilter5[j]==1) ){
				hHF_P_5Th[j]->Fill(hiHF);
				hTrk_P_5Th[j]->Fill(hiNtracks);
				hPix_P_5Th[j]->Fill(hiNpix);
			}
			if ( (HLT_HIZeroBias_v1==1) && (pprimaryVertexFilter && pclusterCompatibilityFilter && CoincFilter5[j])==0 ){
				hHF_F_5Th[j]->Fill(hiHF);
				hTrk_F_5Th[j]->Fill(hiNtracks);
				hPix_F_5Th[j]->Fill(hiNpix);
			}
		}
	}//// End Loop

	outFile->cd();
	//// Write 
	for (int j=0;j<9;j++)
	{
		//ThX
		hHF_Pixels_P[j]->Write(); 
		hHF_Ntracks_P[j]->Write();
		hPix_Trk_P[j]->Write();

		hHF_Pixels_F[j]->Write(); 
		hHF_Ntracks_F[j]->Write();
		hPix_Trk_F[j]->Write();

		hHF_P[j]->Write(); 
		hTrk_P[j]->Write();
		hPix_P[j]->Write();
		hHF_F[j]->Write();
		hTrk_F[j]->Write();
		hPix_F[j]->Write();

		//2ThX
		hHF_Pixels_P_2Th[j]->Write(); 
		hHF_Ntracks_P_2Th[j]->Write();
		hPix_Trk_P_2Th[j]->Write();

		hHF_Pixels_F_2Th[j]->Write(); 
		hHF_Ntracks_F_2Th[j]->Write();
		hPix_Trk_F_2Th[j]->Write();

		hHF_P_2Th[j]->Write(); 
		hHF_F_2Th[j]->Write();
		hTrk_P_2Th[j]->Write();
		hPix_P_2Th[j]->Write();
		hTrk_F_2Th[j]->Write();
		hPix_F_2Th[j]->Write();

		//3ThX
		hHF_Pixels_P_3Th[j]->Write(); 
		hHF_Ntracks_P_3Th[j]->Write();
		hPix_Trk_P_3Th[j]->Write();

		hHF_Pixels_F_3Th[j]->Write(); 
		hHF_Ntracks_F_3Th[j]->Write();
		hPix_Trk_F_3Th[j]->Write();

		hHF_P_3Th[j]->Write(); 
		hHF_F_3Th[j]->Write();
		hTrk_P_3Th[j]->Write();
		hPix_P_3Th[j]->Write();
		hTrk_F_3Th[j]->Write();
		hPix_F_3Th[j]->Write();

		//4ThX
		hHF_Pixels_P_4Th[j]->Write(); 
		hHF_Ntracks_P_4Th[j]->Write();
		hPix_Trk_P_4Th[j]->Write();

		hHF_Pixels_F_4Th[j]->Write(); 
		hHF_Ntracks_F_4Th[j]->Write();
		hPix_Trk_F_4Th[j]->Write();

		hHF_P_4Th[j]->Write(); 
		hHF_F_4Th[j]->Write();
		hTrk_P_4Th[j]->Write();
		hPix_P_4Th[j]->Write();
		hTrk_F_4Th[j]->Write();
		hPix_F_4Th[j]->Write();

		//5ThX
		hHF_Pixels_P_5Th[j]->Write(); 
		hHF_Ntracks_P_5Th[j]->Write();
		hPix_Trk_P_5Th[j]->Write();

		hHF_Pixels_F_5Th[j]->Write(); 
		hHF_Ntracks_F_5Th[j]->Write();
		hPix_Trk_F_5Th[j]->Write();

		hHF_P_5Th[j]->Write(); 
		hHF_F_5Th[j]->Write();
		hTrk_P_5Th[j]->Write();
		hPix_P_5Th[j]->Write();
		hTrk_F_5Th[j]->Write();
		hPix_F_5Th[j]->Write();
	}
	outFile->Close();
	t -> Stop();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());
}
