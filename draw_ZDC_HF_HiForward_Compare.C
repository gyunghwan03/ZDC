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

void draw_ZDC_HF_HiForward_Compare(int nevt=-1) 
{
	TStopwatch *t = new TStopwatch();
	t->Start();
	using namespace std;

	TString DATE = "210726";
	gSystem->mkdir(Form("figs/%s", DATE.Data()),kTRUE);

	gStyle->SetOptStat(00000);
	gROOT->ForceStyle();

	//TFile *fin = new TFile("./roots/HiForward_ForwardPD_Merged_all_*.root","read");
	//TTree *skimtree = (TTree*) fin -> Get("skimtree");
	
	TString fname = "./roots/HiForestAOD_merged.root";
	TString fdata = "./roots/HiForward_ForwardPD_Merged_all_*.root";
//    TString fname = "/eos/cms/store/group/phys_heavyions/bdiab/Centrality/HIForest_MC_HydjetDrum5F2019_2018PbPb/HiForestAOD_merged.root";
	TChain *evttree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *hlttree = new TChain("hltanalysis/HltTree");
	TChain *skimtree = new TChain("skimanalysis/HltTree");
//	TChain *skimtree = new TChain("skimtree");
	skimtree->Add(fname.Data());
	evttree->Add(fname.Data());
	hlttree->Add(fname.Data());

	TChain *dataTree = new TChain("skimtree");
	dataTree->Add(fdata.Data());

	//	int nolines=nolines_region1+nolines_region2+nolines_region3+nolines_region4;
	int nPixBins_P, nZDCBins_P, nTrkBins_P, nHFBins_P;
	int nPixBins_F, nZDCBins_F, nTrkBins_F, nHFBins_F;
	double npix_max_P, zdc_max_P, ntrk_max_P, hf_max_P, npix_min_P, zdc_min_P, ntrk_min_P, hf_min_P;
	nPixBins_P=500;
	nZDCBins_P=100;
	nTrkBins_P=50;
	nHFBins_P =200;
	npix_max_P=500.0;
	zdc_max_P=450000.0;  
	ntrk_max_P=50;
	hf_max_P=50;

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

	TBranch         *b_pprimaryVertexFilter;
	TBranch         *b_phfCoincFilterTh1;
	TBranch         *b_phfCoincFilterTh1p5;
	TBranch         *b_phfCoincFilterTh2;
	TBranch         *b_phfCoincFilterTh2p5;
	TBranch         *b_phfCoincFilterTh3;
	TBranch         *b_phfCoincFilterTh3p5;
	TBranch         *b_phfCoincFilterTh4;
	TBranch         *b_phfCoincFilterTh4p5;
	TBranch         *b_phfCoincFilterTh5;
	TBranch         *b_phfCoincFilter2Th1;
	TBranch         *b_phfCoincFilter2Th2;
	TBranch         *b_phfCoincFilter2Th3;
	TBranch         *b_phfCoincFilter2Th4;
	TBranch         *b_phfCoincFilter2Th5;
	TBranch         *b_phfCoincFilter2Th1p5;
	TBranch         *b_phfCoincFilter2Th2p5;
	TBranch         *b_phfCoincFilter2Th3p5;
	TBranch         *b_phfCoincFilter2Th4p5;
	TBranch         *b_phfCoincFilter2Th5p5;
	TBranch         *b_phfCoincFilter3Th1;
	TBranch         *b_phfCoincFilter3Th2;
	TBranch         *b_phfCoincFilter3Th3;
	TBranch         *b_phfCoincFilter3Th4;
	TBranch         *b_phfCoincFilter3Th5;
	TBranch         *b_phfCoincFilter3Th1p5;
	TBranch         *b_phfCoincFilter3Th2p5;
	TBranch         *b_phfCoincFilter3Th3p5;
	TBranch         *b_phfCoincFilter3Th4p5;
	TBranch         *b_phfCoincFilter3Th5p5;
	TBranch         *b_phfCoincFilter4Th1;
	TBranch         *b_phfCoincFilter4Th2;
	TBranch         *b_phfCoincFilter4Th3;
	TBranch         *b_phfCoincFilter4Th4;
	TBranch         *b_phfCoincFilter4Th5;
	TBranch         *b_phfCoincFilter4Th1p5;
	TBranch         *b_phfCoincFilter4Th2p5;
	TBranch         *b_phfCoincFilter4Th3p5;
	TBranch         *b_phfCoincFilter4Th4p5;
	TBranch         *b_phfCoincFilter4Th5p5;
	TBranch         *b_phfCoincFilter5Th1;
	TBranch         *b_phfCoincFilter5Th2;
	TBranch         *b_phfCoincFilter5Th3;
	TBranch         *b_phfCoincFilter5Th4;
	TBranch         *b_phfCoincFilter5Th5;
	TBranch         *b_phfCoincFilter5Th1p5;
	TBranch         *b_phfCoincFilter5Th2p5;
	TBranch         *b_phfCoincFilter5Th3p5;
	TBranch         *b_phfCoincFilter5Th4p5;
	TBranch         *b_phfCoincFilter5Th5p5;
	TBranch         *b_pclusterCompatibilityFilter;

	TBranch        *b_evt;
    TBranch        *b_hiBin;
    TBranch        *b_hiNpix;
    TBranch        *b_hiNtracks;
    TBranch        *b_hiHF;

	TBranch         *b_HLT_HIZeroBias_v1;
    hlttree->SetBranchAddress("HLT_HIZeroBias_v1", &HLT_HIZeroBias_v1, &b_HLT_HIZeroBias_v1);

	evttree->SetBranchAddress("evt",&evt,&b_evt);
	evttree->SetBranchAddress("hiNpix",&hiNpix,&b_hiNpix);
	evttree->SetBranchAddress("hiNtracks",&hiNtracks,&b_hiNtracks);
	evttree->SetBranchAddress("hiBin",&hiBin,&b_hiBin);
	evttree->SetBranchAddress("hiHF",&hiHF,&b_hiHF);

	skimtree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter,&b_pprimaryVertexFilter);
	skimtree->SetBranchAddress("phfCoincFilterTh1"   ,&phfCoincFilterTh1	,&b_phfCoincFilterTh1);
	skimtree->SetBranchAddress("phfCoincFilterTh2"   ,&phfCoincFilterTh2	,&b_phfCoincFilterTh2);
	skimtree->SetBranchAddress("phfCoincFilterTh3"   ,&phfCoincFilterTh3	,&b_phfCoincFilterTh3);
	skimtree->SetBranchAddress("phfCoincFilterTh4"   ,&phfCoincFilterTh4	,&b_phfCoincFilterTh4);
	skimtree->SetBranchAddress("phfCoincFilterTh5"   ,&phfCoincFilterTh5	,&b_phfCoincFilterTh5);
	skimtree->SetBranchAddress("phfCoincFilterTh1p5" ,&phfCoincFilterTh1p5	,&b_phfCoincFilterTh1p5);
	skimtree->SetBranchAddress("phfCoincFilterTh2p5" ,&phfCoincFilterTh2p5	,&b_phfCoincFilterTh2p5);
	skimtree->SetBranchAddress("phfCoincFilterTh3p5" ,&phfCoincFilterTh3p5	,&b_phfCoincFilterTh3p5);
	skimtree->SetBranchAddress("phfCoincFilterTh4p5" ,&phfCoincFilterTh4p5	,&b_phfCoincFilterTh4p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th1"  ,&phfCoincFilter2Th1	,&b_phfCoincFilter2Th1);
	skimtree->SetBranchAddress("phfCoincFilter2Th2"  ,&phfCoincFilter2Th2	,&b_phfCoincFilter2Th2);
	skimtree->SetBranchAddress("phfCoincFilter2Th3"  ,&phfCoincFilter2Th3	,&b_phfCoincFilter2Th3);
	skimtree->SetBranchAddress("phfCoincFilter2Th4"  ,&phfCoincFilter2Th4	,&b_phfCoincFilter2Th4);
	skimtree->SetBranchAddress("phfCoincFilter2Th5"  ,&phfCoincFilter2Th5	,&b_phfCoincFilter2Th5);
	skimtree->SetBranchAddress("phfCoincFilter2Th1p5",&phfCoincFilter2Th1p5	,&b_phfCoincFilter2Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th2p5",&phfCoincFilter2Th2p5	,&b_phfCoincFilter2Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th3p5",&phfCoincFilter2Th3p5	,&b_phfCoincFilter2Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter2Th4p5",&phfCoincFilter2Th4p5	,&b_phfCoincFilter2Th4p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th1"  ,&phfCoincFilter3Th1	,&b_phfCoincFilter3Th1);
	skimtree->SetBranchAddress("phfCoincFilter3Th2"  ,&phfCoincFilter3Th2	,&b_phfCoincFilter3Th2);
	skimtree->SetBranchAddress("phfCoincFilter3Th3"  ,&phfCoincFilter3Th3	,&b_phfCoincFilter3Th3);
	skimtree->SetBranchAddress("phfCoincFilter3Th4"  ,&phfCoincFilter3Th4	,&b_phfCoincFilter3Th4);
	skimtree->SetBranchAddress("phfCoincFilter3Th5"  ,&phfCoincFilter3Th5	,&b_phfCoincFilter3Th5);
	skimtree->SetBranchAddress("phfCoincFilter3Th1p5",&phfCoincFilter3Th1p5	,&b_phfCoincFilter3Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th2p5",&phfCoincFilter3Th2p5	,&b_phfCoincFilter3Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th3p5",&phfCoincFilter3Th3p5	,&b_phfCoincFilter3Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter3Th4p5",&phfCoincFilter3Th4p5	,&b_phfCoincFilter3Th4p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th1"  ,&phfCoincFilter4Th1	,&b_phfCoincFilter4Th1);
	skimtree->SetBranchAddress("phfCoincFilter4Th2"  ,&phfCoincFilter4Th2	,&b_phfCoincFilter4Th2);
	skimtree->SetBranchAddress("phfCoincFilter4Th3"  ,&phfCoincFilter4Th3	,&b_phfCoincFilter4Th3);
	skimtree->SetBranchAddress("phfCoincFilter4Th4"  ,&phfCoincFilter4Th4	,&b_phfCoincFilter4Th4);
	skimtree->SetBranchAddress("phfCoincFilter4Th5"  ,&phfCoincFilter4Th5	,&b_phfCoincFilter4Th5);
	skimtree->SetBranchAddress("phfCoincFilter4Th1p5",&phfCoincFilter4Th1p5	,&b_phfCoincFilter4Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th2p5",&phfCoincFilter4Th2p5	,&b_phfCoincFilter4Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th3p5",&phfCoincFilter4Th3p5	,&b_phfCoincFilter4Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter4Th4p5",&phfCoincFilter4Th4p5	,&b_phfCoincFilter4Th4p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th1"  ,&phfCoincFilter5Th1	,&b_phfCoincFilter5Th1);
	skimtree->SetBranchAddress("phfCoincFilter5Th2"  ,&phfCoincFilter5Th2	,&b_phfCoincFilter5Th2);
	skimtree->SetBranchAddress("phfCoincFilter5Th3"  ,&phfCoincFilter5Th3	,&b_phfCoincFilter5Th3);
	skimtree->SetBranchAddress("phfCoincFilter5Th4"  ,&phfCoincFilter5Th4	,&b_phfCoincFilter5Th4);
	skimtree->SetBranchAddress("phfCoincFilter5Th5"  ,&phfCoincFilter5Th5	,&b_phfCoincFilter5Th5);
	skimtree->SetBranchAddress("phfCoincFilter5Th1p5",&phfCoincFilter5Th1p5	,&b_phfCoincFilter5Th1p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th2p5",&phfCoincFilter5Th2p5	,&b_phfCoincFilter5Th2p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th3p5",&phfCoincFilter5Th3p5	,&b_phfCoincFilter5Th3p5);
	skimtree->SetBranchAddress("phfCoincFilter5Th4p5",&phfCoincFilter5Th4p5	,&b_phfCoincFilter5Th4p5);
	skimtree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter,&b_pclusterCompatibilityFilter);


	dataTree->SetBranchAddress("evt",&evt);
    dataTree->SetBranchAddress("hiNpix",&hiNpix);
    dataTree->SetBranchAddress("hiNtracks",&hiNtracks);
    dataTree->SetBranchAddress("hiBin",&hiBin);
    dataTree->SetBranchAddress("hiHF",&hiHF);
    dataTree->SetBranchAddress("phfCoincFilterTh1"  ,&phfCoincFilterTh1);
    dataTree->SetBranchAddress("phfCoincFilterTh2"  ,&phfCoincFilterTh2);
    dataTree->SetBranchAddress("phfCoincFilterTh3"  ,&phfCoincFilterTh3);
    dataTree->SetBranchAddress("phfCoincFilterTh4"  ,&phfCoincFilterTh4);
    dataTree->SetBranchAddress("phfCoincFilterTh5"  ,&phfCoincFilterTh5);
    dataTree->SetBranchAddress("phfCoincFilterTh1p5",&phfCoincFilterTh1p5);
    dataTree->SetBranchAddress("phfCoincFilterTh2p5",&phfCoincFilterTh2p5);
    dataTree->SetBranchAddress("phfCoincFilterTh3p5",&phfCoincFilterTh3p5);
    dataTree->SetBranchAddress("phfCoincFilterTh4p5",&phfCoincFilterTh4p5);
    dataTree->SetBranchAddress("phfCoincFilter2Th1"  ,&phfCoincFilter2Th1);
    dataTree->SetBranchAddress("phfCoincFilter2Th2"  ,&phfCoincFilter2Th2);
    dataTree->SetBranchAddress("phfCoincFilter2Th3"  ,&phfCoincFilter2Th3);
    dataTree->SetBranchAddress("phfCoincFilter2Th4"  ,&phfCoincFilter2Th4);
    dataTree->SetBranchAddress("phfCoincFilter2Th5"  ,&phfCoincFilter2Th5);
    dataTree->SetBranchAddress("phfCoincFilter2Th1p5",&phfCoincFilter2Th1p5);
    dataTree->SetBranchAddress("phfCoincFilter2Th2p5",&phfCoincFilter2Th2p5);
    dataTree->SetBranchAddress("phfCoincFilter2Th3p5",&phfCoincFilter2Th3p5);
    dataTree->SetBranchAddress("phfCoincFilter2Th4p5",&phfCoincFilter2Th4p5);
    dataTree->SetBranchAddress("phfCoincFilter3Th1"  ,&phfCoincFilter3Th1);
    dataTree->SetBranchAddress("phfCoincFilter3Th2"  ,&phfCoincFilter3Th2);
    dataTree->SetBranchAddress("phfCoincFilter3Th3"  ,&phfCoincFilter3Th3);
    dataTree->SetBranchAddress("phfCoincFilter3Th4"  ,&phfCoincFilter3Th4);
    dataTree->SetBranchAddress("phfCoincFilter3Th5"  ,&phfCoincFilter3Th5);
    dataTree->SetBranchAddress("phfCoincFilter3Th1p5",&phfCoincFilter3Th1p5);
    dataTree->SetBranchAddress("phfCoincFilter3Th2p5",&phfCoincFilter3Th2p5);
    dataTree->SetBranchAddress("phfCoincFilter3Th3p5",&phfCoincFilter3Th3p5);
    dataTree->SetBranchAddress("phfCoincFilter3Th4p5",&phfCoincFilter3Th4p5);
    dataTree->SetBranchAddress("phfCoincFilter4Th1"  ,&phfCoincFilter4Th1);
    dataTree->SetBranchAddress("phfCoincFilter4Th2"  ,&phfCoincFilter4Th2);
    dataTree->SetBranchAddress("phfCoincFilter4Th3"  ,&phfCoincFilter4Th3);
    dataTree->SetBranchAddress("phfCoincFilter4Th4"  ,&phfCoincFilter4Th4);
    dataTree->SetBranchAddress("phfCoincFilter4Th5"  ,&phfCoincFilter4Th5);
    dataTree->SetBranchAddress("phfCoincFilter4Th1p5",&phfCoincFilter4Th1p5);
	//dataTree->SetBranchAddress("phfCoincFilter4Th2p5",&phfCoincFilter4Th2p5);
	dataTree->SetBranchAddress("phfCoincFilter4Th3p5",&phfCoincFilter4Th3p5);
    dataTree->SetBranchAddress("phfCoincFilter4Th4p5",&phfCoincFilter4Th4p5);
    dataTree->SetBranchAddress("phfCoincFilter5Th1"  ,&phfCoincFilter5Th1);
    dataTree->SetBranchAddress("phfCoincFilter5Th2"  ,&phfCoincFilter5Th2);
    dataTree->SetBranchAddress("phfCoincFilter5Th3"  ,&phfCoincFilter5Th3);
    dataTree->SetBranchAddress("phfCoincFilter5Th4"  ,&phfCoincFilter5Th4);
    dataTree->SetBranchAddress("phfCoincFilter5Th5"  ,&phfCoincFilter5Th5);
    dataTree->SetBranchAddress("phfCoincFilter5Th1p5",&phfCoincFilter5Th1p5);
    dataTree->SetBranchAddress("phfCoincFilter5Th2p5",&phfCoincFilter5Th2p5);
    dataTree->SetBranchAddress("phfCoincFilter5Th3p5",&phfCoincFilter5Th3p5);
    dataTree->SetBranchAddress("phfCoincFilter5Th4p5",&phfCoincFilter5Th4p5);
    dataTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);


	TString filter[5][9] = { {"Th1","Th1p5","Th2","Th2p5","Th3","Th3p5","Th4","Th4p5","Th5"},
		{"2Th1","2Th1p5","2Th2","2Th2p5","2Th3","2Th3p5","2Th4","2Th4p5","2Th5"},
		{"3Th1","3Th1p5","3Th2","3Th2p5","3Th3","3Th3p5","3Th4","3Th4p5","3Th5"},
		{"4Th1","4Th1p5","4Th2","4Th2p5","4Th3","4Th3p5","4Th4","4Th4p5","4Th5"},
		{"5Th1","5Th1p5","5Th2","5Th2p5","5Th3","5Th3p5","5Th4","5Th4p5","5Th5"}
	};
	TString s[5][9];
	for(int i=0;i<5;i++) {
		for (int j=0;j<9;j++) {	s[i][j]="phfCoincFilter"; }
	}
	//cout << "fName : " << fName.Data() << endl;

	TH2F* hHF_Pixels_P[5][9];
	TH2F* hHF_Ntracks_P[5][9];
	TH2F* hPix_Trk_P[5][9];

	TH2F* hHF_Pixels_F[5][9];
	TH2F* hHF_Ntracks_F[5][9];
	TH2F* hPix_Trk_F[5][9];

	TH1F *hHF_P[5][9];
	TH1F *hTrk_P[5][9];
	TH1F *hPix_P[5][9];
	TH1F *hHF_F[5][9];
	TH1F *hTrk_F[5][9];
	TH1F *hPix_F[5][9];

	TH1F *hHF_P_Hydjet = new TH1F("hHF_P_Hydjet","",100,0,7000);
	TH1F *hHF_P_data = new TH1F("hHF_P_data","",100,0,7000);

	for(int i=0;i<5;i++){
		for(int j=0;j<9;j++){
			// ThX
			hHF_Pixels_P[i][j]  = new TH2F(Form("hHF_Pixels_P_%dTh_%d",i+1,j),";No. of pixels;HF",nPixBins_P,npix_min_P,npix_max_P,200,hf_min_P,hf_max_P);
			hHF_Ntracks_P[i][j] = new TH2F(Form("hHF_Tracks_P_%dTh_%d",i+1,j),";No. of tracks;HF",nTrkBins_P,ntrk_min_P,ntrk_max_P,600,hf_min_P,hf_max_P);
			hPix_Trk_P[i][j]    = new TH2F(Form("hPix_Trk_P_%dTh_%d",i+1,j),";No. of pixels;No. of tracks",250,npix_min_F,npix_max_P,nTrkBins_P,ntrk_min_F,ntrk_max_P);

			hHF_Pixels_F[i][j]  = new TH2F(Form("hHF_Pixels_F_%dTh_%d",i+1,j),";No. of pixels;HF",nPixBins_F,npix_min_F,npix_max_F,nHFBins_F,hf_min_F,hf_max_F);
			hHF_Ntracks_F[i][j] = new TH2F(Form("hHF_Tracks_F_%dTh_%d",i+1,j),";No. of tracks;HF",nTrkBins_F,ntrk_min_F,ntrk_max_F,nHFBins_F,hf_min_F,hf_max_F);
			hPix_Trk_F[i][j]    = new TH2F(Form("hPix_Trk_F_%dTh_%d",i+1,j),";No. of pixels;No. of tracks",nPixBins_F,npix_min_F,npix_max_F,nTrkBins_F,ntrk_min_F,ntrk_max_F);

			hHF_P[i][j]  = new TH1F(Form("hHF_P_%dTh_%d",i+1,j),"",100,0,6000);
			hHF_F[i][j]  = new TH1F(Form("hHF_F_%dTh_%d",i+1,j),"",100,0,hf_max_P);
			hTrk_P[i][j] = new TH1F(Form("hTrk_P_%dTh_%d",i+1,j),"",100,0,100);
			hPix_P[i][j] = new TH1F(Form("hPix_P_%dTh_%d",i+1,j),"",100,0,700);
			hTrk_F[i][j] = new TH1F(Form("hTrk_F_%dTh_%d",i+1,j),"",100,0,100);
			hPix_F[i][j] = new TH1F(Form("hPix_F_%dTh_%d",i+1,j),"",100,0,300);
		}
	}

	if (nevt==-1) nevt = skimtree->GetEntries();

	// Start Loop
	cout << "###########  Hydjet Loop  ###########" << endl;
	for (int iev=0;iev<nevt;iev++) {
		if(iev%1000000==0) cout << ">>>>> EVENT " << iev << " / " << skimtree->GetEntries() <<  " ("<<(int)(100.*iev/skimtree->GetEntries()) << "%)" << endl;
		skimtree->GetEntry(iev);
		hlttree->GetEntry(iev);
		evttree->GetEntry(iev);
		//if (HLT_HIZeroBias_v1==0) continue;
		//if (pprimaryVertexFilter==0) continue;
		//if (pprimaryVertexFilter==1 && pclusterCompatibilityFilter==1) hHF_P_Hydjet->Fill(hiHF);
		hHF_P_Hydjet->Fill(hiHF);
		/*if (hiBin>140) {
			int CoincFilter[5][9] = {
				{phfCoincFilterTh1, phfCoincFilterTh1p5, phfCoincFilterTh2, phfCoincFilterTh2p5,phfCoincFilterTh3, phfCoincFilterTh3p5, phfCoincFilterTh4, phfCoincFilterTh4p5, phfCoincFilterTh5},
				{phfCoincFilter2Th1, phfCoincFilter2Th1p5, phfCoincFilter2Th2, phfCoincFilter2Th2p5,phfCoincFilter2Th3, phfCoincFilter2Th3p5, phfCoincFilter2Th4, phfCoincFilter2Th4p5, phfCoincFilter2Th5},
				{phfCoincFilter3Th1, phfCoincFilter3Th1p5, phfCoincFilter3Th2, phfCoincFilter3Th2p5,phfCoincFilter3Th3, phfCoincFilter3Th3p5, phfCoincFilter3Th4, phfCoincFilter3Th4p5, phfCoincFilter3Th5},
				{phfCoincFilter4Th1, phfCoincFilter4Th1p5, phfCoincFilter4Th2, phfCoincFilter4Th2p5,phfCoincFilter4Th3, phfCoincFilter4Th3p5, phfCoincFilter4Th4, phfCoincFilter4Th4p5, phfCoincFilter4Th5},
				{phfCoincFilter5Th1, phfCoincFilter5Th1p5, phfCoincFilter5Th2, phfCoincFilter5Th2p5,phfCoincFilter5Th3, phfCoincFilter5Th3p5, phfCoincFilter5Th4, phfCoincFilter5Th4p5, phfCoincFilter5Th5}
			};
			for (int i=0;i<5;i++){
				for (int j=0;j<9;j++){
					if ( (pclusterCompatibilityFilter==1 && CoincFilter[i][j]==1) ){
						hHF_Pixels_P[i][j]->Fill(hiNpix,hiHF);
						hHF_Ntracks_P[i][j]->Fill(hiNtracks,hiHF);
						hPix_Trk_P[i][j]->Fill(hiNpix,hiNtracks);
						hHF_P[i][j]->Fill(hiHF);
						hTrk_P[i][j]->Fill(hiNtracks);
						hPix_P[i][j]->Fill(hiNpix);
					}
					if ( (pclusterCompatibilityFilter==0 && CoincFilter[i][j]==0) ){
						hHF_Pixels_F[i][j]->Fill(hiNpix,hiHF);
						hHF_Ntracks_F[i][j]->Fill(hiNtracks,hiHF);
						hPix_Trk_F[i][j]->Fill(hiNpix,hiNtracks);
						hHF_F[i][j]->Fill(hiHF);
						hTrk_F[i][j]->Fill(hiNtracks);
						hPix_F[i][j]->Fill(hiNpix);                
					}
				}
			}
		}*/
	}//// End Loop

	cout << "###########  DATA Loop  ###########" << endl;
	// Start data Loop
	for (int iev=0;iev<nevt;iev++) {
		if(iev%1000000==0) cout << ">>>>> EVENT " << iev << " / " << dataTree->GetEntries() <<  " ("<<(int)(100.*iev/dataTree->GetEntries()) << "%)" << endl;
		dataTree->GetEntry(iev);
		//if ( pclusterCompatibilityFilter==1 ) hHF_P_data->Fill(hiHF);
		hHF_P_data->Fill(hiHF);
		/*if (hiBiVn>140) {
			int CoincFilter[5][9] = {
				{phfCoincFilterTh1, phfCoincFilterTh1p5, phfCoincFilterTh2, phfCoincFilterTh2p5,phfCoincFilterTh3, phfCoincFilterTh3p5, phfCoincFilterTh4, phfCoincFilterTh4p5, phfCoincFilterTh5},
				{phfCoincFilter2Th1, phfCoincFilter2Th1p5, phfCoincFilter2Th2, phfCoincFilter2Th2p5,phfCoincFilter2Th3, phfCoincFilter2Th3p5, phfCoincFilter2Th4, phfCoincFilter2Th4p5, phfCoincFilter2Th5},
				{phfCoincFilter3Th1, phfCoincFilter3Th1p5, phfCoincFilter3Th2, phfCoincFilter3Th2p5,phfCoincFilter3Th3, phfCoincFilter3Th3p5, phfCoincFilter3Th4, phfCoincFilter3Th4p5, phfCoincFilter3Th5},
				{phfCoincFilter4Th1, phfCoincFilter4Th1p5, phfCoincFilter4Th2, phfCoincFilter4Th2p5,phfCoincFilter4Th3, phfCoincFilter4Th3p5, phfCoincFilter4Th4, phfCoincFilter4Th4p5, phfCoincFilter4Th5},
				{phfCoincFilter5Th1, phfCoincFilter5Th1p5, phfCoincFilter5Th2, phfCoincFilter5Th2p5,phfCoincFilter5Th3, phfCoincFilter5Th3p5, phfCoincFilter5Th4, phfCoincFilter5Th4p5, phfCoincFilter5Th5}
			};
			for (int i=0;i<5;i++){
				for (int j=0;j<9;j++){
					if ( (pclusterCompatibilityFilter==1 && CoincFilter[i][j]==1) ){
						//hHF_Pixels_P[i][j]->Fill(hiNpix,hiHF);
						//hHF_Ntracks_P[i][j]->Fill(hiNtracks,hiHF);
						//hPix_Trk_P[i][j]->Fill(hiNpix,hiNtracks);
						hHF_P_data[i][j]->Fill(hiHF);
						hTrk_P_data[i][j]->Fill(hiNtracks);
						hPix_P_data[i][j]->Fill(hiNpix);
					}
					if ( (pclusterCompatibilityFilter==0 && CoincFilter[i][j]==0) ){
						hHF_Pixels_F[i][j]->Fill(hiNpix,hiHF);
						hHF_Ntracks_F[i][j]->Fill(hiNtracks,hiHF);
						hPix_Trk_F[i][j]->Fill(hiNpix,hiNtracks);
						hHF_F[i][j]->Fill(hiHF);
						hTrk_F[i][j]->Fill(hiNtracks);
						hPix_F[i][j]->Fill(hiNpix);                
					}
				}
			}
		}*/
	}//// End Loop

	//hHF_P_data->Scale(factor/hHF_P_data->GetEntries());
	//Double_t factor=(hHF_P_data->GetEntries())/(hHF_P_Hydjet->GetEntries());
	//cout << "data Entries : " << hHF_P_data->GetEntries() << " Hydjet Entrie : " << hHF_P_Hydjet->GetEntries() << " factor : " << factor << endl;
	Double_t factor = ((double) dataTree->GetEntries())/( (double) skimtree->GetEntries());
	cout << "factor : " << factor << endl;
	hHF_P_Hydjet->Scale(factor/hHF_P_Hydjet->GetEntries());
	TCanvas *c1 = new TCanvas("c1","",700,700);
	c1->cd();
	hHF_P_Hydjet->Draw();

	hHF_P_data->Scale(1/hHF_P_data->GetEntries());


	//////////////    START Drawing     //////////////
	TCanvas *cHF_compare = new TCanvas("cHF_compare","",700,700);
	TLegend *leg1 = new TLegend(0.58,0.55,0.8,0.8);
	leg1->SetFillColor(0);
    leg1->SetFillStyle(4000);
    leg1->SetBorderSize(0);
    leg1->SetMargin(0.2);
    leg1->SetTextSize(0.029);
    leg1->SetTextFont(42);
	cHF_compare->cd();
	cHF_compare->SetLogy();
	hHF_P_Hydjet->GetXaxis()->SetTitle("hiHF");
	hHF_P_Hydjet->GetYaxis()->SetTitle("Normalized Distribution");
	hHF_P_Hydjet->GetXaxis()->SetRangeUser(0,7000);
	hHF_P_Hydjet->GetYaxis()->SetRangeUser(1e-5,10);
	hHF_P_Hydjet->SetLineColor(kRed+2);
	hHF_P_Hydjet->SetMarkerColor(kRed+2);
	hHF_P_Hydjet->SetMarkerStyle(24);
	hHF_P_Hydjet->Draw();

	hHF_P_data->GetXaxis()->SetTitle("hiHF");
	hHF_P_data->GetYaxis()->SetTitle("Normalized Distribution");
	hHF_P_data->GetXaxis()->SetRangeUser(0,7000);
	hHF_P_data->GetYaxis()->SetRangeUser(1e-5,10);
	hHF_P_data->SetLineColor(kBlue+2);
	hHF_P_data->SetMarkerColor(kBlue+2);
	hHF_P_data->SetMarkerStyle(26);
	hHF_P_data->Draw("same");

	leg1->AddEntry(hHF_P_data,"Data");
	leg1->AddEntry(hHF_P_Hydjet,"Hydjet MC");
	leg1->Draw("same");




	/*
	TString title[5][9];
	TCanvas* cHF1D_P[5];
	TLegend *leg1_P = new TLegend(0.58,0.55,0.8,0.8);
	leg1_P->SetFillColor(0);
	leg1_P->SetFillStyle(4000);
	leg1_P->SetBorderSize(0);
	leg1_P->SetMargin(0.2);
	leg1_P->SetTextSize(0.029);
	leg1_P->SetTextFont(42);
	for(int i=0;i<5;i++){
		cHF1D_P[i] = new TCanvas(Form("cHF1D_P_%d",i),"",700,700);
		cHF1D_P[i]->cd();
		cHF1D_P[i]->SetLogy();
		for(int j=0;j<9;j++)
		{
			title[i][j]=s[i][j]+filter[i][j];
			hHF_P[i][j]->SetLineColor(j+1);
			hHF_P[i][j]->GetXaxis()->SetTitle("HF");
			hHF_P[i][j]->GetYaxis()->SetTitle("Entries");
			hHF_P[i][j]->GetYaxis()->SetRangeUser(1e+2,1e+5);
			hHF_P[i][j]->GetYaxis()->CenterTitle();
			hHF_P[i][j]->GetXaxis()->CenterTitle();
			hHF_P[i][j]->GetYaxis()->SetTitleOffset(1.3);
			hHF_P[i][j]->SetTitle("Passed Condition HF");
			hHF_P[i][j]->Draw("same");
			leg1_P->AddEntry(hHF_P[i][j],Form("%s",title[i][j].Data()));
			leg1_P->Draw("same");
		}
		cHF1D_P[i]->SaveAs(Form("figs/%s/1D_HF_HiForward_%dTh_Passed.png",DATE.Data(),i+1));
		leg1_P->Clear();
	}

	TCanvas* cTrk1D_P[5];
	TLegend *leg2_P = new TLegend(0.58,0.55,0.8,0.8);
	leg2_P->SetFillColor(0);
	leg2_P->SetFillStyle(4000);
	leg2_P->SetBorderSize(0);
	leg2_P->SetMargin(0.2);
	leg2_P->SetTextSize(0.029);
	leg2_P->SetTextFont(42);
	for(int i=0;i<5;i++)
	{
		cTrk1D_P[i]= new TCanvas(Form("cTrk1D_P_%d",i),"",700,700);
		cTrk1D_P[i]->cd();
		cTrk1D_P[i]->SetLogy();
		for(int j=0;j<9;j++)
		{
			title[i][j]=s[i][j]+filter[i][j];
			hTrk_P[i][j]->SetLineColor(j+1);
			hTrk_P[i][j]->GetXaxis()->SetTitle("Num. of Trackers");
			hTrk_P[i][j]->GetYaxis()->SetTitle("Entries");
			hTrk_P[i][j]->GetYaxis()->SetRangeUser(1e+2,1e+5);
			hTrk_P[i][j]->GetYaxis()->CenterTitle();
			hTrk_P[i][j]->GetXaxis()->CenterTitle();
			hTrk_P[i][j]->GetYaxis()->SetTitleOffset(1.3);
			hTrk_P[i][j]->SetTitle("Passed Condition Trk");
			hTrk_P[i][j]->Draw("same");
			leg2_P->AddEntry(hTrk_P[i][j],Form("%s",title[i][j].Data()));
			leg2_P->Draw("same");
		}
		cTrk1D_P[i]->SaveAs(Form("figs/%s/1D_Trk_HiForward_%dTh_Passed.png",DATE.Data(),i+1));
		leg2_P->Clear();
	}

	TCanvas* cPix1D_P[5];
	TLegend *leg3_P = new TLegend(0.58,0.55,0.8,0.8);
	leg3_P->SetFillColor(0);
	leg3_P->SetFillStyle(4000);
	leg3_P->SetBorderSize(0);
	leg3_P->SetMargin(0.2);
	leg3_P->SetTextSize(0.029);
	leg3_P->SetTextFont(42);
	for(int i=0;i<5;i++)
	{
		cPix1D_P[i]= new TCanvas(Form("cPix1D_P_%d",i),"",700,700);
		cPix1D_P[i]->cd();
		cPix1D_P[i]->SetLogy();
		for(int j=0;j<9;j++)
		{
			title[i][j]=s[i][j]+filter[i][j];
			hPix_P[i][j]->SetLineColor(j+1);
			hPix_P[i][j]->GetXaxis()->SetTitle("Num. of Pixels");
			hPix_P[i][j]->GetYaxis()->SetTitle("Entries");
			hPix_P[i][j]->GetYaxis()->SetRangeUser(1e+2,1e+5);
			hPix_P[i][j]->GetYaxis()->CenterTitle();
			hPix_P[i][j]->GetXaxis()->CenterTitle();
			hPix_P[i][j]->GetYaxis()->SetTitleOffset(1.3);
			hPix_P[i][j]->SetTitle("Passed Condition Pix");
			hPix_P[i][j]->Draw("same");
			leg3_P->AddEntry(hPix_P[i][j],Form("%s",title[i][j].Data()));
			leg3_P->Draw("same");
		}
		cPix1D_P[i]->SaveAs(Form("figs/%s/1D_Pix_HiForward_%dTh_Passed.png",DATE.Data(),i+1));
		leg3_P->Clear();
	}

	TCanvas* cHF_Pixels_P[5][9];
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<9;j++)
		{
			cHF_Pixels_P[i][j] = new TCanvas(Form("cHF_Pixels_P_%dTh%d",i+1,j),"",700,700);
			cHF_Pixels_P[i][j]->cd();
			//			cHF_Pixels_P[i][j]->SetLogz();
			title[i][j]=s[i][j]+filter[i][j];
			hHF_Pixels_P[i][j]->SetLineColor(j+1);
			hHF_Pixels_P[i][j]->GetYaxis()->CenterTitle();
			hHF_Pixels_P[i][j]->GetXaxis()->CenterTitle();
			hHF_Pixels_P[i][j]->GetYaxis()->SetTitleOffset(1.3);
			hHF_Pixels_P[i][j]->GetZaxis()->SetLabelSize(0.03);
			hHF_Pixels_P[i][j]->SetTitle(Form("HF vs Number of Pixels Passed Cond. (%s)",title[i][j].Data()));
			hHF_Pixels_P[i][j]->Draw("colz");
			cHF_Pixels_P[i][j]->SaveAs(Form("figs/%s/HF_Pixels_HiForward_%s_Passed.png",DATE.Data(),filter[i][j].Data()));
		}
	}

	TCanvas* cHF_Ntracks_P[5][9];
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<9;j++)
		{
			cHF_Ntracks_P[i][j] = new TCanvas(Form("cHF_Ntracks_P_%dTh%d",i+1,j),"",700,700);
			cHF_Ntracks_P[i][j]->cd();
			//			cHF_Ntracks_P[i][j]->SetLogz();
			title[i][j]=s[i][j]+filter[i][j];
			hHF_Ntracks_P[i][j]->SetLineColor(j+1);
			hHF_Ntracks_P[i][j]->GetYaxis()->CenterTitle();
			hHF_Ntracks_P[i][j]->GetXaxis()->CenterTitle();
			hHF_Ntracks_P[i][j]->GetYaxis()->SetTitleOffset(1.3);
			hHF_Ntracks_P[i][j]->GetZaxis()->SetLabelSize(0.03);
			hHF_Ntracks_P[i][j]->SetTitle(Form("HF vs Number of Ntracks Passed Cond. (%s)",title[i][j].Data()));
			hHF_Ntracks_P[i][j]->Draw("colz");
			cHF_Ntracks_P[i][j]->SaveAs(Form("figs/%s/HF_Ntracks_HiForward_%s_Passed.png",DATE.Data(),filter[i][j].Data()));
		}
	}

	TCanvas* cPix_Trk_P[5][9];
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<9;j++)
		{
			cPix_Trk_P[i][j] = new TCanvas(Form("cPix_Trk_P_%dTh%d",i+1,j),"",700,700);
			cPix_Trk_P[i][j]->cd();
			//			cPix_Trk_P[i][j]->SetLogz();
			title[i][j]=s[i][j]+filter[i][j];
			hPix_Trk_P[i][j]->SetLineColor(j+1);
			hPix_Trk_P[i][j]->GetYaxis()->CenterTitle();
			hPix_Trk_P[i][j]->GetXaxis()->CenterTitle();
			hPix_Trk_P[i][j]->GetYaxis()->SetTitleOffset(1.3);
			hPix_Trk_P[i][j]->GetZaxis()->SetLabelSize(0.03);
			hPix_Trk_P[i][j]->SetTitle(Form("Pixels vs Tracks Passed Cond. (%s)",title[i][j].Data()));
			hPix_Trk_P[i][j]->Draw("colz");
			cPix_Trk_P[i][j]->SaveAs(Form("figs/%s/Pix_Trk_HiForward_%s_Passed.png",DATE.Data(),filter[i][j].Data()));
		}
	}
	/*
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

	   TCanvas* cHF1D_F[5];
	   TLegend *leg1_F = new TLegend(0.58,0.55,0.8,0.8);
	   leg1_F->SetFillColor(0);
	   leg1_F->SetFillStyle(4000);
	   leg1_F->SetBorderSize(0);
	   leg1_F->SetMargin(0.2);
	   leg1_F->SetTextSize(0.029);
	   leg1_F->SetTextFont(42);
	   for(int i=0;i<5;i++){
	   cHF1D_F[i]= new TCanvas(Form("cHF1D_F_%d",i),"",700,700);
	   cHF1D_F[i]->cd();
	   cHF1D_F[i]->SetLogy();
	   for(int j=0;j<9;j++)
	   {
	   title[i][j]=s[i][j]+filter[i][j];
	   hHF_F[i][j]->SetLineColor(j+1);
	   hHF_F[i][j]->GetXaxis()->SetTitle("HF");
	   hHF_F[i][j]->GetYaxis()->SetTitle("Entries");
	   hHF_F[i][j]->GetYaxis()->SetRangeUser(1e+2,1e+7);
	   hHF_F[i][j]->GetYaxis()->CenterTitle();
	hHF_F[i][j]->GetXaxis()->CenterTitle();
	hHF_F[i][j]->GetYaxis()->SetTitleOffset(1.3);
	hHF_F[i][j]->SetTitle("Failed Condition HF");
	hHF_F[i][j]->Draw("same");
	leg1_F->AddEntry(hHF_F[i][j],Form("%s",title[i][j].Data()));
}
leg1_F->Draw("same");
cHF1D_F[i]->SaveAs(Form("figs/%s/1D_HF_HiForward_%dTh_Falied.png",DATE.Data(),i+1));
}

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
*/
	t -> Stop();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());
}
