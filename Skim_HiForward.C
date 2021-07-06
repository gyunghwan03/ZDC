#include <ctime>
#include <TChain.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TFile.h>
#include <TStopwatch.h>
static const long MAXTREESIZE = 10000000000;

void Skim_HiForward()
{

	using namespace std;
	TStopwatch *t = new TStopwatch;
	t -> Start();

	TString fname1 = "roots/HiForestAOD_220.root";
	TString fname2 = "roots/HiForestAOD_684.root";

	TChain *evtTree = new TChain("hiEvtAnalyzer/HiTree");
    TChain *hltTree = new TChain("hltanalysis/HltTree");
    TChain *skimTree = new TChain("skimanalysis/HltTree");

	evtTree->Add(fname1.Data()); evtTree->Add(fname2.Data());
	hltTree->Add(fname1.Data()); hltTree->Add(fname2.Data());
	skimTree->Add(fname1.Data()); skimTree->Add(fname2.Data());

	TFile *outFile = new TFile("HiForward.root","recreate");

	ULong64_t		evt;
	Int_t			hiBin;
	Int_t			hiNpix;
	Int_t			hiNtracks;
	Float_t			hiHF;
	TBranch        *b_evt;
	TBranch        *b_hiBin;
	TBranch        *b_hiNpix;
	TBranch        *b_hiNtracks;
	TBranch        *b_hiHF;

    evtTree->SetBranchAddress("evt",&evt,&b_evt);
    evtTree->SetBranchAddress("hiNpix",&hiNpix,&b_hiNpix);
    evtTree->SetBranchAddress("hiNtracks",&hiNtracks,&b_hiNtracks);
    evtTree->SetBranchAddress("hiBin",&hiBin,&b_hiBin);
    evtTree->SetBranchAddress("hiHF",&hiHF,&b_hiHF);
	
	Int_t		   pprimaryVertexFilter;
	Int_t          phfCoincFilterTh1;
	Int_t          phfCoincFilterTh1p5;
	Int_t          phfCoincFilterTh2;
	Int_t          phfCoincFilterTh2p5;
	Int_t          phfCoincFilterTh3;
	Int_t          phfCoincFilterTh3p5;
	Int_t          phfCoincFilterTh4;
	Int_t          phfCoincFilterTh4p5;
	Int_t          phfCoincFilterTh5;
	Int_t          phfCoincFilter2Th1;
	Int_t          phfCoincFilter2Th2;
	Int_t          phfCoincFilter2Th3;
	Int_t          phfCoincFilter2Th4;
	Int_t          phfCoincFilter2Th5;
	Int_t          phfCoincFilter2Th1p5;
	Int_t          phfCoincFilter2Th2p5;
	Int_t          phfCoincFilter2Th3p5;
	Int_t          phfCoincFilter2Th4p5;
	Int_t          phfCoincFilter2Th5p5;
	Int_t          phfCoincFilter3Th1;
	Int_t          phfCoincFilter3Th2;
	Int_t          phfCoincFilter3Th3;
	Int_t          phfCoincFilter3Th4;
	Int_t          phfCoincFilter3Th5;
	Int_t          phfCoincFilter3Th1p5;
	Int_t          phfCoincFilter3Th2p5;
	Int_t          phfCoincFilter3Th3p5;
	Int_t          phfCoincFilter3Th4p5;
	Int_t          phfCoincFilter3Th5p5;
	Int_t          phfCoincFilter4Th1;
	Int_t          phfCoincFilter4Th2;
	Int_t          phfCoincFilter4Th3;
	Int_t          phfCoincFilter4Th4;
	Int_t          phfCoincFilter4Th5;
	Int_t          phfCoincFilter4Th1p5;
//	Int_t          phfCoincFilter4Th2p5;
	Int_t          phfCoincFilter4Th3p5;
	Int_t          phfCoincFilter4Th4p5;
	Int_t          phfCoincFilter4Th5p5;
	Int_t          phfCoincFilter5Th1;
	Int_t          phfCoincFilter5Th2;
	Int_t          phfCoincFilter5Th3;
	Int_t          phfCoincFilter5Th4;
	Int_t          phfCoincFilter5Th5;
	Int_t          phfCoincFilter5Th1p5;
	Int_t          phfCoincFilter5Th2p5;
	Int_t          phfCoincFilter5Th3p5;
	Int_t          phfCoincFilter5Th4p5;
	Int_t          phfCoincFilter5Th5p5;
	Int_t          pclusterCompatibilityFilter;
	TBranch			*b_pprimaryVertexFilter;
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
//	TBranch         *b_phfCoincFilter4Th2p5;
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

	skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter,&b_pprimaryVertexFilter);
    skimTree->SetBranchAddress("phfCoincFilterTh1",&phfCoincFilterTh1,&b_phfCoincFilterTh1);
    skimTree->SetBranchAddress("phfCoincFilterTh2",&phfCoincFilterTh2,&b_phfCoincFilterTh2);
    skimTree->SetBranchAddress("phfCoincFilterTh3",&phfCoincFilterTh3,&b_phfCoincFilterTh3);
    skimTree->SetBranchAddress("phfCoincFilterTh4",&phfCoincFilterTh4,&b_phfCoincFilterTh4);
    skimTree->SetBranchAddress("phfCoincFilterTh5",&phfCoincFilterTh5,&b_phfCoincFilterTh5);
    skimTree->SetBranchAddress("phfCoincFilterTh1p5",&phfCoincFilterTh1p5,&b_phfCoincFilterTh1p5);
    skimTree->SetBranchAddress("phfCoincFilterTh2p5",&phfCoincFilterTh2p5,&b_phfCoincFilterTh2p5);
    skimTree->SetBranchAddress("phfCoincFilterTh3p5",&phfCoincFilterTh3p5,&b_phfCoincFilterTh3p5);
    skimTree->SetBranchAddress("phfCoincFilterTh4p5",&phfCoincFilterTh4p5,&b_phfCoincFilterTh4p5);
    skimTree->SetBranchAddress("phfCoincFilter2Th1",&phfCoincFilter2Th1,&b_phfCoincFilter2Th1);
    skimTree->SetBranchAddress("phfCoincFilter2Th2",&phfCoincFilter2Th2,&b_phfCoincFilter2Th2);
    skimTree->SetBranchAddress("phfCoincFilter2Th3",&phfCoincFilter2Th3,&b_phfCoincFilter2Th3);
    skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4,&b_phfCoincFilter2Th4);
    skimTree->SetBranchAddress("phfCoincFilter2Th5",&phfCoincFilter2Th5,&b_phfCoincFilter2Th5);
    skimTree->SetBranchAddress("phfCoincFilter2Th1p5",&phfCoincFilter2Th1p5,&b_phfCoincFilter2Th1p5);
    skimTree->SetBranchAddress("phfCoincFilter2Th2p5",&phfCoincFilter2Th2p5,&b_phfCoincFilter2Th2p5);
    skimTree->SetBranchAddress("phfCoincFilter2Th3p5",&phfCoincFilter2Th3p5,&b_phfCoincFilter2Th3p5);
    skimTree->SetBranchAddress("phfCoincFilter2Th4p5",&phfCoincFilter2Th4p5,&b_phfCoincFilter2Th4p5);
    skimTree->SetBranchAddress("phfCoincFilter3Th1",&phfCoincFilter3Th1,&b_phfCoincFilter3Th1);
    skimTree->SetBranchAddress("phfCoincFilter3Th2",&phfCoincFilter3Th2,&b_phfCoincFilter3Th2);
    skimTree->SetBranchAddress("phfCoincFilter3Th3",&phfCoincFilter3Th3,&b_phfCoincFilter3Th3);
    skimTree->SetBranchAddress("phfCoincFilter3Th4",&phfCoincFilter3Th4,&b_phfCoincFilter3Th4);
    skimTree->SetBranchAddress("phfCoincFilter3Th5",&phfCoincFilter3Th5,&b_phfCoincFilter3Th5);
    skimTree->SetBranchAddress("phfCoincFilter3Th1p5",&phfCoincFilter3Th1p5,&b_phfCoincFilter3Th1p5);
    skimTree->SetBranchAddress("phfCoincFilter3Th2p5",&phfCoincFilter3Th2p5,&b_phfCoincFilter3Th2p5);
    skimTree->SetBranchAddress("phfCoincFilter3Th3p5",&phfCoincFilter3Th3p5,&b_phfCoincFilter3Th3p5);
    skimTree->SetBranchAddress("phfCoincFilter3Th4p5",&phfCoincFilter3Th4p5,&b_phfCoincFilter3Th4p5);
    skimTree->SetBranchAddress("phfCoincFilter4Th1",&phfCoincFilter4Th1,&b_phfCoincFilter4Th1);
    skimTree->SetBranchAddress("phfCoincFilter4Th2",&phfCoincFilter4Th2,&b_phfCoincFilter4Th2);
    skimTree->SetBranchAddress("phfCoincFilter4Th3",&phfCoincFilter4Th3,&b_phfCoincFilter4Th3);
    skimTree->SetBranchAddress("phfCoincFilter4Th4",&phfCoincFilter4Th4,&b_phfCoincFilter4Th4);
    skimTree->SetBranchAddress("phfCoincFilter4Th5",&phfCoincFilter4Th5,&b_phfCoincFilter4Th5);
    skimTree->SetBranchAddress("phfCoincFilter4Th1p5",&phfCoincFilter4Th1p5,&b_phfCoincFilter4Th1p5);
    //skimTree->SetBranchAddress("phfCoincFilter4Th2p5",&phfCoincFilter4Th2p5,&b_phfCoincFilter4Th2p5);
    skimTree->SetBranchAddress("phfCoincFilter4Th3p5",&phfCoincFilter4Th3p5,&b_phfCoincFilter4Th3p5);
    skimTree->SetBranchAddress("phfCoincFilter4Th4p5",&phfCoincFilter4Th4p5,&b_phfCoincFilter4Th4p5);
    skimTree->SetBranchAddress("phfCoincFilter5Th1",&phfCoincFilter5Th1,&b_phfCoincFilter5Th1);
    skimTree->SetBranchAddress("phfCoincFilter5Th2",&phfCoincFilter5Th2,&b_phfCoincFilter5Th2);
    skimTree->SetBranchAddress("phfCoincFilter5Th3",&phfCoincFilter5Th3,&b_phfCoincFilter5Th3);
    skimTree->SetBranchAddress("phfCoincFilter5Th4",&phfCoincFilter5Th4,&b_phfCoincFilter5Th4);
    skimTree->SetBranchAddress("phfCoincFilter5Th5",&phfCoincFilter5Th5,&b_phfCoincFilter5Th5);
    skimTree->SetBranchAddress("phfCoincFilter5Th1p5",&phfCoincFilter5Th1p5,&b_phfCoincFilter5Th1p5);
    skimTree->SetBranchAddress("phfCoincFilter5Th2p5",&phfCoincFilter5Th2p5,&b_phfCoincFilter5Th2p5);
    skimTree->SetBranchAddress("phfCoincFilter5Th3p5",&phfCoincFilter5Th3p5,&b_phfCoincFilter5Th3p5);
    skimTree->SetBranchAddress("phfCoincFilter5Th4p5",&phfCoincFilter5Th4p5,&b_phfCoincFilter5Th4p5);
    skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter,&b_pclusterCompatibilityFilter);

	Int_t          HLT_HIZeroBias_v1;
	TBranch 		*b_HLT_HIZeroBias_v1;

 	hltTree->SetBranchAddress("HLT_HIZeroBias_v1", &HLT_HIZeroBias_v1, &b_HLT_HIZeroBias_v1);

	TTree *evttree = new TTree("evttree","");
    evttree->Branch("evt",&evt,"evt/I");
    evttree->Branch("hiNpix",&hiNpix,"hiNpix/I");
    evttree->Branch("hiNtracks",&hiNtracks,"hiNtracks/I");
    evttree->Branch("hiBin",&hiBin,"hiBin/I");
    evttree->Branch("hiHF",&hiHF,"hiHF/F");

	TTree *hlttree = new TTree("hlttree","");
 	hlttree->Branch("HLT_HIZeroBias_v1", &HLT_HIZeroBias_v1, "HLT_HIZeroBias_v1/I");

	TTree *skimtree = new TTree("skimtree","");
	skimtree->Branch("pprimaryVertexFilter",&pprimaryVertexFilter,"pprimaryVertexFilter/I");
    skimtree->Branch("phfCoincFilterTh1",&phfCoincFilterTh1,"phfCoincFilterTh1/I");
    skimtree->Branch("phfCoincFilterTh2",&phfCoincFilterTh2,"phfCoincFilterTh2/I");
    skimtree->Branch("phfCoincFilterTh3",&phfCoincFilterTh3,"phfCoincFilterTh3/I");
    skimtree->Branch("phfCoincFilterTh4",&phfCoincFilterTh4,"phfCoincFilterTh4/I");
    skimtree->Branch("phfCoincFilterTh5",&phfCoincFilterTh5,"phfCoincFilterTh5/I");
    skimtree->Branch("phfCoincFilter2Th1",&phfCoincFilter2Th1,"phfCoincFilter2Th1/I");
    skimtree->Branch("phfCoincFilter2Th2",&phfCoincFilter2Th2,"phfCoincFilter2Th2/I");
    skimtree->Branch("phfCoincFilter2Th3",&phfCoincFilter2Th3,"phfCoincFilter2Th3/I");
    skimtree->Branch("phfCoincFilter2Th4",&phfCoincFilter2Th4,"phfCoincFilter2Th4/I");
    skimtree->Branch("phfCoincFilter2Th5",&phfCoincFilter2Th5,"phfCoincFilter2Th5/I");
    skimtree->Branch("phfCoincFilter3Th1",&phfCoincFilter3Th1,"phfCoincFilter3Th1/I");
    skimtree->Branch("phfCoincFilter3Th2",&phfCoincFilter3Th2,"phfCoincFilter3Th2/I");
    skimtree->Branch("phfCoincFilter3Th3",&phfCoincFilter3Th3,"phfCoincFilter3Th3/I");
    skimtree->Branch("phfCoincFilter3Th4",&phfCoincFilter3Th4,"phfCoincFilter3Th4/I");
    skimtree->Branch("phfCoincFilter3Th5",&phfCoincFilter3Th5,"phfCoincFilter3Th5/I");
    skimtree->Branch("phfCoincFilter4Th1",&phfCoincFilter4Th1,"phfCoincFilter4Th1/I");
    skimtree->Branch("phfCoincFilter4Th2",&phfCoincFilter4Th2,"phfCoincFilter4Th2/I");
    skimtree->Branch("phfCoincFilter4Th3",&phfCoincFilter4Th3,"phfCoincFilter4Th3/I");
    skimtree->Branch("phfCoincFilter4Th4",&phfCoincFilter4Th4,"phfCoincFilter4Th4/I");
    skimtree->Branch("phfCoincFilter4Th5",&phfCoincFilter4Th5,"phfCoincFilter4Th5/I");
    skimtree->Branch("phfCoincFilter5Th1",&phfCoincFilter5Th1,"phfCoincFilter5Th1/I");
    skimtree->Branch("phfCoincFilter5Th2",&phfCoincFilter5Th2,"phfCoincFilter5Th2/I");
    skimtree->Branch("phfCoincFilter5Th3",&phfCoincFilter5Th3,"phfCoincFilter5Th3/I");
    skimtree->Branch("phfCoincFilter5Th4",&phfCoincFilter5Th4,"phfCoincFilter5Th4/I");
    skimtree->Branch("phfCoincFilter5Th5",&phfCoincFilter5Th5,"phfCoincFilter5Th5/I");
    skimtree->Branch("phfCoincFilterTh1p5",&phfCoincFilterTh1p5,"phfCoincFilterTh1p5/I");
    skimtree->Branch("phfCoincFilterTh2p5",&phfCoincFilterTh2p5,"phfCoincFilterTh2p5/I");
    skimtree->Branch("phfCoincFilterTh3p5",&phfCoincFilterTh3p5,"phfCoincFilterTh3p5/I");
    skimtree->Branch("phfCoincFilterTh4p5",&phfCoincFilterTh4p5,"phfCoincFilterTh4p5/I");
    skimtree->Branch("phfCoincFilter2Th1p5",&phfCoincFilter2Th1p5,"phfCoincFilter2Th1p5/I");
    skimtree->Branch("phfCoincFilter2Th2p5",&phfCoincFilter2Th2p5,"phfCoincFilter2Th2p5/I");
    skimtree->Branch("phfCoincFilter2Th3p5",&phfCoincFilter2Th3p5,"phfCoincFilter2Th3p5/I");
    skimtree->Branch("phfCoincFilter2Th4p5",&phfCoincFilter2Th4p5,"phfCoincFilter2Th4p5/I");
    skimtree->Branch("phfCoincFilter3Th1p5",&phfCoincFilter3Th1p5,"phfCoincFilter3Th1p5/I");
    skimtree->Branch("phfCoincFilter3Th2p5",&phfCoincFilter3Th2p5,"phfCoincFilter3Th2p5/I");
    skimtree->Branch("phfCoincFilter3Th3p5",&phfCoincFilter3Th3p5,"phfCoincFilter3Th3p5/I");
    skimtree->Branch("phfCoincFilter3Th4p5",&phfCoincFilter3Th4p5,"phfCoincFilter3Th4p5/I");
    skimtree->Branch("phfCoincFilter4Th1p5",&phfCoincFilter4Th1p5,"phfCoincFilter4Th1p5/I");
    skimtree->Branch("phfCoincFilter4Th3p5",&phfCoincFilter4Th3p5,"phfCoincFilter4Th3p5/I");
    skimtree->Branch("phfCoincFilter4Th4p5",&phfCoincFilter4Th4p5,"phfCoincFilter4Th4p5/I");
    skimtree->Branch("phfCoincFilter5Th1p5",&phfCoincFilter5Th1p5,"phfCoincFilter5Th1p5/I");
    skimtree->Branch("phfCoincFilter5Th2p5",&phfCoincFilter5Th2p5,"phfCoincFilter5Th2p5/I");
    skimtree->Branch("phfCoincFilter5Th3p5",&phfCoincFilter5Th3p5,"phfCoincFilter5Th3p5/I");
    skimtree->Branch("phfCoincFilter5Th4p5",&phfCoincFilter5Th4p5,"phfCoincFilter5Th4p5/I");
    skimtree->Branch("pclusterCompatibilityFilter",&pclusterCompatibilityFilter,"pclusterCompatibilityFilter/I");
    //skimtree->SetBranchAddress("phfCoincFilter4Th2p5",&phfCoincFilter4Th2p5,&b_phfCoincFilter4Th2p5);

	int nevt1 = evtTree->GetEntries();
	int nevt2 = hltTree->GetEntries();
	int nevt3 = skimTree->GetEntries();

	cout << "nevt = " << nevt1 << ", " << nevt2 << ", " << nevt3 << endl;
	for(int iev=0; iev<nevt1 ; ++iev)     {
		if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << evtTree->GetEntries() <<  " ("<<(int)(100.*iev/evtTree->GetEntries()) << "%)" << endl;
		evtTree->GetEntry(iev);
		hltTree->GetEntry(iev);
		skimTree->GetEntry(iev);

		evttree->Fill();
		hlttree->Fill();
		skimtree->Fill();
	}
	cout << " " << endl;
	cout << "Writing Files ... " << endl;

	outFile->cd();
	evttree->Write();
	hlttree->Write();
	skimtree->Write();

	outFile->Close();

	t -> Stop();
	printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());

	cout << "Merge is DONE" << endl;
}
