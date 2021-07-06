#include "TFile.h"
#include "TH1.h"

#include "iostream"

void test()
{
	using namespace std;

	TFile *f1 = new TFile("./roots/hist_HiForward_0000.root","read");
	TFile *f2 = new TFile("./roots/hist_HiForward_0001.root","read");
	TFile *f3 = new TFile("./roots/hist_HiForward_0002.root","read");
	TFile *f4 = new TFile("./roots/hist_HiForward_0003.root","read");
	TFile *f5 = new TFile("./roots/hist_HiForward_0004.root","read");
	TFile *f6 = new TFile("./roots/hist_HiForward_0005.root","read");
	TFile *f7 = new TFile("./roots/hist_HiForward_0006.root","read");
	TFile *f8 = new TFile("./roots/hist_HiForward_0007.root","read");
	TFile *f9 = new TFile("./roots/hist_HiForward_0008.root","read");

	TFile *f = new TFile("./roots/hist_HiForward.root","read");

	TH1F *h1 = (TH1F*)f1->Get("hHF_P_8");
	TH1F *h2 = (TH1F*)f2->Get("hHF_P_8");
	TH1F *h3 = (TH1F*)f3->Get("hHF_P_8");
	TH1F *h4 = (TH1F*)f4->Get("hHF_P_8");
	TH1F *h5 = (TH1F*)f5->Get("hHF_P_8");
	TH1F *h6 = (TH1F*)f6->Get("hHF_P_8");
	TH1F *h7 = (TH1F*)f7->Get("hHF_P_8");
	TH1F *h8 = (TH1F*)f8->Get("hHF_P_8");
	TH1F *h9 = (TH1F*)f9->Get("hHF_P_8");

	TH1F *h = (TH1F*)f->Get("hHF_P_8");

	double sum = ( h1->GetEntries()+h2->GetEntries()+h3->GetEntries()+h4->GetEntries()+h5->GetEntries()+h6->GetEntries()+h7->GetEntries()+h8->GetEntries() );

	cout << "h1 : " << h1->GetEntries() << endl;
	cout << "h_i : " << sum << endl;
	cout << "h : " << h->GetEntries() << endl;
}
