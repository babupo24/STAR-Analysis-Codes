#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <THStack.h>
#include <TH2.h>

using namespace std;

void draw_signal_bkg_stacked()
{

    // get histgram root file
    TFile *fin = new TFile("hist4ResAndUnfolding.root", "R");
    if (!fin)
        cout << "Input file doesn't exist" << endl;

    // get signal histograms
    TH1D *hsigJP0 = (TH1D *)fin->Get("hMRecJP0");
    hsigJP0->SetFillColor(kGreen);
    TH1D *hsigJP1 = (TH1D *)fin->Get("hMRecJP1");
    hsigJP1->SetFillColor(kGreen);
    TH1D *hsigJP2 = (TH1D *)fin->Get("hMRecJP2");
    hsigJP2->SetFillColor(kGreen);
    // get signal+bkg histograms
    TH1D *hsigJP0bk = (TH1D *)fin->Get("hMRecJP0All");
    hsigJP0bk->SetFillColor(4);
    TH1D *hsigJP1bk = (TH1D *)fin->Get("hMRecJP1All");
    hsigJP1bk->SetFillColor(4);
    TH1D *hsigJP2bk = (TH1D *)fin->Get("hMRecJP2All");
    hsigJP2bk->SetFillColor(4);

    // get background histograms
    TH1D *hJP0bk1 = (TH1D *)fin->Get("hMRecJP0bk1");
    TH1D *hJP0bk2 = (TH1D *)fin->Get("hMRecJP0bk2");
    TH1D *hJP1bk1 = (TH1D *)fin->Get("hMRecJP1bk1");
    TH1D *hJP1bk2 = (TH1D *)fin->Get("hMRecJP1bk2");
    TH1D *hJP2bk1 = (TH1D *)fin->Get("hMRecJP2bk1");
    TH1D *hJP2bk2 = (TH1D *)fin->Get("hMRecJP2bk2");

    // combine backgrounds
    TH1D *hJP0bk = (TH1D *)fin->Get("hMRecJP0bk1");
    hJP0bk->Add(hJP0bk2);
    hJP0bk->SetFillColor(kRed);
    TH1D *hJP1bk = (TH1D *)fin->Get("hMRecJP1bk1");
    hJP1bk->Add(hJP1bk2);
    hJP1bk->SetFillColor(kRed);
    TH1D *hJP2bk = (TH1D *)fin->Get("hMRecJP2bk1");
    hJP2bk->Add(hJP2bk2);
    hJP2bk->SetFillColor(kRed);

    TH1D *htestJP0All = (TH1D *)hsigJP0->Clone();
    htestJP0All->Add(hJP0bk);
    htestJP0All->SetFillColor(5);

    THStack *hstackJP0 = new THStack("hstackJP0", "JP0 signal and background");
    hstackJP0->Add(htestJP0All);
    hstackJP0->Add(hsigJP0bk);
    hstackJP0->Add(hsigJP0);
    hstackJP0->Add(hJP0bk);

    TCanvas *cst = new TCanvas("cst", "", 500, 500);
    cst->cd();
    cst->SetGrid(0, 0);
    cst->SetLogy();
    hstackJP0->Draw("hist nostack");
    cst->Update();
}
