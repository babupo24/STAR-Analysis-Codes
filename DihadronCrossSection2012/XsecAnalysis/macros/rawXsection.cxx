#include <iostream>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;

void rawXsection()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);

    // const int nBins = 11;
    // Double_t nBinsEdges[nBins + 1] = {0.24, 0.34, 0.44, 0.56, 0.68, 0.83, 1.0, 1.25, 1.6, 2.0, 2.75, 4.0};

    const int nBins = 16;
    // Double_t nBinsEdges[nBins + 1] = {0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.58, 0.64, 0.72, 0.82, 0.95, 1.10, 1.40, 1.90, 2.5, 3.1, 4.0};
    Double_t nBinsEdges[nBins + 1] = {0.28, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};
    TFile *fin = new TFile("hist_Xsec.root", "r");

    Double_t avgM[nBins];
    Double_t avgMrms[nBins];
    Double_t statErr[nBins];
    Double_t yields[nBins];

    TH1D *hperBin[nBins];

    for (int i = 0; i < nBins; i++)
    {
        hperBin[i] = (TH1D *)fin->Get(Form("hxecMinAvgBin_%i", i));
        avgM[i] = hperBin[i]->GetMean();
        avgMrms[i] = hperBin[i]->GetRMS();
        yields[i] = hperBin[i]->GetEntries();
        statErr[i] = sqrt(yields[i]);
    }
    TH1D *hxsecM = (TH1D *)fin->Get("hxsecMinv");
    hxsecM->Sumw2();

    Double_t lumJPtotal = 26.64; // 26.64 total;
    Double_t dipiYields[nBins] = {0};
    Double_t sigmaStat[nBins] = {0};
    Double_t xsecRaw[nBins] = {0};
    Double_t xsecSigma[nBins] = {0};
    for (int i = 0; i < nBins; i++)
    {
        dipiYields[i] = hxsecM->GetBinContent(i);
        sigmaStat[i] = hxsecM->GetBinError(i);
        xsecRaw[i] = (Double_t)dipiYields[i] / lumJPtotal;
    }
    TCanvas *cxsec = new TCanvas("cxsec", "", 600, 500);
    cxsec->cd();
    cxsec->SetGrid(0, 0);
    cxsec->SetFillColor(18);
    cxsec->SetLeftMargin(0.12);
    cxsec->SetBottomMargin(0.12);
    cxsec->SetLogy();

    hxsecM->Scale(1. / lumJPtotal); // N / L
    hxsecM->SetLineColor(2);
    hxsecM->SetLineWidth(2);
    hxsecM->GetXaxis()->SetLimits(0.0, 4.50);
    hxsecM->GetXaxis()->SetNdivisions(510);
    hxsecM->GetXaxis()->SetTitle("M^{#pi^{+}#pi^{-}}_{inv} (GeV/c^{2})");
    hxsecM->GetYaxis()->SetTitle("(d#sigma^{#pi^{+}#pi^{-}}/dM^{#pi^{+}#pi^{-}})_{raw}, L_{int} = 26.64 pb^{-1}");
    hxsecM->GetYaxis()->CenterTitle();
    hxsecM->GetYaxis()->SetTitleOffset(1.2);
    hxsecM->GetXaxis()->CenterTitle();

    hxsecM->SetMarkerStyle(20);
    hxsecM->SetMarkerColor(2);
    hxsecM->Draw("hist E");

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.03);
    tex.DrawLatex(2.6, 1e6, "#color[4]{-1 < #eta^{#pi^{+}#pi^{-}} < 1}");
    tex.DrawLatex(2.6, 6e5, "#color[4]{-1 < n#sigma_{#pi} < 2}");
    tex.DrawLatex(2.6, 3e5, "#color[4]{#sigma^{#pi^{+}#pi^{-}} = #frac{N^{#pi^{+}#pi^{-}}}{L_{int}}}, L_{int} = 26.64 pb^{-1}");
    tex.DrawLatex(2.6, 1.6e5, "#color[4]{Triggers: JP0, JP1, JP2}");
    tex.DrawLatex(2.6, 1e5, "#color[2]{#minus Bin Width}");
    tex.DrawLatex(2.6, 6e4, "#color[2]{#bullet  Bin Center}");
    tex.DrawLatex(2.6, 4e4, "#color[2]{#void8  Bin Average}");

    gPad->Update();

    TLine *l[nBins];
    for (int i = 0; i < nBins; i++)
    {
        // l[i] = new TLine(avgM[i], 1.2e6, avgM[i], 1.55e6);
        l[i] = new TLine(avgM[i], gPad->GetUymax() - 20, avgM[i], gPad->GetUymax());
        l[i]->SetLineColor(2);
        l[i]->SetLineWidth(2);
        l[i]->Draw();
        gPad->Update();
    }

    cxsec->Update();
    cxsec->SaveAs("rawXsection.pdf");

    // create histogram for  the average invariant mass per bin
}
