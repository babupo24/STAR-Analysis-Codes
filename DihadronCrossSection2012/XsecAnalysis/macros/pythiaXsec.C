#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iomanip>

using namespace std;

void pythiaXsec()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    // gStyle->SetLabelFont(13, "xyz");
    gStyle->SetLegendBorderSize(0);

    const int nBins = 13;
    Double_t nBinsEdges[nBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};
    // binning for cross-section money plot. Add bins at underflow and overflow
    const int nBinsU = 15;
    Double_t nBinsEdgesU[nBinsU + 1] = {0.05, 0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0, 4.25};

    // measured cross section
    double xsec_mes[13] = {7.10334e+09, 8.71813e+09, 7.59875e+09, 2.75747e+09, 7.40509e+08, 1.41126e+08, 3.77774e+07, 8.4353e+06, 1.3966e+06, 359442, 76575.7, 7387.86, 616.115};
    double xsec_staterr_mes[13] = {48659.9, 53907.7, 50328.1, 30317.6, 15711, 6858.72, 3548.59, 1676.83, 682.301, 346.141, 159.766, 49.6248, 14.3308};

    // partonic pt parameters
    double partPtRange[nBins + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., -1};

    Double_t crossSection[nBins] = {9.0012e+00, 1.46253e+00, 3.54466e-01, 1.51622e-01, 2.49062e-02, 5.84527e-03, 2.30158e-03, 3.42755e-04, 4.57002e-05, 9.72535e-06, 4.69889e-07, 2.69202e-08, 1.43453e-09}; // in mb

    Double_t numEvents[nBins] = {3318626, 3301413, 3291662, 3280010, 3282543, 3275693, 3276437, 3276795, 3272804, 2179660, 2183230, 1091927, 1090857}; // number of events in file after processing MuDsts

    Double_t ptBinLumiMb[13] = {1.67958e+09, 3.20524e+08, 8.07797e+07, 3.51516e+07, 5.76973e+06, 1.35694e+06, 534174, 79541.3, 10618.3, 3392.93, 163.664, 18.7475, 1.0};

    double numGen[13] = {3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 2458337, 2459112, 1229556, 1229556};

    double numTried[13] = {17681975, 14266371, 15245667, 17945677, 19050998, 19385968, 22015921, 25182153, 30278209, 22632389, 23503610, 12736750, 12487743};

    // pythia luminosity for partonic pT bins

    double ptBinLumi[13] = {0};
    for (int i = 0; i < 13; i++)
    {
        ptBinLumi[i] = (numEvents[i] / (crossSection[i] * 1e9));
    }

    TFile *frawX = new TFile("pythiaXsecAll.root", "R");
    // TFile *frawX = new TFile("PythiaXsec/test.root", "R");

    // get raw cross section in each pt bins
    TH1D *hrawX[13];
    TH1D *hevProcess[13];
    double lum[13] = {0};
    double evProcessed[13] = {0};

    ofstream foutpyth;
    foutpyth.open("fpythia_events_lum.tex");
    foutpyth << "\\hline \\hline" << endl;
    foutpyth << "$p_T$ Bin (GeV/c) & $N_{p_T}^{tried}$ & $N_{p_T}^{gen}$ & $N_{p_T}^{pro}$ & $\\sigma_{p_T}$ (pb) & L ($pb^{-1}$) \\\\" << endl;

    for (int i = 0; i < 13; i++)
    {
        hrawX[i] = (TH1D *)frawX->Get(Form("hpythiaMinv_%g_%g", partPtRange[i], partPtRange[i + 1]));
        hevProcess[i] = (TH1D *)frawX->Get(Form("hevProcessed_%g_%g", partPtRange[i], partPtRange[i + 1]));
        evProcessed[i] = hevProcess[i]->Integral();
        lum[i] = evProcessed[i] / (crossSection[i] * 1e9);
        hrawX[i]->Scale(1. / lum[i]);
        foutpyth << partPtRange[i] << "$-$" << partPtRange[i + 1] << "&" << setprecision(3) << numTried[i] << "&" << setprecision(3) << numGen[i] << "&" << setprecision(3) << evProcessed[i] << "&" << setprecision(3) << crossSection[i] * 1e9 << "&" << setprecision(3) << lum[i] << "\\\\" << endl;

        for (int j = 1; j <= hrawX[i]->GetNbinsX(); j++)
        {
            // double
        }
    }
    foutpyth << "\\hline \\hline" << endl;

    TH1D *hpythiaXsec = (TH1D *)hrawX[0]->Clone();
    for (int i = 1; i < 13; i++)
    {
        hpythiaXsec->Add(hrawX[i]);
    }

    TH1D *hpythiaXsecW = (TH1D *)hpythiaXsec->Clone();
    hpythiaXsecW->Reset();
    TH1D *hmesXsec = (TH1D *)hpythiaXsec->Clone();
    hmesXsec->Reset();
    for (int i = 1; i <= 13; i++)
    {
        double nn = hpythiaXsec->GetBinContent(i);
        double bw = hpythiaXsec->GetBinWidth(i);
        double ne = hpythiaXsec->GetBinError(i);
        hpythiaXsecW->SetBinContent(i, nn / bw);
        hpythiaXsecW->SetBinError(i, ne);

        hmesXsec->SetBinContent(i, xsec_mes[i - 1]);
        hmesXsec->SetBinError(i, xsec_staterr_mes[i - 1]);
    }

    // hpythiaXsecW->SetMaximum(5e12);
    hpythiaXsecW->SetName("PYTHIA");

    TCanvas *c0 = new TCanvas("c0", "", 500, 400);
    c0->cd();
    hpythiaXsecW->Draw("hist E ");
    hmesXsec->SetLineColor(2);
    hmesXsec->SetName("Measured");
    hmesXsec->Draw("hist E SAME");
    c0->SetLogy();
    c0->BuildLegend();
    // c0->SaveAs("ResultsMinConeCut/pythia_mes_xsec.pdf");
    // c0->Delete();
    TCanvas *c1 = new TCanvas("c1", "", 500, 400);
    c1->cd();
    TH1D *hratio = (TH1D *)hmesXsec->Clone();
    hratio->Add(hpythiaXsecW, -1);
    hratio->Divide(hpythiaXsecW);
    hratio->GetYaxis()->SetRangeUser(-2.5, 2.5);
    hratio->Draw();
    c1->Update();
}