#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <cmath>
#include <iomanip>

using namespace std;

void drawCanvas(int nbins, TH2D *hist[13], const char *species, const char *cname, double binEdges);

void drawnSigma()
{

    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.5);

    const int xnBins = 13;
    Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

    TFile *fin_true = new TFile("hist4CombBkg.root", "R");
    TFile *fin_det = new TFile("../PID/hist4CombBkg.root", "R");
    // TFile *fin_data = new TFile("../EmbeddingComparison/Data/histData.root", "R");
    TFile *fin_data = new TFile("histData.root", "R");

    TH2D *hPiPi_det[13];
    TH2D *hPiPi_data[13];

    TH2D *hPiPi[13];
    TH2D *hPiP[13];
    TH2D *hPiK[13];
    TH2D *hPiE[13];
    TH2D *hPP[13];
    TH2D *hPK[13];
    TH2D *hPE[13];
    TH2D *hKK[13];
    TH2D *hKE[13];

    for (int i = 0; i < 13; i++)
    {
        hPiPi_data[i] = (TH2D *)fin_data->Get(Form("hnsigmaPionPiPi_bin%i", i));

        hPiPi_det[i] = (TH2D *)fin_det->Get(Form("hnsigmaPionPiPi_bin%i", i));

        hPiPi[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPiPi_bin%i", i));
        hPiP[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPiP_bin%i", i));
        hPiK[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPiK_bin%i", i));
        hPiE[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPiE_bin%i", i));
        hPP[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPP_bin%i", i));
        hPK[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPK_bin%i", i));
        hPE[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionPE_bin%i", i));
        hKK[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionKK_bin%i", i));
        hKE[i] = (TH2D *)fin_true->Get(Form("hnsigmaPionKE_bin%i", i));
    }

    drawCanvas(13, hPiPi, "PiPi", "cPiPiTrue", xnBinsEdges);
    drawCanvas(13, hPiP, "PiP", "cPiPTrue", xnBinsEdges);
    drawCanvas(13, hPiK, "PiK", "cPiKTrue", xnBinsEdges);
    drawCanvas(13, hPiE, "PiE", "cPiETrue", xnBinsEdges);
    drawCanvas(13, hPP, "PP", "cPPTrue", xnBinsEdges);
    drawCanvas(13, hPK, "PK", "cPKTrue", xnBinsEdges);
    drawCanvas(13, hPE, "PE", "cPETrue", xnBinsEdges);
    drawCanvas(13, hPE, "KK", "cKKTrue", xnBinsEdges);
    drawCanvas(13, hPE, "KE", "cKETrue", xnBinsEdges);
    drawCanvas(13, hPiPi_data, "PiPi", "cPiPiData", xnBinsEdges);
    drawCanvas(13, hPiPi_det, "PiPi", "cPiPiDet", xnBinsEdges);
}

void drawCanvas(int nbins, TH2D *hist[13], const char *species, const char *cname, double binEdges[14])
{
    TCanvas *cc = new TCanvas(Form("%s", cname), "", 900, 600);
    cc->Divide(5, 3);

    for (int i = 0; i < nbins; i++)
    {
        if (i < 13)
        {
            cc->cd(i + 1);

            gPad->SetGrid(0, 0);
            gPad->SetRightMargin(0.15);

            gPad->SetBottomMargin(0.15);
            gPad->SetLeftMargin(0.15);
            gPad->SetLogz();

            hist[i]->SetTitle(Form("M_{inv}::[%g - %g) GeV/c^{2}", binEdges[i], binEdges[i + 1]));
            hist[i]->SetTitleSize(0.08);
            hist[i]->GetYaxis()->SetTitle("n#sigma#pi (+)");
            hist[i]->GetYaxis()->SetTitleSize(0.08);
            hist[i]->GetYaxis()->CenterTitle();
            hist[i]->GetYaxis()->SetLabelSize(0.06);

            hist[i]->GetZaxis()->SetLabelSize(0.06);

            hist[i]->GetXaxis()->SetTitle("n#sigma#pi (-)");
            hist[i]->GetXaxis()->SetTitleSize(0.10);
            hist[i]->GetXaxis()->CenterTitle();
            hist[i]->GetXaxis()->SetTitleOffset(0.6);

            hist[i]->GetXaxis()->SetLabelSize(0.06);
            hist[i]->Draw("colz");
        }

        if (i >= 13)

            continue;
    }
    cc->Update();
    cc->SaveAs(Form("Plots/nSigma_%s.pdf", cname));
}
