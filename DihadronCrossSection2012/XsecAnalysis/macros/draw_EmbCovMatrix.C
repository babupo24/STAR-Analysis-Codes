#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace std;
void draw_hist(TH2D *hist);
void draw_EmbCovMatrix()
{

    TFile *fin = new TFile("../BinResolutionAndUnfolding/UnfoldingResultsMinConeCut/unfoldingOutput.root", "r");
    TH2D *hcov_jp0 = (TH2D *)fin->Get("hCovMatrixJP0tot");
    TH2D *hcov_jp1 = (TH2D *)fin->Get("hCovMatrixJP1tot");
    TH2D *hcov_jp2 = (TH2D *)fin->Get("hCovMatrixJP2tot");

    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    TCanvas *can = new TCanvas("can", "", 900, 350);
    can->Divide(3, 1);
    can->cd(1);
    draw_hist(hcov_jp0);
    can->cd(2);
    draw_hist(hcov_jp1);
    can->cd(3);
    draw_hist(hcov_jp2);
    can->Update();
    can->SaveAs("ResultsMinConeCut/unfold_cov_matrix.pdf");
}
void draw_hist(TH2D *hist)
{
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogz();
    hist->GetXaxis()->SetTitle("Minv (det) (GeV/c^{2})");
    hist->GetYaxis()->SetTitle("Minv (gen) (GeV/c^{2})");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->Draw("colz");
}