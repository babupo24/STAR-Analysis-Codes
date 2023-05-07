#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TMath.h>

using namespace std;

void datembSys()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);

    TFile *fin = new TFile("ResultsMinConeCut/fxsecOutput.root", "R");
    TH1D *hdat = (TH1D *)fin->Get("hdatComb");
    TH1D *hemb = (TH1D *)fin->Get("hembComb");

    TCanvas *can = new TCanvas("can", "", 800, 400);
    can->Divide(2, 1);
    can->cd(1);
    gPad->SetLogy();
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    TH1D *hdatc = (TH1D *)hdat->Clone();
    hdatc->SetName("Data");
    hdatc->SetLineColor(2);
    hdatc->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
    hdatc->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
    hdatc->GetYaxis()->SetTitleOffset(1.4);
    hdatc->Draw("hist E");
    TH1D *hembc = (TH1D *)hemb->Clone();
    hembc->SetName("Embed");
    hembc->SetLineColor(4);
    hembc->Draw("hist E same");
    TLegend *leg = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg->AddEntry(hdatc, "Data", "le");
    leg->AddEntry(hembc, "Embed", "le");
    leg->Draw();

    can->cd(2);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    TH1D *hratio = (TH1D *)hemb->Clone();
    hratio->SetName("Ratio");
    hratio->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hratio->GetYaxis()->SetTitle("(Emb - Dat)/Dat");
    hratio->GetYaxis()->SetTitleOffset(1.4);
    hratio->Add(hdat, -1);
    hratio->Divide(hdat);
    hratio->GetYaxis()->SetRangeUser(-0.41, 0.41);
    hratio->SetLineColor(2);
    hratio->Draw("E1");
    gPad->Update();
    TLine *lin = new TLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);
    lin->SetLineStyle(2);
    lin->Draw();

    gPad->Update();
    can->SaveAs("ResultsMinConeCut/plot_dataembsys.pdf");
    can->SaveAs("ResultsMinConeCut/plot_dataembsys.png");
}