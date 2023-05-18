#include <iostream>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iomanip>
#include <cmath>

using namespace std;

void printNum(ofstream &fname, TH1D *hist, const char *type);

void get_fraction_JP1()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetLegendBorderSize(0);

    TFile *fin = new TFile("hist4CombBkg.root", "r");

    const TH1D *htrueKE_JP1 = (TH1D *)fin->Get("htrueKaonElectronPairJP1");
    const TH1D *htrueKK_JP1 = (TH1D *)fin->Get("htrueKaonPairJP1");
    const TH1D *htruePP_JP1 = (TH1D *)fin->Get("htrueProtonPairJP1");
    const TH1D *htruePE_JP1 = (TH1D *)fin->Get("htrueProtonElectronPairJP1");
    const TH1D *htruePK_JP1 = (TH1D *)fin->Get("htrueProtonKaonPairJP1");
    const TH1D *htruePiE_JP1 = (TH1D *)fin->Get("htruePionElectronPairJP1");
    const TH1D *htruePiK_JP1 = (TH1D *)fin->Get("htruePionKaonPairJP1");
    const TH1D *htruePiP_JP1 = (TH1D *)fin->Get("htruePionProtonPairJP1");
    const TH1D *htruePi_JP1 = (TH1D *)fin->Get("htruePionPairJP1"); // det and true
    const TH1D *hrecPi_JP1 = (TH1D *)fin->Get("hrecPionPairJP1");   // only detector pid

    TH1D *hdeno = (TH1D *)hrecPi_JP1->Clone();
    TH1D *hfPiPi = (TH1D *)htruePi_JP1->Clone();
    hfPiPi->Divide(hdeno);
    TH1D *hfPiP = (TH1D *)htruePiP_JP1->Clone();
    hfPiP->Divide(hdeno);
    TH1D *hfPiK = (TH1D *)htruePiK_JP1->Clone();
    hfPiK->Divide(hdeno);
    TH1D *hfPiE = (TH1D *)htruePiE_JP1->Clone();
    hfPiE->Divide(hdeno);
    TH1D *hfPE = (TH1D *)htruePE_JP1->Clone();
    hfPE->Divide(hdeno);
    TH1D *hfPK = (TH1D *)htruePK_JP1->Clone();
    hfPK->Divide(hdeno);
    TH1D *hfPP = (TH1D *)htruePP_JP1->Clone();
    hfPP->Divide(hdeno);
    TH1D *hfKK = (TH1D *)htrueKK_JP1->Clone();
    hfKK->Divide(hdeno);
    TH1D *hfKE = (TH1D *)htrueKE_JP1->Clone();
    hfKE->Divide(hdeno);

    TCanvas *can = new TCanvas("can", "", 500, 450);
    can->cd();
    TH1D *cPiPi = (TH1D *)hfPiPi->Clone();
    cPiPi->SetMarkerStyle(8);
    cPiPi->SetMarkerColor(1);
    cPiPi->SetTitle("Trig::JP1");
    cPiPi->GetXaxis()->SetTitle("M_{inv}");
    cPiPi->GetYaxis()->SetTitle("fraction");
    cPiPi->GetYaxis()->SetRangeUser(0, 1.10);
    cPiPi->Draw("E");
    TH1D *cPiP = (TH1D *)hfPiP->Clone();
    cPiP->SetMarkerStyle(8);
    cPiP->SetMarkerColor(2);
    cPiP->SetLineColor(2);
    cPiP->Draw("E SAME");
    TH1D *cPiK = (TH1D *)hfPiK->Clone();
    cPiK->SetMarkerStyle(8);
    cPiK->SetMarkerColor(3);
    cPiK->SetLineColor(3);
    cPiK->Draw("E SAME");
    TH1D *cPiE = (TH1D *)hfPiE->Clone();
    cPiE->SetMarkerStyle(8);
    cPiE->SetMarkerColor(4);
    cPiE->SetLineColor(4);
    cPiE->Draw("E SAME");
    TH1D *cPE = (TH1D *)hfPE->Clone();
    cPE->SetMarkerStyle(8);
    cPE->SetMarkerColor(5);
    cPE->SetLineColor(5);
    cPE->Draw("E SAME");
    TH1D *cPK = (TH1D *)hfPK->Clone();
    cPK->SetMarkerStyle(8);
    cPK->SetMarkerColor(6);
    cPK->SetLineColor(6);
    cPK->Draw("E SAME");
    TH1D *cPP = (TH1D *)hfPP->Clone();
    cPP->SetMarkerStyle(8);
    cPP->SetMarkerColor(7);
    cPP->SetLineColor(7);
    cPP->Draw("E SAME");
    TH1D *cKK = (TH1D *)hfKK->Clone();
    cKK->SetMarkerStyle(8);
    cKK->SetMarkerColor(8);
    cKK->SetLineColor(8);
    cKK->Draw("E SAME");
    TH1D *cKE = (TH1D *)hfKE->Clone();
    cKE->SetMarkerStyle(8);
    cKE->SetMarkerColor(9);
    cKE->SetLineColor(9);
    cKE->Draw("E SAME");
    can->Update();

    TLegend *leg = new TLegend(0.35, 0.35, 0.7, 0.65);
    leg->SetNColumns(2);
    leg->AddEntry(cPiPi, " #pi#pi", "lep");
    leg->AddEntry(cPiP, " #pi p", "lep");
    leg->AddEntry(cPiK, " #pi K", "lep");
    leg->AddEntry(cPiE, " #pi e", "lep");
    leg->AddEntry(cPE, " pe", "lep");
    leg->AddEntry(cPK, " pK", "lep");
    leg->AddEntry(cPP, " pp", "lep");
    leg->AddEntry(cKK, " KK", "lep");
    leg->AddEntry(cKE, " KE", "lep");
    leg->SetTextSize(0.03);
    leg->Draw();
    TLatex tex;
    tex.SetTextSize(0.035);
    tex.DrawLatex(0.5, 0.25, "#font[22]{fraction = #frac{N_{#pi#pi(det) && xy(true)}}{N_{#pi#pi(det)}}, x,y = particle species}");
    can->Update();
    can->SaveAs("Plots/fractions_JP1.pdf");

    ofstream fout;
    fout.open("Plots/fractions_JP1.txt");
    printNum(fout, hfPiPi, "PiPi");
    printNum(fout, hfPiP, "PiP");
    printNum(fout, hfPiK, "PiK");
    printNum(fout, hfPiE, "PiE");
    printNum(fout, hfPP, "PP");
    printNum(fout, hfPK, "PK");
    printNum(fout, hfPE, "PE");
    printNum(fout, hfPE, "PE");
    printNum(fout, hfKK, "KK");
    printNum(fout, hfKE, "KE");
}

void printNum(ofstream &fname, TH1D *hist, const char *type)
{
    fname << Form("double fraction_%s[13]={", type);
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double nf = hist->GetBinContent(i);
        if (i == hist->GetNbinsX())
            fname << nf << "};" << endl;
        else
            fname << nf << ", ";
    }
}
