#include <iostream>
#include <TFile.h>
#include <cmath>
#include <TH1.h>

using namespace std;

void getTriggerBias()
{

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    // open parton level root file
    // open particle level root file
    TFile *fin = new TFile("hist4TrigEffGoodRuns_v2.root");
    //TFile *fin = new TFile("testEmbedHist.root");

    // get parton level hists
    TH1D *hQuarkPar = (TH1D *)fin->Get("hQuarksPar");
    TH1D *hGluonPar = (TH1D *)fin->Get("hGluonsPar");
    TH1D *hPartonPar = (TH1D *)fin->Get("hPartonsPar");
    // get detector level hists
    TH1D *hQuarkDet = (TH1D *)fin->Get("hQuarksDet");
    TH1D *hGluonDet = (TH1D *)fin->Get("hGluonsDet");
    TH1D *hPartonDet = (TH1D *)fin->Get("hDettonsDet");

    // quark fraction at particle level
    TH1D *hQ_frac_par = (TH1D *)hQuarkPar->Clone();
    hQ_frac_par->Divide(hPartonPar);
    TH1D *hG_frac_par = (TH1D *)hGluonPar->Clone();
    hG_frac_par->Divide(hPartonPar);

    // quark fraction at particle level for jp triggers
    TH1D *hQ_frac_det = (TH1D *)hQuarkDet->Clone();
    hQ_frac_det->Divide(hPartonDet);
    TH1D *hG_frac_det = (TH1D *)hGluonDet->Clone();
    hG_frac_det->Divide(hPartonDet);

    // trigger bias
    TH1D *hQBias = (TH1D *)hQ_frac_det->Clone();
    hQBias->Divide(hQ_frac_par);

    // draw ratios for parton level
    TCanvas *crp = new TCanvas("cratios", "", 900, 450);
    crp->Divide(2, 1);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.20);
    crp->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hQ_frac_par->SetLineColor(2);
    hQ_frac_par->SetLineWidth(2);
    hQ_frac_par->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hQ_frac_par->GetYaxis()->SetTitle("fractions");
    hQ_frac_par->GetYaxis()->SetTitleOffset(1.5);
    hQ_frac_par->GetYaxis()->SetRangeUser(0.0, 1.0);
    hQ_frac_par->Draw("HIST E");

    hG_frac_par->SetLineColor(4);
    hG_frac_par->SetLineWidth(2);
    hG_frac_par->SetLineStyle(1);
    hG_frac_par->Draw("HIST E SAME");

    hG_frac_det->SetLineColor(4);
    hG_frac_det->SetLineWidth(2);
    hG_frac_det->SetLineStyle(2);
    hG_frac_det->Draw("HIST E SAME");

    hQ_frac_det->SetLineColor(2);
    hQ_frac_det->SetLineWidth(2);
    hQ_frac_det->SetLineStyle(2);
    hQ_frac_det->Draw("HIST E SAME");

    TLegend *lg = new TLegend(0.3, 0.75, 0.8, 0.90);
    lg->SetNColumns(2);
    lg->AddEntry(hQ_frac_par, " Quark (particle)", "l");
    lg->AddEntry(hQ_frac_det, " Quark (detector)", "l");
    lg->AddEntry(hG_frac_par, " Gluon (particle)", "l");
    lg->AddEntry(hG_frac_det, " Gluon (detector)", "l");
    lg->Draw();
    gPad->Update();
    crp->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    hQBias->SetLineColor(2);
    hQBias->SetMarkerStyle(8);
    hQBias->SetMarkerColor(2);
    hQBias->SetLineWidth(2);
    hQBias->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hQBias->GetYaxis()->SetRangeUser(0.5, 1.5);
    hQBias->GetYaxis()->SetTitle("Trigger Bias");
    hQBias->GetYaxis()->SetTitleOffset(1.5);
    hQBias->Draw("E");
    hQBias->Fit("pol0");

    crp->Update();
    crp->SaveAs("triggerBias.pdf");
}
