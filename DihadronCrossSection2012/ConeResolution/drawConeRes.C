#include <iostream>
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

using namespace std;

void draw(int numBins, TH1D **hist, const char *hname);

void drawConeRes()
{
    // ROOT->Reset("Default");
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetLegendBorderSize(0);
    TFile *ifile = new TFile("hist4ConeRes.root");
    assert(ifile);

    TH1D *hc7[12];
    TH1D *hc6[12];
    TH1D *hc5[12];
    TH1D *hc4[12];
    TH1D *hc3[12];

    TH1D *hcres = (TH1D *)ifile->Get("hConeResInt");

    for (int i = 0; i < 12; i++)
    {
        hc7[i] = (TH1D *)ifile->Get(Form("hConeResLt7%i", i));
        hc6[i] = (TH1D *)ifile->Get(Form("hConeResLt6%i", i));
        hc5[i] = (TH1D *)ifile->Get(Form("hConeResLt5%i", i));
        hc4[i] = (TH1D *)ifile->Get(Form("hConeResLt4%i", i));
        hc3[i] = (TH1D *)ifile->Get(Form("hConeResLt3%i", i));
    }
    TH1D *hc7all;
    TH1D *hc6all;
    TH1D *hc5all;
    TH1D *hc4all;
    TH1D *hc3all;
    for (int i = 0; i < 12; i++)
    {
        if (i == 0)
        {
            hc7all = (TH1D *)hc7[i]->Clone();
            hc6all = (TH1D *)hc6[i]->Clone();
            hc5all = (TH1D *)hc5[i]->Clone();
            hc4all = (TH1D *)hc4[i]->Clone();
            hc3all = (TH1D *)hc3[i]->Clone();
        }
        else
        {
            hc7all->Add(hc7[i]);
            hc6all->Add(hc6[i]);
            hc5all->Add(hc5[i]);
            hc4all->Add(hc4[i]);
            hc3all->Add(hc3[i]);
        }
    }

    TCanvas *canres = new TCanvas("canres", "", 600, 450);
    canres->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogz();
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    hcres->GetXaxis()->SetTitle("Cone (Part - Det.)");
    hcres->GetYaxis()->SetTitle("Cone (particle) ");
    hcres->GetXaxis()->SetTitle("Cone (detector) ");
    hcres->Draw("colz");
    canres->Update();
    canres->SaveAs("Plots/coneresolution2d.pdf");

    TCanvas *canall = new TCanvas("canall", "", 600, 450);
    canall->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogy();
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    hc7all->GetXaxis()->SetTitle("Cone (Part - Det.)");
    hc7all->GetYaxis()->SetTitle("Dipion Yields");
    hc7all->SetLineColor(1);
    hc7all->Draw();
    hc7all->Fit("gaus");
    hc7all->GetFunction("gaus")->SetLineColorAlpha(1, 0.5);

    gPad->Update();
    hc6all->SetLineColorAlpha(2, 0.8);
    hc6all->Draw(" same");
    hc6all->Fit("gaus");
    hc6all->GetFunction("gaus")->SetLineColorAlpha(2, 0.5);
    gPad->Update();
    hc5all->SetLineColorAlpha(3, 0.8);
    hc5all->Draw(" same");
    hc5all->Fit("gaus");
    hc5all->GetFunction("gaus")->SetLineColorAlpha(3, 0.5);
    gPad->Update();
    hc4all->SetLineColorAlpha(4, 0.8);
    hc4all->Draw("same");
    hc4all->Fit("gaus");
    hc4all->GetFunction("gaus")->SetLineColorAlpha(4, 0.5);
    gPad->Update();
    hc3all->SetLineColorAlpha(6, 0.8);
    hc3all->Draw(" same");
    hc3all->Fit("gaus");
    hc3all->GetFunction("gaus")->SetLineColorAlpha(6, 0.5);
    gPad->Update();

    TLegend *lg1 = new TLegend(0.15, 0.5, 0.35, 0.85);
    lg1->AddEntry(hc7all, Form("Cone < 0.7, #sigma = %4.4g", hc7all->GetFunction("gaus")->GetParameter(2)), " l");
    lg1->AddEntry(hc6all, Form("Cone < 0.6, #sigma = %4.4g", hc6all->GetFunction("gaus")->GetParameter(2)), " l");
    lg1->AddEntry(hc5all, Form("Cone < 0.5, #sigma = %4.4g", hc5all->GetFunction("gaus")->GetParameter(2)), " l");
    lg1->AddEntry(hc4all, Form("Cone < 0.4, #sigma = %4.4g", hc4all->GetFunction("gaus")->GetParameter(2)), " l");
    lg1->AddEntry(hc3all, Form("Cone < 0.3, #sigma = %4.4g", hc3all->GetFunction("gaus")->GetParameter(2)), " l");
    lg1->SetTextSize(0.03);
    lg1->Draw();
    canall->Update();

    canall->SaveAs("./Plots/resolution_combined.pdf");
    delete canall;
    delete canres;
    // draw(12, hc7, "ConeLt7");
    // draw(12, hc6, "ConeLt6");
    // draw(12, hc5, "ConeLt5");
    // draw(12, hc4, "ConeLt4");
    // draw(12, hc3, "ConeLt3");

    delete ifile;

} // main

void draw(int numBins, TH1D **hist, const char *hname)
{
    TCanvas *can = new TCanvas("can", "", 600, 450);
    can->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogy();
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    for (int i = 0; i < numBins; i++)
    {
        if (i == 0)
        {
            hist[0]->GetXaxis()->SetTitle("Cone (Part. - Det.)");
            hist[0]->GetYaxis()->SetTitle("Dipion Yields)");
            hist[0]->SetLineColor(1);
            hist[0]->Draw("hist");
        }
        else
        {
            hist[i]->SetLineColor(i + 1);
            hist[i]->Draw("hist same");
        }
    }
    can->Update();
    can->SaveAs(Form("resolution_%s.pdf", hname));

    delete can;
}
