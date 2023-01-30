
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

using namespace std;
void getBinError(double, double, double, double);
void setBinError(TH1D * /*ratio*/, TH1D * /*num*/, TH1D * /*den*/); // calculate and set error when histograms are divided

void drawMatchingFraction()
{
    gROOT->Reset();
    gStyle->SetOptDate(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);

    TFile *fembed = new TFile("hist4ResAndUnfolding.root");
    TFile *fdata = new TFile("histData4unfolding.root");
    if (!fembed || !fdata)
        cout << "File doesn't exist." << endl;

    // histograms for the matching fraction
    // denominator (gen&rec + rec)
    TH1D *hrecJP0_All = (TH1D *)fembed->Get("hMRecJP0All");
    hrecJP0_All->Rebin(2.0);
    TH1D *hrecJP1_All = (TH1D *)fembed->Get("hMRecJP1All");
    hrecJP1_All->Rebin(2.0);
    TH1D *hrecJP2_All = (TH1D *)fembed->Get("hMRecJP2All");
    hrecJP2_All->Rebin(2.0);
    // numerator (gen&rec)
    TH1D *hrecJP0 = (TH1D *)fembed->Get("hMRecJP0");
    hrecJP0->Rebin(2.0);
    TH1D *hrecJP1 = (TH1D *)fembed->Get("hMRecJP1");
    hrecJP1->Rebin(2.0);
    TH1D *hrecJP2 = (TH1D *)fembed->Get("hMRecJP2");
    hrecJP2->Rebin(2.0);

    // matching ratio
    TH1D *hmratioJP0 = (TH1D *)hrecJP0->Clone();
    hmratioJP0->Divide(hrecJP0_All);
    TH1D *hmratioJP1 = (TH1D *)hrecJP1->Clone();
    hmratioJP1->Divide(hrecJP1_All);
    TH1D *hmratioJP2 = (TH1D *)hrecJP2->Clone();
    hmratioJP2->Divide(hrecJP2_All);

    // get data input and perform the background correction multiplying with the matching ratio bin-by-bin;
    TH1D *hMDatJP0 = (TH1D *)fdata->Get("hMDatJP0"); // Data input for unfolding
    hMDatJP0->Rebin(2.0);
    TH1D *hMDatJP0_Corr = (TH1D *)hMDatJP0->Clone(); // Data input for unfolding
    hMDatJP0_Corr->Multiply(hmratioJP0);

    TH1D *hMDatJP1 = (TH1D *)fdata->Get("hMDatJP1"); // Data input for unfolding
    hMDatJP1->Rebin(2.0);
    TH1D *hMDatJP1_Corr = (TH1D *)hMDatJP1->Clone(); // Data input for unfolding
    hMDatJP1_Corr->Multiply(hmratioJP1);
    TH1D *hMDatJP2 = (TH1D *)fdata->Get("hMDatJP2"); // Data input for unfolding
    hMDatJP2->Rebin(2.0);
    TH1D *hMDatJP2_Corr = (TH1D *)hMDatJP2->Clone(); // Data input for unfolding
    hMDatJP2_Corr->Multiply(hmratioJP2);

    TCanvas *cmratio = new TCanvas("mratio", "", 900, 900);
    cmratio->DivideSquare(9, 0.001, 0.001);
    cmratio->cd(1);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.06);
    gPad->SetLogy();
    hrecJP0_All->SetLineColor(2);
    hrecJP0_All->SetLineWidth(2);
    hrecJP0_All->SetLineStyle(1);
    hrecJP0_All->SetTitle("JP0: Signal & Background");
    hrecJP0_All->Draw("hist E");
    hrecJP0->SetLineColor(4);
    hrecJP0->SetLineWidth(2);
    hrecJP0->SetLineStyle(2);
    hrecJP0->Draw("hist E same");

    TLegend *leg11 = new TLegend(0.4, 0.75, 0.85, 0.9);
    leg11->AddEntry(hrecJP0_All, "Embed. Signal + Bkg", "l");
    leg11->AddEntry(hrecJP0, "Embed. Signal", "l");
    leg11->Draw();
    leg11->SetTextSize(0.04);
    gPad->Update();

    cmratio->cd(2);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.06);
    gPad->SetBottomMargin(0.06);

    setBinError(hmratioJP0, hrecJP0, hrecJP0_All);

    hmratioJP0->SetLineColor(4);
    hmratioJP0->SetLineWidth(2);
    hmratioJP0->SetTitle("Matching Fraction: Rec(Gen)/ Rec(Gen+Bkg)");
    hmratioJP0->SetTitle("JP0: Matching Fraction");
    hmratioJP0->Draw("hist E");

    TLegend *leg12 = new TLegend(0.1, 0.75, 0.3, 0.9);
    leg12->AddEntry("", "#frac{Signal}{Signal+Bkg}", "");
    leg12->SetTextSize(0.04);
    leg12->Draw();

    gPad->Update();

    cmratio->cd(3);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.06);
    gPad->SetBottomMargin(0.06);
    gPad->SetLogy();
    hMDatJP0->SetLineWidth(2);
    hMDatJP0->SetLineColor(2);
    hMDatJP0->SetTitle("JP0: Data Background Correction ");
    hMDatJP0->Draw("hist E");
    hMDatJP0_Corr->SetLineColor(4);
    hMDatJP0_Corr->SetLineWidth(2);
    hMDatJP0_Corr->SetLineStyle(2);
    hMDatJP0_Corr->Draw("hist E same");

    TLegend *leg13 = new TLegend(0.4, 0.75, 0.85, 0.9);
    leg13->AddEntry(hMDatJP0, "Data (Signal + Bkg)", "l");
    leg13->AddEntry(hMDatJP0_Corr, "Data Corrected", "l");
    leg13->SetTextSize(0.04);
    leg13->Draw();
    gPad->Update();
    // draw JP1
    cmratio->cd(4);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.06);
    gPad->SetLogy();
    hrecJP1_All->SetLineColor(2);
    hrecJP1_All->SetLineWidth(2);
    hrecJP1_All->SetLineStyle(1);
    hrecJP1_All->SetTitle("JP1: Signal & Background");
    hrecJP1_All->Draw("hist E");
    hrecJP1->SetLineColor(4);
    hrecJP1->SetLineWidth(2);
    hrecJP1->SetLineStyle(2);
    hrecJP1->Draw("hist E same");

    TLegend *leg21 = new TLegend(0.4, 0.75, 0.85, 0.9);
    leg21->AddEntry(hrecJP1_All, "Embed. Signal + Bkg", "l");
    leg21->AddEntry(hrecJP1, "Embed. Signal", "l");
    leg21->Draw();
    leg21->SetTextSize(0.04);
    gPad->Update();

    cmratio->cd(5);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.06);
    gPad->SetBottomMargin(0.06);

    setBinError(hmratioJP1, hrecJP1, hrecJP1_All);

    hmratioJP1->SetLineColor(4);
    hmratioJP1->SetLineWidth(2);
    hmratioJP1->SetTitle("Matching Fraction: Rec(Gen)/ Rec(Gen+Bkg)");
    hmratioJP1->SetTitle("JP1: Matching Fraction");
    hmratioJP1->Draw("hist E");

    TLegend *leg22 = new TLegend(0.1, 0.75, 0.3, 0.9);
    leg22->AddEntry("", "#frac{Signal}{Signal+Bkg}", "");
    leg22->SetTextSize(0.04);
    leg22->Draw();

    gPad->Update();

    cmratio->cd(6);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.06);
    gPad->SetBottomMargin(0.06);
    gPad->SetLogy();
    hMDatJP1->SetLineWidth(2);
    hMDatJP1->SetLineColor(2);
    hMDatJP1->SetTitle("JP1: Data Background Correction ");
    hMDatJP1->Draw("hist E");
    hMDatJP1_Corr->SetLineColor(4);
    hMDatJP1_Corr->SetLineWidth(2);
    hMDatJP1_Corr->SetLineStyle(2);
    hMDatJP1_Corr->Draw("hist E same");

    TLegend *leg23 = new TLegend(0.4, 0.75, 0.85, 0.9);
    leg23->AddEntry(hMDatJP1, "Data (Signal + Bkg)", "l");
    leg23->AddEntry(hMDatJP1_Corr, "Data Corrected", "l");
    leg23->SetTextSize(0.04);
    leg23->Draw();
    gPad->Update();

    // draw JP2
    cmratio->cd(7);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.06);
    gPad->SetLogy();
    hrecJP2_All->SetLineColor(2);
    hrecJP2_All->SetLineWidth(2);
    hrecJP2_All->SetLineStyle(1);
    hrecJP2_All->SetTitle("JP2: Signal & Background");
    hrecJP2_All->Draw("hist E");
    hrecJP2->SetLineColor(4);
    hrecJP2->SetLineWidth(2);
    hrecJP2->SetLineStyle(2);
    hrecJP2->Draw("hist E same");

    TLegend *leg31 = new TLegend(0.4, 0.75, 0.85, 0.9);
    leg31->AddEntry(hrecJP2_All, "Embed. Signal + Bkg", "l");
    leg31->AddEntry(hrecJP2, "Embed. Signal", "l");
    leg31->Draw();
    leg31->SetTextSize(0.04);
    gPad->Update();

    cmratio->cd(8);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.06);
    gPad->SetBottomMargin(0.06);

    setBinError(hmratioJP2, hrecJP2, hrecJP2_All);

    hmratioJP2->SetLineColor(4);
    hmratioJP2->SetLineWidth(2);
    hmratioJP2->SetTitle("Matching Fraction: Rec(Gen)/ Rec(Gen+Bkg)");
    hmratioJP2->SetTitle("JP2: Matching Fraction");
    hmratioJP2->Draw("hist E");

    TLegend *leg32 = new TLegend(0.1, 0.75, 0.3, 0.9);
    leg32->AddEntry("", "#frac{Signal}{Signal+Bkg}", "");
    leg32->SetTextSize(0.04);
    leg32->Draw();

    gPad->Update();

    cmratio->cd(9);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.06);
    gPad->SetBottomMargin(0.06);
    gPad->SetLogy();
    hMDatJP2->SetLineWidth(2);
    hMDatJP2->SetLineColor(2);
    hMDatJP2->SetTitle("JP2: Data Background Correction ");
    hMDatJP2->Draw("hist E");
    hMDatJP2_Corr->SetLineColor(4);
    hMDatJP2_Corr->SetLineWidth(2);
    hMDatJP2_Corr->SetLineStyle(2);
    hMDatJP2_Corr->Draw("hist E same");

    TLegend *leg33 = new TLegend(0.4, 0.75, 0.85, 0.9);
    leg33->AddEntry(hMDatJP2, "Data (Signal + Bkg)", "l");
    leg33->AddEntry(hMDatJP2_Corr, "Data Corrected", "l");
    leg33->SetTextSize(0.04);
    leg33->Draw();
    gPad->Update();

    cmratio->Update();
    cmratio->SaveAs("./UnfoldingResults/matchingFractionAngbkgCorrection.pdf");
}

double getBinError(double numBinCt, double numBinErr, double denBinCt, double denBinErr)
{
    return (numBinCt / denBinCt) * sqrt(pow(numBinErr / numBinCt, 2) + pow(denBinErr / denBinCt, 2));
}

void setBinError(TH1D *hratio, TH1D *hnum, TH1D *hden)
{
    for (int i = 1; i < hnum->GetNbinsX(); i++)
    {
        double binError = 0;
        if (hnum->GetBinContent(i) == 0 && hden->GetBinContent(i) == 0)
            continue;
        binError = getBinError(hnum->GetBinContent(i), hnum->GetBinError(i), hden->GetBinContent(i), hden->GetBinError(i));
        hratio->SetBinError(i, binError);
    }
}