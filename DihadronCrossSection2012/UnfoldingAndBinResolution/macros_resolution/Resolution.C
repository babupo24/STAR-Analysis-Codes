// synopsis
// read histogram files from the data and embedding. Get the invariant mass distribution for the pythia and geant JP0, which has the best match so far. Subtract them. The subtracted distribution should be a gaussian type, which gives the idea of reslution for binning.
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include <cmath>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph.h"

using namespace std;

void Resolution()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetLegendBorderSize(0);

    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.25);

    const int nBins = 12;
    Double_t nBinsEdges[nBins + 1] = {0.26, 0.56, 0.86, 1.16, 1.46, 1.76, 2.06, 2.36, 2.66, 2.96, 3.26, 3.56, 4.0};
    Double_t binCenter[12];

    TFile *fRes = new TFile("test_hist4ResAndUnfolding.root", "R");
    if (!fRes)
    {
        cout << "File not found" << endl;
    }
    Double_t norm = 1.0;
    Double_t scale[12];

    TH1D *hist[12];
    for (int i = 0; i < 12; i++)
    {

        binCenter[i] = (nBinsEdges[i] + nBinsEdges[i + 1]) / 2.0;
        hist[i] = (TH1D *)fRes->Get(Form("hMGenRec%i", i));
        hist[i]->SetName("M_{inv, Rec.}");
        hist[i]->GetXaxis()->SetLabelSize(0.08);
        hist[i]->GetXaxis()->SetTitleSize(0.08);
        hist[i]->GetYaxis()->SetLabelSize(0.08);
        hist[i]->GetYaxis()->SetTitleOffset(0.06);
        hist[i]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
        // hist[i]->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
        //  hist[i]->GetYaxis()->SetTitleSize(0.06);
    }
    TCanvas *can = new TCanvas("can", "", 1200, 900);
    can->Divide(3, 4, 0, 0);

    TF1 *fgaus[12];
    Double_t mean[12];
    Double_t width[12];
    Double_t xerr[12] = {0.0};
    TLegend *leg[12];
    for (int i = 0; i < 12; i++)
    {
        can->cd(i + 1);
        gPad->SetLogy();
        hist[i]->GetXaxis()->SetTitle("M_{inv, gen}^{#pi^{+}#pi^{-}} - M_{inv, rec}^{#pi^{+}#pi^{-}}");
        hist[i]->Draw();
        fgaus[i] = new TF1("fgaus", "gaus", -.5, 0.5);
        hist[i]->Fit(fgaus[i], "R");
        mean[i] = fgaus[i]->GetParameter(1);
        width[i] = fgaus[i]->GetParameter(2);

        leg[i] = new TLegend(0.2, 0.8, 0.4, 0.9);
        leg[i]->AddEntry("", Form("#color[2]{M: %3.2f - %3.2f}", nBinsEdges[i], nBinsEdges[i + 1]), "");
        leg[i]->SetTextSize(0.07);
        leg[i]->Draw();

        gPad->Update();
    }
    can->Update();
    can->SaveAs("Resolution_fit.pdf");

    TGraphErrors *gres = new TGraphErrors(12, binCenter, mean, xerr, width);
    TCanvas *cres = new TCanvas("cres", "", 500, 500);
    cres->cd();
    gres->SetTitle("");
    gres->SetMarkerStyle(8);
    // gres->SetFillStyle(1001);
    gres->SetMarkerColor(2);
    // gres->SetFillColor(3);
    gres->GetXaxis()->SetTitle("M_{inv}");
    gres->GetYaxis()->SetTitle("Resolution");
    gres->GetYaxis()->SetTitleOffset(1.2);
    gres->GetYaxis()->SetRangeUser(-0.15, 0.15);
    gres->Draw("ALP");
    gPad->Update();

    gPad->Update();
    TLine *lup = new TLine(0.06, fabs(mean[11]) + width[11], 4.1, fabs(mean[11]) + width[11]);
    lup->SetLineStyle(2);
    lup->SetLineWidth(2);
    lup->SetLineColor(2);
    lup->Draw();

    TLine *l0 = new TLine(0.06, 0, 4.1, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineWidth(2);
    l0->SetLineColor(1);
    l0->Draw();

    TLine *ldn = new TLine(0.06, (-1) * (fabs(mean[11]) + width[11]), 4.1, (-1) * (fabs(mean[11]) + width[11]));
    ldn->SetLineStyle(2);
    ldn->SetLineWidth(2);
    ldn->SetLineColor(2);
    ldn->Draw();

    TLatex tex;
    tex.SetTextSize(0.03);
    tex.DrawLatex(1.5, -0.09, Form("Max. Resolution: #pm %4.4f", fabs(mean[11]) + width[11]));
    tex.DrawLatex(1.5, 0.13, "#color[2]{#bullet}  Gaussian Mean");
    tex.DrawLatex(1.5, 0.11, "| Gaussian Width");
    tex.DrawLatex(1.5, 0.09, "#color[2]{---} Absolute resolution");

    cres->Update();
    cres->SaveAs("Resolution.pdf");
    cout << "*********************************" << endl;
    cout << "Double_t mean[12] = {";
    for (int i = 0; i < 12; i++)
    {
        if (i == 11)
        {
            cout << mean[i] << "};" << endl;
        }
        else
        {
            cout << mean[i] << ",";
        }
    }
    cout << "Double_t width[12] = {";
    for (int i = 0; i < 12; i++)
    {
        if (i == 11)
        {
            cout << width[i] << "};" << endl;
        }
        else
        {
            cout << width[i] << ",";
        }
    }

    cout << "*********************************" << endl;
}
