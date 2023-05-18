// get the average invariant mass and pt of dipion in cross section bins

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void getAvgQt()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    TFile *fhist = new TFile("hist_Xsec_MinConeCut.root");
    // TFile *fhist = new TFile("test.root");
    if (!fhist->IsOpen())
        break;

    const int nbins = 13;
    Double_t nBinsEdges[nbins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};
    TH1D *hM[nbins] = {NULL};
    Double_t avgM[nbins] = {0};
    TH1D *hPt[nbins] = {NULL};
    Double_t avgPt[nbins] = {0};
    TH1D *hEta[nbins] = {NULL};
    Double_t avgEta[nbins] = {0};

    TH1D *hMc3[8] = {NULL};
    Double_t avgMc3[8] = {0};
    TH1D *hPtc3[8] = {NULL};
    Double_t avgPtc3[8] = {0};
    TH1D *hEtac3[8] = {NULL};
    Double_t avgEtac3[8] = {0};

    TH1D *hjp0 = (TH1D *)fhist->Get("hMDatJP0");
    TH1D *hjp1 = (TH1D *)fhist->Get("hMDatJP1");
    TH1D *hjp2 = (TH1D *)fhist->Get("hMDatJP2");
    TCanvas *can = new TCanvas("can", "", 500, 450);
    can->cd();
    can->SetGrid(0, 0);
    can->SetLeftMargin(0.1);
    can->SetBottomMargin(0.15);
    can->SetLogy();
    hjp0->SetLineColor(2);
    hjp0->SetLineWidth(2);
    hjp0->SetMaximum(5e7);
    hjp0->GetXaxis()->SetRangeUser(0, 4);
    hjp0->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}} GeV/c^{2}");
    hjp0->GetYaxis()->SetTitle("Raw #pi^{+}#pi^{-} Yields");
    hjp0->Draw("hist E");

    hjp1->SetLineColor(4);
    hjp1->SetLineWidth(2);
    hjp1->Draw("hist E same");
    hjp2->SetLineColor(6);
    hjp2->SetLineWidth(2);
    hjp2->Draw("hist E same");

    TLegend *leg = new TLegend(0.31, 0.7, 0.7, 0.88);
    leg->AddEntry(hjp0, Form("#color[2]{JP0, N^{#pi^{+}#pi^{-}}_{ev} = %3.2g, <M> = %3.2g GeV/c^{2}}", hjp0->GetEntries(), hjp0->GetMean()), "l");
    leg->AddEntry(hjp1, Form("#color[4]{JP1, N^{#pi^{+}#pi^{-}}_{ev} = %3.2g, <M> = %3.2g GeV/c^{2}}", hjp1->GetEntries(), hjp1->GetMean()), "l");
    leg->AddEntry(hjp2, Form("#color[6]{JP2, N^{#pi^{+}#pi^{-}}_{ev} = %3.2g, <M> = %3.2g GeV/c^{2}}", hjp2->GetEntries(), hjp2->GetMean()), "l");
    leg->SetTextSize(0.0325);
    leg->Draw();
    can->Update();
    can->SaveAs("ResultsMinConeCut/rawyields.pdf");

    ofstream ftab;
    ftab.open("ResultsMinConeCut/table_avgMpT.tex");
    ofstream ftabc3;
    ftabc3.open("ResultsMinConeCut/table_avgMpTc3.tex");
    ofstream fout;
    fout.open("ResultsMinConeCut/avgMpT.txt");
    fout << "#Average dipion invariant mass, p_T, and eta in cross section bins" << endl;
    fout << "----------------------------------------------------------------------------------" << endl;
    fout << "Bin Id \t M_{inv} Range \t <M_{inv}>(GeV/c^2) \t <p_{T}>(GeV/c) \t <Eta> " << endl;
    fout << "----------------------------------------------------------------------------------" << endl;

    ftab << "Bin Id & $M_{inv}$-Bin (GeV/$c^2$) & BW & $<M_{inv}>$(GeV/$c^2)$ & $<p_{T}>$(GeV/c) & $<\\eta> $\\\\ \n \\hline " << endl;
    ftabc3 << "Bin Id & M_{inv} Bin (GeV/c^2) & <M_{inv}>(GeV/c^2) & <p_{T}>(GeV/c) & <Eta> \\\\" << endl;
    for (int i = 0; i < nbins; i++)
    {
        hM[i] = (TH1D *)fhist->Get(Form("hMinv_%i", i));
        avgM[i] = hM[i]->GetMean();

        hPt[i] = (TH1D *)fhist->Get(Form("hPt_%i", i));
        avgPt[i] = hPt[i]->GetMean();

        hEta[i] = (TH1D *)fhist->Get(Form("hEta_%i", i));
        avgEta[i] = hEta[i]->GetMean();

        if (i < 8)
        {
            hMc3[i] = (TH1D *)fhist->Get(Form("hMinv_C3%i", i));
            avgMc3[i] = hMc3[i]->GetMean();

            hPtc3[i] = (TH1D *)fhist->Get(Form("hPt_C3%i", i));
            avgPtc3[i] = hPtc3[i]->GetMean();

            hEtac3[i] = (TH1D *)fhist->Get(Form("hEta_C3%i", i));
            avgEtac3[i] = hEtac3[i]->GetMean();

            ftabc3 << i + 1 << " & " << setprecision(3) << nBinsEdges[i] << " - " << nBinsEdges[i + 1] << " &   " << setprecision(3) << avgMc3[i] << " & " << setprecision(3) << avgPtc3[i] << "& " << setprecision(3) << avgEtac3[i] << "\\\\" << endl;
        }

        fout << i + 1 << " \t " << nBinsEdges[i] << " - " << nBinsEdges[i + 1] << " \t   " << avgM[i] << " \t \t   " << avgPt[i] << " \t      " << avgEta[i] << endl;
        ftab << i + 1 << " & " << setprecision(3) << nBinsEdges[i] << " - " << nBinsEdges[i + 1] << " &   " << setprecision(3) << fabs(nBinsEdges[i] - nBinsEdges[i + 1]) << " & " << setprecision(3) << avgM[i] << " & " << setprecision(3) << avgPt[i] << "& " << setprecision(3) << avgEta[i] << "\\\\" << endl;
    }
    ftab << "\\hline \\hline " << endl;
    fout << "----------------------------------------------------------------------------------" << endl;
    fout.close();
}
