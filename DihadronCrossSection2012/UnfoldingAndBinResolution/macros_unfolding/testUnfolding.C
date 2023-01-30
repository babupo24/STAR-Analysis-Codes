#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

using namespace std;

void testUnfolding()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetLegendBorderSize(0);

    const double low_x = 0.27;
    const double high_x = 4.0;

    const int xnBins = 16;
    Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.};
    const int xnBins2 = 18;
    Double_t xnBinsEdges2[xnBins2 + 1] = {0.0, 0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4., 4.2};

    TFile *fembed = new TFile("hist4ResAndUnfolding.root", "R");
    if (!fembed)
        cout << "input file doesn't exist!" << endl;

    // histogram input
    TH1D *hinput = (TH1D *)fembed->Get("hMRecTestInput");
    hinput->Rebin(2.5);
    // migration matrix
    TH2D *hmig_mtx = (TH2D *)fembed->Get("hMGenVsRecTestM");
    hmig_mtx->RebinX(2.5);

    int nDet = hmig_mtx->GetNbinsX();
    int nGen = hmig_mtx->GetNbinsY();

    // rec control plot
    TH1D *hrec_ctr1 = (TH1D *)fembed->Get("hMRecTestControlPlot1st");
    TH1D *hrec_ctr2 = (TH1D *)fembed->Get("hMRecTestControlPlot2nd");
    // gen control plot
    TH1D *hgen_ctr1 = (TH1D *)fembed->Get("hMGenTestControlPlot1st");
    TH1D *hgen_ctr2 = (TH1D *)fembed->Get("hMGenTestControlPlot2nd");

    // set up unfolding
    //      regularize curvature
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    // preserve the area
    TUnfoldDensity::EConstraint constraintMode = TUnfoldDensity::kEConstraintArea;
    // bin content is divided by the bin width
    TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
    // set up unfolding
    TUnfoldDensity unfold(hmig_mtx, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);

    if (unfold.SetInput(hinput) >= 10000)
    {
        cout << "Unfolding result may be wrong....\n";
    }

    // do the unfolding here
    Double_t tauMin = 0.0;
    Double_t tauMax = 0.0;
    Int_t nScan = 50;
    Int_t iBest;
    TSpline *logTauX, *logTauY;
    TGraph *lCurve;
    // this method scans the parameter tau and finds the kink in the L curve
    // finally, the unfolding is done for the "best" choice of tau
    iBest = unfold.ScanLcurve(nScan, tauMin, tauMax, &lCurve, &logTauX, &logTauY);

    std::cout << " In tau=" << unfold.GetTau() << "\n";
    std::cout << "chi**2=" << unfold.GetChi2A() << "+" << unfold.GetChi2L()
              << " / " << unfold.GetNdf() << "\n";

    // save point corresponding to the kink in the L curve as TGraph
    Double_t t[1], x[1], y[1];
    logTauX->GetKnot(iBest, t[0], x[0]);
    logTauY->GetKnot(iBest, t[0], y[0]);

    TGraph *bestLcurve = new TGraph(1, x, y);

    // get outputs
    TH1D *histMunfold = new TH1D("Unfolded", "Unfolded ", xnBins, xnBinsEdges);
    unfold.GetOutput(histMunfold);
    // get correlation matrix
    unfold.GetRhoIJtotal("hRhoIJ");

    // get probablity matrix
    TUnfold::EHistMap histmap = TUnfold::kHistMapOutputVert;
    TH2D *hProbMatrix = new TH2D("hProbMatrix", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
    unfold.GetProbabilityMatrix(hProbMatrix, histmap);

    TCanvas *c1 = new TCanvas("c1", "", 600, 550);
    c1->cd();
    c1->SetGrid(0, 0);
    c1->SetLogy();

    // histMunfold->SetLineColor(2);
    // histMunfold->Draw("hist E");

    hrec_ctr1->SetName("Data Rec.");
    hrec_ctr1->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hrec_ctr1->SetLineColor(2);
    hrec_ctr1->SetMarkerStyle(8);
    hrec_ctr1->SetMarkerColor(2);
    hrec_ctr1->SetLineWidth(2);
    hrec_ctr1->SetLineStyle(1);
    hrec_ctr1->Draw("hist E");

    hrec_ctr2->SetName("MC Rec.");
    hrec_ctr2->SetLineColor(6);
    hrec_ctr2->SetLineWidth(2);
    hrec_ctr2->SetMarkerStyle(24);
    hrec_ctr2->SetMarkerColor(6);
    hrec_ctr2->SetLineStyle(3);
    hrec_ctr2->Draw("same hist E");

    hgen_ctr1->SetName("Data Gen.");
    hgen_ctr1->SetLineColor(4);
    hgen_ctr1->SetLineWidth(2);
    hgen_ctr1->SetLineStyle(4);
    hgen_ctr1->Draw("same hist E");

    hgen_ctr2->SetName("MC Gen.");
    hgen_ctr2->SetLineColor(7);
    hgen_ctr2->SetLineWidth(2);
    hgen_ctr2->SetLineStyle(6);
    hgen_ctr2->Draw("same hist E");

    c1->BuildLegend();
    c1->Update();
    c1->SaveAs("test_controlplots.pdf");

    // draw l-curve
    TCanvas *clcurve = new TCanvas("clcurve", "", 500, 500);
    clcurve->cd();
    lCurve->Draw("AL");
    lCurve->SetTitle(" L-Curve");
    bestLcurve->SetMarkerColor(kRed);
    bestLcurve->Draw("*");
    clcurve->SaveAs("test_lcurve.pdf");

    TCanvas *c2 = new TCanvas("c2", "", 500, 500);
    c2->cd();
    gPad->SetLogz();
    hmig_mtx->SetTitle("Migration Matrix");
    hmig_mtx->GetXaxis()->SetTitle("Reconstructed");
    hmig_mtx->GetYaxis()->SetTitle("Generated");
    hmig_mtx->Draw("colz");
    c2->Update();
    c2->SaveAs("test_migrationmtx.pdf");

    TCanvas *c3 = new TCanvas("c3", "", 600, 700);
    c3->cd();
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    pad1->SetRightMargin(0.12);
    pad1->SetBottomMargin(0.0);
    hgen_ctr1->SetLineStyle(1);
    hgen_ctr1->SetMarkerStyle(20);
    hgen_ctr1->SetLineColor(1);
    hgen_ctr1->SetLineWidth(1);
    hgen_ctr1->GetXaxis()->SetLabelSize(0);
    hgen_ctr1->SetName("Generated");
    hgen_ctr1->Draw("hist E");

    histMunfold->SetLineColor(2);
    histMunfold->SetLineWidth(2);
    histMunfold->SetLineStyle(2);
    histMunfold->SetName("Unfolded");
    histMunfold->Draw("hist E same");

    pad1->BuildLegend();
    // TLegend *leg1 = new TLegend(0.7, 0.85, 0.9, 0.85);
    // leg1->AddEntry(hgen_ctr1, "Generated", "l");
    // leg1->AddEntry(histMunfold, "Unfolded", "l");
    // leg1->Draw();

    c3->Update();
    c3->cd();
    TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 0.98, 0.3);
    pad2->Draw();
    pad2->cd();
    pad2->SetGrid(2, 2);
    pad2->SetBottomMargin(0.22);
    pad2->SetTopMargin(0.01);
    TH1D *r_rec = (TH1D *)histMunfold->Clone();
    r_rec->Divide(hgen_ctr1);
    r_rec->SetTitle("");
    r_rec->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    r_rec->GetXaxis()->SetTitleSize(0.08);
    r_rec->GetXaxis()->SetLabelSize(0.08);
    r_rec->GetYaxis()->SetLabelSize(0.08);
    r_rec->GetYaxis()->SetTitle("#frac{Gen}{Unf.}");
    r_rec->GetYaxis()->SetTitleSize(0.08);
    r_rec->GetYaxis()->SetTitleOffset(0.5);
    r_rec->GetYaxis()->CenterTitle();
    r_rec->GetYaxis()->SetRangeUser(0.0, 2.0);
    r_rec->SetMarkerStyle(20);
    r_rec->SetLineColor(2);
    r_rec->SetLineWidth(2);
    r_rec->Draw("hist E");

    pad2->Update();

    c3->Update();
    c3->SaveAs("test_unfoldresult.pdf");

    /*
    TCanvas *c5 = new TCanvas("c5", "", 600, 700);
                c5->cd();
                TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
                pad1->Draw();
                pad1->cd();
                pad1->SetLogy();
                pad1->SetRightMargin(0.12);
                pad1->SetBottomMargin(0.0);
                hgen_ctr1->SetLineStyle(1);
                hgen_ctr1->SetLineColor(1);
                hgen_ctr1->SetLineWidth(1);
                hgen_ctr1->GetXaxis()->SetLabelSize(0);
                hgen_ctr1->Draw("hist E");

                hgen_ctr2->SetLineColor(2);
                hgen_ctr2->SetLineWidth(2);
                hgen_ctr2->SetLineStyle(2);
                hgen_ctr2->Draw("hist E same");

                c5->Update();
                c5->cd();
                TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 0.98, 0.3);
                pad2->Draw();
                pad2->cd();
                pad2->SetGrid(2, 2);
                pad2->SetBottomMargin(0.2);
                pad2->SetTopMargin(0.01);
                TH1D *r_rec = (TH1D *)hgen_ctr1->Clone();
                r_rec->Divide(hgen_ctr2);
                r_rec->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
                r_rec->GetXaxis()->SetTitleSize(0.08);
                r_rec->GetXaxis()->SetLabelSize(0.08);
                r_rec->GetYaxis()->SetLabelSize(0.08);
                r_rec->GetYaxis()->SetRangeUser(0.0, 2.0);
                r_rec->SetLineColor(2);
                r_rec->SetLineWidth(2);
                r_rec->Draw("hist E");

                pad2->Update();

                c5->Update();
            */
}