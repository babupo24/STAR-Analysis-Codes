#include <iostream>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TFile.h"

using namespace std;

void unfold_subtract_bkg()
{

	TH1::SetDefaultSumw2();
	gStyle->SetPaintTextFormat("4.3f");
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	// number of bins in histograms
	// const int nDet = 40;
	// const int nGen = 16;

	// open file to write unfolding outputs
	TFile *fout = new TFile("UnfoldingResults/unfoldingOutput_bkgsubtracted.root", "recreate");

	Double_t presJP0 = 0.00707; //(1. / 141.35);
	Double_t presJP1 = 0.3978;	//(1. / 2.514);
	Double_t presJP2 = 1.0;

	const double low_x = 0.27;
	const double high_x = 4.0;

	const int xnBins = 16;
	Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.};
	const int xnBins2 = 18;
	Double_t xnBinsEdges2[xnBins2 + 1] = {0.0, 0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4., 4.2};
	TFile *fdata = new TFile("Histograms/histData4Unfolding.root", "R");
	// data control plots in x-section bins.
	TH1D *hxsecMJP0Ctr = (TH1D *)fdata->Get("hxsecMJP0Ctr");
	TH1D *hxsecMJP1Ctr = (TH1D *)fdata->Get("hxsecMJP1Ctr");
	TH1D *hxsecMJP2Ctr = (TH1D *)fdata->Get("hxsecMJP2Ctr");

	TFile *fembed = new TFile("Histograms/hist4ResAndUnfolding.root", "R");
	// generated control plots
	TH1D *hxsecMGenCtr = (TH1D *)fembed->Get("hxsecMinvGen");
	TH1D *hxsecMJP0GenCtr = (TH1D *)fembed->Get("hMGenJP0");
	TH1D *hxsecMJP1GenCtr = (TH1D *)fembed->Get("hMGenJP1");
	TH1D *hxsecMJP2GenCtr = (TH1D *)fembed->Get("hMGenJP2");

	// response matrix for  |eta|<0.5
	TH2D *hMGenVsRecJP0 = (TH2D *)fembed->Get("hMGenVsRecJP0"); // migration matrix
	// hMGenVsRecJP0->RebinX(3.125);
	TH2D *hMGenVsRecJP1 = (TH2D *)fembed->Get("hMGenVsRecJP1"); // migration matrix
	// hMGenVsRecJP1->RebinX(3.125);
	TH2D *hMGenVsRecJP2 = (TH2D *)fembed->Get("hMGenVsRecJP2"); // migration matrix
	// hMGenVsRecJP2->RebinX(3.125);

	// background
	TH1D *hMRecJP0bk1 = (TH1D *)fembed->Get("hMRecJP0bk1");
	TH1D *hMRecJP1bk1 = (TH1D *)fembed->Get("hMRecJP1bk1");
	TH1D *hMRecJP2bk1 = (TH1D *)fembed->Get("hMRecJP2bk1");

	TH1D *hMRecJP0bk2 = (TH1D *)fembed->Get("hMRecJP0bk2");
	TH1D *hMRecJP1bk2 = (TH1D *)fembed->Get("hMRecJP1bk2");
	TH1D *hMRecJP2bk2 = (TH1D *)fembed->Get("hMRecJP2bk2");

	int nDet = hMGenVsRecJP0->GetNbinsX();
	int nGen = hMGenVsRecJP0->GetNbinsY();
	// TH1D *hMDatJP0 = (TH1D *)fdata->Get("hxsecMinvJP0"); // Data input for unfolding
	TH1D *hMDatJP0 = (TH1D *)fdata->Get("hMDatJP0"); // Data input for unfolding
	// hMDatJP0->Rebin(3.125);
	//   TH1D *hMDatJP1 = (TH1D *)fdata->Get("hxsecMinvJP1"); // Data input for unfolding
	TH1D *hMDatJP1 = (TH1D *)fdata->Get("hMDatJP1"); // Data input for unfolding
	// hMDatJP1->Rebin(3.125);
	//   TH1D *hMDatJP2 = (TH1D *)fdata->Get("hxsecMinvJP2"); // Data input for unfolding
	TH1D *hMDatJP2 = (TH1D *)fdata->Get("hMDatJP2"); // Data input for unfolding
	// hMDatJP2->Rebin(3.125);
	//      regularize curvature
	TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
	// preserve the area
	TUnfoldDensity::EConstraint constraintMode = TUnfoldDensity::kEConstraintArea;
	// bin content is divided by the bin width
	TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
	// set up unfolding
	TUnfoldDensity unfoldJP0(hMGenVsRecJP0, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP1(hMGenVsRecJP1, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP2(hMGenVsRecJP2, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	// pass data hist as an input
	if (unfoldJP0.SetInput(hMDatJP0) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	unfoldJP1.SetInput(hMDatJP1); // pass data hist as an input
	if (unfoldJP1.SetInput(hMDatJP1) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP2.SetInput(hMDatJP2) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}

	// subtract background
	// unfoldJP0.SubtractBackground(hMRecJP0bk1, "backgroundJP0_1", 1.0 /*background normalization scale*/, 0.0 /*background notmalization error*/);
	// unfoldJP0.SubtractBackground(hMRecJP0bk2, "backgroundJP0_2", 1.0 /*background normalization scale*/, 0.0 /*background notmalization error*/);
	// unfoldJP1.SubtractBackground(hMRecJP1bk1, "backgroundJP1_1", 1.0 /*background normalization scale*/, 0.0 /*background notmalization error*/);
	// unfoldJP1.SubtractBackground(hMRecJP1bk2, "backgroundJP1_2", 1.0 /*background normalization scale*/, 0.0 /*background notmalization error*/);
	// unfoldJP2.SubtractBackground(hMRecJP2bk1, "backgroundJP2_1", 1.0 /*background normalization scale*/, 0.0 /*background notmalization error*/);
	// unfoldJP2.SubtractBackground(hMRecJP2bk2, "backgroundJP2_2", 1.0 /*background normalization scale*/, 0.0 /*background notmalization error*/);

	// do the unfolding here
	Double_t tauMin = 0.0;
	Double_t tauMax = 0.0;
	Int_t nScan = 30;
	Int_t iBestJP0, iBestJP1, iBestJP2;
	TSpline *logTauXJP0, *logTauYJP0, *logTauXJP1, *logTauYJP1, *logTauXJP2, *logTauYJP2;
	TGraph *lCurveJP0, *lCurveJP1, *lCurveJP2;
	// this method scans the parameter tau and finds the kink in the L curve
	// finally, the unfolding is done for the "best" choice of tau
	iBestJP0 = unfoldJP0.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP0, &logTauXJP0, &logTauYJP0);
	iBestJP1 = unfoldJP1.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP1, &logTauXJP1, &logTauYJP1);
	iBestJP2 = unfoldJP2.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP2, &logTauXJP2, &logTauYJP2);
	std::cout << " In tau=" << unfoldJP2.GetTau() << "\n";
	std::cout << "chi**2=" << unfoldJP2.GetChi2A() << "+" << unfoldJP2.GetChi2L()
			  << " / " << unfoldJP2.GetNdf() << "\n";

	// save point corresponding to the kink in the L curve as TGraph
	Double_t t_JP0[1], x_JP0[1], y_JP0[1], t_JP1[1], x_JP1[1], y_JP1[1], t_JP2[1], x_JP2[1], y_JP2[1];
	logTauXJP0->GetKnot(iBestJP0, t_JP0[0], x_JP0[0]);
	logTauYJP0->GetKnot(iBestJP1, t_JP0[0], y_JP0[0]);
	logTauXJP1->GetKnot(iBestJP1, t_JP1[0], x_JP1[0]);
	logTauYJP1->GetKnot(iBestJP1, t_JP1[0], y_JP1[0]);
	logTauXJP2->GetKnot(iBestJP2, t_JP2[0], x_JP2[0]);
	logTauYJP2->GetKnot(iBestJP2, t_JP2[0], y_JP2[0]);

	TGraph *bestLcurveJP0 = new TGraph(1, x_JP0, y_JP0);
	TGraph *bestLcurveJP1 = new TGraph(1, x_JP1, y_JP1);
	TGraph *bestLcurveJP2 = new TGraph(1, x_JP2, y_JP2);

	//============================================================
	// extract unfolding results into histograms

	// set up a bin map, excluding underflow and overflow bins
	// the binMap relates the the output of the unfolding to the final
	// histogram bins
	Int_t *binMap = new Int_t[nGen + 2]; //
	for (Int_t i = 1; i <= nGen; i++)
		binMap[i] = i;
	binMap[0] = 0;
	binMap[nGen + 1] = 0;
	// Bin mapping is not used in this case
	//-----------
	// get outputs
	TH1D *histMunfoldJP0 = new TH1D("UnfoldedJP0", " ", xnBins, xnBinsEdges);
	unfoldJP0.GetOutput(histMunfoldJP0);
	TH1D *histMunfoldJP1 = new TH1D("UnfoldedJP1", " ", xnBins, xnBinsEdges);
	unfoldJP1.GetOutput(histMunfoldJP1);
	TH1D *histMunfoldJP2 = new TH1D("UnfoldedJP2", " ", xnBins, xnBinsEdges);
	unfoldJP2.GetOutput(histMunfoldJP2);
	// get correlation matrix
	unfoldJP0.GetRhoIJtotal("hRhoIJJP0");
	unfoldJP1.GetRhoIJtotal("hRhoIJJP1");
	unfoldJP2.GetRhoIJtotal("hRhoIJJP2");
	// get probablity matrix
	TUnfold::EHistMap histmap = TUnfold::kHistMapOutputVert;
	TH2D *hProbMatrixJP0 = new TH2D("hProbMatrixJP0", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP0.GetProbabilityMatrix(hProbMatrixJP0, histmap);
	TH2D *hProbMatrixJP1 = new TH2D("hProbMatrixJP1", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP1.GetProbabilityMatrix(hProbMatrixJP1, histmap);
	TH2D *hProbMatrixJP2 = new TH2D("hProbMatrixJP2", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP2.GetProbabilityMatrix(hProbMatrixJP2, histmap);

	delete[] binMap;
	binMap = 0;
	// draw inputs
	TCanvas *cinput = new TCanvas("cinput", "", 900, 300);
	cinput->Divide(3, 1);
	cinput->cd(1);
	gPad->SetLogy();
	// back+signal
	TH1D *hMRecJP0All_c = (TH1D *)hMRecJP0All->Clone();
	hMRecJP0All_c->SetTitle("JP0, #color[2]{Signal}, #color[3]{Background}");
	hMRecJP0All_c->Rebin(2.5);
	hMRecJP0All_c->Scale(1. / (hMRecJP0All_c->Integral()));
	hMRecJP0All_c->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
	hMRecJP0All_c->Draw("hist E");
	// signal
	TH1D *hMRecJP0_c = (TH1D *)hMGenVsRecJP0->ProjectionX();
	hMRecJP0_c->Rebin(2.5);
	hMRecJP0_c->Scale(1. / (hMRecJP0_c->Integral()));
	hMRecJP0_c->SetLineColor(2);
	hMRecJP0_c->Draw("hist E SAME");

	// bkg
	TH1D *hMRecJP0bk1_c = (TH1D *)hMRecJP0bk1->Clone();
	hMRecJP0bk1_c->Add(hMRecJP0bk2);
	hMRecJP0bk1_c->Rebin(2.5);
	hMRecJP0bk1_c->Scale(1. / (hMRecJP0bk1_c->Integral()));
	hMRecJP0bk1_c->SetLineColor(3);
	hMRecJP0bk1_c->Draw("hist E SAME");

	cinput->cd(2);
	gPad->SetLogy();
	// back+signal
	TH1D *hMRecJP1All_c = (TH1D *)hMRecJP1All->Clone();
	hMRecJP1All_c->SetTitle("JP1, #color[2]{Signal}, #color[3]{Background}");
	hMRecJP1All_c->Rebin(2.5);
	hMRecJP1All_c->Scale(1. / (hMRecJP1All_c->Integral()));
	hMRecJP1All_c->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
	hMRecJP1All_c->Draw("hist E");
	// signal
	TH1D *hMRecJP1_c = (TH1D *)hMGenVsRecJP1->ProjectionX();
	hMRecJP1_c->Rebin(2.5);
	hMRecJP1_c->Scale(1. / (hMRecJP1_c->Integral()));
	hMRecJP1_c->SetLineColor(2);
	hMRecJP1_c->Draw("hist E SAME");

	// bkg
	TH1D *hMRecJP1bk1_c = (TH1D *)hMRecJP1bk1->Clone();
	hMRecJP1bk1_c->Add(hMRecJP1bk2);
	hMRecJP1bk1_c->Rebin(2.5);
	hMRecJP1bk1_c->Scale(1. / (hMRecJP1bk1_c->Integral()));
	hMRecJP1bk1_c->SetLineColor(3);
	hMRecJP1bk1_c->Draw("hist E SAME");

	cinput->cd(3);
	gPad->SetLogy();
	// back+signal
	TH1D *hMRecJP2All_c = (TH1D *)hMRecJP2All->Clone();
	hMRecJP2All_c->SetTitle("JP2, #color[2]{Signal}, #color[3]{Background}");
	hMRecJP2All_c->Rebin(2.5);
	hMRecJP2All_c->Scale(1. / (hMRecJP2All_c->Integral()));
	hMRecJP2All_c->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
	hMRecJP2All_c->Draw("hist E");
	// signal
	TH1D *hMRecJP2_c = (TH1D *)hMGenVsRecJP2->ProjectionX();
	hMRecJP2_c->Rebin(2.5);
	hMRecJP2_c->Scale(1. / (hMRecJP2_c->Integral()));
	hMRecJP2_c->SetLineColor(2);
	hMRecJP2_c->Draw("hist E SAME");

	// bkg
	TH1D *hMRecJP2bk1_c = (TH1D *)hMRecJP2bk1->Clone();
	hMRecJP2bk1_c->Add(hMRecJP2bk2);
	hMRecJP2bk1_c->Rebin(2.5);
	hMRecJP2bk1_c->Scale(1. / (hMRecJP2bk1_c->Integral()));
	hMRecJP2bk1_c->SetLineColor(3);
	hMRecJP2bk1_c->Draw("hist E SAME");

	cinput->Update();
	cinput->SaveAs("UnfoldingResults/UnfoldingSignalBkg.pdf");

	/*	// draw response matrix In
		TCanvas *cresIn = new TCanvas("cresIn", "", 900, 900);
		cresIn->Divide(2, 2);
		cresIn->cd(1);
		hMGenVsRecJP0->Draw("colz");
		hMGenVsRecJP0->SetTitle(" Migration Matrix JP0");
		hMGenVsRecJP0->GetXaxis()->SetTitle("M_{inv}(Recons.)");
		hMGenVsRecJP0->GetYaxis()->SetTitle("M_{inv}(Generated)");
		gPad->SetLogz();
		cresIn->cd(2);
		gPad->SetLogz();
		hMGenVsRecJP1->Draw("colz");
		hMGenVsRecJP1->SetTitle(" Migration Matrix JP1");
		hMGenVsRecJP1->GetXaxis()->SetTitle("M_{inv}(Recons.)");
		hMGenVsRecJP1->GetYaxis()->SetTitle("M_{inv}(Generated)");
		gPad->SetLogz();
		// hProbMatrixIn->Draw("colz ");
		cresIn->cd(3);
		hMGenVsRecJP2->Draw("colz");
		hMGenVsRecJP2->SetTitle(" Migration Matrix JP2");
		hMGenVsRecJP2->GetXaxis()->SetTitle("M_{inv}(Recons.)");
		hMGenVsRecJP2->GetYaxis()->SetTitle("M_{inv}(Generated)");
		gPad->SetLogz();

		cresIn->cd(4);

		hxsecMJP0Ctr->SetLineColor(2);
		hxsecMJP0Ctr->SetLineStyle(2);
		hxsecMJP0Ctr->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
		hxsecMJP0Ctr->SetTitle("Unfolded Distribution");
		hxsecMJP0Ctr->SetStats(0);
		hxsecMJP0Ctr->Draw("hist E");
		hxsecMJP0Ctr->SetMaximum(1e8);

		histMunfoldJP0->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
		histMunfoldJP0->SetTitle("Unfolded Distribution");
		histMunfoldJP0->SetStats(0);
		histMunfoldJP0->SetLineColor(2);
		histMunfoldJP0->SetMaximum(1e8);
		histMunfoldJP0->Draw("hist E");

		// histMunfoldJP0->Draw("hist same E");

		histMunfoldJP1->Draw("hist E same");
		histMunfoldJP1->SetStats(0);
		histMunfoldJP1->SetLineColor(3);

		// hxsecMJP1Ctr->SetLineColor(3);
		// hxsecMJP1Ctr->SetLineStyle(2);
		// hxsecMJP1Ctr->Draw("hist same E");

		histMunfoldJP2->Draw("hist E same");
		histMunfoldJP2->SetStats(0);
		histMunfoldJP2->SetLineColor(4);

		// hxsecMJP2Ctr->SetLineColor(4);
		// hxsecMJP2Ctr->SetLineStyle(2);
		// hxsecMJP2Ctr->Draw("hist same E");

		TLegend *legIn = new TLegend(0.7, 0.5, 0.85, 0.85);
		legIn->AddEntry(histMunfoldJP0, "JP0 Unfolded", "lp");
		// legIn->AddEntry(hxsecMJP0Ctr, "JP0 Input", "lp");
		legIn->AddEntry(histMunfoldJP1, "JP1 Unfolded", "lp");
		// legIn->AddEntry(hxsecMJP1Ctr, "JP1 Input", "lp");
		legIn->AddEntry(histMunfoldJP2, "JP2 Ufolded", "lp");
		// legIn->AddEntry(hxsecMJP2Ctr, "JP2 Input", "lp");
		legIn->Draw();
		gPad->SetLogy();

		cresIn->Update();
		// cresIn->SaveAs("UnfoldingResults/unfoldingMigrationMatrixnResults.pdf");

		// draw l-curve
		TCanvas *clcurve = new TCanvas("clcurve", "", 900, 900);
		clcurve->Divide(2, 2);
		clcurve->cd(1);
		lCurveJP0->Draw("AL");
		lCurveJP0->SetTitle(" L-Curve JP0");
		bestLcurveJP0->SetMarkerColor(kRed);
		bestLcurveJP0->Draw("*");
		clcurve->cd(2);
		lCurveJP1->Draw("AL");
		lCurveJP1->SetTitle(" L-Curve JP1");
		bestLcurveJP1->SetMarkerColor(kRed);
		bestLcurveJP1->Draw("*");
		clcurve->cd(3);
		lCurveJP2->Draw("AL");
		lCurveJP2->SetTitle(" L-Curve JP2");
		bestLcurveJP2->SetMarkerColor(kRed);
		bestLcurveJP2->Draw("*");
		clcurve->Update();
		clcurve->SaveAs("UnfoldingResults/lcurve.pdf");
		// draw correlation matrix
		TCanvas *crhoij = new TCanvas("rhoij", "rho", 900, 900);
		crhoij->Divide(2, 2);
		crhoij->cd(1);
		hRhoIJJP0->Draw("colz ");
		hRhoIJJP0->SetTitle("Correlation Matrix JP0 (#rho_{ij})");
		hRhoIJJP0->GetXaxis()->SetTitle("M-Reconstructed");
		hRhoIJJP0->GetYaxis()->SetTitle("M-Generated");
		crhoij->cd(2);
		hRhoIJJP1->Draw("colz ");
		hRhoIJJP1->SetTitle("Correlation Matrix JP1 (#rho_{ij})");
		hRhoIJJP1->GetXaxis()->SetTitle("M-Reconstructed");
		hRhoIJJP1->GetYaxis()->SetTitle("M-Generated");
		crhoij->cd(3);
		hRhoIJJP2->Draw("colz ");
		hRhoIJJP2->SetTitle("Correlation Matrix JP2 (#rho_{ij})");
		hRhoIJJP2->GetXaxis()->SetTitle("M-Reconstructed");
		hRhoIJJP2->GetYaxis()->SetTitle("M-Generated");
		crhoij->SaveAs("UnfoldingResults/correlationMatrix.pdf");

		// draw probablity matrix
		TCanvas *cprobmtx = new TCanvas("probmtx", "rho", 900, 900);
		cprobmtx->Divide(2, 2);
		cprobmtx->cd(1);
		gPad->SetLogz();
		hProbMatrixJP0->Draw("colz");
		// hProbMatrixJP0->GetZaxis()->SetRangeUser(1e-10, 1);
		hProbMatrixJP0->GetZaxis()->SetRangeUser(1e-10, 1);
		hProbMatrixJP0->SetTitle("Probablity Matrix JP0)");
		hProbMatrixJP0->GetXaxis()->SetTitle("M-Reconstructed");
		hProbMatrixJP0->GetYaxis()->SetTitle("M-Generated");
		cprobmtx->cd(2);
		gPad->SetLogz();
		hProbMatrixJP1->Draw("colz");
		hProbMatrixJP1->GetZaxis()->SetRangeUser(1e-10, 1);
		hProbMatrixJP1->SetTitle("Probablity Matrix JP1)");
		hProbMatrixJP1->GetXaxis()->SetTitle("M-Reconstructed");
		hProbMatrixJP1->GetYaxis()->SetTitle("M-Generated");
		cprobmtx->cd(3);
		gPad->SetLogz();
		hProbMatrixJP2->Draw("colz");
		hProbMatrixJP2->GetZaxis()->SetRangeUser(1e-10, 1);
		hProbMatrixJP2->SetTitle("Probablity Matrix JP2)");
		hProbMatrixJP2->GetXaxis()->SetTitle("M-Reconstructed");
		hProbMatrixJP2->GetYaxis()->SetTitle("M-Generated");
		cprobmtx->Update();
		cprobmtx->SaveAs("UnfoldingResults/probablityMatrix.pdf");
		// cout << "JP0  = " << histMunfoldJP0->Integral() << " " << hMDatJP0->Integral() << endl;
		// cout << "JP1  = " << histMunfoldJP1->Integral() << " " << hMDatJP1->Integral() << endl;
		// cout << "JP2  = " << histMunfoldJP2->Integral() << " " << hMDatJP2->Integral() << endl;

		// copy unfolded output with wider x-axis range (for cosmetic purpose)
		TH1D *dataJP0 = new TH1D("dataJP0 ", "", xnBins2, xnBinsEdges2);
		TH1D *dataJP1 = new TH1D("dataJP1 ", "", xnBins2, xnBinsEdges2);
		TH1D *dataJP2 = new TH1D("dataJP2 ", "", xnBins2, xnBinsEdges2);

		TH1D *genJP0 = new TH1D("genJP0 ", "", xnBins2, xnBinsEdges2);
		TH1D *genJP1 = new TH1D("genJP1 ", "", xnBins2, xnBinsEdges2);
		TH1D *genJP2 = new TH1D("genJP2 ", "", xnBins2, xnBinsEdges2);

		TH1D *unfoldedJP0 = new TH1D("unfoldedJP0 ", "", xnBins2, xnBinsEdges2);
		TH1D *unfoldedJP1 = new TH1D("unfoldedJP1 ", "", xnBins2, xnBinsEdges2);
		TH1D *unfoldedJP2 = new TH1D("unfoldedJP2 ", "", xnBins2, xnBinsEdges2);

		TH1D *hgen = new TH1D("generatedM", "", xnBins2, xnBinsEdges2);
		TH1D *hgenJP0 = new TH1D("hgenJP0", "", xnBins2, xnBinsEdges2);
		TH1D *hgenJP1 = new TH1D("hgenJP1", "", xnBins2, xnBinsEdges2);
		TH1D *hgenJP2 = new TH1D("hgenJP2", "", xnBins2, xnBinsEdges2);

		for (int i = 0; i < xnBins2; i++)
		{
			if (i == 0 || i == 17)
				continue;
			Double_t dct0 = hxsecMJP0Ctr->GetBinContent(i);
			Double_t dct_er0 = hxsecMJP0Ctr->GetBinError(i);
			if (dct0 != 0)
			{
				dataJP0->SetBinContent(i + 1, dct0);
				dataJP0->SetBinError(i + 1, dct_er0);
			}
			Double_t ct0 = histMunfoldJP0->GetBinContent(i);
			Double_t ct_er0 = histMunfoldJP0->GetBinError(i);
			if (ct0 != 0)
			{
				unfoldedJP0->SetBinContent(i + 1, ct0);
				unfoldedJP0->SetBinError(i + 1, ct_er0);
			}
		}
		for (int i = 0; i < xnBins2; i++)
		{
			if (i == 0 || i == 17)
				continue;
			Double_t dct1 = hxsecMJP1Ctr->GetBinContent(i);
			Double_t dct_er1 = hxsecMJP1Ctr->GetBinError(i);
			if (dct1 != 0)
			{
				dataJP1->SetBinContent(i + 1, dct1);
				dataJP1->SetBinError(i + 1, dct_er1);
			}
			Double_t ct1 = histMunfoldJP1->GetBinContent(i);
			Double_t ct_er1 = histMunfoldJP1->GetBinError(i);
			if (ct1 != 0)
			{
				unfoldedJP1->SetBinContent(i + 1, ct1);
				unfoldedJP1->SetBinError(i + 1, ct_er1);
			}
		}
		for (int i = 0; i < xnBins2; i++)
		{
			if (i == 0 || i == 17)
				continue;
			Double_t dct2 = hxsecMJP2Ctr->GetBinContent(i);
			Double_t dct_er2 = hxsecMJP2Ctr->GetBinError(i);
			if (dct2 != 0)
			{
				dataJP2->SetBinContent(i + 1, dct2);
				dataJP2->SetBinError(i + 1, dct_er2);
			}

			Double_t ct2 = histMunfoldJP2->GetBinContent(i);
			Double_t ct_er2 = histMunfoldJP2->GetBinError(i);
			if (ct2 != 0)
			{
				unfoldedJP2->SetBinContent(i + 1, ct2);
				unfoldedJP2->SetBinError(i + 1, ct_er2);
			}
		}
		for (int i = 0; i < xnBins2; i++)
		{
			if (i == 0 || i == 17)
				continue;
			Double_t ct = hxsecMGenCtr->GetBinContent(i);
			Double_t ct_er = hxsecMGenCtr->GetBinError(i);
			if (ct != 0)
			{
				hgen->SetBinContent(i + 1, ct);
				hgen->SetBinError(i + 1, ct_er);
			}
		}
		for (int i = 0; i < xnBins2; i++)
		{
			if (i == 0 || i == 17)
				continue;
			Double_t ct0 = hxsecMJP0GenCtr->GetBinContent(i);
			Double_t ct_er0 = hxsecMJP0GenCtr->GetBinError(i);
			Double_t ct1 = hxsecMJP1GenCtr->GetBinContent(i);
			Double_t ct_er1 = hxsecMJP1GenCtr->GetBinError(i);
			Double_t ct2 = hxsecMJP2GenCtr->GetBinContent(i);
			Double_t ct_er2 = hxsecMJP2GenCtr->GetBinError(i);
			if (ct0 != 0 && ct1 != 0 && ct1 != 0)
			{
				// cout << "i = " << i << " content = " << ct0 << endl;
				hgenJP0->SetBinContent(i + 1, ct0);
				hgenJP0->SetBinError(i + 1, ct_er0);
				hgenJP1->SetBinContent(i + 1, ct1);
				hgenJP1->SetBinError(i + 1, ct_er1);
				hgenJP2->SetBinContent(i + 1, ct2);
				hgenJP2->SetBinError(i + 1, ct_er2);
			}
		}

		TCanvas *can = new TCanvas("can", "", 600, 700);
		can->cd();
		TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
		pad1->SetBottomMargin(0);
		pad1->SetLeftMargin(0.15);
		pad1->Draw();
		pad1->cd();

		unfoldedJP0->SetMaximum(1e8);
		unfoldedJP0->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}} GeV/c^{2}");
		unfoldedJP0->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
		unfoldedJP0->GetYaxis()->SetTitleOffset(1.2);
		unfoldedJP0->SetTitle("");
		unfoldedJP0->SetLineColor(2);
		unfoldedJP0->SetLineWidth(2);
		unfoldedJP0->Draw("hist E");

		dataJP0->SetLineColor(2);
		dataJP0->SetMarkerStyle(20);
		dataJP0->SetMarkerColor(2);
		dataJP0->SetLineWidth(1.5);
		dataJP0->SetLineStyle(2);
		dataJP0->SetLineWidth(2);
		dataJP0->Draw("hist P E same");

		unfoldedJP1->SetLineColor(4);
		unfoldedJP1->SetLineWidth(2);
		unfoldedJP1->Draw("hist same E");
		dataJP1->SetLineColor(4);
		dataJP1->SetMarkerColor(4);
		dataJP1->SetMarkerStyle(20);
		dataJP1->SetLineWidth(1.5);
		dataJP1->SetLineStyle(2);
		dataJP1->SetLineWidth(2);
		dataJP1->Draw("hist same E");

		unfoldedJP2->SetLineColor(6);
		unfoldedJP2->SetLineWidth(2);
		unfoldedJP2->Draw("hist same E");
		dataJP2->SetLineColor(6);
		dataJP2->SetMarkerColor(6);
		dataJP2->SetMarkerStyle(20);
		dataJP2->SetLineWidth(1.5);
		dataJP2->SetLineStyle(2);
		dataJP2->SetLineWidth(2);
		dataJP2->Draw("hist same E");

		TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg->AddEntry(unfoldedJP0, "JP0 Unfolded", "lp");
		leg->AddEntry(dataJP0, "JP0 Input", "lp");
		leg->AddEntry(unfoldedJP1, "JP1 Unfolded", "lp");
		leg->AddEntry(dataJP1, "JP1 Input", "lp");
		leg->AddEntry(unfoldedJP2, "JP2 Ufolded", "lp");
		leg->AddEntry(dataJP2, "JP2 Input", "lp");
		leg->Draw();
		gPad->SetLogy();

		can->cd();
		TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
		pad2->SetTopMargin(0);
		pad2->SetLeftMargin(0.15);
		pad2->SetBottomMargin(0.3);
		pad2->Draw();
		pad2->cd();
		gPad->SetGrid(1, 1);
		TH1D *hr_JP0 = (TH1D *)unfoldedJP0->Clone();
		hr_JP0->Add(dataJP0, -1);
		hr_JP0->Divide(dataJP0);
		hr_JP0->GetXaxis()->SetTitleSize(0.12);
		hr_JP0->GetXaxis()->SetLabelSize(0.12);
		hr_JP0->GetYaxis()->SetTitle("(Un. - In.)/ In.");
		hr_JP0->GetYaxis()->SetTitleOffset(0.4);
		hr_JP0->GetYaxis()->SetTitleSize(0.12);
		hr_JP0->GetYaxis()->SetLabelSize(0.12);

		hr_JP0->Draw("hist E");
		TH1D *hr_JP1 = (TH1D *)unfoldedJP1->Clone();
		hr_JP1->Add(dataJP1, -1);
		hr_JP1->Divide(dataJP1);
		hr_JP1->Draw("hist E same");
		TH1D *hr_JP2 = (TH1D *)unfoldedJP2->Clone();
		hr_JP2->Add(dataJP2, -1);
		hr_JP2->Divide(dataJP2);
		hr_JP2->Draw("hist E same");

		can->Update();
		can->SaveAs("UnfoldingResults/unfoldingResModefied.pdf");

		// draw unfolded result with the control plot
		TCanvas *canc = new TCanvas("canc", "", 600, 500);
		canc->cd();
		gPad->SetBottomMargin(0.15);
		pad1->SetLeftMargin(0.15);
		unfoldedJP0->SetMaximum(1e8);
		unfoldedJP0->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}} GeV/c^{2}");
		unfoldedJP0->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
		unfoldedJP0->GetYaxis()->SetTitleOffset(1.2);
		unfoldedJP0->SetTitle("");
		unfoldedJP0->SetLineColor(2);
		unfoldedJP0->SetLineWidth(2);
		unfoldedJP0->DrawNormalized("hist E");

		hgenJP0->SetLineColor(2);
		hgenJP0->SetMarkerStyle(20);
		hgenJP0->SetMarkerColor(2);
		hgenJP0->SetLineWidth(1.5);
		hgenJP0->SetLineStyle(2);
		hgenJP0->SetLineWidth(2);
		hgenJP0->DrawNormalized("hist P E same");

		unfoldedJP1->SetLineColor(4);
		unfoldedJP1->SetLineWidth(2);
		unfoldedJP1->DrawNormalized("hist same E");
		hgenJP1->SetLineColor(4);
		hgenJP1->SetMarkerColor(4);
		hgenJP1->SetMarkerStyle(20);
		hgenJP1->SetLineWidth(1.5);
		hgenJP1->SetLineStyle(2);
		hgenJP1->SetLineWidth(2);
		hgenJP1->DrawNormalized("hist same E");

		unfoldedJP2->SetLineColor(6);
		unfoldedJP2->SetLineWidth(2);
		unfoldedJP2->DrawNormalized("hist same E");
		hgenJP2->SetLineColor(6);
		hgenJP2->SetMarkerColor(6);
		hgenJP2->SetMarkerStyle(20);
		hgenJP2->SetLineWidth(1.5);
		hgenJP2->SetLineStyle(2);
		hgenJP2->SetLineWidth(2);
		hgenJP2->DrawNormalized("hist same E");

		TLegend *leg2 = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg2->AddEntry(unfoldedJP0, "JP0 Unfolded", "lp");
		leg2->AddEntry(hgenJP0, "JP0 Gen", "lp");
		leg2->AddEntry(unfoldedJP1, "JP1 Unfolded", "lp");
		leg2->AddEntry(hgenJP1, "JP1 Gen", "lp");
		leg2->AddEntry(unfoldedJP2, "JP2 Ufolded", "lp");
		leg2->AddEntry(hgenJP2, "JP2 Gen", "lp");
		leg2->Draw();
		gPad->SetLogy();
		canc->Update();
		canc->SaveAs("./UnfoldingResults/Unfolding_Generated_comparison.pdf");
		//// with control plot
		// write outputs on file
		fout->cd();
		hMGenVsRecJP0->Write();
		hMGenVsRecJP1->Write();
		hMGenVsRecJP2->Write();

		hProbMatrixJP0->Write();
		hProbMatrixJP1->Write();
		hProbMatrixJP2->Write();

		bestLcurveJP0->Write();
		bestLcurveJP1->Write();
		bestLcurveJP2->Write();

		cresIn->Write();
		crhoij->Write();
		clcurve->Write();
		canc->Write();

		unfoldedJP0->Write();
		unfoldedJP1->Write();
		unfoldedJP2->Write();

		dataJP0->Write();
		dataJP1->Write();
		dataJP2->Write();

		hgen->Write();

		fout->Close();
	*/
}
