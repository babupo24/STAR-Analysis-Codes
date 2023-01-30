/*
-> dihadron cross section unfolding approach
-> Same event and selection cuts are applied for the embedding and data
Response matrix:
	>> 2D histogram with the truth level events on y-axis and detector level events along x-axis. Number of y-bins (16 bins) is same as the x-section binning (unfolding bins) and x-axis binning is choosen uniform for best unfolding (should be greater than y-axis binning).
	>> ONLY detector level events that has corresponding generated pair goes in to the response matrix.
	>> If a pair has a tracks (single or both) that is not coming from the generated track, the pair is considered as a backgrond and goes for the matching ratio. The matching ratio is the fraction of reconstructed events associated to the generated events. Matching fraction is multiplied to the data input bin-by-bin, which is similar to the background correction.
Input:
	>> 1D histogram from the real data. Number of bins should be same as in the x-bins in the unfolding response matrix.

Unfolding:
	>> Unfolding is done using TUnfoldDensity with the L-Scan curve method.
*/
/***********UPDATES***********/

/***********UPDATES**********/
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

void tunfolddensity_trgCombined()
{

	TH1::SetDefaultSumw2();
	gStyle->SetPaintTextFormat("4.3f");
	gStyle->SetOptStat(0);
	// number of bins in histograms
	// const int nDet = 40;
	// const int nGen = 16;

	// open file to write unfolding outputs
	TFile *fout = new TFile("UnfoldingResults/unfoldingOutput.root", "recreate");

	Double_t presJP0 = 0.00707; //(1. / 141.35);
	Double_t presJP1 = 0.3978;	//(1. / 2.514);
	Double_t presJP2 = 1.0;

	const double low_x = 0.27;
	const double high_x = 4.0;

	const int xnBins = 16;
	Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};
	const int xnBins2 = 18;
	Double_t xnBinsEdges2[xnBins2 + 1] = {0.0, 0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.20, 4.0, 4.2};
	// input histogram files
	TFile *fembed = new TFile("Histograms/hist4ResAndUnfolding341.root", "R");
	// TFile *fembed = new TFile("Histograms/hist4ResAndUnfolding.root", "R");
	TFile *fdata = new TFile("Histograms/histData4unfolding.root", "R");

	//  get data input
	const TH1D *hMDatJP0 = (TH1D *)fdata->Get("hMDatJP0"); // Data input for unfolding
	const TH1D *hMDatJP1 = (TH1D *)fdata->Get("hMDatJP1"); // Data input for unfolding
	const TH1D *hMDatJP2 = (TH1D *)fdata->Get("hMDatJP2"); // Data input for unfolding

	// combine all trigger for the data input
	TH1D *hdata = (TH1D *)hMDatJP0->Clone();
	hdata->Add(hMDatJP1);
	hdata->Add(hMDatJP2);

	//  get embedding bkg+signal
	const TH1D *hMEmbJP0 = (TH1D *)fembed->Get("hMRecJP0All");
	const TH1D *hMEmbJP1 = (TH1D *)fembed->Get("hMRecJP1All");
	const TH1D *hMEmbJP2 = (TH1D *)fembed->Get("hMRecJP2All");
	// combine all trigger for the embedding
	TH1D *hembed = (TH1D *)hMEmbJP0->Clone();
	hembed->Scale(presJP0);
	hembed->Add(hMEmbJP1, presJP1);
	hembed->Add(hMEmbJP2, presJP2);

	// response matrix
	TH2D *hMGenVsRecJP0 = (TH2D *)fembed->Get("hMGenVsRecJP0NoUF"); // migration matrix
	hMGenVsRecJP0->RebinX(2.0);
	TH2D *hMGenVsRecJP1 = (TH2D *)fembed->Get("hMGenVsRecJP1NoUF"); // migration matrix
	hMGenVsRecJP1->RebinX(2.0);
	TH2D *hMGenVsRecJP2 = (TH2D *)fembed->Get("hMGenVsRecJP2NoUF"); // migration matrix
	hMGenVsRecJP2->RebinX(2.0);

	TH2D *hMGenVsRec = (TH2D *)hMGenVsRecJP0->Clone();
	hMGenVsRec->Scale(presJP0);
	hMGenVsRec->Add(hMGenVsRecJP1, presJP1);
	hMGenVsRec->Add(hMGenVsRecJP2, presJP2);

	int nDet = hMGenVsRecJP0->GetNbinsX();
	int nGen = hMGenVsRecJP0->GetNbinsY();

	// Unfolding with the matching ratio

	Double_t int_data = hdata->Integral(); // will be used to normalize the bkg from embedding

	Double_t int_embed = hembed->Integral(); // will be used to normalize the bkg from embedding

	// unfolding with the background subtraction--------
	// get backgrounds from embedding
	TH1D *hBkJP0 = (TH1D *)fembed->Get("hMRecJP0bk1");
	TH1D *hBk2JP0 = (TH1D *)fembed->Get("hMRecJP0bk2");
	hBkJP0->Add(hBk2JP0, 1);
	hBkJP0->Rebin(2.0);
	TH1D *hBkJP1 = (TH1D *)fembed->Get("hMRecJP1bk1");
	TH1D *hBk2JP1 = (TH1D *)fembed->Get("hMRecJP1bk2");
	hBkJP1->Add(hBk2JP1, 1);
	hBkJP1->Rebin(2.0);
	TH1D *hBkJP2 = (TH1D *)fembed->Get("hMRecJP2bk1");
	TH1D *hBk2JP2 = (TH1D *)fembed->Get("hMRecJP2bk2");
	hBkJP2->Add(hBk2JP2, 1);
	hBkJP2->Rebin(2.0);

	TH1D *hbkg = (TH1D *)hBkJP0->Clone();
	hbkg->Scale(presJP0);
	hbkg->Add(hBkJP1, presJP1);
	hbkg->Add(hBkJP2, presJP2);

	hbkg->Scale(1. / int_embed);
	hbkg->Scale(int_data); // data <=> embedding
	// data - background
	TH1D *hdata_bkg = (TH1D *)hdata->Clone();
	hdata_bkg->Add(hbkg, -1);

	//  regularize curvature
	TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
	// preserve the area
	TUnfoldDensity::EConstraint constraintMode = TUnfoldDensity::kEConstraintArea;
	// bin content is divided by the bin width
	TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
	// set up unfolding
	TUnfoldDensity unfoldJP(hMGenVsRec, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);

	// pass background subtracted data  as an input

	if (unfoldJP.SetInput(hdata_bkg) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}

	// do the unfolding here
	Double_t tauMin = 0.0;
	Double_t tauMax = 0.0;
	Int_t nScan = 50;
	Int_t iBestJP;
	TSpline *logTauXJP, *logTauYJP, ;
	TGraph *lCurveJP;
	// this method scans the parameter tau and finds the kink in the L curve
	// finally, the unfolding is done for the "best" choice of tau
	iBestJP = unfoldJP.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP, &logTauXJP, &logTauYJP);
	std::cout << "JP In tau=" << unfoldJP.GetTau() << "\n";
	std::cout << "chi**2=" << unfoldJP.GetChi2A() << "+" << unfoldJP.GetChi2L()
			  << " / " << unfoldJP.GetNdf() << "\n";

	// save point corresponding to the kink in the L curve as TGraph
	Double_t t_JP[1], x_JP[1], y_JP[1];
	logTauXJP->GetKnot(iBestJP, t_JP[0], x_JP[0]);
	logTauYJP->GetKnot(iBestJP, t_JP[0], y_JP[0]);

	TGraph *bestLcurveJP = new TGraph(1, x_JP, y_JP);

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
	TH1D *histMunfoldJP = new TH1D("UnfoldedJP", " ", xnBins, xnBinsEdges);
	unfoldJP.GetOutput(histMunfoldJP);
	// get correlation matrix
	unfoldJP.GetRhoIJtotal("hRhoIJJP");
	// get probablity matrix
	TUnfold::EHistMap histmap = TUnfold::kHistMapOutputVert;
	TH2D *hProbMatrixJP = new TH2D("hProbMatrixInJP", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP.GetProbabilityMatrix(hProbMatrixJP, histmap);

	TCanvas *can = new TCanvas("can", "", 900, 500);
	can->Divide(2, 1);
	can->cd(1);
	gPad->SetGrid(0, 0);
	hMGenVsRec->Draw("colz");
	gPad->SetLogz();
	can->cd(2);
	gPad->SetGrid(0, 0);
	histMunfoldJP->Draw("hist E");
	gPad->SetLogy();
	can->Update();
}
