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
void drawHist2x2(TH2D *hmig, TGraph *lcurve, TGraph *knot, TH2D *hbinCorr, TH1D *hunfold, TH1D *hinput, TH1D *hdat_bkg);
void drawHist2x2E(TH2D *hmig, TGraph *lcurve, TGraph *knot, TH2D *hbinCorr, TH1D *hunfold);
void drawOutput3x1(TH1D *hinput0, TH1D *hbkg0, TH1D *hdat_bkg0, TH1D *hinput1, TH1D *hbkg1, TH1D *hdat_bkg1, TH1D *hinput2, TH1D *hbkg2, TH1D *hdat_bkg2);
void tunfolddensity()
{

	TH1::SetDefaultSumw2();
	gStyle->SetPaintTextFormat("4.3f");
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	// number of bins in histograms
	// const int nDet = 40;
	// const int nGen = 16;

	// input histogram files
	TFile *fembed = new TFile("Histograms/histEmbed4Unfolding.root", "R");
	// TFile *fembed = new TFile("Histograms/hist4ResAndUnfolding.root", "R");
	TFile *fdata = new TFile("Histograms/histData4Unfolding.root", "R");

	// open file to write unfolding outputs
	TFile *fout = new TFile("UnfoldingResults/unfoldingOutput.root", "recreate");

	Double_t presJP0 = 0.00707; //(1. / 141.35);
	Double_t presJP1 = 0.3978;	//(1. / 2.514);
	Double_t presJP2 = 1.0;

	const double low_x = 0.27;
	const double high_x = 4.0;

	// const int xnBins = 16;
	// Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};
	const int xnBins2 = 15;
	Double_t xnBinsEdges2[xnBins2 + 1] = {0.0, 0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0, 4.2};

	const int xnBins = 13;
	Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

	// data control plots in x-section bins.
	TH1D *hxsecMJP0Ctr = (TH1D *)fdata->Get("hxsecMJP0Ctr");
	TH1D *hxsecMJP1Ctr = (TH1D *)fdata->Get("hxsecMJP1Ctr");
	TH1D *hxsecMJP2Ctr = (TH1D *)fdata->Get("hxsecMJP2Ctr");

	//  get data input
	const TH1D *hMDatJP0 = (TH1D *)fdata->Get("hMDatJP0");
	const TH1D *hMDatJP1 = (TH1D *)fdata->Get("hMDatJP1");
	const TH1D *hMDatJP2 = (TH1D *)fdata->Get("hMDatJP2");

	const TH1D *hMDatJP0Gt = (TH1D *)fdata->Get("hMDatJP0Gt");
	const TH1D *hMDatJP1Gt = (TH1D *)fdata->Get("hMDatJP1Gt");
	const TH1D *hMDatJP2Gt = (TH1D *)fdata->Get("hMDatJP2Gt");

	const TH1D *hMDatJP0Lt = (TH1D *)fdata->Get("hMDatJP0Lt");
	const TH1D *hMDatJP1Lt = (TH1D *)fdata->Get("hMDatJP1Lt");
	const TH1D *hMDatJP2Lt = (TH1D *)fdata->Get("hMDatJP2Lt");

	//  get embedding bkg+signal
	const TH1D *hMEmbJP0 = (TH1D *)fembed->Get("hMRecJP0All");
	const TH1D *hMEmbJP1 = (TH1D *)fembed->Get("hMRecJP1All");
	const TH1D *hMEmbJP2 = (TH1D *)fembed->Get("hMRecJP2All");

	const TH1D *hMEmbJP0Lt = (TH1D *)fembed->Get("hMRecJP0AllLt");
	const TH1D *hMEmbJP1Lt = (TH1D *)fembed->Get("hMRecJP1AllLt");
	const TH1D *hMEmbJP2Lt = (TH1D *)fembed->Get("hMRecJP2AllLt");

	const TH1D *hMEmbJP0Gt = (TH1D *)fembed->Get("hMRecJP0AllGt");
	const TH1D *hMEmbJP1Gt = (TH1D *)fembed->Get("hMRecJP1AllGt");
	const TH1D *hMEmbJP2Gt = (TH1D *)fembed->Get("hMRecJP2AllGt");

	// generated control plots
	TH1D *hxsecMGenCtr = (TH1D *)fembed->Get("hxsecMinvGen");
	TH1D *hxsecMJP0GenCtr = (TH1D *)fembed->Get("hMGenJP0");
	TH1D *hxsecMJP1GenCtr = (TH1D *)fembed->Get("hMGenJP1");
	TH1D *hxsecMJP2GenCtr = (TH1D *)fembed->Get("hMGenJP2");

	// response matrix
	TH2D *hMGenVsRecJP0 = (TH2D *)fembed->Get("hMGenVsRecJP0NoUF"); // migration matrix
	TH2D *hMGenVsRecJP0Gt = (TH2D *)fembed->Get("hMGenVsRecJP0Gt"); // migration matrix
	TH2D *hMGenVsRecJP0Lt = (TH2D *)fembed->Get("hMGenVsRecJP0Lt"); // migration matrix
	TH2D *hMGenVsRecJP1 = (TH2D *)fembed->Get("hMGenVsRecJP1NoUF"); // migration matrix
	TH2D *hMGenVsRecJP1Gt = (TH2D *)fembed->Get("hMGenVsRecJP1Gt"); // migration matrix
	TH2D *hMGenVsRecJP1Lt = (TH2D *)fembed->Get("hMGenVsRecJP1Lt"); // migration matrix
	TH2D *hMGenVsRecJP2 = (TH2D *)fembed->Get("hMGenVsRecJP2NoUF"); // migration matrix
	TH2D *hMGenVsRecJP2Gt = (TH2D *)fembed->Get("hMGenVsRecJP2Gt"); // migration matrix
	TH2D *hMGenVsRecJP2Lt = (TH2D *)fembed->Get("hMGenVsRecJP2Lt"); // migration matrix

	// histograms for the matching fraction
	// denominator fine binning (gen&rec + rec)
	TH1D *hrecJP0_All = (TH1D *)fembed->Get("hMRecJP0All");
	TH1D *hrecJP1_All = (TH1D *)fembed->Get("hMRecJP1All");
	TH1D *hrecJP2_All = (TH1D *)fembed->Get("hMRecJP2All");
	// numerator fine binning (gen&rec)
	TH1D *hrecJP0 = (TH1D *)fembed->Get("hMRecJP0");
	TH1D *hrecJP1 = (TH1D *)fembed->Get("hMRecJP1");
	TH1D *hrecJP2 = (TH1D *)fembed->Get("hMRecJP2");

	//**********   for control plots  *********************
	// denominator x-sec binning (gen&rec + rec)
	// histograms for the matching fraction
	// denominator (gen&rec + rec)
	TH1D *hrecJP0_Allx = (TH1D *)fembed->Get("hMRecJP0Allx");
	TH1D *hrecJP1_Allx = (TH1D *)fembed->Get("hMRecJP1Allx");
	TH1D *hrecJP2_Allx = (TH1D *)fembed->Get("hMRecJP2Allx");
	// numerator (gen&rec)
	TH1D *hrecJP0x = (TH1D *)fembed->Get("hMRecJP0x");
	TH1D *hrecJP1x = (TH1D *)fembed->Get("hMRecJP1x");
	TH1D *hrecJP2x = (TH1D *)fembed->Get("hMRecJP2x");
	// backgrounds in cross-section bins
	TH1D *hbkJP0x = (TH1D *)fembed->Get("hMRecJP0bkx");
	TH1D *hbkJP1x = (TH1D *)fembed->Get("hMRecJP1bkx");
	TH1D *hbkJP2x = (TH1D *)fembed->Get("hMRecJP2bkx");

	Double_t int_dat_JP0 = hMDatJP0->Integral();	 // will be used to normalize the bkg from embedding
	Double_t int_dat_JP0Gt = hMDatJP0Gt->Integral(); // will be used to normalize the bkg from embedding
	Double_t int_dat_JP0Lt = hMDatJP0Lt->Integral(); // will be used to normalize the bkg from embedding
	Double_t int_dat_JP1 = hMDatJP1->Integral();
	Double_t int_dat_JP1Gt = hMDatJP1Gt->Integral();
	Double_t int_dat_JP1Lt = hMDatJP1Lt->Integral();
	Double_t int_dat_JP2 = hMDatJP2->Integral();
	Double_t int_dat_JP2Gt = hMDatJP2Gt->Integral();
	Double_t int_dat_JP2Lt = hMDatJP2Lt->Integral();

	Double_t int_emb_JP0 = hMEmbJP0->Integral();	 // will be used to normalize the bkg from embedding
	Double_t int_emb_JP0Gt = hMEmbJP0Gt->Integral(); // will be used to normalize the bkg from embedding
	Double_t int_emb_JP0Lt = hMEmbJP0Lt->Integral(); // will be used to normalize the bkg from embedding
	Double_t int_emb_JP1 = hMEmbJP1->Integral();
	Double_t int_emb_JP1Gt = hMEmbJP1Gt->Integral();
	Double_t int_emb_JP1Lt = hMEmbJP1Lt->Integral();
	Double_t int_emb_JP2 = hMEmbJP2->Integral();
	Double_t int_emb_JP2Gt = hMEmbJP2Gt->Integral();
	Double_t int_emb_JP2Lt = hMEmbJP2Lt->Integral();

	int nDet = hMGenVsRecJP0->GetNbinsX();
	int nGen = hMGenVsRecJP0->GetNbinsY();

	// bkg subtraction in x-section bins for control plots

	TH1D *hMDatJP0x = (TH1D *)fdata->Get("hxsecMJP0Ctr"); // Data input for unfolding
	TH1D *hMDatJP0_Corrx = (TH1D *)hMDatJP0x->Clone();	  // Data input for unfolding
	hMDatJP0_Corrx->Add(hbkJP0x, (-1) * (int_dat_JP0 / int_emb_JP0));

	TH1D *hMDatJP1x = (TH1D *)fdata->Get("hxsecMJP1Ctr"); // Data input for unfolding
	TH1D *hMDatJP1_Corrx = (TH1D *)hMDatJP1x->Clone();	  // Data input for unfolding
	hMDatJP1_Corrx->Add(hbkJP1x, (-1) * (int_dat_JP1 / int_emb_JP1));

	TH1D *hMDatJP2x = (TH1D *)fdata->Get("hxsecMJP2Ctr"); // Data input for unfolding
	TH1D *hMDatJP2_Corrx = (TH1D *)hMDatJP2x->Clone();	  // Data input for unfolding
	hMDatJP2_Corrx->Add(hbkJP2x, (-1) * (int_dat_JP2 / int_emb_JP2));

	//----------------------------------
	// unfolding with the background subtraction--------
	// get backgrounds from embedding
	TH1D *hBkJP0 = (TH1D *)fembed->Get("hMRecJP0bk1");
	TH1D *hBk2JP0 = (TH1D *)fembed->Get("hMRecJP0bk2");
	hBkJP0->Add(hBk2JP0);
	hBkJP0->Scale((Double_t)int_dat_JP0 / int_emb_JP0);

	TH1D *hBkJP1 = (TH1D *)fembed->Get("hMRecJP1bk1");
	TH1D *hBk2JP1 = (TH1D *)fembed->Get("hMRecJP1bk2");
	hBkJP1->Add(hBk2JP1);
	hBkJP1->Scale((Double_t)int_dat_JP1 / int_emb_JP1);

	TH1D *hBkJP2 = (TH1D *)fembed->Get("hMRecJP2bk1");
	TH1D *hBk2JP2 = (TH1D *)fembed->Get("hMRecJP2bk2");
	hBkJP2->Add(hBk2JP2);
	hBkJP2->Scale((Double_t)int_dat_JP2 / int_emb_JP2);

	TH1D *hBkJP0Gt = (TH1D *)fembed->Get("hMRecJP0bkGt");
	hBkJP0Gt->Scale((Double_t)int_dat_JP0Gt / int_emb_JP0Gt);
	TH1D *hBkJP0Lt = (TH1D *)fembed->Get("hMRecJP0bkLt");
	hBkJP0Lt->Scale((Double_t)int_dat_JP0Lt / int_emb_JP0Lt);

	TH1D *hBkJP1Gt = (TH1D *)fembed->Get("hMRecJP1bkGt");
	hBkJP1Gt->Scale((Double_t)int_dat_JP1Gt / int_emb_JP1Gt);
	TH1D *hBkJP1Lt = (TH1D *)fembed->Get("hMRecJP1bkLt");
	hBkJP1Lt->Scale((Double_t)int_dat_JP1Lt / int_emb_JP1Lt);

	TH1D *hBkJP2Gt = (TH1D *)fembed->Get("hMRecJP2bkGt");
	hBkJP2Gt->Scale((Double_t)int_dat_JP2Gt / int_emb_JP2Gt);
	TH1D *hBkJP2Lt = (TH1D *)fembed->Get("hMRecJP2bkLt");
	hBkJP2Lt->Scale((Double_t)int_dat_JP2Lt / int_emb_JP2Lt);

	// data - background
	TH1D *hdata_bkgJP0 = (TH1D *)hMDatJP0->Clone();
	hdata_bkgJP0->Rebin(1.0);
	hdata_bkgJP0->Add(hBkJP0, -1);
	TH1D *hdata_bkgJP1 = (TH1D *)hMDatJP1->Clone();
	hdata_bkgJP1->Rebin(1.0);
	hdata_bkgJP1->Add(hBkJP1, -1);
	TH1D *hdata_bkgJP2 = (TH1D *)hMDatJP2->Clone();
	hdata_bkgJP2->Rebin(1.0);
	hdata_bkgJP2->Add(hBkJP2, -1);

	// input histograms before bkg subtraction
	TH1D *hinputJP0 = (TH1D *)hMDatJP0->Clone();
	TH1D *hinputJP0Gt = (TH1D *)hMDatJP0Gt->Clone();
	TH1D *hinputJP0Lt = (TH1D *)hMDatJP0Lt->Clone();

	TH1D *hinputJP1 = (TH1D *)hMDatJP1->Clone();
	TH1D *hinputJP1Gt = (TH1D *)hMDatJP1Gt->Clone();
	TH1D *hinputJP1Lt = (TH1D *)hMDatJP1Lt->Clone();

	TH1D *hinputJP2 = (TH1D *)hMDatJP2->Clone();
	TH1D *hinputJP2Gt = (TH1D *)hMDatJP2Gt->Clone();
	TH1D *hinputJP2Lt = (TH1D *)hMDatJP2Lt->Clone();

	//  regularize curvature
	TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
	// preserve the area
	TUnfoldDensity::EConstraint constraintMode = TUnfoldDensity::kEConstraintArea;
	// bin content is divided by the bin width
	TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
	// set up unfolding
	TUnfoldDensity unfoldJP0(hMGenVsRecJP0, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP0Gt(hMGenVsRecJP0Gt, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP0Lt(hMGenVsRecJP0Lt, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);

	TUnfoldDensity unfoldJP1(hMGenVsRecJP1, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP1Gt(hMGenVsRecJP1Gt, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP1Lt(hMGenVsRecJP1Lt, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);

	TUnfoldDensity unfoldJP2(hMGenVsRecJP2, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP2Gt(hMGenVsRecJP2Gt, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldJP2Lt(hMGenVsRecJP2Lt, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);

	// pass background subtracted data  as an input
	/*
	if (unfoldJP0.SetInput(hdata_bkgJP0) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP1.SetInput(hdata_bkgJP1) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP2.SetInput(hdata_bkgJP2) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	*/
	if (unfoldJP0.SetInput(hinputJP0) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP0Lt.SetInput(hinputJP0Lt) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP0Gt.SetInput(hinputJP0Gt) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}

	if (unfoldJP1.SetInput(hinputJP1) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP1Gt.SetInput(hinputJP1Gt) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP1Lt.SetInput(hinputJP1Lt) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}

	if (unfoldJP2.SetInput(hinputJP2) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP2Gt.SetInput(hinputJP2Gt) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}
	if (unfoldJP2Lt.SetInput(hinputJP2Lt) >= 10000)
	{
		cout << "Unfolding result may be wrong....\n";
	}

	// subtract background
	// unfoldJP0.SubtractBackground(hBkJP0, "backgroundJP0", (Double_t)(int_dat_JP0 / int_emb_JP0) , 0.0 );
	// unfoldJP1.SubtractBackground(hBkJP1, "backgroundJP1", (Double_t)(int_dat_JP1 / int_emb_JP1) , 0.0 );
	// unfoldJP2.SubtractBackground(hBkJP2, "backgroundJP2", (Double_t)(int_dat_JP2 / int_emb_JP2) , 0.0 );

	unfoldJP0.SubtractBackground(hBkJP0, "backgroundJP0", 1.0, 0.0); // backgrounds are already normalized
	unfoldJP0Gt.SubtractBackground(hBkJP0Gt, "backgroundJP0Gt", 1.0, 0.0);
	unfoldJP0Lt.SubtractBackground(hBkJP0Lt, "backgroundJP0Lt", 1.0, 0.0);

	unfoldJP1.SubtractBackground(hBkJP1, "backgroundJP1", 1.0, 0.0);
	unfoldJP1Gt.SubtractBackground(hBkJP1Gt, "backgroundJP1Gt", 1.0, 0.0);
	unfoldJP1Lt.SubtractBackground(hBkJP1Lt, "backgroundJP1Lt", 1.0, 0.0);

	unfoldJP2.SubtractBackground(hBkJP2, "backgroundJP2", 1.0, 0.0);
	unfoldJP2Gt.SubtractBackground(hBkJP2Gt, "backgroundJP2Gt", 1.0, 0.0);
	unfoldJP2Lt.SubtractBackground(hBkJP2Lt, "backgroundJP2Lt", 1.0, 0.0);

	// do the unfolding here
	Double_t tauMin = 0.0;
	Double_t tauMax = 0.0;
	Int_t nScan = 30;
	Int_t iBestJP0, iBestJP1, iBestJP2;
	Int_t iBestJP0Gt, iBestJP1Gt, iBestJP2Gt;
	Int_t iBestJP0Lt, iBestJP1Lt, iBestJP2Lt;

	TSpline *logTauXJP0, *logTauYJP0, *logTauXJP1, *logTauYJP1, *logTauXJP2, *logTauYJP2;
	TSpline *logTauXJP0Gt, *logTauYJP0Gt, *logTauXJP1Gt, *logTauYJP1Gt, *logTauXJP2Gt, *logTauYJP2Gt;
	TSpline *logTauXJP0Lt, *logTauYJP0Lt, *logTauXJP1Lt, *logTauYJP1Lt, *logTauXJP2Lt, *logTauYJP2Lt;

	TGraph *lCurveJP0, *lCurveJP1, *lCurveJP2;
	TGraph *lCurveJP0Gt, *lCurveJP1Gt, *lCurveJP2Gt;
	TGraph *lCurveJP0Lt, *lCurveJP1Lt, *lCurveJP2Lt;

	// this method scans the parameter tau and finds the kink in the L curve
	// finally, the unfolding is done for the "best" choice of tau
	iBestJP0 = unfoldJP0.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP0, &logTauXJP0, &logTauYJP0);
	iBestJP0Gt = unfoldJP0Gt.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP0Gt, &logTauXJP0Gt, &logTauYJP0Gt);
	iBestJP0Lt = unfoldJP0Lt.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP0Lt, &logTauXJP0Lt, &logTauYJP0Lt);

	iBestJP1 = unfoldJP1.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP1, &logTauXJP1, &logTauYJP1);
	iBestJP1Gt = unfoldJP1Gt.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP1Gt, &logTauXJP1Gt, &logTauYJP1Gt);
	iBestJP1Lt = unfoldJP1Lt.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP1Lt, &logTauXJP1Lt, &logTauYJP1Lt);

	iBestJP2 = unfoldJP2.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP2, &logTauXJP2, &logTauYJP2);
	iBestJP2Gt = unfoldJP2Gt.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP2Gt, &logTauXJP2Gt, &logTauYJP2Gt);
	iBestJP2Lt = unfoldJP2Lt.ScanLcurve(nScan, tauMin, tauMax, &lCurveJP2Lt, &logTauXJP2Lt, &logTauYJP2Lt);

	std::cout << "JP0 In tau=" << unfoldJP0.GetTau() << "\n";
	std::cout << "chi**2=" << unfoldJP0.GetChi2A() << "+" << unfoldJP0.GetChi2L()
			  << " / " << unfoldJP0.GetNdf() << "\n";
	std::cout << "JP1 In tau=" << unfoldJP1.GetTau() << "\n";
	std::cout << "chi**2=" << unfoldJP1.GetChi2A() << "+" << unfoldJP1.GetChi2L()
			  << " / " << unfoldJP1.GetNdf() << "\n";
	std::cout << "JP2 In tau=" << unfoldJP2.GetTau() << "\n";
	std::cout << "chi**2=" << unfoldJP2.GetChi2A() << "+" << unfoldJP2.GetChi2L()
			  << " / " << unfoldJP2.GetNdf() << "\n";

	// save point corresponding to the kink in the L curve as TGraph
	Double_t t_JP0[1], x_JP0[1], y_JP0[1], t_JP1[1], x_JP1[1], y_JP1[1], t_JP2[1], x_JP2[1], y_JP2[1];
	Double_t t_JP0Gt[1], x_JP0Gt[1], y_JP0Gt[1], t_JP1Gt[1], x_JP1Gt[1], y_JP1Gt[1], t_JP2Gt[1], x_JP2Gt[1], y_JP2Gt[1];
	Double_t t_JP0Lt[1], x_JP0Lt[1], y_JP0Lt[1], t_JP1Lt[1], x_JP1Lt[1], y_JP1Lt[1], t_JP2Lt[1], x_JP2Lt[1], y_JP2Lt[1];

	logTauXJP0->GetKnot(iBestJP0, t_JP0[0], x_JP0[0]);
	logTauXJP0Gt->GetKnot(iBestJP0Gt, t_JP0Gt[0], x_JP0Gt[0]);
	logTauXJP0Lt->GetKnot(iBestJP0Lt, t_JP0Lt[0], x_JP0Lt[0]);

	logTauYJP0->GetKnot(iBestJP0, t_JP0[0], y_JP0[0]);
	logTauYJP0Gt->GetKnot(iBestJP0Gt, t_JP0Gt[0], y_JP0Gt[0]);
	logTauYJP0Lt->GetKnot(iBestJP0Lt, t_JP0Lt[0], y_JP0Lt[0]);

	logTauXJP1->GetKnot(iBestJP1, t_JP1[0], x_JP1[0]);
	logTauXJP1Gt->GetKnot(iBestJP1Gt, t_JP1Gt[0], x_JP1Gt[0]);
	logTauXJP1Lt->GetKnot(iBestJP1Lt, t_JP1Lt[0], x_JP1Lt[0]);

	logTauYJP1->GetKnot(iBestJP1, t_JP1[0], y_JP1[0]);
	logTauYJP1Gt->GetKnot(iBestJP1Gt, t_JP1Gt[0], y_JP1Gt[0]);
	logTauYJP1Lt->GetKnot(iBestJP1Lt, t_JP1Lt[0], y_JP1Lt[0]);

	logTauXJP2->GetKnot(iBestJP2, t_JP2[0], x_JP2[0]);
	logTauXJP2Gt->GetKnot(iBestJP2Gt, t_JP2Gt[0], x_JP2Gt[0]);
	logTauXJP2Lt->GetKnot(iBestJP2Lt, t_JP2Lt[0], x_JP2Lt[0]);

	logTauYJP2->GetKnot(iBestJP2, t_JP2[0], y_JP2[0]);
	logTauYJP2Gt->GetKnot(iBestJP2Gt, t_JP2Gt[0], y_JP2Gt[0]);
	logTauYJP2Lt->GetKnot(iBestJP2Lt, t_JP2Lt[0], y_JP2Lt[0]);

	TGraph *knotJP0 = new TGraph(1, x_JP0, y_JP0);
	TGraph *knotJP0Gt = new TGraph(1, x_JP0Gt, y_JP0Gt);
	TGraph *knotJP0Lt = new TGraph(1, x_JP0Lt, y_JP0Lt);

	TGraph *knotJP1 = new TGraph(1, x_JP1, y_JP1);
	TGraph *knotJP1Gt = new TGraph(1, x_JP1Gt, y_JP1Gt);
	TGraph *knotJP1Lt = new TGraph(1, x_JP1Lt, y_JP1Lt);

	TGraph *knotJP2 = new TGraph(1, x_JP2, y_JP2);
	TGraph *knotJP2Gt = new TGraph(1, x_JP2Gt, y_JP2Gt);
	TGraph *knotJP2Lt = new TGraph(1, x_JP2Lt, y_JP2Lt);

	//============================================================
	// extract unfolding results into histograms

	// set up a bin map, excluding underflow and overflow bins
	// the binMap relates the the output of the unfolding to the final
	// histogram bins
	Int_t *binMap = new Int_t[nGen + 2];
	for (Int_t i = 1; i <= nGen; i++)
		binMap[i] = i;
	binMap[0] = 0;
	binMap[nGen + 1] = 0;
	// Bin mapping is not used in this case
	//-----------
	// get outputs for eta integrated
	// get result for different number of iterations

	TH1D *histMunfoldJP0 = new TH1D("UnfoldedJP0", " ", xnBins, xnBinsEdges);
	unfoldJP0.GetOutput(histMunfoldJP0);
	TH1D *histMunfoldJP1 = new TH1D("UnfoldedJP1", " ", xnBins, xnBinsEdges);
	unfoldJP1.GetOutput(histMunfoldJP1);
	TH1D *histMunfoldJP2 = new TH1D("UnfoldedJP2", " ", xnBins, xnBinsEdges);
	unfoldJP2.GetOutput(histMunfoldJP2);
	// get outputs for eta > 0
	TH1D *histMunfoldJP0Gt = new TH1D("UnfoldedJP0Gt", " ", xnBins, xnBinsEdges);
	unfoldJP0Gt.GetOutput(histMunfoldJP0Gt);
	TH1D *histMunfoldJP1Gt = new TH1D("UnfoldedJP1Gt", " ", xnBins, xnBinsEdges);
	unfoldJP1Gt.GetOutput(histMunfoldJP1Gt);
	TH1D *histMunfoldJP2Gt = new TH1D("UnfoldedJP2Gt", " ", xnBins, xnBinsEdges);
	unfoldJP2Gt.GetOutput(histMunfoldJP2Gt);
	// get outputs for eta < 0
	TH1D *histMunfoldJP0Lt = new TH1D("UnfoldedJP0Lt", " ", xnBins, xnBinsEdges);
	unfoldJP0Lt.GetOutput(histMunfoldJP0Lt);
	TH1D *histMunfoldJP1Lt = new TH1D("UnfoldedJP1Lt", " ", xnBins, xnBinsEdges);
	unfoldJP1Lt.GetOutput(histMunfoldJP1Lt);
	TH1D *histMunfoldJP2Lt = new TH1D("UnfoldedJP2Lt", " ", xnBins, xnBinsEdges);
	unfoldJP2Lt.GetOutput(histMunfoldJP2Lt);

	// get correlation matrix
	// eta integrated
	unfoldJP0.GetRhoIJtotal("hRhoIJJP0");
	unfoldJP1.GetRhoIJtotal("hRhoIJJP1");
	unfoldJP2.GetRhoIJtotal("hRhoIJJP2");
	// eta > 0
	unfoldJP0Gt.GetRhoIJtotal("hRhoIJJP0Gt");
	unfoldJP1Gt.GetRhoIJtotal("hRhoIJJP1Gt");
	unfoldJP2Gt.GetRhoIJtotal("hRhoIJJP2Gt");
	// eta < 0
	unfoldJP0Lt.GetRhoIJtotal("hRhoIJJP0Lt");
	unfoldJP1Lt.GetRhoIJtotal("hRhoIJJP1Lt");
	unfoldJP2Lt.GetRhoIJtotal("hRhoIJJP2Lt");

	// get probablity matrix
	TUnfold::EHistMap histmap = TUnfold::kHistMapOutputVert;
	// eta integrated
	TH2D *hProbMatrixJP0 = new TH2D("hProbMatrixInJP0", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP0.GetProbabilityMatrix(hProbMatrixJP0, histmap);
	TH2D *hProbMatrixJP1 = new TH2D("hProbMatrixInJP1", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP1.GetProbabilityMatrix(hProbMatrixJP1, histmap);
	TH2D *hProbMatrixJP2 = new TH2D("hProbMatrixInJP2", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP2.GetProbabilityMatrix(hProbMatrixJP2, histmap);
	// eta>0
	TH2D *hProbMatrixJP0Gt = new TH2D("hProbMatrixInJP0Gt", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP0Gt.GetProbabilityMatrix(hProbMatrixJP0Gt, histmap);
	TH2D *hProbMatrixJP1Gt = new TH2D("hProbMatrixInJP1Gt", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP1Gt.GetProbabilityMatrix(hProbMatrixJP1Gt, histmap);
	TH2D *hProbMatrixJP2Gt = new TH2D("hProbMatrixInJP2Gt", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP2Gt.GetProbabilityMatrix(hProbMatrixJP2Gt, histmap);
	// eta<0
	TH2D *hProbMatrixJP0Lt = new TH2D("hProbMatrixInJP0Lt", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP0Lt.GetProbabilityMatrix(hProbMatrixJP0Lt, histmap);
	TH2D *hProbMatrixJP1Lt = new TH2D("hProbMatrixInJP1Lt", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP1Lt.GetProbabilityMatrix(hProbMatrixJP1Lt, histmap);
	TH2D *hProbMatrixJP2Lt = new TH2D("hProbMatrixInJP2Lt", " Probability Matrix", nDet, low_x, high_x, nGen, xnBinsEdges);
	unfoldJP2Lt.GetProbabilityMatrix(hProbMatrixJP2Lt, histmap);

	delete[] binMap;
	binMap = 0;

	cout << "JP2 Yields data = " << hdata_bkgJP2->Integral() << " unfolded = " << histMunfoldJP2->Integral() << endl;
	cout << "JP1 Yields data = " << hdata_bkgJP1->Integral() << " unfolded = " << histMunfoldJP1->Integral() << endl;
	cout << "JP0 Yields data = " << hdata_bkgJP0->Integral() << " unfolded = " << histMunfoldJP0->Integral() << endl;

	// draw input and background
	drawOutput3x1(hMDatJP0, hBkJP0, hdata_bkgJP0, hMDatJP1, hBkJP1, hdata_bkgJP1, hMDatJP2, hBkJP2, hdata_bkgJP2);
	can3x1->SaveAs("UnfoldingResults/input_bkg_3x1.pdf");
	// draw 2x2
	// eta integrated
	drawHist2x2(hMGenVsRecJP0, lCurveJP0, knotJP0, hRhoIJJP0, histMunfoldJP0, hMDatJP0x, hMDatJP0_Corrx);
	can->SaveAs("UnfoldingResults/c2x2_JP0.pdf");
	drawHist2x2(hMGenVsRecJP1, lCurveJP1, knotJP1, hRhoIJJP1, histMunfoldJP1, hMDatJP1x, hMDatJP1_Corrx);
	can->SaveAs("UnfoldingResults/c2x2_JP1.pdf");
	drawHist2x2(hMGenVsRecJP2, lCurveJP2, knotJP2, hRhoIJJP2, histMunfoldJP2, hMDatJP2x, hMDatJP2_Corrx);
	can->SaveAs("UnfoldingResults/c2x2_JP2.pdf");
	// eta > 0
	drawHist2x2E(hMGenVsRecJP0Gt, lCurveJP0Gt, knotJP0Gt, hRhoIJJP0Gt, histMunfoldJP0Gt);
	can->SaveAs("UnfoldingResults/c2x2_JP0Gt.pdf");
	drawHist2x2E(hMGenVsRecJP1Gt, lCurveJP1Gt, knotJP1Gt, hRhoIJJP1Gt, histMunfoldJP1Gt);
	can->SaveAs("UnfoldingResults/c2x2_JP1Gt.pdf");
	drawHist2x2E(hMGenVsRecJP2Gt, lCurveJP2Gt, knotJP2Gt, hRhoIJJP2Gt, histMunfoldJP2Gt);
	can->SaveAs("UnfoldingResults/c2x2_JP2Gt.pdf");
	// eta > 0
	drawHist2x2E(hMGenVsRecJP0Lt, lCurveJP0Lt, knotJP0Lt, hRhoIJJP0Lt, histMunfoldJP0Lt);
	can->SaveAs("UnfoldingResults/c2x2_JP0Lt.pdf");
	drawHist2x2E(hMGenVsRecJP1Lt, lCurveJP1Lt, knotJP1Lt, hRhoIJJP1Lt, histMunfoldJP1Lt);
	can->SaveAs("UnfoldingResults/c2x2_JP1Lt.pdf");
	drawHist2x2E(hMGenVsRecJP2Lt, lCurveJP2Lt, knotJP2Lt, hRhoIJJP2Lt, histMunfoldJP2Lt);
	can->SaveAs("UnfoldingResults/c2x2_JP2Lt.pdf");

	fout->cd();
	// write input histograms
	hinputJP0->Write();
	hBkJP0->Write(); // from embedding already scaled to the data
	hinputJP1->Write();
	hBkJP1->Write();
	hinputJP2->Write();
	hBkJP2->Write();
	hinputJP0Gt->Write();
	hBkJP0Gt->Write();
	hinputJP1Gt->Write();
	hBkJP1Gt->Write();
	hinputJP2Gt->Write();
	hBkJP2Gt->Write();
	hinputJP0Lt->Write();
	hBkJP0Lt->Write();
	hinputJP1Lt->Write();
	hBkJP1Lt->Write();
	hinputJP2Lt->Write();
	hBkJP2Lt->Write();
	// write unfolding output
	// unfolding distribution
	histMunfoldJP0->Write();
	histMunfoldJP0Gt->Write();
	histMunfoldJP0Lt->Write();
	histMunfoldJP1->Write();
	histMunfoldJP1Gt->Write();
	histMunfoldJP1Lt->Write();
	histMunfoldJP2->Write();
	histMunfoldJP2Gt->Write();
	histMunfoldJP2Lt->Write();
	// migration matrix
	hMGenVsRecJP0->Write();
	hMGenVsRecJP0Gt->Write();
	hMGenVsRecJP0Lt->Write();
	hMGenVsRecJP1->Write();
	hMGenVsRecJP1Gt->Write();
	hMGenVsRecJP1Lt->Write();
	hMGenVsRecJP2->Write();
	hMGenVsRecJP2Gt->Write();
	hMGenVsRecJP2Lt->Write();

	// bin correlation
	hRhoIJJP0->Write();
	hRhoIJJP0Gt->Write();
	hRhoIJJP0Lt->Write();
	hRhoIJJP1->Write();
	hRhoIJJP1Gt->Write();
	hRhoIJJP1Lt->Write();
	hRhoIJJP2->Write();
	hRhoIJJP2Gt->Write();
	hRhoIJJP2Lt->Write();
	fout->Close();
}

void drawHist2x2(TH2D *hmig, TGraph *lcurve, TGraph *knot, TH2D *hbinCorr, TH1D *hunfold, TH1D *hinput, TH1D *hdat_bkg)
{
	can = new TCanvas("can", "c2x2", 900, 900);
	can->Divide(2, 2);
	can->cd(1);
	// draw migration matrix
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	hmig->SetTitle("Migration Matrix");
	hmig->GetXaxis()->SetTitle("M_{rec} GeV/c^{2}");
	hmig->GetYaxis()->SetTitle("M_{gen} GeV/c^{2}");
	hmig->GetYaxis()->SetTitleOffset(1.4);
	hmig->Draw("colz");
	gPad->SetLogz();
	can->cd(2);
	// draw lcurve
	gPad->SetLeftMargin(0.15);
	lcurve->SetTitle("L-Curve");
	lcurve->GetXaxis()->SetTitle("log(#it{L_{1}})");
	lcurve->GetYaxis()->SetTitle("log(#it{L_{2}/#tau^{2}})");
	lcurve->GetYaxis()->SetTitleOffset(1.4);
	lcurve->GetYaxis()->SetNdivisions(510);
	lcurve->Draw();
	knot->Draw("* same");
	knot->SetMarkerColor(2);
	can->cd(3);
	// draw bin correlation
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	hbinCorr->SetTitle("Bin Correlation Matrix");
	hbinCorr->GetXaxis()->SetTitle("M_{rec} GeV/c^{2}");
	hbinCorr->GetYaxis()->SetTitle("M_{gen} GeV/c^{2}");
	hbinCorr->GetYaxis()->SetTitleOffset(1.4);
	hbinCorr->Draw("colz");
	can->cd(4);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	// draw unfolding result
	hunfold->SetTitle("Unfolding Result");
	hunfold->SetLineColor(2);
	hunfold->GetXaxis()->SetTitle("M_{gen} GeV/c^{2}");
	hunfold->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
	hunfold->GetYaxis()->SetTitleOffset(1.4);
	hunfold->Draw("hist E");

	hinput->SetLineColor(3);
	hinput->Draw("same hist E");

	hdat_bkg->SetLineColor(4);
	hdat_bkg->SetLineStyle(2);
	hdat_bkg->Draw("same hist E");

	TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.85);
	leg->AddEntry(hinput, "Raw Data", "l");
	leg->AddEntry(hdat_bkg, "Raw Data - Bkg", "l");
	leg->AddEntry(hunfold, "Unfolded", "l");
	leg->Draw();

	can->Update();
}
void drawHist2x2E(TH2D *hmig, TGraph *lcurve, TGraph *knot, TH2D *hbinCorr, TH1D *hunfold)
{
	can = new TCanvas("can", "c2x2", 900, 900);
	can->Divide(2, 2);
	can->cd(1);
	// draw migration matrix
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	hmig->SetTitle("Migration Matrix");
	hmig->GetXaxis()->SetTitle("M_{rec} GeV/c^{2}");
	hmig->GetYaxis()->SetTitle("M_{gen} GeV/c^{2}");
	hmig->GetYaxis()->SetTitleOffset(1.4);
	hmig->Draw("colz");
	gPad->SetLogz();
	can->cd(2);
	// draw lcurve
	gPad->SetLeftMargin(0.15);
	lcurve->SetTitle("L-Curve");
	lcurve->GetXaxis()->SetTitle("log(#it{L_{1}})");
	lcurve->GetYaxis()->SetTitle("log(#it{L_{2}/#tau^{2}})");
	lcurve->GetYaxis()->SetTitleOffset(1.4);
	lcurve->GetYaxis()->SetNdivisions(510);
	lcurve->Draw();
	knot->Draw("* same");
	knot->SetMarkerColor(2);
	can->cd(3);
	// draw bin correlation
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	hbinCorr->SetTitle("Bin Correlation Matrix");
	hbinCorr->GetXaxis()->SetTitle("M_{rec} GeV/c^{2}");
	hbinCorr->GetYaxis()->SetTitle("M_{gen} GeV/c^{2}");
	hbinCorr->GetYaxis()->SetTitleOffset(1.4);
	hbinCorr->Draw("colz");
	can->cd(4);
	gPad->SetLeftMargin(0.15);
	// draw unfolding result
	hunfold->SetTitle("Unfolding Result");
	hunfold->SetLineColor(2);
	hunfold->GetXaxis()->SetTitle("M_{gen} GeV/c^{2}");
	hunfold->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
	hunfold->GetYaxis()->SetTitleOffset(1.4);
	hunfold->Draw("hist E");

	gPad->SetLogy();
	can->Update();
}

void drawOutput3x1(TH1D *hinput0, TH1D *hbkg0, TH1D *hdat_bkg0, TH1D *hinput1, TH1D *hbkg1, TH1D *hdat_bkg1, TH1D *hinput2, TH1D *hbkg2, TH1D *hdat_bkg2)
{
	TCanvas *c3x1 = new TCanvas("can3x1", "", 1200, 400);
	can3x1->Divide(3, 1);

	can3x1->cd(1);
	gPad->SetLogy();
	gPad->SetRightMargin(0.15);
	hinput0->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
	hinput0->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
	hinput0->GetYaxis()->SetTitleOffset(1.3);
	hinput0->SetTitle("JP0");
	hinput0->SetLineColor(2);
	hinput0->SetFillColorAlpha(2, 0.25);
	hinput0->Draw("hist E");

	hbkg0->SetLineColor(4);
	hbkg0->SetFillColorAlpha(4, 0.25);
	hbkg0->Draw("hist E same");

	hdat_bkg0->SetLineColor(3);
	hdat_bkg0->Draw("hist E same");

	can3x1->cd(2);
	gPad->SetLogy();
	gPad->SetRightMargin(0.15);
	hinput1->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
	hinput1->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
	hinput1->GetYaxis()->SetTitleOffset(1.3);
	hinput1->SetTitle("JP1");
	hinput1->SetLineColor(2);
	hinput1->SetFillColorAlpha(2, 0.25);
	hinput1->Draw("hist E");

	hbkg1->SetLineColor(4);
	hbkg1->SetFillColorAlpha(4, 0.25);
	hbkg1->Draw("hist E same");

	hdat_bkg1->SetLineColor(3);
	hdat_bkg1->Draw("hist E same");

	TLegend *lgd1 = new TLegend(0.50, 0.65, 0.85, 0.9);
	lgd1->AddEntry(hinput1, " Raw Data", "f");
	lgd1->AddEntry(hbkg1, " Bkg. (from Embed.)", "f");
	lgd1->AddEntry(hdat_bkg1, " Corrected Data", "l");
	lgd1->Draw();
	gPad->Update();

	can3x1->cd(3);
	gPad->SetLogy();
	gPad->SetRightMargin(0.15);
	hinput2->GetYaxis()->SetTitle("#pi^{+}#pi^{-} Yields");
	hinput2->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
	hinput2->GetYaxis()->SetTitleOffset(1.3);
	hinput2->SetTitle("JP2");
	hinput2->SetLineColor(2);
	hinput2->SetFillColorAlpha(2, 0.25);
	hinput2->Draw("hist E");
	gPad->Update();
	hbkg2->SetLineColor(4);
	hbkg2->SetFillColorAlpha(4, 0.25);
	hbkg2->Draw("hist E same");
	gPad->Update();

	hdat_bkg2->SetLineColor(3);
	hdat_bkg2->Draw("hist E same");

	can3x1->Update();
}
