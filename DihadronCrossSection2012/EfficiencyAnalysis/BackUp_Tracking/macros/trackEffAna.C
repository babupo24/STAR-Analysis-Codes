#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TGraphErrors.h"
using namespace std;
// Helper functions
// function to draw efficiency histigrams for jp0, jp1 and jp2
void drawEff(TH1D *, TH1D *, TH1D *, TH1D *, TH1D *, TH1D *, const char *);
void drawHist(TH1D *, TH1D *, TH1D *, TH1D *, TH1D *, TH1D *, const char *);
// function to calculate efficiency using histograms
Double_t getEff(TH1D *, TH1D *);
// function to calculate efficiency error (stat)
Double_t getEffErr(TH1D *, TH1D *);
// function to calculate err of ratio (a/b)
Double_t getStatError_Ratio(Double_t, Double_t, Double_t, Double_t);
Double_t getStatError_Product(Double_t, Double_t, Double_t, Double_t);

void drawEffGr(TGraphErrors *gr_m, TGraphErrors *gr_c, TGraphErrors *gr_pt, const char *, const char *);
void drawEffGrX(TGraphErrors *gr_m, TGraphErrors *gr_c, TGraphErrors *gr_pt, const char *, const char *);

// main function for efficiency analysis
void trackEffAna()
{
	gStyle->Reset("Default");
	gStyle->SetOptDate(0);
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	const int size = 20;
	// gStyle->SetTitleSize(size, "xyz");
	//  gStyle->SetLabelSize(size, "xyz");
	gStyle->SetTextSize(size);

	TH1::SetDefaultSumw2();
	TFile *fhist = new TFile("hist4EfficiencyParticle.root", "R");
	if (!fhist)
	{
		cout << "file doesn't exist! ...break" << endl;
		break;
	}

	// efficiencies and errors
	Double_t avgM[13] = {0.318715, 0.401093, 0.518795, 0.667178, 0.832331, 1.03511, 1.2368, 1.45406, 1.7209, 2.02024, 2.34893, 2.81033, 3.45647};
	// efficiency in x-section bins with the same quality cuts as in the cross section
	Double_t effx_jp0_Mpos[13] = {0};
	Double_t efferrx_jp0_Mpos[13] = {0};
	Double_t effx_jp0_Mneg[13] = {0};
	Double_t efferrx_jp0_Mneg[13] = {0};
	Double_t effx_jp1_Mpos[13] = {0};
	Double_t efferrx_jp1_Mpos[13] = {0};
	Double_t effx_jp1_Mneg[13] = {0};
	Double_t efferrx_jp1_Mneg[13] = {0};
	Double_t effx_jp2_Mpos[13] = {0};
	Double_t efferrx_jp2_Mpos[13] = {0};
	Double_t effx_jp2_Mneg[13] = {0};
	Double_t efferrx_jp2_Mneg[13] = {0};

	Double_t pair_trkeffx_jp0[13] = {0};
	Double_t pair_trkefferrx_jp0[13] = {0};
	Double_t pair_trkeffx_jp1[13] = {0};
	Double_t pair_trkefferrx_jp1[13] = {0};
	Double_t pair_trkeffx_jp2[13] = {0};
	Double_t pair_trkefferrx_jp2[13] = {0};

	Double_t eff_jp0_Mpos[13] = {0};
	Double_t efferr_jp0_Mpos[13] = {0};
	Double_t eff_jp0_Mneg[13] = {0};
	Double_t efferr_jp0_Mneg[13] = {0};
	Double_t eff_jp1_Mpos[13] = {0};
	Double_t efferr_jp1_Mpos[13] = {0};
	Double_t eff_jp1_Mneg[13] = {0};
	Double_t efferr_jp1_Mneg[13] = {0};
	Double_t eff_jp2_Mpos[13] = {0};
	Double_t efferr_jp2_Mpos[13] = {0};
	Double_t eff_jp2_Mneg[13] = {0};
	Double_t efferr_jp2_Mneg[13] = {0};

	Double_t pair_trkeff_jp0[13] = {0};
	Double_t pair_trkefferr_jp0[13] = {0};
	Double_t pair_trkeff_jp1[13] = {0};
	Double_t pair_trkefferr_jp1[13] = {0};
	Double_t pair_trkeff_jp2[13] = {0};
	Double_t pair_trkefferr_jp2[13] = {0};

	Double_t eff_all[13] = {0};
	Double_t efferr_all[13] = {0};
	Double_t eff_all2[13] = {0}; // squared for pair
	Double_t efferr_all2[13] = {0};

	const int nxbins = 13;

	// histograms for eff in x-sec bins
	//  jp0
	TH1D *hGenRecJP0PosMx[nxbins];
	TH1D *hGenJP0PosMx[nxbins];
	TH1D *hGenRecJP0NegMx[nxbins];
	TH1D *hGenJP0NegMx[nxbins];

	TH1D *hGenRecJP1PosMx[nxbins];
	TH1D *hGenJP1PosMx[nxbins];
	TH1D *hGenRecJP1NegMx[nxbins];
	TH1D *hGenJP1NegMx[nxbins];

	TH1D *hGenRecJP2PosMx[nxbins];
	TH1D *hGenJP2PosMx[nxbins];
	TH1D *hGenRecJP2NegMx[nxbins];
	TH1D *hGenJP2NegMx[nxbins];

	TH1D *hGenRecAll[nxbins];
	TH1D *hGenAll[nxbins];

	TH1D *heff_all[13];

	// with the acceptance effect at the track level (loose cuts at gen level)
	TH1D *hGenRecJP0Pos[nxbins];
	TH1D *hGenJP0Pos[nxbins];
	TH1D *hGenRecJP0Neg[nxbins];
	TH1D *hGenJP0Neg[nxbins];

	TH1D *hGenRecJP1Pos[nxbins];
	TH1D *hGenJP1Pos[nxbins];
	TH1D *hGenRecJP1Neg[nxbins];
	TH1D *hGenJP1Neg[nxbins];

	TH1D *hGenRecJP2Pos[nxbins];
	TH1D *hGenJP2Pos[nxbins];
	TH1D *hGenRecJP2Neg[nxbins];
	TH1D *hGenJP2Neg[nxbins];

	// pair level tracking efficiency
	// with the acceptance effect at the pair level (loose cuts at gen level)
	TH1D *hpairGenJP0[13];
	TH1D *hpairRecJP0[13];
	TH1D *hpairEffJP0[13];

	TH1D *hpairGenJP1[13];
	TH1D *hpairRecJP1[13];
	TH1D *hpairEffJP1[13];

	TH1D *hpairGenJP2[13];
	TH1D *hpairRecJP2[13];
	TH1D *hpairEffJP2[13];

	// same cuts at at the gen level
	TH1D *hpairGenJP0X[13];
	TH1D *hpairRecJP0X[13];
	TH1D *hpairEffJP0X[13];

	TH1D *hpairGenJP1X[13];
	TH1D *hpairRecJP1X[13];
	TH1D *hpairEffJP1X[13];

	TH1D *hpairGenJP2X[13];
	TH1D *hpairRecJP2X[13];
	TH1D *hpairEffJP2X[13];

	Double_t pairEffJP0[13] = {0}, pairEffJP1[13] = {0}, pairEffJP2[13] = {0};
	Double_t pairEffErrJP0[13] = {0}, pairEffErrJP1[13] = {0}, pairEffErrJP2[13] = {0};

	Double_t pairEffJP0X[13] = {0}, pairEffJP1X[13] = {0}, pairEffJP2X[13] = {0};
	Double_t pairEffErrJP0X[13] = {0}, pairEffErrJP1X[13] = {0}, pairEffErrJP2X[13] = {0};

	//------

	// prescale values for combining triggers
	Double_t presJP0 = 0.00707; //(1. / 141.35);
	Double_t presJP1 = 0.3978;	//(1. / 2.514);
	Double_t presJP2 = 1.0;

	// hist for eff in cross section mass bins
	for (int i = 0; i < nxbins; i++)
	{
		hGenRecAll[i] = (TH1D *)fhist->Get(Form("hGenRecpT_MxAll_%i", i));
		hGenAll[i] = (TH1D *)fhist->Get(Form("hGenpT_MxAll_%i", i));
		// jp0
		hGenRecJP0PosMx[i] = (TH1D *)fhist->Get(Form("hGenRecpT_Mx_JP0Pos%i", i));
		hGenJP0PosMx[i] = (TH1D *)fhist->Get(Form("hGenpT_Mx_JP0Pos%i", i));
		hGenRecJP0NegMx[i] = (TH1D *)fhist->Get(Form("hGenRecpT_Mx_JP0Neg%i", i));
		hGenJP0NegMx[i] = (TH1D *)fhist->Get(Form("hGenpT_Mx_JP0Neg%i", i));

		hGenRecJP0Pos[i] = (TH1D *)fhist->Get(Form("hGenRecpT_JP0Pos%i", i));
		hGenJP0Pos[i] = (TH1D *)fhist->Get(Form("hGenpT_JP0Pos%i", i));
		hGenRecJP0Neg[i] = (TH1D *)fhist->Get(Form("hGenRecpT_JP0Neg%i", i));
		hGenJP0Neg[i] = (TH1D *)fhist->Get(Form("hGenpT_JP0Neg%i", i));
		// jp1
		hGenRecJP1PosMx[i] = (TH1D *)fhist->Get(Form("hGenRecpT_Mx_JP1Pos%i", i));
		hGenJP1PosMx[i] = (TH1D *)fhist->Get(Form("hGenpT_Mx_JP1Pos%i", i));
		hGenRecJP1NegMx[i] = (TH1D *)fhist->Get(Form("hGenRecpT_Mx_JP1Neg%i", i));
		hGenJP1NegMx[i] = (TH1D *)fhist->Get(Form("hGenpT_Mx_JP1Neg%i", i));

		hGenRecJP1Pos[i] = (TH1D *)fhist->Get(Form("hGenRecpT_JP1Pos%i", i));
		hGenJP1Pos[i] = (TH1D *)fhist->Get(Form("hGenpT_JP1Pos%i", i));
		hGenRecJP1Neg[i] = (TH1D *)fhist->Get(Form("hGenRecpT_JP1Neg%i", i));
		hGenJP1Neg[i] = (TH1D *)fhist->Get(Form("hGenpT_JP1Neg%i", i));

		// jp2
		hGenRecJP2PosMx[i] = (TH1D *)fhist->Get(Form("hGenRecpT_Mx_JP2Pos%i", i));
		hGenJP2PosMx[i] = (TH1D *)fhist->Get(Form("hGenpT_Mx_JP2Pos%i", i));
		hGenRecJP2NegMx[i] = (TH1D *)fhist->Get(Form("hGenRecpT_Mx_JP2Neg%i", i));
		hGenJP2NegMx[i] = (TH1D *)fhist->Get(Form("hGenpT_Mx_JP2Neg%i", i));

		hGenRecJP2Pos[i] = (TH1D *)fhist->Get(Form("hGenRecpT_JP2Pos%i", i));
		hGenJP2Pos[i] = (TH1D *)fhist->Get(Form("hGenpT_JP2Pos%i", i));
		hGenRecJP2Neg[i] = (TH1D *)fhist->Get(Form("hGenRecpT_JP2Neg%i", i));
		hGenJP2Neg[i] = (TH1D *)fhist->Get(Form("hGenpT_JP2Neg%i", i));

		// pair level histograms
		hpairGenJP0X[i] = (TH1D *)fhist->Get(Form("hpairGenJP0X_%i", i));
		hpairRecJP0X[i] = (TH1D *)fhist->Get(Form("hpairRecJP0X_%i", i));

		hpairGenJP1X[i] = (TH1D *)fhist->Get(Form("hpairGenJP1X_%i", i));
		hpairRecJP1X[i] = (TH1D *)fhist->Get(Form("hpairRecJP1X_%i", i));

		hpairGenJP2X[i] = (TH1D *)fhist->Get(Form("hpairGenJP2X_%i", i));
		hpairRecJP2X[i] = (TH1D *)fhist->Get(Form("hpairRecJP2X_%i", i));

		hpairGenJP0[i] = (TH1D *)fhist->Get(Form("hpairGenJP0_%i", i));
		hpairRecJP0[i] = (TH1D *)fhist->Get(Form("hpairRecJP0_%i", i));

		hpairGenJP1[i] = (TH1D *)fhist->Get(Form("hpairGenJP1_%i", i));
		hpairRecJP1[i] = (TH1D *)fhist->Get(Form("hpairRecJP1_%i", i));

		hpairGenJP2[i] = (TH1D *)fhist->Get(Form("hpairGenJP2_%i", i));
		hpairRecJP2[i] = (TH1D *)fhist->Get(Form("hpairRecJP2_%i", i));
	}

	// eff in xsec bins
	for (Int_t i = 0; i < nxbins; i++)
	{

		// pair level efficiency
		hpairEffJP0[i] = (TH1D *)hpairRecJP0[i]->Clone();
		hpairEffJP0[i]->Divide(hpairGenJP0[i]);
		pairEffJP0[i] = getEff(hpairRecJP0[i], hpairGenJP0[i]);
		pairEffErrJP0[i] = getEffErr(hpairRecJP0[i], hpairGenJP0[i]);

		hpairEffJP1[i] = (TH1D *)hpairRecJP1[i]->Clone();
		hpairEffJP1[i]->Divide(hpairGenJP1[i]);
		pairEffJP1[i] = getEff(hpairRecJP1[i], hpairGenJP1[i]);
		pairEffErrJP1[i] = getEffErr(hpairRecJP1[i], hpairGenJP1[i]);

		hpairEffJP2[i] = (TH1D *)hpairRecJP2[i]->Clone();
		hpairEffJP2[i]->Divide(hpairGenJP2[i]);
		pairEffJP2[i] = getEff(hpairRecJP2[i], hpairGenJP2[i]);
		pairEffErrJP2[i] = getEffErr(hpairRecJP2[i], hpairGenJP2[i]);

		hpairEffJP0X[i] = (TH1D *)hpairRecJP0X[i]->Clone();
		hpairEffJP0X[i]->Divide(hpairGenJP0X[i]);
		pairEffJP0X[i] = getEff(hpairRecJP0X[i], hpairGenJP0X[i]);
		pairEffErrJP0X[i] = getEffErr(hpairRecJP0X[i], hpairGenJP0X[i]);

		hpairEffJP1X[i] = (TH1D *)hpairRecJP1X[i]->Clone();
		hpairEffJP1X[i]->Divide(hpairGenJP1X[i]);
		pairEffJP1X[i] = getEff(hpairRecJP1X[i], hpairGenJP1X[i]);
		pairEffErrJP1X[i] = getEffErr(hpairRecJP1X[i], hpairGenJP1X[i]);

		hpairEffJP2X[i] = (TH1D *)hpairRecJP2X[i]->Clone();
		hpairEffJP2X[i]->Divide(hpairGenJP2X[i]);
		pairEffJP2X[i] = getEff(hpairRecJP2X[i], hpairGenJP2X[i]);
		pairEffErrJP2X[i] = getEffErr(hpairRecJP2X[i], hpairGenJP2X[i]);

		eff_all[i] = getEff(hGenRecAll[i], hGenAll[i]);
		efferr_all[i] = getEffErr(hGenRecAll[i], hGenAll[i]);

		eff_all2[i] = pow(eff_all[i], 2);
		efferr_all2[i] = getStatError_Product(eff_all[i], efferr_all[i], eff_all[i], efferr_all[i]);

		// jp0
		//+ve pion
		effx_jp0_Mpos[i] = getEff(hGenRecJP0PosMx[i], hGenJP0PosMx[i]);
		efferrx_jp0_Mpos[i] = getEffErr(hGenRecJP0PosMx[i], hGenJP0PosMx[i]);
		//-ve pion
		effx_jp0_Mneg[i] = getEff(hGenRecJP0NegMx[i], hGenJP0NegMx[i]);
		efferrx_jp0_Mneg[i] = getEffErr(hGenRecJP0NegMx[i], hGenJP0NegMx[i]);
		// jp1
		//+ve pion
		effx_jp1_Mpos[i] = getEff(hGenRecJP1PosMx[i], hGenJP1PosMx[i]);
		efferrx_jp1_Mpos[i] = getEffErr(hGenRecJP1PosMx[i], hGenJP1PosMx[i]);
		//-ve pion
		effx_jp1_Mneg[i] = getEff(hGenRecJP1NegMx[i], hGenJP1NegMx[i]);
		efferrx_jp1_Mneg[i] = getEffErr(hGenRecJP1NegMx[i], hGenJP1NegMx[i]);
		// jp2
		//+ve pion
		effx_jp2_Mpos[i] = getEff(hGenRecJP2PosMx[i], hGenJP2PosMx[i]);
		efferrx_jp2_Mpos[i] = getEffErr(hGenRecJP2PosMx[i], hGenJP2PosMx[i]);
		//-ve pion
		effx_jp2_Mneg[i] = getEff(hGenRecJP2NegMx[i], hGenJP2NegMx[i]);
		efferrx_jp2_Mneg[i] = getEffErr(hGenRecJP2NegMx[i], hGenJP2NegMx[i]);
		// product +ve and -ve track efficiency in mass bins
		// jp0
		pair_trkeffx_jp0[i] = (Double_t)(effx_jp0_Mneg[i] * effx_jp0_Mpos[i]);
		pair_trkefferrx_jp0[i] = getStatError_Product(effx_jp0_Mneg[i], efferrx_jp0_Mneg[i], effx_jp0_Mpos[i], efferrx_jp0_Mpos[i]);
		// jp1
		pair_trkeffx_jp1[i] = (Double_t)(effx_jp1_Mneg[i] * effx_jp1_Mpos[i]);
		pair_trkefferrx_jp1[i] = getStatError_Product(effx_jp1_Mneg[i], efferrx_jp1_Mneg[i], effx_jp1_Mpos[i], efferrx_jp1_Mpos[i]);
		// jp2
		pair_trkeffx_jp2[i] = (Double_t)(effx_jp2_Mneg[i] * effx_jp2_Mpos[i]);
		pair_trkefferrx_jp2[i] = getStatError_Product(effx_jp2_Mneg[i], efferrx_jp2_Mneg[i], effx_jp2_Mpos[i], efferrx_jp2_Mpos[i]);
		// jp0
		//+ve pion
		eff_jp0_Mpos[i] = getEff(hGenRecJP0Pos[i], hGenJP0Pos[i]);
		efferr_jp0_Mpos[i] = getEffErr(hGenRecJP0Pos[i], hGenJP0Pos[i]);
		//-ve pion
		eff_jp0_Mneg[i] = getEff(hGenRecJP0Neg[i], hGenJP0Neg[i]);
		efferr_jp0_Mneg[i] = getEffErr(hGenRecJP0Neg[i], hGenJP0Neg[i]);
		// jp1
		//+ve pion
		eff_jp1_Mpos[i] = getEff(hGenRecJP1Pos[i], hGenJP1Pos[i]);
		efferr_jp1_Mpos[i] = getEffErr(hGenRecJP1Pos[i], hGenJP1Pos[i]);
		//-ve pion
		eff_jp1_Mneg[i] = getEff(hGenRecJP1Neg[i], hGenJP1Neg[i]);
		efferr_jp1_Mneg[i] = getEffErr(hGenRecJP1Neg[i], hGenJP1Neg[i]);
		// jp2
		//+ve pion
		eff_jp2_Mpos[i] = getEff(hGenRecJP2Pos[i], hGenJP2Pos[i]);
		efferr_jp2_Mpos[i] = getEffErr(hGenRecJP2Pos[i], hGenJP2Pos[i]);
		//-ve pion
		eff_jp2_Mneg[i] = getEff(hGenRecJP2Neg[i], hGenJP2Neg[i]);
		efferr_jp2_Mneg[i] = getEffErr(hGenRecJP2Neg[i], hGenJP2Neg[i]);
		// product +ve and -ve track efficiency in mass bins
		// jp0
		pair_trkeff_jp0[i] = (Double_t)(eff_jp0_Mneg[i] * eff_jp0_Mpos[i]);
		pair_trkefferr_jp0[i] = getStatError_Product(eff_jp0_Mneg[i], efferr_jp0_Mneg[i], eff_jp0_Mpos[i], efferr_jp0_Mpos[i]);
		// jp1
		pair_trkeff_jp1[i] = (Double_t)(eff_jp1_Mneg[i] * eff_jp1_Mpos[i]);
		pair_trkefferr_jp1[i] = getStatError_Product(eff_jp1_Mneg[i], efferr_jp1_Mneg[i], eff_jp1_Mpos[i], efferr_jp1_Mpos[i]);
		// jp2
		pair_trkeff_jp2[i] = (Double_t)(eff_jp2_Mneg[i] * eff_jp2_Mpos[i]);
		pair_trkefferr_jp2[i] = getStatError_Product(eff_jp2_Mneg[i], efferr_jp2_Mneg[i], eff_jp2_Mpos[i], efferr_jp2_Mpos[i]);
	}

	// print values

	const char *val_name = "trkeff_jp0";
	printArrayValues(13, pair_trkeffx_jp0, val_name);
	const char *val_name = "trkefferr_jp0";
	printArrayValues(13, pair_trkefferrx_jp0, val_name);
	const char *val_name = "trkeff_jp1";
	printArrayValues(13, pair_trkeffx_jp1, val_name);
	const char *val_name = "trkefferr_jp1";
	printArrayValues(13, pair_trkefferrx_jp1, val_name);
	const char *val_name = "trkeff_jp2";
	printArrayValues(13, pair_trkeffx_jp2, val_name);
	const char *val_name = "trkefferr_jp2";
	printArrayValues(13, pair_trkefferrx_jp2, val_name);

	/*
		// pair level efficiency
		TH1D *h_peff0[13];
		TH1D *h_peff1[13];
		TH1D *h_peff2[13];
		TCanvas *cpair44 = new TCanvas("cpair44", "", 900, 700);
		cpair44->Divide(4, 4);
		for (int i = 0; i < 13; i++)
		{
			cpair44->cd(i + 1);
			h_peff0[i] = (TH1D *)hpairEffJP0[i]->Clone();
			h_peff0[i]->SetLineColor(2);
			h_peff0[i]->SetLineWidth(2);
			h_peff0[i]->GetYaxis()->SetTitle("Efficiency");
			h_peff0[i]->GetYaxis()->SetLabelSize(0.08);
			h_peff0[i]->GetYaxis()->SetRangeUser(0.0, 1.10);
			h_peff0[i]->GetXaxis()->SetLabelSize(0.08);
			h_peff0[i]->SetTitle(Form("Eff. vs p_{T} Mass Bin %i", i + 1));
			h_peff0[i]->SetTitleSize(0.12);
			h_peff0[i]->Draw("hist E");

			h_peff1[i] = (TH1D *)hpairEffJP1[i]->Clone();
			h_peff1[i]->SetLineColor(3);
			h_peff1[i]->SetLineWidth(2);
			h_peff1[i]->SetLineStyle(2);
			h_peff1[i]->Draw("same");

			h_peff2[i] = (TH1D *)hpairEffJP2[i]->Clone();
			h_peff2[i]->SetLineColor(4);
			h_peff2[i]->SetLineWidth(2);
			h_peff2[i]->SetLineStyle(4);
			h_peff2[i]->Draw("same");

			gPad->Update();
		}
		cpair44->Update();
		cpair44->SaveAs("Plots/pairEffvsM_mbins.pdf");

			// pi+- combined eff per mass bin
			TH1D *c_heff_all[13];
			TCanvas *c4x4 = new TCanvas("c4x4", "", 900, 700);
			c4x4->Divide(4, 4);
			for (int i = 0; i < 13; i++)
			{
				c4x4->cd(i + 1);
				// drawHist(hGenRecJP2PosMx[i], hGenJP2PosMx[i], hGenRecJP1PosMx[i], hGenJP1PosMx[i], hGenRecJP0PosMx[i], hGenJP0PosMx[i], Form("#pi^{+} Bin %i", i + 1));
				c_heff_all[i] = (TH1D *)heff_all[i]->Clone();
				c_heff_all[i]->SetLineColor(2);
				c_heff_all[i]->SetLineWidth(2);
				// c_heff_all[i]->GetXaxis()->SetTitle("p_{T,particle} GeV/c");
				c_heff_all[i]->GetYaxis()->SetTitle("Efficiency");
				c_heff_all[i]->GetYaxis()->SetLabelSize(0.08);
				c_heff_all[i]->GetYaxis()->SetRangeUser(0.0, 1.10);
				c_heff_all[i]->GetXaxis()->SetLabelSize(0.08);
				c_heff_all[i]->SetTitle(Form("Eff. vs p_{T} Mass Bin %i", i + 1));
				c_heff_all[i]->Draw("hist E");
				c_heff_all[i]->SetTitleSize(0.12);
				gPad->Update();
			}
			c4x4->Update();
			c4x4->SaveAs("Plots/effvspt_mbins.pdf");
			*/
	// TGraphErrors
	Double_t xbinerrx[13] = {0};

	// in 13 bins
	TGraphErrors *grx_all = new TGraphErrors(13, avgM, eff_all, xbinerrx, efferr_all);
	TGraphErrors *grx_all2 = new TGraphErrors(13, avgM, eff_all2, xbinerrx, efferr_all2);

	TGraphErrors *grx_jp0Mpos = new TGraphErrors(13, avgM, effx_jp0_Mpos, xbinerrx, efferrx_jp0_Mpos);
	TGraphErrors *grx_jp0Mneg = new TGraphErrors(13, avgM, effx_jp0_Mneg, xbinerrx, efferrx_jp0_Mneg);
	TGraphErrors *grxPair_jp0M = new TGraphErrors(13, avgM, pair_trkeffx_jp0, xbinerrx, pair_trkefferrx_jp0);

	TGraphErrors *grx_jp1Mpos = new TGraphErrors(13, avgM, effx_jp1_Mpos, xbinerrx, efferrx_jp1_Mpos);
	TGraphErrors *grx_jp1Mneg = new TGraphErrors(13, avgM, effx_jp1_Mneg, xbinerrx, efferrx_jp1_Mneg);
	TGraphErrors *grxPair_jp1M = new TGraphErrors(13, avgM, pair_trkeffx_jp1, xbinerrx, pair_trkefferrx_jp1);

	TGraphErrors *grx_jp2Mpos = new TGraphErrors(13, avgM, effx_jp2_Mpos, xbinerrx, efferrx_jp2_Mpos);
	TGraphErrors *grx_jp2Mneg = new TGraphErrors(13, avgM, effx_jp2_Mneg, xbinerrx, efferrx_jp2_Mneg);
	TGraphErrors *grxPair_jp2M = new TGraphErrors(13, avgM, pair_trkeffx_jp2, xbinerrx, pair_trkefferrx_jp2);

	TGraphErrors *grx_pairEff_jp2 = new TGraphErrors(13, avgM, pairEffJP2X, xbinerrx, pairEffErrJP2X);
	TGraphErrors *grx_pairEff_jp1 = new TGraphErrors(13, avgM, pairEffJP1X, xbinerrx, pairEffErrJP1X);
	TGraphErrors *grx_pairEff_jp0 = new TGraphErrors(13, avgM, pairEffJP0X, xbinerrx, pairEffErrJP0X);

	// with the acceptance effect
	TGraphErrors *gr_jp0Mpos = new TGraphErrors(13, avgM, eff_jp0_Mpos, xbinerrx, efferr_jp0_Mpos);
	TGraphErrors *gr_jp0Mneg = new TGraphErrors(13, avgM, eff_jp0_Mneg, xbinerrx, efferr_jp0_Mneg);
	TGraphErrors *gr_Pair_jp0M = new TGraphErrors(13, avgM, pair_trkeff_jp0, xbinerrx, pair_trkefferr_jp0);

	TGraphErrors *gr_jp1Mpos = new TGraphErrors(13, avgM, eff_jp1_Mpos, xbinerrx, efferr_jp1_Mpos);
	TGraphErrors *gr_jp1Mneg = new TGraphErrors(13, avgM, eff_jp1_Mneg, xbinerrx, efferr_jp1_Mneg);
	TGraphErrors *gr_Pair_jp1M = new TGraphErrors(13, avgM, pair_trkeff_jp1, xbinerrx, pair_trkefferr_jp1);

	TGraphErrors *gr_jp2Mpos = new TGraphErrors(13, avgM, eff_jp2_Mpos, xbinerrx, efferr_jp2_Mpos);
	TGraphErrors *gr_jp2Mneg = new TGraphErrors(13, avgM, eff_jp2_Mneg, xbinerrx, efferr_jp2_Mneg);
	TGraphErrors *gr_Pair_jp2M = new TGraphErrors(13, avgM, pair_trkeff_jp2, xbinerrx, pair_trkefferr_jp2);

	TGraphErrors *gr_pairEff_jp2 = new TGraphErrors(13, avgM, pairEffJP2, xbinerrx, pairEffErrJP2);
	TGraphErrors *gr_pairEff_jp1 = new TGraphErrors(13, avgM, pairEffJP1, xbinerrx, pairEffErrJP1);
	TGraphErrors *gr_pairEff_jp0 = new TGraphErrors(13, avgM, pairEffJP0, xbinerrx, pairEffErrJP0);

	// draw eff vs bin id in 13 xsec bins
	TCanvas *cx_gr = new TCanvas("cx_gr", "", 700, 700);
	cx_gr->Divide(2, 2);
	cx_gr->cd(1);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(grx_jp0Mpos, grx_jp1Mpos, grx_jp2Mpos, " #pi^{+}", "Mass ");
	TLegend *lgx = new TLegend(0.3, 0.25, 0.7, 0.3);
	lgx->SetNColumns(3);
	lgx->AddEntry(grx_jp0Mpos, " JP0", " lp");
	lgx->AddEntry(grx_jp1Mpos, " JP1", " lp");
	lgx->AddEntry(grx_jp2Mpos, " JP2", " lp");
	lgx->Draw();
	gPad->SetGrid(0, 0);
	gPad->Update();
	cx_gr->cd(2);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(grx_jp0Mneg, grx_jp1Mneg, grx_jp2Mneg, " #pi^{-}", "Mass ");

	cx_gr->cd(3);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(grxPair_jp0M, grxPair_jp1M, grxPair_jp2M, "#pi^{+} x #pi^{-}", "Mass ");
	cx_gr->Update();
	cx_gr->SaveAs("Plots/graph_xbin_chargeseparated.pdf");

	cx_gr->cd(4);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	grx_all->SetMarkerStyle(20);
	drawEffGrX(grx_pairEff_jp0, grx_pairEff_jp1, grx_pairEff_jp2, " #pi^{+}#pi^{-} Efficiency", " Mass ");
	cx_gr->Update();
	cx_gr->SaveAs("Plots/graph_xbin_chargetrigcombined.pdf");

	// draw eff vs mass with acceptance
	TCanvas *c_gr = new TCanvas("c_gr", "", 700, 700);
	c_gr->Divide(2, 2);
	c_gr->cd(1);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGr(gr_jp0Mpos, gr_jp1Mpos, gr_jp2Mpos, " #pi^{+}", "Mass ");
	TLegend *lgx = new TLegend(0.3, 0.25, 0.7, 0.3);
	lgx->SetNColumns(3);
	lgx->AddEntry(gr_jp0Mpos, " JP0", " lp");
	lgx->AddEntry(gr_jp1Mpos, " JP1", " lp");
	lgx->AddEntry(gr_jp2Mpos, " JP2", " lp");
	lgx->Draw();
	gPad->SetGrid(0, 0);
	gPad->Update();
	c_gr->cd(2);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGr(gr_jp0Mneg, gr_jp1Mneg, gr_jp2Mneg, " #pi^{-}", "Mass ");

	c_gr->cd(3);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGr(gr_Pair_jp0M, gr_Pair_jp1M, gr_Pair_jp2M, "#pi^{+} x #pi^{-}", "Mass ");
	c_gr->Update();
	c_gr->SaveAs("Plots/graph_withacceptance.pdf");

	c_gr->cd(4);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);

	drawEffGr(gr_pairEff_jp0, gr_pairEff_jp1, gr_pairEff_jp2, " #pi^{+}#pi^{-} Efficiency", " Mass ");
	c_gr->Update();
	c_gr->SaveAs("Plots/graph_xbin_withacceptance.pdf");

} // main

Double_t getEff(TH1D *hnum, TH1D *hden)
{
	Double_t err_num = 0;
	Double_t err_den = 0;
	Double_t int_num = 0;
	Double_t int_den = 0;
	Double_t eff = 0;
	Double_t efferr = 0;

	int_num = hnum->IntegralAndError(1, hnum->GetNbinsX(), err_num, " ");
	int_den = hden->IntegralAndError(1, hden->GetNbinsX(), err_den, " ");
	eff = (Double_t)int_num / int_den;
	efferr = getStatError_Binary(int_num, int_den);
	return eff;
}
Double_t getEffErr(TH1D *hnum, TH1D *hden)
{
	Double_t err_num = 0;
	Double_t err_den = 0;
	Double_t int_num = 0;
	Double_t int_den = 0;
	Double_t eff = 0;
	Double_t efferr = 0;

	int_num = hnum->IntegralAndError(1, hnum->GetNbinsX(), err_num, " ");
	int_den = hden->IntegralAndError(1, hden->GetNbinsX(), err_den, " ");
	eff = (Double_t)int_num / int_den;
	efferr = getStatError_Binary(int_num, int_den);
	return efferr;
}

void drawEff(TH1D *hnum2, TH1D *hden2, TH1D *hnum1, TH1D *hden1, TH1D *hnum0, TH1D *hden0, const char *cname)
{

	TCanvas *can = new TCanvas(Form("can_%s", cname), "", 600, 500);
	can->cd();
	gPad->SetLeftMargin(.15);
	gPad->SetBottomMargin(.15);
	TH1D *heff2 = (TH1D *)hnum2->Clone();
	heff2->Divide(hden2);
	// heff2->SetTitle("|#eta_{thrown}|<1.5 && |#eta_{detector}|<1.5 && track QA cuts");
	heff2->SetTitle(Form("%s", cname));
	heff2->GetXaxis()->SetTitle("#eta_{thrown}");
	heff2->GetYaxis()->SetTitle("Tracking Efficiency");
	heff2->GetYaxis()->SetRangeUser(0., 1.2);
	heff2->SetLineColor(6);
	heff2->SetLineWidth(1.5);
	heff2->Draw("hist E");

	TH1D *heff1 = (TH1D *)hnum1->Clone();
	heff1->Divide(hden1);
	heff1->SetLineColor(4);
	heff1->SetLineWidth(1.5);
	heff1->Draw("hist E same");

	TH1D *heff0 = (TH1D *)hnum0->Clone();
	heff0->Divide(hden0);
	heff0->SetLineColor(2);
	heff0->SetLineWidth(1.5);
	heff0->Draw("hist E same");

	TLegend *leg = new TLegend(0.75, 0.7, 0.9, 0.9);
	leg->AddEntry(heff2, " JP2 ", "l");
	leg->AddEntry(heff1, " JP1 ", "l");
	leg->AddEntry(heff0, " JP0 ", "l");
	leg->Draw();

	// jp2 efficiency
	Double_t errjp2_num = 0;
	Double_t errjp2_den = 0;
	Double_t intjp2_num = 0;
	Double_t intjp2_den = 0;
	Double_t effjp2 = 0;
	Double_t efferrjp2 = 0;

	intjp2_num = hnum2->IntegralAndError(1, 50, errjp2_num, " ");
	intjp2_den = hden2->IntegralAndError(1, 50, errjp2_den, " ");
	effjp2 = (Double_t)intjp2_num / intjp2_den;
	efferrjp2 = getStatError_Ratio(intjp2_num, errjp2_num, intjp2_den, errjp2_den);
	// jp1 efficiency
	Double_t errjp1_num = 0;
	Double_t errjp1_den = 0;
	Double_t intjp1_num = 0;
	Double_t intjp1_den = 0;
	Double_t effjp1 = 0;
	Double_t efferrjp1 = 0;

	intjp1_num = hnum1->IntegralAndError(1, 50, errjp1_num, " ");
	intjp1_den = hden1->IntegralAndError(1, 50, errjp1_den, " ");
	effjp1 = (Double_t)intjp1_num / intjp1_den;
	efferrjp1 = getStatError_Ratio(intjp1_num, errjp1_num, intjp1_den, errjp1_den);
	// jp0 efficiency
	Double_t errjp0_num = 0;
	Double_t errjp0_den = 0;
	Double_t intjp0_num = 0;
	Double_t intjp0_den = 0;
	Double_t effjp0 = 0;
	Double_t efferrjp0 = 0;

	intjp0_num = hnum1->IntegralAndError(1, 50, errjp0_num, " ");
	intjp0_den = hden1->IntegralAndError(1, 50, errjp0_den, " ");
	effjp0 = (Double_t)intjp0_num / intjp0_den;
	efferrjp0 = getStatError_Ratio(intjp0_num, errjp0_num, intjp0_den, errjp0_den);

	TLatex text;
	text.SetTextSize(0.04);
	text.SetTextAlign(13);
	text.DrawLatex(-1.1, 0.7, Form("#color[6]{Eff. JP2 = %5.4g #pm %5.4g}", effjp2, efferrjp2));
	text.DrawLatex(-1.1, 0.6, Form("#color[4]{Eff. JP1 = %5.4g #pm %5.4g}", effjp1, efferrjp1));
	text.DrawLatex(-1.1, 0.5, Form("#color[2]{Eff. JP0 = %5.4g #pm %5.4g}", effjp0, efferrjp0));
	can->Update();
	// can->SaveAs(Form("./EffResults/%s.pdf", cname));
}
void drawHist(TH1D *hnum2, TH1D *hden2, TH1D *hnum1, TH1D *hden1, TH1D *hnum0, TH1D *hden0, const char *hname)
{

	TH1D *heff2 = (TH1D *)hnum2->Clone();
	heff2->Divide(hden2);
	// heff2->SetTitle("|#eta_{thrown}|<1.5 && |#eta_{detector}|<1.5 && track QA cuts");
	heff2->SetTitle(Form("%s", hname));
	heff2->GetXaxis()->SetTitle("#eta_{thrown}");
	heff2->GetYaxis()->SetTitle("Tracking Efficiency");
	heff2->GetYaxis()->SetRangeUser(0., 1.2);
	heff2->SetLineColor(6);
	heff2->SetLineWidth(1.5);
	heff2->Draw("hist E");

	TH1D *heff1 = (TH1D *)hnum1->Clone();
	heff1->Divide(hden1);
	heff1->SetLineColor(4);
	heff1->SetLineWidth(1.5);
	heff1->Draw("hist E same");

	TH1D *heff0 = (TH1D *)hnum0->Clone();
	heff0->Divide(hden0);
	heff0->SetLineColor(2);
	heff0->SetLineWidth(1.5);
	heff0->Draw("hist E same");

	TLegend *leg = new TLegend(0.75, 0.7, 0.9, 0.9);
	leg->AddEntry(heff2, " JP2 ", "l");
	leg->AddEntry(heff1, " JP1 ", "l");
	leg->AddEntry(heff0, " JP0 ", "l");
	leg->Draw();

	// jp2 efficiency
	Double_t errjp2_num = 0;
	Double_t errjp2_den = 0;
	Double_t intjp2_num = 0;
	Double_t intjp2_den = 0;
	Double_t effjp2 = 0;
	Double_t efferrjp2 = 0;

	intjp2_num = hnum2->IntegralAndError(1, 50, errjp2_num, " ");
	intjp2_den = hden2->IntegralAndError(1, 50, errjp2_den, " ");
	effjp2 = (Double_t)intjp2_num / intjp2_den;
	efferrjp2 = getStatError_Ratio(intjp2_num, errjp2_num, intjp2_den, errjp2_den);
	// jp1 efficiency
	Double_t errjp1_num = 0;
	Double_t errjp1_den = 0;
	Double_t intjp1_num = 0;
	Double_t intjp1_den = 0;
	Double_t effjp1 = 0;
	Double_t efferrjp1 = 0;

	intjp1_num = hnum1->IntegralAndError(1, 50, errjp1_num, " ");
	intjp1_den = hden1->IntegralAndError(1, 50, errjp1_den, " ");
	effjp1 = (Double_t)intjp1_num / intjp1_den;
	efferrjp1 = getStatError_Ratio(intjp1_num, errjp1_num, intjp1_den, errjp1_den);
	// jp0 efficiency
	Double_t errjp0_num = 0;
	Double_t errjp0_den = 0;
	Double_t intjp0_num = 0;
	Double_t intjp0_den = 0;
	Double_t effjp0 = 0;
	Double_t efferrjp0 = 0;

	intjp0_num = hnum1->IntegralAndError(1, 50, errjp0_num, " ");
	intjp0_den = hden1->IntegralAndError(1, 50, errjp0_den, " ");
	effjp0 = (Double_t)intjp0_num / intjp0_den;
	efferrjp0 = getStatError_Ratio(intjp0_num, errjp0_num, intjp0_den, errjp0_den);

	TLatex text;
	text.SetTextSize(0.04);
	text.SetTextAlign(13);
	text.DrawLatex(-1.1, 0.7, Form("#color[6]{Eff. JP2 = %5.4g #pm %5.4g}", effjp2, efferrjp2));
	text.DrawLatex(-1.1, 0.6, Form("#color[4]{Eff. JP1 = %5.4g #pm %5.4g}", effjp1, efferrjp1));
	text.DrawLatex(-1.1, 0.5, Form("#color[2]{Eff. JP0 = %5.4g #pm %5.4g}", effjp0, efferrjp0));
}

void drawEffGr(TGraphErrors *gr_m, TGraphErrors *gr_c, TGraphErrors *gr_pt, const char *gname, const char *nbin)
{
	gr_m->SetMarkerStyle(20);
	gr_m->SetMarkerColor(1);
	gr_m->SetLineColor(1);
	gr_m->SetTitle(Form("%s", gname));
	gr_m->GetXaxis()->SetTitle("<M_{inv}>");
	gr_m->GetXaxis()->SetNdivisions(505);
	gr_m->GetYaxis()->SetTitle("Trk. Eff (%)");
	gr_m->GetYaxis()->CenterTitle();
	gr_m->GetYaxis()->SetTitleOffset(1.5);
	gr_m->GetYaxis()->SetRangeUser(0.25, 1.0);
	gr_m->Draw("AP");

	gr_c->SetMarkerStyle(25);
	gr_c->SetMarkerColor(2);
	gr_c->SetLineColor(2);
	gr_c->Draw("P SAME");

	gr_pt->SetMarkerStyle(22);
	gr_pt->SetMarkerColor(4);
	gr_pt->SetLineColor(4);
	gr_pt->Draw("P SAME");
}
void drawEffGrX(TGraphErrors *gr_m, TGraphErrors *gr_c, TGraphErrors *gr_pt, const char *gname, const char *nbin)
{
	gr_m->SetMarkerStyle(20);
	gr_m->SetMarkerColor(1);
	gr_m->SetLineColor(1);
	gr_m->SetTitle(Form("%s", gname));
	gr_m->GetXaxis()->SetTitle("<M_{inv}>");
	gr_m->GetXaxis()->SetNdivisions(505);
	gr_m->GetYaxis()->SetTitle("Trk. Eff (%)");
	gr_m->GetYaxis()->CenterTitle();
	gr_m->GetYaxis()->SetTitleOffset(1.5);
	gr_m->GetYaxis()->SetRangeUser(0.5, 1.05);
	gr_m->Draw("AP");

	gr_c->SetMarkerStyle(25);
	gr_c->SetMarkerColor(2);
	gr_c->SetLineColor(2);
	gr_c->Draw("P SAME");

	gr_pt->SetMarkerStyle(22);
	gr_pt->SetMarkerColor(4);
	gr_pt->SetLineColor(4);
	gr_pt->Draw("P SAME");
}

Double_t getStatError_Ratio(Double_t num, Double_t numerr, Double_t den, Double_t denerr)
{
	Double_t statErrorRatio = 0;
	statErrorRatio = (num / den) * sqrt(pow(numerr / num, 2) + pow(denerr / den, 2));
	return statErrorRatio;
}
Double_t getStatError_Product(Double_t a, Double_t aerr, Double_t b, Double_t berr)
{
	Double_t statErrorProduct = 0;
	statErrorProduct = sqrt(pow(b * aerr, 2) + pow(a * berr, 2));
	return statErrorProduct;
}
double getStatError_Binary(double n_num, double n_den)
{
	double binaryError = 0;
	binaryError = n_num / n_den * 1 / n_den * sqrt(n_num * (1 - n_num / n_den)); // binary error propagation
	return binaryError;
}

void printArrayValues(int size, double val[], const char *valName)
{
	cout << "double " << valName << "[" << size << "] = { ";
	for (int i = 0; i < size; i++)
	{
		if (i == size - 1)
			cout << val[i] << "};" << endl;
		else
			cout << val[i] << ",";
	}
}
