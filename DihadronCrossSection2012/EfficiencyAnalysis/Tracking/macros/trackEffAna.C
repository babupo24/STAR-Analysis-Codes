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
Double_t getEff(TH1D *, TH1D *, TH1D *);
// function to calculate efficiency error (stat)
Double_t getEffErr(TH1D *, TH1D *);
// function to calculate err of ratio (a/b)
Double_t getStatError_Ratio(Double_t, Double_t, Double_t, Double_t);
Double_t getStatError_Product(Double_t, Double_t, Double_t, Double_t);

void drawEffDiffC(TH1D *h0, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5);
void drawEffGr(TH1D *gr_m, TH1D *gr_c, TH1D *gr_pt, const char *, const char *);
void drawEffGrX(TH1D *gr_m, TH1D *gr_c, TH1D *gr_pt, const char *, const char *);

// main function for efficiency analysis
void trackEffAna()
{
	// gStyle->Reset("Default");
	gStyle->SetOptDate(0);
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	const int size = 20;
	// gStyle->SetTitleSize(size, "xyz");
	//  gStyle->SetLabelSize(size, "xyz");
	gStyle->SetTextSize(size);

	TH1::SetDefaultSumw2();
	TFile *fhist = new TFile("hist4Efficiency.root", "R");
	if (!fhist)
	{
		cout << "file doesn't exist! ...break" << endl;
		break;
	}

	// efficiencies and errors
	Double_t avgM[13] = {0.318715, 0.401093, 0.518795, 0.667178, 0.832331, 1.03511, 1.2368, 1.45406, 1.7209, 2.02024, 2.34893, 2.81033, 3.45647};
	// efficiency in x-section bins with the same quality cuts as in the cross section
	Double_t effx_jp0_pos[13] = {0};
	Double_t efferrx_jp0_pos[13] = {0};
	Double_t effx_jp0_neg[13] = {0};
	Double_t efferrx_jp0_neg[13] = {0};
	Double_t effx_jp1_pos[13] = {0};
	Double_t efferrx_jp1_pos[13] = {0};
	Double_t effx_jp1_neg[13] = {0};
	Double_t efferrx_jp1_neg[13] = {0};
	Double_t effx_jp2_pos[13] = {0};
	Double_t efferrx_jp2_pos[13] = {0};
	Double_t effx_jp2_neg[13] = {0};
	Double_t efferrx_jp2_neg[13] = {0};

	Double_t pair_trkeffx_jp0[13] = {0};
	Double_t pair_trkefferrx_jp0[13] = {0};
	Double_t pair_trkeffx_jp1[13] = {0};
	Double_t pair_trkefferrx_jp1[13] = {0};
	Double_t pair_trkeffx_jp2[13] = {0};
	Double_t pair_trkefferrx_jp2[13] = {0};

	Double_t eff_jp0_pos[13] = {0};
	Double_t efferr_jp0_pos[13] = {0};
	Double_t eff_jp0_neg[13] = {0};
	Double_t efferr_jp0_neg[13] = {0};
	Double_t eff_jp1_pos[13] = {0};
	Double_t efferr_jp1_pos[13] = {0};
	Double_t eff_jp1_neg[13] = {0};
	Double_t efferr_jp1_neg[13] = {0};
	Double_t eff_jp2_pos[13] = {0};
	Double_t efferr_jp2_pos[13] = {0};
	Double_t eff_jp2_neg[13] = {0};
	Double_t efferr_jp2_neg[13] = {0};

	Double_t pair_trkeff_jp0[13] = {0};
	Double_t pair_trkefferr_jp0[13] = {0};
	Double_t pair_trkeff_jp1[13] = {0};
	Double_t pair_trkefferr_jp1[13] = {0};
	Double_t pair_trkeff_jp2[13] = {0};
	Double_t pair_trkefferr_jp2[13] = {0};

	Double_t pairEffJP0[13] = {0}, pairEffJP1[13] = {0}, pairEffJP2[13] = {0};
	Double_t pairEffErrJP0[13] = {0}, pairEffErrJP1[13] = {0}, pairEffErrJP2[13] = {0};

	Double_t pairEffJP0X[13] = {0}, pairEffJP1X[13] = {0}, pairEffJP2X[13] = {0};
	Double_t pairEffErrJP0X[13] = {0}, pairEffErrJP1X[13] = {0}, pairEffErrJP2X[13] = {0};

	const int nxbins = 13;

	// histograms for eff in x-sec bins
	const TH1D *hGenRecJP0Pos = (TH1D *)fhist->Get("hGenRecJP0Pos");
	const TH1D *hGenJP0Pos = (TH1D *)fhist->Get("hGenJP0Pos");
	const TH1D *hGenRecJP0Neg = (TH1D *)fhist->Get("hGenRecJP0Neg");
	const TH1D *hGenJP0Neg = (TH1D *)fhist->Get("hGenJP0Neg");

	const TH1D *hGenRecJP1Pos = (TH1D *)fhist->Get("hGenRecJP1Pos");
	const TH1D *hGenJP1Pos = (TH1D *)fhist->Get("hGenJP1Pos");
	const TH1D *hGenRecJP1Neg = (TH1D *)fhist->Get("hGenRecJP1Neg");
	const TH1D *hGenJP1Neg = (TH1D *)fhist->Get("hGenJP1Neg");

	const TH1D *hGenRecJP2Pos = (TH1D *)fhist->Get("hGenRecJP2Pos");
	const TH1D *hGenJP2Pos = (TH1D *)fhist->Get("hGenJP2Pos");
	const TH1D *hGenRecJP2Neg = (TH1D *)fhist->Get("hGenRecJP2Neg");
	const TH1D *hGenJP2Neg = (TH1D *)fhist->Get("hGenJP2Neg");

	const TH1D *hGenRecJP0PosX = (TH1D *)fhist->Get("hGenRecJP0PosX");
	const TH1D *hGenJP0PosX = (TH1D *)fhist->Get("hGenJP0PosX");
	const TH1D *hGenRecJP0NegX = (TH1D *)fhist->Get("hGenRecJP0NegX");
	const TH1D *hGenJP0NegX = (TH1D *)fhist->Get("hGenJP0NegX");

	const TH1D *hGenRecJP1PosX = (TH1D *)fhist->Get("hGenRecJP1PosX");
	const TH1D *hGenJP1PosX = (TH1D *)fhist->Get("hGenJP1PosX");
	const TH1D *hGenRecJP1NegX = (TH1D *)fhist->Get("hGenRecJP1NegX");
	const TH1D *hGenJP1NegX = (TH1D *)fhist->Get("hGenJP1NegX");

	const TH1D *hGenRecJP2PosX = (TH1D *)fhist->Get("hGenRecJP2PosX");
	const TH1D *hGenJP2PosX = (TH1D *)fhist->Get("hGenJP2PosX");
	const TH1D *hGenRecJP2NegX = (TH1D *)fhist->Get("hGenRecJP2NegX");
	const TH1D *hGenJP2NegX = (TH1D *)fhist->Get("hGenJP2NegX");

	// pair level tracking efficiency
	// with the acceptance effect at the pair level (loose cuts at gen level)
	const TH1D *hpairGenJP0 = (TH1D *)fhist->Get("hpairGenJP0");
	const TH1D *hpairGenRecJP0 = (TH1D *)fhist->Get("hpairGenRecJP0");
	const TH1D *hpairGenJP1 = (TH1D *)fhist->Get("hpairGenJP1");
	const TH1D *hpairGenRecJP1 = (TH1D *)fhist->Get("hpairGenRecJP1");
	const TH1D *hpairGenJP2 = (TH1D *)fhist->Get("hpairGenJP2");
	const TH1D *hpairGenRecJP2 = (TH1D *)fhist->Get("hpairGenRecJP2");

	// same cuts at at the gen level
	const TH1D *hpairGenJP0X = (TH1D *)fhist->Get("hpairGenJP0X");
	const TH1D *hpairGenRecJP0X = (TH1D *)fhist->Get("hpairGenRecJP0X");
	const TH1D *hpairGenJP1X = (TH1D *)fhist->Get("hpairGenJP1X");
	const TH1D *hpairGenRecJP1X = (TH1D *)fhist->Get("hpairGenRecJP1X");
	const TH1D *hpairGenJP2X = (TH1D *)fhist->Get("hpairGenJP2X");
	const TH1D *hpairGenRecJP2X = (TH1D *)fhist->Get("hpairGenRecJP2X");

	const TH1D *hpairGenRecJP0C = (TH1D *)fhist->Get("hGenRecJP0C");
	const TH1D *hpairGenRecJP1C = (TH1D *)fhist->Get("hGenRecJP1C");
	const TH1D *hpairGenRecJP2C = (TH1D *)fhist->Get("hGenRecJP2C");

	// pair level reconstructed with the different minimum cone cut
	const TH1D *hpairGenRecJP0XC2 = (TH1D *)fhist->Get("hGenRecJP0XC2");
	const TH1D *hpairGenRecJP0XC3 = (TH1D *)fhist->Get("hGenRecJP0XC3");
	const TH1D *hpairGenRecJP0XC4 = (TH1D *)fhist->Get("hGenRecJP0XC4");
	const TH1D *hpairGenRecJP0XC5 = (TH1D *)fhist->Get("hGenRecJP0XC5");

	const TH1D *hpairGenRecJP1XC2 = (TH1D *)fhist->Get("hGenRecJP1XC2");
	const TH1D *hpairGenRecJP1XC3 = (TH1D *)fhist->Get("hGenRecJP1XC3");
	const TH1D *hpairGenRecJP1XC4 = (TH1D *)fhist->Get("hGenRecJP1XC4");
	const TH1D *hpairGenRecJP1XC5 = (TH1D *)fhist->Get("hGenRecJP1XC5");

	const TH1D *hpairGenRecJP2XC2 = (TH1D *)fhist->Get("hGenRecJP2XC2");
	const TH1D *hpairGenRecJP2XC3 = (TH1D *)fhist->Get("hGenRecJP2XC3");
	const TH1D *hpairGenRecJP2XC4 = (TH1D *)fhist->Get("hGenRecJP2XC4");
	const TH1D *hpairGenRecJP2XC5 = (TH1D *)fhist->Get("hGenRecJP2XC5");

	// track level efficiency
	TH1D *heff_pos_jp0 = (TH1D *)hGenRecJP0Pos->Clone();
	getEff(heff_pos_jp0, hGenRecJP0Pos, hGenJP0Pos);
	TH1D *heff_neg_jp0 = (TH1D *)hGenRecJP0Neg->Clone();
	getEff(heff_neg_jp0, hGenRecJP0Neg, hGenJP0Neg);

	TH1D *heff_pos_jp1 = (TH1D *)hGenRecJP1Pos->Clone();
	getEff(heff_pos_jp1, hGenRecJP1Pos, hGenJP1Pos);
	TH1D *heff_neg_jp1 = (TH1D *)hGenRecJP1Neg->Clone();
	getEff(heff_neg_jp1, hGenRecJP1Neg, hGenJP1Neg);

	TH1D *heff_pos_jp2 = (TH1D *)hGenRecJP2Pos->Clone();
	getEff(heff_pos_jp2, hGenRecJP2Pos, hGenJP2Pos);
	TH1D *heff_neg_jp2 = (TH1D *)hGenRecJP2Neg->Clone();
	getEff(heff_neg_jp2, hGenRecJP2Neg, hGenJP2Neg);

	TH1D *heff_pos_jp0x = (TH1D *)hGenRecJP0PosX->Clone();
	getEff(heff_pos_jp0x, hGenRecJP0PosX, hGenJP0PosX);
	TH1D *heff_neg_jp0x = (TH1D *)hGenRecJP0NegX->Clone();
	getEff(heff_neg_jp0x, hGenRecJP0NegX, hGenJP0NegX);

	TH1D *heff_pos_jp1x = (TH1D *)hGenRecJP1PosX->Clone();
	getEff(heff_pos_jp1x, hGenRecJP1PosX, hGenJP1PosX);
	TH1D *heff_neg_jp1x = (TH1D *)hGenRecJP1NegX->Clone();
	getEff(heff_neg_jp1x, hGenRecJP1NegX, hGenJP1NegX);

	TH1D *heff_pos_jp2x = (TH1D *)hGenRecJP2PosX->Clone();
	getEff(heff_pos_jp2x, hGenRecJP2PosX, hGenJP2PosX);
	TH1D *heff_neg_jp2x = (TH1D *)hGenRecJP2NegX->Clone();
	getEff(heff_neg_jp2x, hGenRecJP2NegX, hGenJP2NegX);

	// product of pi+ and pi-
	TH1D *heff_product_jp0x = (TH1D *)heff_pos_jp0x->Clone();
	heff_product_jp0x->Multiply(heff_neg_jp0x);
	TH1D *heff_product_jp1x = (TH1D *)heff_pos_jp1x->Clone();
	heff_product_jp1x->Multiply(heff_neg_jp1x);
	TH1D *heff_product_jp2x = (TH1D *)heff_pos_jp2x->Clone();
	heff_product_jp2x->Multiply(heff_neg_jp2x);

	// pair level efficiency
	TH1D *heff_pair_jp0 = (TH1D *)hpairGenRecJP0->Clone();
	getEff(heff_pair_jp0, hpairGenRecJP0, hpairGenJP0);
	TH1D *heff_pair_jp1 = (TH1D *)hpairGenRecJP1->Clone();
	getEff(heff_pair_jp1, hpairGenRecJP1, hpairGenJP1);
	TH1D *heff_pair_jp2 = (TH1D *)hpairGenRecJP2->Clone();
	getEff(heff_pair_jp2, hpairGenRecJP2, hpairGenJP2);

	TH1D *heff_pair_jp0x = (TH1D *)hpairGenRecJP0X->Clone();
	getEff(heff_pair_jp0x, hpairGenRecJP0X, hpairGenJP0X);
	TH1D *heff_pair_jp1x = (TH1D *)hpairGenRecJP1X->Clone();
	getEff(heff_pair_jp1x, hpairGenRecJP1X, hpairGenJP1X);
	TH1D *heff_pair_jp2x = (TH1D *)hpairGenRecJP2X->Clone();
	getEff(heff_pair_jp2x, hpairGenRecJP2X, hpairGenJP2X);

	// pair level efficiency with the minimum cone cut at the reconstructed level
	TH1D *heff_pair_jp0xc2 = (TH1D *)hpairGenRecJP0XC2->Clone();
	getEff(heff_pair_jp0xc2, hpairGenRecJP0XC2, hpairGenJP0X);
	TH1D *heff_pair_jp0xc3 = (TH1D *)hpairGenRecJP0XC3->Clone();
	getEff(heff_pair_jp0xc3, hpairGenRecJP0XC3, hpairGenJP0X);
	TH1D *heff_pair_jp0xc4 = (TH1D *)hpairGenRecJP0XC4->Clone();
	getEff(heff_pair_jp0xc4, hpairGenRecJP0XC4, hpairGenJP0X);
	TH1D *heff_pair_jp0xc5 = (TH1D *)hpairGenRecJP0XC5->Clone();
	getEff(heff_pair_jp0xc5, hpairGenRecJP0XC5, hpairGenJP0X);

	TH1D *heff_pair_jp1xc2 = (TH1D *)hpairGenRecJP1XC2->Clone();
	getEff(heff_pair_jp1xc2, hpairGenRecJP1XC2, hpairGenJP1X);
	TH1D *heff_pair_jp1xc3 = (TH1D *)hpairGenRecJP1XC3->Clone();
	getEff(heff_pair_jp1xc3, hpairGenRecJP1XC3, hpairGenJP1X);
	TH1D *heff_pair_jp1xc4 = (TH1D *)hpairGenRecJP1XC4->Clone();
	getEff(heff_pair_jp1xc4, hpairGenRecJP1XC4, hpairGenJP1X);
	TH1D *heff_pair_jp1xc5 = (TH1D *)hpairGenRecJP1XC5->Clone();
	getEff(heff_pair_jp1xc5, hpairGenRecJP1XC5, hpairGenJP0X);

	TH1D *heff_pair_jp2xc2 = (TH1D *)hpairGenRecJP2XC2->Clone();
	getEff(heff_pair_jp2xc2, hpairGenRecJP2XC2, hpairGenJP2X);
	TH1D *heff_pair_jp2xc3 = (TH1D *)hpairGenRecJP2XC3->Clone();
	getEff(heff_pair_jp2xc3, hpairGenRecJP2XC3, hpairGenJP2X);
	TH1D *heff_pair_jp2xc4 = (TH1D *)hpairGenRecJP2XC4->Clone();
	getEff(heff_pair_jp2xc4, hpairGenRecJP2XC4, hpairGenJP2X);
	TH1D *heff_pair_jp2xc5 = (TH1D *)hpairGenRecJP2XC5->Clone();
	getEff(heff_pair_jp2xc5, hpairGenRecJP2XC5, hpairGenJP2X);

	// output efficiencies
	ofstream fout;
	fout.open("Plots/tracking_eff_MinConeCut.txt");
	fout << "double eff_track_jp0[13]={";
	for (int i = 1; i <= heff_pair_jp0xc2->GetNbinsX(); i++)
	{
		if (i == heff_pair_jp0xc2->GetNbinsX())
			fout << heff_pair_jp0xc2->GetBinContent(i) << "};" << endl;
		else
			fout << heff_pair_jp0xc2->GetBinContent(i) << ",";
	}
	fout << "double efferr_track_jp0[13]={";
	for (int i = 1; i <= heff_pair_jp0xc2->GetNbinsX(); i++)
	{
		if (i == heff_pair_jp0xc2->GetNbinsX())
			fout << heff_pair_jp0xc2->GetBinError(i) << "};" << endl;
		else
			fout << heff_pair_jp0xc2->GetBinError(i) << ",";
	}
	fout << "double eff_track_jp1[13]={";
	for (int i = 1; i <= heff_pair_jp1xc2->GetNbinsX(); i++)
	{
		if (i == heff_pair_jp1xc2->GetNbinsX())
			fout << heff_pair_jp1xc2->GetBinContent(i) << "};" << endl;
		else
			fout << heff_pair_jp1xc2->GetBinContent(i) << ",";
	}
	fout << "double efferr_track_jp1[13]={";
	for (int i = 1; i <= heff_pair_jp1xc2->GetNbinsX(); i++)
	{
		if (i == heff_pair_jp1xc2->GetNbinsX())
			fout << heff_pair_jp1xc2->GetBinError(i) << "};" << endl;
		else
			fout << heff_pair_jp1xc2->GetBinError(i) << ",";
	}
	fout << "double eff_track_jp2[13]={";
	for (int i = 1; i <= heff_pair_jp2xc2->GetNbinsX(); i++)
	{
		if (i == heff_pair_jp2xc2->GetNbinsX())
			fout << heff_pair_jp2xc2->GetBinContent(i) << "};" << endl;
		else
			fout << heff_pair_jp2xc2->GetBinContent(i) << ",";
	}
	fout << "double efferr_track_jp2[13]={";
	for (int i = 1; i <= heff_pair_jp2xc2->GetNbinsX(); i++)
	{
		if (i == heff_pair_jp2xc2->GetNbinsX())
			fout << heff_pair_jp2xc2->GetBinError(i) << "};" << endl;
		else
			fout << heff_pair_jp2xc2->GetBinError(i) << ",";
	}

	// draw eff vs bin id in 13 xsec bins
	TCanvas *cx_gr = new TCanvas("cx_gr", "", 700, 700);
	cx_gr->Divide(2, 2);
	cx_gr->cd(1);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(heff_pos_jp0x, heff_pos_jp1x, heff_pos_jp2x, " #pi^{+}", "Mass ");
	TLegend *lgx = new TLegend(0.3, 0.25, 0.7, 0.3);
	lgx->SetNColumns(3);
	lgx->AddEntry(heff_pos_jp0x, " JP0", " lep");
	lgx->AddEntry(heff_pos_jp1x, " JP1", " lep");
	lgx->AddEntry(heff_pos_jp2x, " JP2", " lep");
	lgx->Draw();
	gPad->SetGrid(0, 0);
	gPad->Update();
	cx_gr->cd(2);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(heff_neg_jp0x, heff_neg_jp1x, heff_neg_jp2x, " #pi^{-}", "Mass ");

	cx_gr->cd(3);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(heff_product_jp0x, heff_product_jp1x, heff_product_jp2x, " #pi^{+} x #pi^{-}", "Mass ");
	cx_gr->Update();
	cx_gr->SaveAs("Plots/graph_xbin_chargeseparated.pdf");

	cx_gr->cd(4);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffGrX(heff_pair_jp0x, heff_pair_jp1x, heff_pair_jp2x, " #pi^{+}#pi^{-} Efficiency", "Mass ");
	cx_gr->Update();
	cx_gr->SaveAs("Plots/graph_xbin_chargetrigcombined.pdf");

	// ratio of product pair to the pair
	TCanvas *c_r = new TCanvas("c_r", "", 500, 450);
	c_r->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	gPad->SetGrid(0, 0);
	TH1D *hr_2 = (TH1D *)heff_pair_jp2x->Clone();
	hr_2->Divide(heff_product_jp2x);
	hr_2->GetYaxis()->SetRangeUser(0.98, 1.02);
	hr_2->GetYaxis()->SetNdivisions(505);
	hr_2->GetYaxis()->SetTitle("eff(#pi^{+}#pi^{-})/(eff(#pi^{+}) x eff(#pi^{-}))");
	hr_2->GetXaxis()->SetTitle("<M_{inv}>");
	hr_2->SetMarkerStyle(22);
	hr_2->SetMarkerColor(4);
	hr_2->SetLineColor(4);
	hr_2->Draw("E");
	TH1D *hr_1 = (TH1D *)heff_pair_jp1x->Clone();
	hr_1->Divide(heff_product_jp1x);
	hr_1->SetMarkerStyle(25);
	hr_1->SetMarkerColor(2);
	hr_1->SetLineColor(2);
	hr_1->Draw("E SAME");
	TH1D *hr_0 = (TH1D *)heff_pair_jp0x->Clone();
	hr_0->Divide(heff_product_jp0x);
	hr_0->SetMarkerStyle(8);
	hr_0->SetMarkerColor(1);
	hr_0->SetLineColor(1);
	hr_0->Draw("E SAME");
	c_r->Update();
	// TLine *l1 = new TLine(c_r->GetUxmin(), 0, c_r->GetUxmax(), 0);
	TLine *l1 = new TLine(c_r->GetUxmin(), 1, c_r->GetUxmax(), 1);
	l1->SetLineColor(3);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	c_r->Update();
	c_r->SaveAs("Plots/eff_ratio.pdf");

	// ratio of pair  efficiency without and with the cone cut
	TCanvas *cc_r = new TCanvas("cc_r", "", 900, 450);
	cc_r->Divide(2, 1);
	cc_r->cd(1);
	gPad->SetGrid(0, 0);
	gPad->SetLeftMargin(0.12);
	drawEffDiffC(heff_pair_jp0x, heff_pair_jp0xc2, heff_pair_jp0xc3, heff_pair_jp0xc4, heff_pair_jp0xc5);
	TLegend *lgx = new TLegend(0.15, 0.25, 0.85, 0.3);
	lgx->SetNColumns(5);
	lgx->AddEntry(heff_pair_jp0x, " No Min cut", " lep");
	lgx->AddEntry(heff_pair_jp0xc2, " >0.02", " lep");
	lgx->AddEntry(heff_pair_jp0xc3, " >0.03", " lep");
	lgx->AddEntry(heff_pair_jp0xc4, " >0.04", " lep");
	lgx->AddEntry(heff_pair_jp0xc5, " >0.05", " lep");
	lgx->SetTextSize(0.035);
	lgx->Draw();
	cc_r->cd(2);
	// gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.12);
	gPad->SetGrid(0, 0);
	TH1D *hr_2 = (TH1D *)heff_pair_jp0x->Clone();
	hr_2->Divide(heff_pair_jp0xc2);
	hr_2->GetYaxis()->SetRangeUser(0.99, 1.05);
	hr_2->GetYaxis()->SetNdivisions(505);
	hr_2->GetYaxis()->SetTitle("Ratio");
	hr_2->GetXaxis()->SetTitle("<M_{inv}>");
	hr_2->SetMarkerStyle(25);
	hr_2->SetTitle("Ratio W/O MinCone to W/ MinCone");
	hr_2->SetMarkerColor(2);
	hr_2->SetLineColor(2);
	hr_2->Draw("E");
	TH1D *hr_3 = (TH1D *)heff_pair_jp0x->Clone();
	hr_3->Divide(heff_pair_jp0xc3);
	hr_3->SetMarkerStyle(22);
	hr_3->SetMarkerColor(4);
	hr_3->SetLineColor(4);
	hr_3->Draw("E SAME");
	TH1D *hr_4 = (TH1D *)heff_pair_jp0x->Clone();
	hr_4->Divide(heff_pair_jp0xc4);
	hr_4->SetMarkerStyle(23);
	hr_4->SetMarkerColor(3);
	hr_4->SetLineColor(3);
	hr_4->Draw("E SAME");

	TH1D *hr_5 = (TH1D *)heff_pair_jp0x->Clone();
	hr_5->Divide(heff_pair_jp0xc5);
	hr_5->SetMarkerStyle(21);
	hr_5->SetMarkerColor(6);
	hr_5->SetLineColor(6);
	hr_5->Draw("E SAME");
	gPad->Update();
	// TLine *l1 = new TLine(cc_r->GetUxmin(), 0, cc_r->GetUxmax(), 0);
	TLine *l1 = new TLine(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
	l1->SetLineColor(3);
	l1->SetLineWidth(2);
	l1->SetLineStyle(2);
	l1->Draw();
	gPad->Update();
	cc_r->Update();
	cc_r->SaveAs("Plots/eff_ratio_NoMinCone_to_minCone.pdf");

} // main

void getEff(TH1D *heff, TH1D *hnum, TH1D *hden)
{
	heff->Reset();
	for (int i = 1; i <= hnum->GetNbinsX(); i++)
	{
		double n1 = hnum->GetBinContent(i);
		double n2 = hden->GetBinContent(i);
		if (n1 == 0 || n2 == 0)
			continue;
		double eff = (Double_t)n1 / n2;
		double efferr = getStatError_Binary(n1, n2);
		heff->SetBinContent(i, eff);
		heff->SetBinError(i, efferr);
	}
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

void drawEffGr(TH1D *gr_m, TH1D *gr_c, TH1D *gr_pt, const char *gname, const char *nbin)
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
void drawEffGrX(TH1D *gr_m, TH1D *gr_c, TH1D *gr_pt, const char *gname, const char *nbin)
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
	gr_m->GetYaxis()->SetRangeUser(0.75, 1.05);
	gr_m->Draw("E");

	gr_c->SetMarkerStyle(25);
	gr_c->SetMarkerColor(2);
	gr_c->SetLineColor(2);
	gr_c->Draw("E SAME");

	gr_pt->SetMarkerStyle(22);
	gr_pt->SetMarkerColor(4);
	gr_pt->SetLineColor(4);
	gr_pt->Draw("E SAME");
}
void drawEffDiffC(TH1D *h0, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5)
{
	h0->SetMarkerStyle(20);
	h0->SetMarkerColor(1);
	h0->SetLineColor(1);
	h0->GetXaxis()->SetTitle("<M_{inv}>");
	h0->GetXaxis()->SetNdivisions(505);
	h0->GetYaxis()->SetTitle("Trk. Eff (%)");
	h0->GetYaxis()->CenterTitle();
	h0->GetYaxis()->SetTitleOffset(1.5);
	h0->GetYaxis()->SetRangeUser(0.78, 1.);
	h0->Draw("E");

	h2->SetMarkerStyle(25);
	h2->SetMarkerColor(2);
	h2->SetLineColor(2);
	h2->Draw("E SAME");

	h3->SetMarkerStyle(22);
	h3->SetMarkerColor(4);
	h3->SetLineColor(4);
	h3->Draw("E SAME");

	h4->SetMarkerStyle(23);
	h4->SetMarkerColor(3);
	h4->SetLineColor(3);
	h4->Draw("E SAME");

	h5->SetMarkerStyle(21);
	h5->SetMarkerColor(6);
	h5->SetLineColor(6);
	h5->Draw("E SAME");
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
