#include <iostream>
#include <TH1.h>
#include <TFile.h>
#include <TMath.h>
#include <cmath>

using namespace std;
void setBinError(TH1D *, TH1D *, TH1D *);
// print array
void print(int size, Double_t arr[]);
void trigEff()
{
	gStyle->SetOptStat(0);
	gStyle->SetLegendBorderSize(0);

	TFile *fbiased = new TFile("hist4EfficiencyParticle.root", "R");
	TFile *funbiased = new TFile("../../TriggerBias/GetUnbiased/histUnbiased.root", "R");
	assert(fbiased && funbiased);

	// in x-sec mass bins
	// get unbiased histograms
	TH1D *hMTrigAll = (TH1D *)funbiased->Get("hUnbiasedMinv");

	// get biased histograms
	TH1D *hMTrigJP0 = (TH1D *)fbiased->Get("hpairGenJP0");
	TH1D *hMTrigJP1 = (TH1D *)fbiased->Get("hpairGenJP1");
	TH1D *hMTrigJP2 = (TH1D *)fbiased->Get("hpairGenJP2");

	// trigger efficiencies
	TH1D *heff_jp0 = (TH1D *)hMTrigJP0->Clone();
	heff_jp0->Divide(hMTrigAll);
	TH1D *heff_jp1 = (TH1D *)hMTrigJP1->Clone();
	heff_jp1->Divide(hMTrigAll);
	TH1D *heff_jp2 = (TH1D *)hMTrigJP2->Clone();
	heff_jp2->Divide(hMTrigAll);

	Double_t eff_trg_jp0[13] = {0};
	Double_t efferr_trg_jp0[13] = {0};
	Double_t eff_trg_jp1[13] = {0};
	Double_t efferr_trg_jp1[13] = {0};
	Double_t eff_trg_jp2[13] = {0};
	Double_t efferr_trg_jp2[13] = {0};

	for (int i = 1; i <= heff_jp0->GetNbinsX(); i++)
	{
		// jp0
		eff_trg_jp0[i - 1] = heff_jp0->GetBinContent(i);
		efferr_trg_jp0[i - 1] = heff_jp0->GetBinError(i);
		// jp1
		eff_trg_jp1[i - 1] = heff_jp1->GetBinContent(i);
		efferr_trg_jp1[i - 1] = heff_jp1->GetBinError(i);
		// jp2
		eff_trg_jp2[i - 1] = heff_jp2->GetBinContent(i);
		efferr_trg_jp2[i - 1] = heff_jp2->GetBinError(i);
		// cout << "Bin = " << i << " eff = " << eff_trg_jp2[i - 1] << "  err = " << efferr_trg_jp2[i - 1] << endl;
	}
	cout << " // jp0 efficiency and error" << endl;
	cout << "Double_t eff_trg_jp0[13]={";
	print(13, eff_trg_jp0);
	cout << "Double_t efferr_trg_jp0[13]={";
	print(13, efferr_trg_jp0);
	cout << " // jp1 efficiency and error" << endl;
	cout << "Double_t eff_trg_jp1[13]={";
	print(13, eff_trg_jp1);
	cout << "Double_t efferr_trg_jp1[13]={";
	print(13, efferr_trg_jp1);
	cout << " // jp2 efficiency and error" << endl;
	cout << "Double_t eff_trg_jp2[13]={";
	print(13, eff_trg_jp2);
	cout << "Double_t efferr_trg_jp2[13]={";
	print(13, efferr_trg_jp2);

	drawEff(hMTrigJP0, hMTrigJP1, hMTrigJP2, hMTrigAll, "M_{inv, true}^{#pi^{+}#pi^{-}}");
}

void drawEff(TH1D *hnumJP0, TH1D *hnumJP1, TH1D *hnumJP2, TH1D *hden, const char *axis_name)
{

	TH1D *heffJP0 = (TH1D *)hnumJP0->Clone();
	heffJP0->Divide(hden);
	heffJP0->SetTitle("|v_{z}|<60 cm, V_{z, rank}>1e6");
	gStyle->SetTitleAlign(23);
	gStyle->SetTitleX(0.5);
	heffJP0->SetLineColor(2);
	heffJP0->SetLineWidth(2);
	heffJP0->SetMarkerStyle(20);
	heffJP0->SetMarkerColor(2);
	heffJP0->GetXaxis()->SetTitle(axis_name);
	heffJP0->GetYaxis()->SetTitle("Trigger Efficiency");
	heffJP0->GetYaxis()->SetTitleOffset(1.2);
	heffJP0->GetYaxis()->SetRangeUser(0., 1.2);
	setBinError(hnumJP0, hden, heffJP0);
	TH1D *heffJP1 = (TH1D *)hnumJP1->Clone();
	heffJP1->Divide(hden);
	heffJP1->SetLineColor(4);
	heffJP1->SetLineWidth(2);
	heffJP1->SetMarkerStyle(20);
	heffJP1->SetMarkerColor(4);
	setBinError(hnumJP1, hden, heffJP1);
	TH1D *heffJP2 = (TH1D *)hnumJP2->Clone();
	heffJP2->Divide(hden);
	heffJP2->SetLineColor(6);
	heffJP2->SetLineWidth(2);
	heffJP2->SetMarkerStyle(20);
	heffJP2->SetMarkerColor(6);
	setBinError(hnumJP2, hden, heffJP2);

	TCanvas *can = new TCanvas(Form("can%s", hden->GetName()), "", 600, 500);
	can->cd();
	can->SetLeftMargin(0.15);
	can->SetBottomMargin(0.15);

	heffJP0->Draw("hist E");
	heffJP1->Draw("hist E SAME");
	heffJP2->Draw("hist E SAME");

	TLegend *leg = new TLegend(0.25, 0.60, 0.45, 0.80);
	leg->AddEntry(heffJP0, "JP0", " lp");
	leg->AddEntry(heffJP1, "JP1", " lp");
	leg->AddEntry(heffJP2, "JP2", " lp");
	leg->Draw();

	can->Update();
	can->SaveAs(Form("Plots/trigEff_%s.pdf", hden->GetName()));
}

void setBinError(TH1D *hnum, TH1D *hden, TH1D *heff)
{
	// binary error propagation formula
	Double_t binError = 0;
	Double_t binCntNum = 0;
	Double_t binCntDen = 0;
	for (int i = 1; i <= hnum->GetNbinsX(); i++)
	{
		binCntNum = hnum->GetBinContent(i);
		binCntDen = hden->GetBinContent(i);
		if (binCntNum == 0 || binCntDen == 0)
			continue;
		binError = (binCntNum / binCntDen) * (1. / binCntDen) * sqrt(binCntNum * (1 - binCntNum / binCntDen));
		heff->SetBinError(i, binError);
	}
}

float getStatError_Binary(float n_num, float n_den)
{
	float binaryError = 0;
	binaryError = n_num / n_den * 1 / n_den * sqrt(n_num * (1 - n_num / n_den)); // binary error propagation
	return binaryError;
}

void print(int size, Double_t arr[])
{
	for (int i = 0; i < size; i++)
	{
		if (i < size - 1)
		{
			cout << arr[i] << ", ";
		}
		else
		{
			cout << arr[i] << "};" << endl;
			;
		}
	}
}
