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

	TFile *fhist = new TFile("hist4TrigEff.root", "R");
	if (!fhist)
		break;
	// in x-sec mass bins
	TH1D *hMTrigAll = (TH1D *)fhist->Get("hMTrigAll");
	TH1D *hMTrigJP0 = (TH1D *)fhist->Get("hMTrigJP0");
	TH1D *hMTrigJP1 = (TH1D *)fhist->Get("hMTrigJP1");
	TH1D *hMTrigJP2 = (TH1D *)fhist->Get("hMTrigJP2");
	// in dipion pt bins
	TH1D *hptTrigAll = (TH1D *)fhist->Get("hptTrigAll");
	TH1D *hptTrigJP0 = (TH1D *)fhist->Get("hptTrigJP0");
	TH1D *hptTrigJP1 = (TH1D *)fhist->Get("hptTrigJP1");
	TH1D *hptTrigJP2 = (TH1D *)fhist->Get("hptTrigJP2");

	// mass vs pt correlation.
	TH1D *hprofJP0 = (TH1D *)fhist->Get("hprofMvspTJP0");
	TH1D *hprofJP1 = (TH1D *)fhist->Get("hprofMvspTJP1");
	TH1D *hprofJP2 = (TH1D *)fhist->Get("hprofMvspTJP2");

	// trigger efficiencies
	TH1D *heff_jp0 = (TH1D *)hMTrigJP0->Clone();
	heff_jp0->Divide(hMTrigAll);
	TH1D *heff_jp1 = (TH1D *)hMTrigJP1->Clone();
	heff_jp1->Divide(hMTrigAll);
	TH1D *heff_jp2 = (TH1D *)hMTrigJP2->Clone();
	heff_jp2->Divide(hMTrigAll);

	TH1D *hpteff_jp0 = (TH1D *)hptTrigJP0->Clone();
	hpteff_jp0->Divide(hptTrigAll);
	TH1D *hpteff_jp1 = (TH1D *)hptTrigJP1->Clone();
	hpteff_jp1->Divide(hptTrigAll);
	TH1D *hpteff_jp2 = (TH1D *)hptTrigJP2->Clone();
	hpteff_jp2->Divide(hptTrigAll);

	Double_t eff_trg_jp0[13] = {0};
	Double_t efferr_trg_jp0[13] = {0};
	Double_t eff_trg_jp1[13] = {0};
	Double_t efferr_trg_jp1[13] = {0};
	Double_t eff_trg_jp2[13] = {0};
	Double_t efferr_trg_jp2[13] = {0};

	Double_t pteff_trg_jp0[9] = {0};
	Double_t ptefferr_trg_jp0[9] = {0};
	Double_t pteff_trg_jp1[9] = {0};
	Double_t ptefferr_trg_jp1[9] = {0};
	Double_t pteff_trg_jp2[9] = {0};
	Double_t ptefferr_trg_jp2[9] = {0};

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
	drawEff(hptTrigJP0, hptTrigJP1, hptTrigJP2, hptTrigAll, "p_{T, true}^{#pi^{+}#pi^{-}}");

	// draw mass vs pt correlation
	TCanvas *c_corr = new TCanvas("c_corr", "", 600, 500);
	c_corr->cd();
	c_corr->SetLeftMargin(0.15);
	c_corr->SetBottomMargin(0.15);
	hprofJP2->SetLineColor(6);
	hprofJP2->SetLineWidth(2);
	hprofJP2->SetMarkerStyle(20);
	hprofJP2->SetMarkerColor(6);
	hprofJP2->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}} GeV/c");
	hprofJP2->GetYaxis()->SetTitle("<M_{inv}^{#pi^{+}#pi^{-}}> GeV/c^{2} ");
	hprofJP2->GetYaxis()->SetTitleOffset(1.2);
	hprofJP2->GetYaxis()->SetRangeUser(0.5, 1.2);
	hprofJP2->Draw("hist E");

	hprofJP1->SetLineColor(4);
	hprofJP1->SetLineWidth(2);
	hprofJP1->SetMarkerStyle(20);
	hprofJP1->SetMarkerColor(4);
	hprofJP1->Draw("hist E SAME");

	hprofJP0->SetLineColor(2);
	hprofJP0->SetLineWidth(2);
	hprofJP0->SetMarkerStyle(20);
	hprofJP0->SetMarkerColor(2);
	hprofJP0->Draw("hist E SAME");

	TLegend *cleg = new TLegend(0.6, 0.2, 0.8, 0.4);
	cleg->AddEntry(hprofJP0, "JP0", " lp");
	cleg->AddEntry(hprofJP1, "JP1", " lp");
	cleg->AddEntry(hprofJP2, "JP2", " lp");
	cleg->Draw();
	c_corr->Update();
	c_corr->SaveAs("tprofile_MvspT.pdf");
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
	can->SaveAs(Form("trigEff_%s.pdf", hden->GetName()));
}

void setBinError(TH1D *hnum, TH1D *hden, TH1D *heff)
{

	Double_t binError = 0;
	Double_t binCntNum = 0;
	Double_t binCntDen = 0;
	for (int i = 1; i <= hnum->GetNbinsX(); i++)
	{
		binCntNum = hnum->GetBinContent(i);
		binCntDen = hden->GetBinContent(i);
		if (binCntNum == 0 || binCntDen == 0)
			continue;
		binError = (binCntNum / binCntDen) * sqrt(pow(hnum->GetBinError(i) / binCntNum, 2) + pow(hden->GetBinError(i) / binCntDen, 2));
		heff->SetBinError(i, binError);
	}
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

double getStatError_Binary(double n_num, double n_den)
{
	double binaryError = 0;
	binaryError = n_num / n_den * 1 / n_den * sqrt(n_num * (1 - n_num / n_den)); // binary error propagation
	return binaryError;
}