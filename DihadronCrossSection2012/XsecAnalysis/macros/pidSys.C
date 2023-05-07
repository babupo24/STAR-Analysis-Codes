// This macro plots the cross-section measured with the response constructed with the nSigmaPion cut and both the nSigmaPion selected True pions and compare the final measured cross-sections. The relative difference is calculated between the two cross-sections and assigned as a systematic uncertainty for the f_loss, which gives the fraction of pion loss due to the

#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <iomanip>
#include <TCanvas.h>

using namespace std;

void pidSys()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    TFile *fdedx = new TFile("ResultsMinConeCut/fxsecOutput.root");
    TH1D *h0 = (TH1D *)fdedx->Get("hxsec_jp0");
    TH1D *h1 = (TH1D *)fdedx->Get("hxsec_jp1");
    TH1D *h2 = (TH1D *)fdedx->Get("hxsec_jp2");
    TH1D *hcomb = (TH1D *)fdedx->Get("hxsec_comb");
    TFile *fdedx_true = new TFile("ResultsMinConeCut/fxsecOutputTrue.root");
    TH1D *h0_tr = (TH1D *)fdedx_true->Get("hxsec_jp0");
    TH1D *h1_tr = (TH1D *)fdedx_true->Get("hxsec_jp1");
    TH1D *h2_tr = (TH1D *)fdedx_true->Get("hxsec_jp2");
    TH1D *hcomb_tr = (TH1D *)fdedx_true->Get("hxsec_comb");

    TCanvas *can = new TCanvas("can", "", 800, 450);
    can->Divide(2, 1);
    can->cd(1);
    TH1D *h0c = (TH1D *)h0->Clone();
    h0c->SetLineColor(1);
    h0c->SetLineWidth(2);
    h0c->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    h0c->GetXaxis()->SetLabelSize(0.04);
    h0c->GetYaxis()->SetTitle("d#sigma/dM");
    h0c->GetYaxis()->SetTitleOffset(1.2);
    h0c->SetLineStyle(1);
    h0c->Draw();
    TH1D *h1c = (TH1D *)h1->Clone();
    h1c->SetLineColor(2);
    h1c->SetLineWidth(2);
    h1c->Draw("same");
    TH1D *h2c = (TH1D *)h2->Clone();
    h2c->SetLineColor(3);
    h2c->SetLineWidth(2);
    h2c->Draw("same");

    TH1D *h0c_tr = (TH1D *)h0_tr->Clone();
    h0c_tr->SetLineColor(1);
    h0c_tr->SetLineWidth(1);
    h0c_tr->SetLineStyle(2);
    h0c_tr->Draw("same");

    TH1D *h1c_tr = (TH1D *)h1_tr->Clone();
    h1c_tr->SetLineColor(2);
    h1c_tr->SetLineStyle(2);
    h1c_tr->SetLineWidth(2);
    h1c_tr->Draw("same");

    TH1D *h2c_tr = (TH1D *)h2_tr->Clone();
    h2c_tr->SetLineColor(3);
    h2c_tr->SetLineStyle(2);
    h2c_tr->SetLineWidth(2);
    h2c_tr->Draw("same");
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h0c, " JP0", " le");
    leg->AddEntry(h1c, " JP1", " le");
    leg->AddEntry(h2c, " JP2", " le");
    leg->Draw();

    gPad->SetLogy();
    gPad->Update();

    can->cd(2);
    gPad->SetLeftMargin(0.12);
    TH1D *hratio0 = (TH1D *)h0->Clone();
    hratio0->Add(h0_tr, -1);
    hratio0->Divide(h0_tr);
    hratio0->SetLineColor(1);
    hratio0->SetLineWidth(2);
    // hratio0->SetTitle("Xsec Ratio (#sigma_{n#sigma_{#pi}} / #sigma_{n#sigma_{#pi} && True})");
    hratio0->SetTitle("");
    hratio0->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hratio0->GetXaxis()->SetLabelSize(0.04);
    hratio0->GetYaxis()->SetTitle("(#sigma_{n#sigma_{#pi}} - #sigma_{n#sigma_{#pi} && True})/ #sigma_{n#sigma_{#pi} && True}");
    hratio0->GetYaxis()->SetTitleOffset(1.4);
    hratio0->GetYaxis()->SetRangeUser(-0.21, 0.21);
    hratio0->Draw();
    TH1D *hratio1 = (TH1D *)h1->Clone();
    hratio1->Add(h1_tr, -1);
    hratio1->Divide(h1_tr);
    hratio1->SetLineColor(2);
    hratio1->SetLineWidth(2);
    hratio1->Draw("same");
    TH1D *hratio2 = (TH1D *)h2->Clone();
    hratio2->Add(h2_tr, -1);
    hratio2->Divide(h2_tr);
    hratio2->SetLineColor(3);
    hratio2->SetLineWidth(2);
    hratio2->Draw("same");
    gPad->Update();
    can->SaveAs("ResultsMinConeCut/true_dedx_sys_Trgs.pdf");
    // compare final  xsection
    TCanvas *canc = new TCanvas("canc", "", 800, 450);
    canc->Divide(2, 1);
    canc->cd(1);
    TH1D *hcombc = (TH1D *)hcomb->Clone();
    hcombc->SetMarkerSize(0);
    hcombc->SetLineColor(1);
    hcombc->SetLineStyle(1);
    hcombc->SetLineWidth(2);
    hcombc->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
    hcombc->GetYaxis()->SetTitle("d#sigma/dM");
    hcombc->GetYaxis()->SetTitleOffset(1.2);
    hcombc->SetMinimum(10);
    hcombc->Draw();
    TH1D *hcombc_tr = (TH1D *)hcomb_tr->Clone();
    hcombc_tr->SetLineColor(2);
    hcombc_tr->SetLineStyle(2);
    hcombc_tr->SetLineWidth(2);
    hcombc_tr->Draw("same");
    gPad->SetLogy();
    TLegend *leg2 = new TLegend(0.4, 0.75, 0.9, 0.9);
    leg2->AddEntry(hcombc, " -1<n#sigma#pi<2", " le");
    leg2->AddEntry(hcombc_tr, " -1<n#sigma#pi<2 && True Pion", " le");
    leg2->SetTextSize(0.03);
    leg2->Draw();
    canc->cd(2);
    gPad->SetLeftMargin(0.12);
    TH1D *hratio = (TH1D *)hcomb->Clone();
    hratio->Add(hcomb_tr, -1);
    hratio->Divide(hcomb_tr);
    hratio->SetMarkerSize(0);
    hratio->SetLineColor(1);
    hratio->SetLineWidth(2);
    hratio->SetLineStyle(1);
    // hratio->SetTitle("Xsec Ratio (dEdx / dEdx && True)");
    hratio->GetXaxis()->SetTitle("M_{inv}");
    hratio->GetYaxis()->SetTitle("(#sigma_{n#sigma_{#pi}} - #sigma_{n#sigma_{#pi} && True})/ #sigma_{n#sigma_{#pi} && True}");
    hratio->GetYaxis()->SetTitleOffset(1.4);
    hratio->GetYaxis()->SetRangeUser(-0.21, 0.21);
    hratio->Draw();
    canc->Update();
    canc->SaveAs("ResultsMinConeCut/true_dedx_sys_comb.pdf");

    cout << "double sys_pid_true0[13]={";
    for (int i = 2; i < h0->GetNbinsX(); i++)
    {
        if (i == h0->GetNbinsX() - 1)
            cout << hratio0->GetBinContent(i) << "};" << endl;
        else
            cout << hratio0->GetBinContent(i) << ", ";
    }
    cout << "double sys_pid_true1[13]={";
    for (int i = 2; i < h1->GetNbinsX(); i++)
    {
        if (i == h1->GetNbinsX() - 1)
            cout << hratio1->GetBinContent(i) << "};" << endl;
        else
            cout << hratio1->GetBinContent(i) << ", ";
    }
    cout << "double sys_pid_true2[13]={";
    for (int i = 2; i < h2->GetNbinsX(); i++)
    {
        if (i == h2->GetNbinsX() - 1)
            cout << hratio2->GetBinContent(i) << "};" << endl;
        else
            cout << hratio2->GetBinContent(i) << ", ";
    }
    cout << "double sys_pid_trueComb[13]={";
    for (int i = 2; i < h2->GetNbinsX(); i++)
    {
        if (i == h2->GetNbinsX() - 1)
            cout << hratio->GetBinContent(i) << "};" << endl;
        else
            cout << hratio->GetBinContent(i) << ", ";
    }
}
