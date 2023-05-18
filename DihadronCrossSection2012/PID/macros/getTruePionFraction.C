#include <iostream>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iomanip>
#include <cmath>

using namespace std;

void drawSys(TCanvas *can, TH1D *hnom, TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, const char *trigname);

void getTruePionFraction()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile *fin = new TFile("hist4CombBkg.root", "r");
    // TFile *fin = new TFile("hist.root", "r");
    assert(fin);
    // hists for true fraciton in with the default dE/dx cut
    TH1D *htrue2 = (TH1D *)fin->Get("htruePionPairJP2");
    TH1D *hrec2 = (TH1D *)fin->Get("hrecPionPairJP2");
    TH1D *htrue1 = (TH1D *)fin->Get("htruePionPairJP1");
    TH1D *hrec1 = (TH1D *)fin->Get("hrecPionPairJP1");
    TH1D *htrue0 = (TH1D *)fin->Get("htruePionPairJP0");
    TH1D *hrec0 = (TH1D *)fin->Get("hrecPionPairJP0");

    TH1D *htrueFraction2 = (TH1D *)htrue2->Clone();
    htrueFraction2->Divide(hrec2);
    TH1D *htrueFraction1 = (TH1D *)htrue1->Clone();
    htrueFraction1->Divide(hrec1);
    TH1D *htrueFraction0 = (TH1D *)htrue0->Clone();
    htrueFraction0->Divide(hrec0);

    TCanvas *c1 = new TCanvas("c1", "", 600, 450);
    c1->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    htrueFraction2->GetXaxis()->SetTitle("M_{inv, true} (GeV/c^{2})");
    htrueFraction2->GetYaxis()->SetTitle("True pair fraction");
    htrueFraction2->SetMarkerStyle(8);
    htrueFraction2->SetMarkerColor(8);
    htrueFraction2->SetLineColor(8);
    htrueFraction2->GetYaxis()->SetRangeUser(0.5, 1.05);
    htrueFraction2->Draw("E");

    htrueFraction1->SetMarkerStyle(8);
    htrueFraction1->SetMarkerColor(2);
    htrueFraction1->SetLineColor(2);
    htrueFraction1->Draw(" SAME E");

    htrueFraction0->SetMarkerStyle(24);
    htrueFraction0->SetMarkerColor(4);
    htrueFraction0->SetLineColor(4);
    htrueFraction0->Draw(" SAME E");

    TLegend *leg1 = new TLegend(0.3, 0.2, 0.7, 0.3);
    leg1->SetNColumns(3);
    leg1->AddEntry(htrueFraction2, " JP2", "lep");
    leg1->AddEntry(htrueFraction1, " JP1", "lep");
    leg1->AddEntry(htrueFraction0, " JP0", "lep");
    leg1->Draw();

    gPad->Update();
    c1->SaveAs("Plots/truePairFraction.pdf");

    // fit last 4 bins
    TF1 *fit2 = new TF1("fit2", "pol2", 1.6, 4);
    fit2->SetLineColor(8);
    fit2->SetLineStyle(2);
    htrueFraction2->Fit("fit2", "R");

    TF1 *fit1 = new TF1("fit1", "pol2", 1.6, 4);
    fit1->SetLineColor(2);
    fit1->SetLineStyle(2);
    htrueFraction1->Fit("fit1", "R");

    TF1 *fit0 = new TF1("fit0", "pol2", 1.6, 4);
    fit0->SetLineColor(4);
    fit0->SetLineStyle(2);
    htrueFraction0->Fit("fit0", "R");
    cout << "JP0 functional form " << fit0->Eval(htrueFraction0->GetBinCenter(11)) << "  " << fit0->Eval(htrueFraction0->GetBinCenter(12)) << "  " << fit0->Eval(htrueFraction0->GetBinCenter(13)) << endl;
    cout << "JP1 functional form " << fit1->Eval(htrueFraction1->GetBinCenter(11)) << "  " << fit1->Eval(htrueFraction1->GetBinCenter(12)) << "  " << fit1->Eval(htrueFraction1->GetBinCenter(13)) << endl;
    cout << "JP2 functional form " << fit2->Eval(htrueFraction2->GetBinCenter(11)) << "  " << fit2->Eval(htrueFraction2->GetBinCenter(12)) << "  " << fit2->Eval(htrueFraction2->GetBinCenter(13)) << endl;

    gPad->Update();
    c1->SaveAs("Plots/truePairFraction_fitted.pdf");

    cout << "double truePairFractionJP2[13]={";
    for (int i = 1; i <= htrueFraction2->GetNbinsX(); i++)
    {
        double nf = htrueFraction2->GetBinContent(i);
        if (nf == 0)
            continue;
        if (i == htrueFraction2->GetNbinsX())
            cout << nf << "};" << endl;
        else
            cout << nf << ",";
    }
    cout << "double truePairFractionJP1[13]={";
    for (int i = 1; i <= htrueFraction1->GetNbinsX(); i++)
    {
        double nf = htrueFraction1->GetBinContent(i);
        if (nf == 0)
            continue;
        if (i == htrueFraction1->GetNbinsX())
            cout << nf << "};" << endl;
        else
            cout << nf << ",";
    }
    cout << "double truePairFractionJP0[13]={";
    for (int i = 1; i <= htrueFraction0->GetNbinsX(); i++)
    {
        double nf = htrueFraction0->GetBinContent(i);
        if (nf == 0)
            continue;
        if (i == htrueFraction0->GetNbinsX())
            cout << nf << "};" << endl;
        else
            cout << nf << ",";
    }

    // systematic studies for the combinatorial background
    TH1D *hrecSys2_0 = (TH1D *)fin->Get("hrecSysJP2_0");
    TH1D *htrueSys2_0 = (TH1D *)fin->Get("htrueSysJP2_0");
    TH1D *hfracSys2_0 = (TH1D *)htrueSys2_0->Clone();
    hfracSys2_0->Divide(hrecSys2_0);

    TH1D *hrecSys2_1 = (TH1D *)fin->Get("hrecSysJP2_1");
    TH1D *htrueSys2_1 = (TH1D *)fin->Get("htrueSysJP2_1");
    TH1D *hfracSys2_1 = (TH1D *)htrueSys2_1->Clone();
    hfracSys2_1->Divide(hrecSys2_1);

    TH1D *hrecSys2_2 = (TH1D *)fin->Get("hrecSysJP2_2");
    TH1D *htrueSys2_2 = (TH1D *)fin->Get("htrueSysJP2_2");
    TH1D *hfracSys2_2 = (TH1D *)htrueSys2_2->Clone();
    hfracSys2_2->Divide(hrecSys2_2);

    TH1D *hrecSys2_3 = (TH1D *)fin->Get("hrecSysJP2_3");
    TH1D *htrueSys2_3 = (TH1D *)fin->Get("htrueSysJP2_3");
    TH1D *hfracSys2_3 = (TH1D *)htrueSys2_3->Clone();
    hfracSys2_3->Divide(hrecSys2_3);

    TH1D *hrecSys1_0 = (TH1D *)fin->Get("hrecSysJP1_0");
    TH1D *htrueSys1_0 = (TH1D *)fin->Get("htrueSysJP1_0");
    TH1D *hfracSys1_0 = (TH1D *)htrueSys1_0->Clone();
    hfracSys1_0->Divide(hrecSys1_0);

    TH1D *hrecSys1_1 = (TH1D *)fin->Get("hrecSysJP1_1");
    TH1D *htrueSys1_1 = (TH1D *)fin->Get("htrueSysJP1_1");
    TH1D *hfracSys1_1 = (TH1D *)htrueSys1_1->Clone();
    hfracSys1_1->Divide(hrecSys1_1);

    TH1D *hrecSys1_2 = (TH1D *)fin->Get("hrecSysJP1_2");
    TH1D *htrueSys1_2 = (TH1D *)fin->Get("htrueSysJP1_2");
    TH1D *hfracSys1_2 = (TH1D *)htrueSys1_2->Clone();
    hfracSys1_2->Divide(hrecSys1_2);

    TH1D *hrecSys1_3 = (TH1D *)fin->Get("hrecSysJP1_3");
    TH1D *htrueSys1_3 = (TH1D *)fin->Get("htrueSysJP1_3");
    TH1D *hfracSys1_3 = (TH1D *)htrueSys1_3->Clone();
    hfracSys1_3->Divide(hrecSys1_3);

    TH1D *hrecSys0_0 = (TH1D *)fin->Get("hrecSysJP0_0");
    TH1D *htrueSys0_0 = (TH1D *)fin->Get("htrueSysJP0_0");
    TH1D *hfracSys0_0 = (TH1D *)htrueSys0_0->Clone();
    hfracSys0_0->Divide(hrecSys0_0);

    TH1D *hrecSys0_1 = (TH1D *)fin->Get("hrecSysJP0_1");
    TH1D *htrueSys0_1 = (TH1D *)fin->Get("htrueSysJP0_1");
    TH1D *hfracSys0_1 = (TH1D *)htrueSys0_1->Clone();
    hfracSys0_1->Divide(hrecSys0_1);

    TH1D *hrecSys0_2 = (TH1D *)fin->Get("hrecSysJP0_2");
    TH1D *htrueSys0_2 = (TH1D *)fin->Get("htrueSysJP0_2");
    TH1D *hfracSys0_2 = (TH1D *)htrueSys0_2->Clone();
    hfracSys0_2->Divide(hrecSys0_2);

    TH1D *hrecSys0_3 = (TH1D *)fin->Get("hrecSysJP0_3");
    TH1D *htrueSys0_3 = (TH1D *)fin->Get("htrueSysJP0_3");
    TH1D *hfracSys0_3 = (TH1D *)htrueSys0_3->Clone();
    hfracSys0_3->Divide(hrecSys0_3);

    // sys values
    double pidSys_jp0[13] = {0};
    double pidSys_jp1[13] = {0};
    double pidSys_jp2[13] = {0};
    for (int i = 1; i <= hfracSys0_3->GetNbinsX(); i++)
    {
        double n3_0 = hfracSys0_3->GetBinContent(i);
        double n2_0 = hfracSys0_2->GetBinContent(i);
        double sys_0 = abs(n3_0 - n2_0);
        pidSys_jp0[i - 1] = sys_0;

        double n3_1 = hfracSys1_3->GetBinContent(i);
        double n2_1 = hfracSys1_2->GetBinContent(i);
        double sys_1 = abs(n3_1 - n2_1);
        pidSys_jp1[i - 1] = sys_1;

        double n3_2 = hfracSys2_3->GetBinContent(i);
        double n2_2 = hfracSys2_2->GetBinContent(i);
        double sys_2 = abs(n3_2 - n2_2);
        pidSys_jp2[i - 1] = sys_2;
    }
    cout << "double pid_sys_jp0[13]={";
    for (int i = 0; i < 13; i++)
    {
        if (i == 12)
            cout << pidSys_jp0[i] << "};" << endl;
        else
            cout << pidSys_jp0[i] << ", ";
    }
    cout << "double pid_sys_jp1[13]={";
    for (int i = 0; i < 13; i++)
    {
        if (i == 12)
            cout << pidSys_jp1[i] << "};" << endl;
        else
            cout << pidSys_jp1[i] << ", ";
    }
    cout << "double pid_sys_jp2[13]={";
    for (int i = 0; i < 13; i++)
    {
        if (i == 12)
            cout << pidSys_jp2[i] << "};" << endl;
        else
            cout << pidSys_jp2[i] << ", ";
    }

    /*
        TCanvas *can0 = new TCanvas("can0", "", 900, 450);
        drawSys(can0, htrueFraction0, hfracSys0_0, hfracSys0_1, hfracSys0_2, hfracSys0_3, "trig::JP0");
        can0->SaveAs("Plots/pidsys_jp0.pdf");
        TCanvas *can1 = new TCanvas("can1", "", 900, 450);
        drawSys(can1, htrueFraction1, hfracSys1_0, hfracSys1_1, hfracSys1_2, hfracSys1_3, "trig::JP1");
        can1->SaveAs("Plots/pidsys_jp1.pdf");
        TCanvas *can2 = new TCanvas("can2", "", 900, 450);
        drawSys(can2, htrueFraction2, hfracSys2_0, hfracSys2_1, hfracSys2_2, hfracSys2_3, "trig::JP2");
        can2->SaveAs("Plots/pidsys_jp2.pdf");
        */
} // main

void drawSys(TCanvas *can, TH1D *hnom, TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, const char *trigname)
{
    // can = new TCanvas("can", "", 900, 450);
    can->Divide(2, 1);
    can->cd(1);
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetGrid(0, 0);

    TH1D *chnom = (TH1D *)hnom->Clone();
    chnom->GetXaxis()->SetTitle("M_{inv, true} (GeV/c^{2})");
    chnom->GetYaxis()->SetTitle("True pair fraction");
    chnom->SetTitle(trigname);
    chnom->SetMarkerStyle(8);
    chnom->SetMarkerColor(8);
    chnom->SetLineColor(8);
    chnom->GetYaxis()->SetRangeUser(0.5, 1.05);
    chnom->Draw("E");

    TH1D *ch0 = (TH1D *)h0->Clone();
    ch0->SetMarkerStyle(8);
    ch0->SetMarkerColor(2);
    ch0->SetLineColor(2);
    ch0->Draw(" SAME E");

    TH1D *ch1 = (TH1D *)h1->Clone();
    ch1->SetMarkerStyle(8);
    ch1->SetMarkerColor(4);
    ch1->SetLineColor(4);
    ch1->Draw(" SAME E");

    TH1D *ch2 = (TH1D *)h2->Clone();
    ch2->SetMarkerStyle(8);
    ch2->SetMarkerColor(1);
    ch2->SetLineColor(1);
    ch2->Draw(" SAME E");

    TH1D *ch3 = (TH1D *)h3->Clone();
    ch3->SetMarkerStyle(8);
    ch3->SetMarkerColor(6);
    ch3->SetLineColor(6);
    ch3->Draw(" SAME E");

    TLegend *leg1 = new TLegend(0.4, 0.2, 0.7, 0.5);
    // leg1->SetNColumns(3);
    leg1->AddEntry(chnom, " -1<n#sigma#pi<2", "lep");
    leg1->AddEntry(ch0, " -0.5<n#sigma#pi<2", "lep");
    leg1->AddEntry(ch1, " -1<n#sigma#pi<1.5", "lep");
    leg1->AddEntry(ch2, " -0.5<n#sigma#pi<1.5", "lep");
    leg1->AddEntry(ch3, " -1.5<n#sigma#pi<2.5", "lep");
    leg1->SetTextSize(0.035);
    leg1->Draw();
    gPad->Update();
    can->cd(2);
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetGrid(0, 0);

    TH1D *chr0 = (TH1D *)hnom->Clone();
    chr0->Divide(h0);
    chr0->GetXaxis()->SetTitle("M_{inv, true} (GeV/c^{2})");
    chr0->GetYaxis()->SetTitle("Default / X");
    chr0->SetMarkerStyle(8);
    chr0->SetMarkerColor(2);
    chr0->SetLineColor(2);
    chr0->GetYaxis()->SetRangeUser(0.5, 1.5);
    chr0->Draw(" E");

    TH1D *chr1 = (TH1D *)hnom->Clone();
    chr1->Divide(h1);
    chr1->SetMarkerStyle(8);
    chr1->SetMarkerColor(4);
    chr1->SetLineColor(4);
    chr1->Draw("SAME E");

    TH1D *chr2 = (TH1D *)hnom->Clone();
    chr2->Divide(h2);
    chr2->SetMarkerStyle(8);
    chr2->SetMarkerColor(1);
    chr2->SetLineColor(1);
    chr2->Draw("SAME E");
    chr2->Fit("pol0");
    chr2->GetFunction("pol0")->SetLineColor(1);
    chr2->GetFunction("pol0")->SetLineStyle(2);

    TH1D *chr3 = (TH1D *)hnom->Clone();
    chr3->Divide(h3);
    chr3->SetMarkerStyle(8);
    chr3->SetMarkerColor(6);
    chr3->SetLineColor(6);
    chr3->Draw("SAME E");
    chr3->Fit("pol0");
    chr3->GetFunction("pol0")->SetLineColor(6);
    chr3->GetFunction("pol0")->SetLineStyle(2);

    TLatex tex;
    tex.SetTextSize(0.03);
    tex.DrawLatex(1.0, 1.4, Form("#color[6]{chi2/NDF = %3.4f / %i}", chr3->GetFunction("pol0")->GetChisquare(), chr3->GetFunction("pol0")->GetNDF()));
    tex.DrawLatex(1.0, 1.35, Form("#color[6]{p0 #pm #delta p0 = %3.4f #pm %3.4f}", chr3->GetFunction("pol0")->GetParameter(0), chr3->GetFunction("pol0")->GetParError(0)));
    tex.DrawLatex(1.0, 0.6, Form("chi2/NDF = %3.4f / %i", chr2->GetFunction("pol0")->GetChisquare(), chr2->GetFunction("pol0")->GetNDF()));
    tex.DrawLatex(1.0, 0.55, Form("p0 #pm #delta p0 = %3.4f #pm %3.4f", chr2->GetFunction("pol0")->GetParameter(0), chr2->GetFunction("pol0")->GetParError(0)));

    gPad->Update();
}