#include <iostream>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iomanip>
#include <cmath>

using namespace std;

void drawSys(TCanvas *can, TH1D *hnom, TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, const char *trigname);

void printNum(ofstream &fname, TH1D *hist, const char *type);

void getTruePionFraction()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile *fintrue = new TFile("hist4CombBkg.root", "r");
    TFile *findet = new TFile("../PID/hist4CombBkg.root", "r");
    ofstream fout;
    fout.open("Plots/fractions_CombBkgCor.txt");

    // hist for mis-identified pion loss
    TH1D *htruePiJP0 = (TH1D *)fintrue->Get("htruePiPairJP0");       // true pion pair
    TH1D *htrueDetPiJP0 = (TH1D *)fintrue->Get("htrueDetPiPairJP0"); // true && det pion pair
    TH1D *htruePiJP1 = (TH1D *)fintrue->Get("htruePiPairJP1");
    TH1D *htrueDetPiJP1 = (TH1D *)fintrue->Get("htrueDetPiPairJP1"); // true && det pion pair
    TH1D *htruePiJP2 = (TH1D *)fintrue->Get("htruePiPairJP2");
    TH1D *htrueDetPiJP2 = (TH1D *)fintrue->Get("htrueDetPiPairJP2"); // true && det pion pair

    // get missed fraction
    TH1D *htrueMissed0 = (TH1D *)htruePiJP0->Clone();
    htrueMissed0->Divide(htrueDetPiJP0);
    TH1D *htrueMissed1 = (TH1D *)htruePiJP1->Clone();
    htrueMissed1->Divide(htrueDetPiJP1);
    TH1D *htrueMissed2 = (TH1D *)htruePiJP2->Clone();
    htrueMissed2->Divide(htrueDetPiJP2);

    // hist for systematic studies
    TH1D *htrueJP2 = (TH1D *)fintrue->Get("htruePiPairJP2"); // true pion pair
    TH1D *hdetJP2 = (TH1D *)findet->Get("hrecPionPairJP2");  // det pion pair
    TH1D *hdetJP2_0 = (TH1D *)findet->Get("hrecSysJP2_0");   // det pion pair
    TH1D *hdetJP2_1 = (TH1D *)findet->Get("hrecSysJP2_1");
    TH1D *hdetJP2_2 = (TH1D *)findet->Get("hrecSysJP2_2");
    TH1D *hdetJP2_3 = (TH1D *)findet->Get("hrecSysJP2_3");

    // get 2D nsigma histograms
    TH2D *hnsigma[13];
    for (int i = 0; i < 13; i++)
    {
        hnsigma[i] = (TH2D *)fintrue->Get(Form("hnsigmaPi_bin%i", i));
    }

    // true fractions
    // missed fractions
    fout << "//Missed fractions-----------------" << endl;
    printNum(fout, htrueMissed2, "missedJP2");
    printNum(fout, htrueMissed1, "missedJP1");
    printNum(fout, htrueMissed0, "missedJP0");

    // draw missed fractions
    TCanvas *c2 = new TCanvas("c2", "", 600, 450);
    c2->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    htrueMissed2->GetXaxis()->SetTitle("M_{inv, true} (GeV/c^{2})");
    htrueMissed2->GetYaxis()->SetTitle("True pair Missed");
    htrueMissed2->SetMarkerStyle(8);
    htrueMissed2->SetMarkerColor(8);
    htrueMissed2->SetLineColor(8);
    htrueMissed2->GetYaxis()->SetRangeUser(1.0, 2.0);
    htrueMissed2->Draw("E");

    htrueMissed1->SetMarkerStyle(8);
    htrueMissed1->SetMarkerColor(2);
    htrueMissed1->SetLineColor(2);
    htrueMissed1->Draw(" SAME E");

    htrueMissed0->SetMarkerStyle(8);
    htrueMissed0->SetMarkerColor(4);
    htrueMissed0->SetLineColor(4);
    htrueMissed0->Draw(" SAME E");

    TLegend *leg1 = new TLegend(0.3, 0.2, 0.7, 0.3);
    leg1->SetNColumns(3);
    leg1->AddEntry(htrueMissed2, " JP2", "lep");
    leg1->AddEntry(htrueMissed1, " JP1", "lep");
    leg1->AddEntry(htrueMissed0, " JP0", "lep");
    leg1->Draw();

    gPad->Update();
    c2->SaveAs("Plots/truePairMissed.pdf");

    // fit last 4 bins
    TF1 *fit2 = new TF1("fit2", "pol0", 1.6, 4);
    fit2->SetLineColor(8);
    fit2->SetLineStyle(2);
    htrueMissed2->Fit("fit2", "R");

    TF1 *fit1 = new TF1("fit1", "pol0", 1.6, 4);
    fit1->SetLineColor(2);
    fit1->SetLineStyle(2);
    htrueMissed1->Fit("fit1", "R");

    TF1 *fit0 = new TF1("fit0", "pol0", 1.6, 4);
    fit0->SetLineColor(4);
    fit0->SetLineStyle(2);
    htrueMissed0->Fit("fit0", "R");

    gPad->Update();
    c2->SaveAs("Plots/truePairMissed_fitted.pdf");
    /*
        // systematic studies for the combinatorial background

        TH1D *hdefault2 = (TH1D *)htrueJP2->Clone();
        hdefault2->Divide(hdetJP2);

        TH1D *hdefault2 = (TH1D *)htrueJP2->Clone();
        hdefault2->Divide(hdetJP2);

        // TH1D *hrecSys2_0 = (TH1D *)fintrue->Get("hrecSysJP2_0");
        // TH1D *htrueSys2_0 = (TH1D *)fintrue->Get("htrueSysJP2_0");
        TH1D *hfracSys2_0 = (TH1D *)htrueJP2->Clone();
        hfracSys2_0->Divide(hdetJP2_0);
        TH1D *hfracSys2_1 = (TH1D *)htrueJP2->Clone();
        hfracSys2_1->Divide(hdetJP2_1);
        TH1D *hfracSys2_2 = (TH1D *)htrueJP2->Clone();
        hfracSys2_2->Divide(hdetJP2_2);
        TH1D *hfracSys2_3 = (TH1D *)htrueJP2->Clone();
        hfracSys2_3->Divide(hdetJP2_3);

        TH1D *hrecSys2_1 = (TH1D *)fintrue->Get("hrecSysJP2_1");
        TH1D *htrueSys2_1 = (TH1D *)fintrue->Get("htrueSysJP2_0");
        TH1D *hfracSys2_1 = (TH1D *)htrueSys2_1->Clone();
        hfracSys2_1->Divide(hrecSys2_1);

        TH1D *hrecSys2_2 = (TH1D *)fintrue->Get("hrecSysJP2_2");
        TH1D *htrueSys2_2 = (TH1D *)fintrue->Get("htrueSysJP2_0");
        TH1D *hfracSys2_2 = (TH1D *)htrueSys2_2->Clone();
        hfracSys2_2->Divide(hrecSys2_2);

        TH1D *hrecSys2_3 = (TH1D *)fintrue->Get("hrecSysJP2_3");
        TH1D *htrueSys2_3 = (TH1D *)fintrue->Get("htrueSysJP2_0");
        TH1D *hfracSys2_3 = (TH1D *)htrueSys2_3->Clone();
        hfracSys2_3->Divide(hrecSys2_3);

        TH1D *hrecSys1_0 = (TH1D *)fintrue->Get("hrecSysJP1_0");
        TH1D *htrueSys1_0 = (TH1D *)fintrue->Get("htrueSysJP1_0");
        TH1D *hfracSys1_0 = (TH1D *)htrueSys1_0->Clone();
        hfracSys1_0->Divide(hrecSys1_0);

        TH1D *hrecSys1_1 = (TH1D *)fintrue->Get("hrecSysJP1_1");
        TH1D *htrueSys1_1 = (TH1D *)fintrue->Get("htrueSysJP1_0");
        TH1D *hfracSys1_1 = (TH1D *)htrueSys1_1->Clone();
        hfracSys1_1->Divide(hrecSys1_1);

        TH1D *hrecSys1_2 = (TH1D *)fintrue->Get("hrecSysJP1_2");
        TH1D *htrueSys1_2 = (TH1D *)fintrue->Get("htrueSysJP1_0");
        TH1D *hfracSys1_2 = (TH1D *)htrueSys1_2->Clone();
        hfracSys1_2->Divide(hrecSys1_2);

        TH1D *hrecSys1_3 = (TH1D *)fintrue->Get("hrecSysJP1_3");
        TH1D *htrueSys1_3 = (TH1D *)fintrue->Get("htrueSysJP1_0");
        TH1D *hfracSys1_3 = (TH1D *)htrueSys1_3->Clone();
        hfracSys1_3->Divide(hrecSys1_3);

        TH1D *hrecSys0_0 = (TH1D *)fintrue->Get("hrecSysJP0_0");
        TH1D *htrueSys0_0 = (TH1D *)fintrue->Get("htrueSysJP0_0");
        TH1D *hfracSys0_0 = (TH1D *)htrueSys0_0->Clone();
        hfracSys0_0->Divide(hrecSys0_0);

        TH1D *hrecSys0_1 = (TH1D *)fintrue->Get("hrecSysJP0_1");
        TH1D *htrueSys0_1 = (TH1D *)fintrue->Get("htrueSysJP0_0");
        TH1D *hfracSys0_1 = (TH1D *)htrueSys0_1->Clone();
        hfracSys0_1->Divide(hrecSys0_1);

        TH1D *hrecSys0_2 = (TH1D *)fintrue->Get("hrecSysJP0_2");
        TH1D *htrueSys0_2 = (TH1D *)fintrue->Get("htrueSysJP0_0");
        TH1D *hfracSys0_2 = (TH1D *)htrueSys0_2->Clone();
        hfracSys0_2->Divide(hrecSys0_2);

        TH1D *hrecSys0_3 = (TH1D *)fintrue->Get("hrecSysJP0_3");
        TH1D *htrueSys0_3 = (TH1D *)fintrue->Get("htrueSysJP0_0");
        TH1D *hfracSys0_3 = (TH1D *)htrueSys0_3->Clone();
        hfracSys0_3->Divide(hrecSys0_3);

        TCanvas *can0 = new TCanvas("can0", "", 900, 450);
        drawSys(can0, htrueMissed0, hfracSys0_0, hfracSys0_1, hfracSys0_2, hfracSys0_3, "trig::JP0");
        can0->SaveAs("Plots/pidsys_jp0.pdf");
        TCanvas *can1 = new TCanvas("can1", "", 900, 450);
        drawSys(can1, htrueMissed1, hfracSys1_0, hfracSys1_1, hfracSys1_2, hfracSys1_3, "trig::JP1");
        can1->SaveAs("Plots/pidsys_jp1.pdf");

        // TCanvas *can2 = new TCanvas("can2", "", 900, 450);
        // drawSys(can2, hdefault2, hfracSys2_0, hfracSys2_1, hfracSys2_2, hfracSys2_3, "trig::JP2");
        // can2->SaveAs("Plots/pidsys_jp2.pdf");
    */
}

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
    chnom->GetYaxis()->SetRangeUser(0.5, 1.6);
    chnom->Draw("E");
    chnom->Fit("pol0");

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

void printNum(ofstream &fname, TH1D *hist, const char *type)
{
    fname << Form("double fraction_%s[13]={", type);
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double nf = hist->GetBinContent(i);
        if (i == hist->GetNbinsX())
            fname << nf << "};" << endl;
        else
            fname << nf << ", ";
    }
}
