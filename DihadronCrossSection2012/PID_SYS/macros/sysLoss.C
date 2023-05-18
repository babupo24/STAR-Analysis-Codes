#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

void sysLoss()
{
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);

    TFile *fin = new TFile("hist4CombBkg.root", "R");
    TH1D *htrue0 = (TH1D *)fin->Get("htruePiPairJP0");
    TH1D *hsys0 = (TH1D *)fin->Get("htrueDetSysCutJP0");
    hsys0->Divide(htrue0);

    TH1D *htrue1 = (TH1D *)fin->Get("htruePiPairJP1");
    TH1D *hsys1 = (TH1D *)fin->Get("htrueDetSysCutJP1");
    hsys1->Divide(htrue1);

    TH1D *htrue2 = (TH1D *)fin->Get("htruePiPairJP2");
    TH1D *hsys2 = (TH1D *)fin->Get("htrueDetSysCutJP2");
    hsys2->Divide(htrue2);

    TCanvas *can = new TCanvas("can", "", 700, 450);
    can->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetGrid(0, 0);
    hsys0->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hsys0->GetYaxis()->SetTitle("#delta_{loss} = #frac{N_{true && |n#sigma#pi|>4}}{N_{true}}");
    hsys0->GetYaxis()->CenterTitle();
    hsys0->GetYaxis()->SetTitleOffset(1.3);
    hsys0->SetTitle("");
    hsys0->SetLineColor(1);
    hsys0->GetYaxis()->SetRangeUser(0, 0.15);
    hsys0->Draw("E");
    hsys1->SetLineColor(2);
    hsys1->Draw("E SAME");
    hsys2->SetLineColor(3);
    hsys2->Draw("E SAME");

    TLegend *leg = new TLegend(0.3, 0.2, 0.6, 0.26);
    leg->SetNColumns(3);
    leg->AddEntry(hsys0, " jp0", "le");
    leg->AddEntry(hsys1, " jp1", "le");
    leg->AddEntry(hsys2, " jp2", "le");
    leg->Draw();

    can->Update();
    can->SaveAs("Plots/sys_loss.pdf");
    can->SaveAs("Plots/sys_loss.png");

    TF1 *fit0 = new TF1("fit0", "pol0", 0.27, 4);
    fit0->SetLineColor(1);
    fit0->SetLineStyle(2);
    hsys0->Fit("fit0", "R");
    TF1 *fit1 = new TF1("fit1", "pol0", 0.27, 4);
    fit1->SetLineColor(2);
    fit1->SetLineStyle(2);
    hsys1->Fit("fit1", "R");
    TF1 *fit2 = new TF1("fit2", "pol0", 0.27, 4);
    fit2->SetLineColor(3);
    fit2->SetLineStyle(2);
    hsys2->Fit("fit2", "R");
    can->Update();

    cout << "double sys_loss_jp0[13]={";
    for (int i = 1; i <= hsys0->GetNbinsX(); i++)
    {
        if (i == 13)
            cout << hsys0->GetBinContent(i) << "};" << endl;
        else
            cout << hsys0->GetBinContent(i) << ",";
    }
    cout << "double sys_loss_jp1[13]={";
    for (int i = 1; i <= hsys1->GetNbinsX(); i++)
    {
        if (i == 13)
            cout << hsys1->GetBinContent(i) << "};" << endl;
        else
            cout << hsys1->GetBinContent(i) << ",";
    }
    cout << "double sys_loss_jp2[13]={";
    for (int i = 1; i <= hsys2->GetNbinsX(); i++)
    {
        if (i == 13)
            cout << hsys2->GetBinContent(i) << "};" << endl;
        else
            cout << hsys2->GetBinContent(i) << ",";
    }
}