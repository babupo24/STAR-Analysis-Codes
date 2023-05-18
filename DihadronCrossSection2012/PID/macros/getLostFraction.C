#include <iostream>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iomanip>

using namespace std;

void getLostFraction()
{

    TFile *fin = new TFile("hist4CombBkg.root", "R");
    TH1D *hrecJP0 = (TH1D *)fin->Get("hrecLostJP0");
    TH1D *htrueJP0 = (TH1D *)fin->Get("hrectrueLostJP0");
    TH1D *hrecJP1 = (TH1D *)fin->Get("hrecLostJP1");
    TH1D *htrueJP1 = (TH1D *)fin->Get("hrectrueLostJP1");
    TH1D *hrecJP2 = (TH1D *)fin->Get("hrecLostJP2");
    TH1D *htrueJP2 = (TH1D *)fin->Get("hrectrueLostJP2");

    TH1D *hrJP0 = (TH1D *)htrueJP0->Clone();
    hrJP0->Divide(hrecJP0);
    hrJP0->SetLineColor(2);

    TH1D *hrJP1 = (TH1D *)htrueJP1->Clone();
    hrJP1->Divide(hrecJP1);
    hrJP1->SetLineColor(3);

    TH1D *hrJP2 = (TH1D *)htrueJP2->Clone();
    hrJP2->Divide(hrecJP2);
    hrJP2->SetLineColor(4);

    hrJP0->Draw();
    hrJP1->Draw("same");
    hrJP2->Draw("same");

    cout << "double fraction_missedJP0[13]={";
    for (int i = 1; i <= hrJP0->GetNbinsX(); i++)
    {
        if (i == 13)
            cout << hrJP0->GetBinContent(i) << "};" << endl;
        else
            cout << hrJP0->GetBinContent(i) << ",";
    }

    cout << "double fraction_missedJP1[13]={";
    for (int i = 1; i <= hrJP1->GetNbinsX(); i++)
    {
        if (i == 13)
            cout << hrJP1->GetBinContent(i) << "};" << endl;
        else
            cout << hrJP1->GetBinContent(i) << ",";
    }

    cout << "double fraction_missedJP2[13]={";
    for (int i = 1; i <= hrJP2->GetNbinsX(); i++)
    {
        if (i == 13)
            cout << hrJP2->GetBinContent(i) << "};" << endl;
        else
            cout << hrJP2->GetBinContent(i) << ",";
    }
}