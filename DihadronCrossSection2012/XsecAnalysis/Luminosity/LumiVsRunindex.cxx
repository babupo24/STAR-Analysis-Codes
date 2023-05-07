#include <iostream>
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TCanvas.h"

using namespace std;

void LumiVsRunindex()
{
    ifstream fJP0;
    fJP0.open("trigJP0.txt");
    ifstream fJP1;
    fJP1.open("trigJP1.txt");
    ifstream fJP2;
    fJP2.open("trigJP2.txt");

    if (!fJP0 || !fJP1 || !fJP2)
    {
        cout << "Some or all files couldn't found!" << endl;
        break;
    }
    /* //check values in the file
    int index, run;
    double lum;
    while (fJP1 >> index >> run >> lum)
    {
        cout << index << "  " << lum << endl;
    }
    */

    TGraph *gJP0 = new TGraph("trigJP0.txt", "%lg %*lg %lg", "");
    TGraph *gJP1 = new TGraph("trigJP1.txt", "%lg %*lg %lg", "");
    TGraph *gJP2 = new TGraph("trigJP2.txt", "%lg %*lg %lg", "");

    TCanvas *can = new TCanvas("can", "", 900, 600);
    can->Divide(1, 3);
    can->cd(1);
    gPad->SetBottomMargin(0.005);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.0);
    gPad->SetGrid(5, 5);

    gJP2->SetLineColor(2);
    gJP2->GetYaxis()->SetRangeUser(0.0, 0.10);
    gJP2->GetYaxis()->SetLabelSize(0.09);
    gJP2->GetYaxis()->SetLabelOffset(0.001);
    gJP2->GetYaxis()->SetNdivisions(510);
    gJP2->GetYaxis()->SetTitle("JP2 Lum. (pb^{-1})");
    gJP2->GetYaxis()->SetTitleSize(0.09);
    gJP2->GetYaxis()->CenterTitle();
    gJP2->GetYaxis()->SetTitleOffset(0.5);

    gJP2->GetXaxis()->SetLimits(0.0, 602);
    gJP2->GetXaxis()->SetLabelSize(0.0);
    gJP2->SetTitle("");
    gJP2->Draw();

    can->cd(2);
    gPad->SetBottomMargin(0.005);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.0);
    gPad->SetGrid(5, 5);

    gJP1->SetLineColor(2);
    gJP1->GetYaxis()->SetRangeUser(0.0, 0.035);
    gJP1->GetYaxis()->SetLabelSize(0.09);
    gJP1->GetYaxis()->SetLabelOffset(0.001);
    gJP1->GetYaxis()->SetNdivisions(510);
    gJP1->GetYaxis()->SetTitle("JP1 Lum. (pb^{-1})");
    gJP1->GetYaxis()->SetTitleSize(0.09);
    gJP1->GetYaxis()->CenterTitle();
    gJP1->GetYaxis()->SetTitleOffset(0.5);

    gJP1->GetXaxis()->SetLimits(0.0, 602);
    gJP1->GetXaxis()->SetLabelSize(0.0);
    gJP1->SetTitle("");
    gJP1->Draw();

    can->cd(3);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.0);
    gPad->SetGrid(5, 5);
    gJP0->SetLineColor(2);

    gJP0->GetYaxis()->SetRangeUser(0.0, 0.001);
    gJP0->GetYaxis()->SetNoExponent(kFALSE);
    gJP0->GetYaxis()->SetLabelSize(0.09);
    gJP0->GetYaxis()->SetLabelOffset(0.001);
    gJP0->GetYaxis()->SetNdivisions(510);
    gJP0->GetYaxis()->SetTitle("JP0 Lum. (pb^{-1})");
    gJP0->GetYaxis()->SetTitleSize(0.09);
    gJP0->GetYaxis()->CenterTitle();
    gJP0->GetYaxis()->SetTitleOffset(0.5);

    gJP0->GetXaxis()->SetLimits(0.0, 602);
    gJP0->GetXaxis()->SetLabelSize(0.09);
    gJP0->GetXaxis()->SetLabelOffset(0.001);
    gJP0->GetXaxis()->SetTitle("Run Index");
    gJP0->GetXaxis()->CenterTitle();
    gJP0->GetXaxis()->SetTitleSize(0.09);
    gJP0->GetXaxis()->SetTitleOffset(0.85);
    gJP0->SetTitle("");
    gJP0->Draw();

    can->Update();
    can->SaveAs("JP0-1-2_lum_per_run.pdf");
}