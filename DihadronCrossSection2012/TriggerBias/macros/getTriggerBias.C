#include <iostream>
#include <TFile.h>
#include <cmath>
#include <TH1.h>

using namespace std;

void drawBias(TH1D *hDetQf, TH1D *hDetGf, TH1D *hParQf, TH1D *hParGf, TH1D *hQBias, const char *trig);

void getTriggerBias()
{

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    // open parton level root file
    // open particle level root file
    TFile *fbiased = new TFile("GetBiased/histBiased.root");
    TFile *funbiased = new TFile("GetUnbiased/histUnbiased.root");

    // get particle level hists
    TH1D *hQuarkPar = (TH1D *)funbiased->Get("hQuarksAll");
    TH1D *hGluonPar = (TH1D *)funbiased->Get("hGluonsAll");
    TH1D *hPartonPar = (TH1D *)funbiased->Get("hPartonsAll");
    // get detector level hists (JP0 || JP1 || JP2)
    TH1D *hQuarkDet = (TH1D *)fbiased->Get("hQuarksAll");
    TH1D *hGluonDet = (TH1D *)fbiased->Get("hGluonsAll");
    TH1D *hPartonDet = (TH1D *)fbiased->Get("hPartonsAll");
    // get detector level hists (JP0 )
    TH1D *hQuarkDetJP0 = (TH1D *)fbiased->Get("hQuarksJP0");
    TH1D *hGluonDetJP0 = (TH1D *)fbiased->Get("hGluonsJP0");
    TH1D *hPartonDetJP0 = (TH1D *)fbiased->Get("hPartonsJP0");
    // get detector level hists (JP1 )
    TH1D *hQuarkDetJP1 = (TH1D *)fbiased->Get("hQuarksJP1");
    TH1D *hGluonDetJP1 = (TH1D *)fbiased->Get("hGluonsJP1");
    TH1D *hPartonDetJP1 = (TH1D *)fbiased->Get("hPartonsJP1");
    // get detector level hists (JP2 )
    TH1D *hQuarkDetJP2 = (TH1D *)fbiased->Get("hQuarksJP2");
    TH1D *hGluonDetJP2 = (TH1D *)fbiased->Get("hGluonsJP2");
    TH1D *hPartonDetJP2 = (TH1D *)fbiased->Get("hPartonsJP2");

    // quark fraction at particle level
    TH1D *hQ_frac_par = (TH1D *)hQuarkPar->Clone();
    hQ_frac_par->Divide(hPartonPar);
    TH1D *hG_frac_par = (TH1D *)hGluonPar->Clone();
    hG_frac_par->Divide(hPartonPar);

    // quark fraction at detector level for jp triggers
    TH1D *hQ_frac_det = (TH1D *)hQuarkDet->Clone();
    hQ_frac_det->Divide(hPartonDet);
    TH1D *hG_frac_det = (TH1D *)hGluonDet->Clone();
    hG_frac_det->Divide(hPartonDet);
    // quark fraction at detector level for jp0
    TH1D *hQ_frac_detJP0 = (TH1D *)hQuarkDetJP0->Clone();
    hQ_frac_detJP0->Divide(hPartonDetJP0);
    // gluon fraction
    TH1D *hG_frac_detJP0 = (TH1D *)hGluonDetJP0->Clone();
    hG_frac_detJP0->Divide(hPartonDetJP0);
    // quark fraction at detector level for jp1
    TH1D *hQ_frac_detJP1 = (TH1D *)hQuarkDetJP1->Clone();
    hQ_frac_detJP1->Divide(hPartonDetJP1);
    // gluon fraction
    TH1D *hG_frac_detJP1 = (TH1D *)hGluonDetJP1->Clone();
    hG_frac_detJP1->Divide(hPartonDetJP1);
    // quark fraction at detector level for jp2
    TH1D *hQ_frac_detJP2 = (TH1D *)hQuarkDetJP2->Clone();
    hQ_frac_detJP2->Divide(hPartonDetJP2);
    // gluon fraction
    TH1D *hG_frac_detJP2 = (TH1D *)hGluonDetJP2->Clone();
    hG_frac_detJP2->Divide(hPartonDetJP2);

    // trigger bias jp triggers
    TH1D *hQBias = (TH1D *)hQ_frac_det->Clone();
    hQBias->Divide(hQ_frac_par);
    // trigger bias JP0
    TH1D *hQBiasJP0 = (TH1D *)hQ_frac_detJP0->Clone();
    hQBiasJP0->Divide(hQ_frac_par);
    // trigger bias JP1
    TH1D *hQBiasJP1 = (TH1D *)hQ_frac_detJP1->Clone();
    hQBiasJP1->Divide(hQ_frac_par);
    // trigger bias JP2
    TH1D *hQBiasJP2 = (TH1D *)hQ_frac_detJP2->Clone();
    hQBiasJP2->Divide(hQ_frac_par);

    drawBias(hQ_frac_det, hG_frac_det, hQ_frac_par, hG_frac_par, hQBias, "All");
    drawBias(hQ_frac_detJP0, hG_frac_detJP0, hQ_frac_par, hG_frac_par, hQBiasJP0, "JP0");
    drawBias(hQ_frac_detJP1, hG_frac_detJP1, hQ_frac_par, hG_frac_par, hQBiasJP1, "JP1");
    drawBias(hQ_frac_detJP2, hG_frac_detJP2, hQ_frac_par, hG_frac_par, hQBiasJP2, "JP2");
    /*
        // draw ratios for parton level
        TCanvas *crp = new TCanvas("cratios", "", 900, 450);
        crp->Divide(2, 1);
        // gPad->SetLeftMargin(0.15);
        // gPad->SetBottomMargin(0.20);
        crp->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        hQ_frac_par->SetLineColor(2);
        hQ_frac_par->SetLineWidth(2);
        hQ_frac_par->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
        hQ_frac_par->GetYaxis()->SetTitle("fractions");
        hQ_frac_par->GetYaxis()->SetTitleOffset(1.5);
        hQ_frac_par->GetYaxis()->SetRangeUser(0.0, 1.0);
        hQ_frac_par->Draw("HIST E");

        hG_frac_par->SetLineColor(4);
        hG_frac_par->SetLineWidth(2);
        hG_frac_par->SetLineStyle(1);
        hG_frac_par->Draw("HIST E SAME");

        hG_frac_det->SetLineColor(4);
        hG_frac_det->SetLineWidth(2);
        hG_frac_det->SetLineStyle(2);
        hG_frac_det->Draw("HIST E SAME");

        hQ_frac_det->SetLineColor(2);
        hQ_frac_det->SetLineWidth(2);
        hQ_frac_det->SetLineStyle(2);
        hQ_frac_det->Draw("HIST E SAME");

        TLegend *lg = new TLegend(0.3, 0.75, 0.8, 0.90);
        lg->SetNColumns(2);
        lg->AddEntry(hQ_frac_par, " Quark (particle)", "l");
        lg->AddEntry(hQ_frac_det, " Quark (detector)", "l");
        lg->AddEntry(hG_frac_par, " Gluon (particle)", "l");
        lg->AddEntry(hG_frac_det, " Gluon (detector)", "l");
        lg->Draw();
        gPad->Update();
        crp->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        hQBias->SetLineColor(2);
        hQBias->SetMarkerStyle(8);
        hQBias->SetMarkerColor(2);
        hQBias->SetLineWidth(2);
        hQBias->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
        hQBias->GetYaxis()->SetRangeUser(0.5, 1.5);
        hQBias->GetYaxis()->SetTitle("Trigger Bias");
        hQBias->GetYaxis()->SetTitleOffset(1.5);
        hQBias->Draw("E");
        hQBias->Fit("pol0");
        hQBias->GetFunction("pol0")->SetLineColor(2);
        hQBias->GetFunction("pol0")->SetLineStyle(2);

        crp->Update();
        crp->SaveAs("triggerBias.pdf");

        cout << "double sys_bias[13]={";
        for (int i = 1; i <= hQBias->GetNbinsX(); i++)
        {
            double bias = hQBias->GetBinContent(i);
            double errbias = hQBias->GetBinError(i);
            if (i == hQBias->GetNbinsX())
                cout << 1 - bias << "};" << endl;
            else
                cout << 1 - bias << ",";
        }
        cout << "double sys_biaserr[13]={";
        for (int i = 1; i <= hQBias->GetNbinsX(); i++)
        {
            double bias = hQBias->GetBinContent(i);
            double errbias = hQBias->GetBinError(i);
            if (i == hQBias->GetNbinsX())
                cout << errbias << "};" << endl;
            else
                cout << errbias << ",";
        }
    */
}

void drawBias(TH1D *hDetQf, TH1D *hDetGf, TH1D *hParQf, TH1D *hParGf, TH1D *hQBias, const char *trig)
{
    TCanvas *can = new TCanvas("can", "", 900, 450);
    can->Divide(2, 1);
    can->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hParQf->SetLineColor(2);
    hParQf->SetLineWidth(2);
    hParQf->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hParQf->GetYaxis()->SetTitle("fractions");
    hParQf->GetYaxis()->SetTitleOffset(1.5);
    hParQf->GetYaxis()->SetRangeUser(0.0, 1.0);
    hParQf->Draw("HIST E");

    hParGf->SetLineColor(4);
    hParGf->SetLineWidth(2);
    hParGf->SetLineStyle(1);
    hParGf->Draw("HIST E SAME");

    hDetGf->SetLineColor(4);
    hDetGf->SetLineWidth(2);
    hDetGf->SetLineStyle(2);
    hDetGf->Draw("HIST E SAME");

    hDetQf->SetLineColor(2);
    hDetQf->SetLineWidth(2);
    hDetQf->SetLineStyle(2);
    hDetQf->Draw("HIST E SAME");

    TLegend *leg = new TLegend(0.3, 0.75, 0.8, 0.90);
    leg->SetNColumns(2);
    leg->AddEntry(hParQf, " Quark (particle)", "l");
    leg->AddEntry(hDetQf, " Quark (detector)", "l");
    leg->AddEntry(hParGf, " Gluon (particle)", "l");
    leg->AddEntry(hDetGf, " Gluon (detector)", "l");
    leg->Draw();
    gPad->Update();
    can->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    hQBias->SetLineColor(2);
    hQBias->SetMarkerStyle(8);
    hQBias->SetMarkerColor(2);
    hQBias->SetLineWidth(2);
    hQBias->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    hQBias->GetYaxis()->SetRangeUser(0.5, 1.5);
    hQBias->GetYaxis()->SetTitle("Trigger Bias");
    hQBias->GetYaxis()->SetTitleOffset(1.5);
    hQBias->Draw("E");
    hQBias->Fit("pol0");
    hQBias->GetFunction("pol0")->SetLineColor(2);
    hQBias->GetFunction("pol0")->SetLineStyle(2);

    can->Update();
    can->SaveAs(Form("triggerBias_%s.pdf", trig));

    cout << "double sys_bias[13]={";
    for (int i = 1; i <= hQBias->GetNbinsX(); i++)
    {
        double bias = hQBias->GetBinContent(i);
        double errbias = hQBias->GetBinError(i);
        if (i == hQBias->GetNbinsX())
            cout << 1 - bias << "};" << endl;
        else
            cout << 1 - bias << ",";
    }
    cout << "double sys_biaserr[13]={";
    for (int i = 1; i <= hQBias->GetNbinsX(); i++)
    {
        double bias = hQBias->GetBinContent(i);
        double errbias = hQBias->GetBinError(i);
        if (i == hQBias->GetNbinsX())
            cout << errbias << "};" << endl;
        else
            cout << errbias << ",";
    }
}
