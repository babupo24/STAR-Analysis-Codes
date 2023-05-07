#ifndef HELPER_FUNCTION_H
#define HELPER_FUNCTION_H

void trigEff(TH1D *hcom, TH1D *hjp0, TH1D *hjp1, TH1D *hjp2)
{

    TH1D *h0 = (TH1D *)hjp0->Clone();
    h0->Add(hcom, -1);
    h0->Divide(hcom);
    TH1D *h1 = (TH1D *)hjp1->Clone();
    h1->Add(hcom, -1);
    h1->Divide(hcom);
    TH1D *h2 = (TH1D *)hjp2->Clone();
    h2->Add(hcom, -1);
    h2->Divide(hcom);

    TCanvas *can = new TCanvas("can", "can", 600, 450);
    can->cd();
    gPad->SetGrid(0, 0);
    h0->SetName(" jp0");
    h0->SetLineColor(3);
    h0->SetLineWidth(2);
    h0->GetYaxis()->SetRangeUser(-0.4, 0.4);
    h0->GetYaxis()->SetTitle("(#sigma_{jp*} #minus #sigma_{comb})/#sigma_{comb}");
    h0->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
    h0->Draw("E");

    h1->SetName(" jp1");
    h1->SetLineColor(4);
    h1->SetLineWidth(2);
    h1->Draw("E SAME");

    h2->SetName(" jp2");
    h2->SetLineColor(2);
    h2->SetLineWidth(2);
    h2->Draw("E SAME");
    can->Update();

    TLegend *leg = new TLegend(0.25, 0.2, 0.7, 0.25);
    leg->SetNColumns(3);
    leg->AddEntry(h0, " jp0", "le");
    leg->AddEntry(h1, " jp1", "le");
    leg->AddEntry(h2, " jp2", "le");
    leg->Draw();

    can->SaveAs("ResultsMinConeCut/sys_trigger_eff.pdf");
    can->SaveAs("ResultsMinConeCut/sys_trigger_eff.png");

    TF1 *f0 = new TF1("f0", "pol0", 2.2, 4.);
    f0->SetLineColor(3);
    f0->SetLineStyle(2);
    h0->Fit("f0", "R");
    TF1 *f1 = new TF1("f1", "pol0", 2.2, 4.);
    f1->SetLineColor(4);
    f1->SetLineStyle(2);
    h1->Fit("f1", "R");
    TF1 *f2 = new TF1("f2", "pol0", 2.2, 4.);
    f2->SetLineColor(2);
    f2->SetLineStyle(2);
    h2->Fit("f2", "R");
    gPad->Update();
    TLine *ln = new TLine(can->GetUxmin(), 0, can->GetUxmax(), 0);
    ln->SetLineStyle(2);
    ln->Draw();
    

    cout << "JP0 = " << f0->Eval(h0->GetBinCenter(11)) << ", " << f0->Eval(h0->GetBinCenter(12)) << ", " << f0->Eval(h0->GetBinCenter(13)) << endl;
    cout << "JP1 = " << f1->Eval(h1->GetBinCenter(11)) << ", " << f1->Eval(h1->GetBinCenter(12)) << ", " << f1->Eval(h1->GetBinCenter(13)) << endl;
    cout << "JP2 = " << f2->Eval(h2->GetBinCenter(11)) << ", " << f2->Eval(h2->GetBinCenter(12)) << ", " << f2->Eval(h2->GetBinCenter(13)) << endl;


    can->Update();
    can->SaveAs("ResultsMinConeCut/sys_trigger_eff_fitted.pdf");
    can->SaveAs("ResultsMinConeCut/sys_trigger_eff_fitted.png");
}

#endif