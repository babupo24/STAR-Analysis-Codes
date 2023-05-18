#ifndef XSEC_THEORY_DATA_H
#define XSEC_THEORY_DATA_H

void drawXsecMoneyPlot(TH1D *hTheory, TH1D *hcom, TH1D *hpyth, TH1D *hsys_err, TH1D *hstat_err, TH1D *htotal_err, TCanvas *can)
{
    can->cd();
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.15);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    gPad->SetGrid(0, 0);

    hpyth->GetXaxis()->SetLabelSize(0);
    hpyth->GetYaxis()->SetTitle("#font[22]{#frac{d#sigma^{#pi^{+}#pi^{-}}}{dM^{#pi^{+}#pi^{-}}}  (#frac{pb}{GeV})}");
    hpyth->GetYaxis()->CenterTitle();
    hpyth->GetYaxis()->SetTickLength(0.02);
    hpyth->GetYaxis()->SetTitleOffset(1.5);
    hpyth->GetXaxis()->CenterTitle();

    TAxis *yAxisP = hpyth->GetYaxis();
    yAxisP->SetLabelFont(22);

    hpyth->SetMaximum(5e10);
    hpyth->SetMinimum(10);
    hpyth->SetMarkerStyle(8);
    hpyth->SetMarkerColor(1);
    hpyth->SetLineColor(1);
    hpyth->SetLineWidth(2);
    hpyth->SetLineStyle(1);
    hpyth->Draw("E1");
    pad1->Update();

    hcom->SetLineColor(2);
    hcom->SetLineWidth(2);
    // hcom->SetMarkerStyle(20);
    hcom->SetMarkerColor(2);
    hcom->Draw("SAME E1");

    TLegend *lx = new TLegend(0.45, 0.53, 0.7, 0.65);
    // lx->SetNColumns(2);
    lx->SetTextSize(0.03);
    lx->AddEntry(hpyth, "#font[22]{ Pythia 6.4.28 @ Perugia 2012,}", "lep");
    lx->AddEntry(" ", "#font[22]{ PARP(90)=0.213, CTEQ6 PDFs}", " ");
    lx->AddEntry(hcom, "#font[22]{ Measured}", "lep");
    lx->Draw();

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.045);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.32, "#font[22]{#color[2]{STAR Preliminary 2012}}");
    tex.SetTextSize(0.033);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.07, "#font[22]{p + p #rightarrow #pi^{+}#pi^{-} + X  at #sqrt{s} = 200 GeV }");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.02, "#font[22]{|#eta^{#pi^{+}#pi^{-}}| < 1, 1 < p_{T}^{#pi^{+}#pi^{-}} < 15 (GeV/c), 0.02 < cone < 0.7}");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.006, "#font[22]{0.27 < M_{inv}^{#pi^{+}#pi^{-}} < 4.0 (GeV/c^{2}),  L_{int} = 26 (pb)^{-1} }");

    TPad *pad1clear = new TPad("pad1clear", "", 0.0, 0.0, 1.0, 1.0);
    pad1clear->SetLeftMargin(0.15);
    pad1clear->SetRightMargin(0.15);
    pad1clear->SetBottomMargin(0);
    pad1clear->SetFillColor(0);
    pad1clear->SetFillStyle(4000);
    pad1clear->SetFrameFillStyle(0);
    pad1clear->SetGrid(0, 0);
    pad1clear->Draw();
    pad1clear->cd();
    pad1clear->SetLogy();
    TH1D *hcom_c = (TH1D *)hcom->Clone();
    hcom_c->SetMaximum(5e10);
    hcom_c->SetMinimum(10);
    hcom_c->GetXaxis()->SetTitle("");
    hcom_c->GetXaxis()->SetLabelSize(0);
    hcom_c->GetYaxis()->SetTitle("");
    hcom_c->GetYaxis()->SetLabelSize(0);
    hcom_c->GetYaxis()->SetTickLength(0.02);
    hcom_c->Draw(" X+ Y+");
    pad1clear->Update();

    can->cd();
    TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.15);
    pad2->Draw();
    pad2->cd();
    pad2->SetGrid(0, 0);

    pad2->Update();

    hsys_err->GetYaxis()->SetRangeUser(-0.6, 0.6);
    hsys_err->GetYaxis()->SetLabelSize(0.11);
    hsys_err->GetYaxis()->SetTitleOffset(0.55);
    // hsys_err->GetYaxis()->SetTitle("#font[22]{#frac{#sigma_{pyth} - #sigma_{mes}}{#sigma_{mes}}}");
    hsys_err->GetYaxis()->SetTitle("#font[22]{(#frac{#delta#sigma}{#sigma})_{mes}}");
    hsys_err->GetYaxis()->SetTitleSize(0.11);
    hsys_err->GetXaxis()->SetLabelOffset(0.02);
    TAxis *xAxis = hsys_err->GetXaxis();
    xAxis->SetLabelFont(22);
    TAxis *yAxis = hsys_err->GetYaxis();
    yAxis->SetLabelFont(22);
    hsys_err->GetXaxis()->SetTitle("#font[22]{M_{inv}^{#pi^{+}#pi^{-}} (GeV/c^{2})}");
    hsys_err->GetXaxis()->SetTitleOffset(1.1);
    hsys_err->GetXaxis()->SetTitleSize(0.11);
    hsys_err->GetXaxis()->CenterTitle();
    hsys_err->GetYaxis()->CenterTitle();
    hsys_err->GetXaxis()->SetLabelSize(0.11);

    hsys_err->GetXaxis()->SetTickLength(0.08);
    hsys_err->GetYaxis()->SetNdivisions(505);
    hsys_err->SetFillColorAlpha(8, 1.0);
    hsys_err->SetFillStyle(1000);
    hsys_err->Draw("E2");

    pad2->Update();
    TLine *l1 = new TLine(pad2->GetUxmin(), 0.0, pad2->GetUxmax(), 0.0);
    l1->SetLineWidth(2);
    l1->SetLineStyle(2);
    l1->Draw();
    pad2->Update();

    hstat_err->SetFillColorAlpha(2, 1);
    hstat_err->Draw("E2 SAME");

    // draw relative difference between pythia and measured
    TH1D *hrdiff = (TH1D *)hpyth->Clone();
    hrdiff->Add(hcom, -1);
    hrdiff->Divide(hcom);
    TH1D *hrdiff_c = (TH1D *)hrdiff->Clone();
    hrdiff_c->Reset();

    for (int i = 1; i <= hrdiff->GetNbinsX(); i++)
    {
        if (hrdiff->GetBinContent(i) == 0)
            continue;
        hrdiff_c->SetBinContent(i, hrdiff->GetBinContent(i));
        hrdiff_c->SetBinError(i, 0.0001);
    }
    // hrdiff_c->SetMarkerColor(2);
    //  hrdiff_c->SetMarkerStyle(8);
    /// hrdiff_c->SetMarkerSize(0.5);
    hrdiff_c->SetLineColor(1);
    hrdiff_c->SetLineWidth(2);
    hrdiff_c->Draw("e1 SAME");

    pad2->Update();
    TGaxis *axis2 = new TGaxis(pad2->GetUxmax(), pad2->GetUymin(), pad2->GetUxmax(), pad2->GetUymax(), -0.6, 0.6, 505, "+L");
    axis2->SetTitle("#font[22]{Ratio}");
    axis2->SetTitleSize(0.11);
    axis2->SetLabelSize(0.11);
    axis2->SetTextFont(22);
    axis2->SetLabelFont(22);
    axis2->SetLabelOffset(0.01);
    axis2->SetTitleOffset(0.4);
    axis2->CenterTitle();
    axis2->Draw();

    TLegend *lg1 = new TLegend(0.25, 0.85, 0.75, 0.95);
    lg1->SetNColumns(3);
    // lg1->AddEntry(hsys_err, "#font[22]{ #delta#sigma_{sys}}", "f");
    lg1->AddEntry(hsys_err, "#font[22]{(#frac{#delta#sigma}{#sigma})_{mes}^{Syst.}}", "f");
    // lg1->AddEntry(hstat_err, "#font[22]{ #delta#sigma_{stat}}", "f");
    lg1->AddEntry(hstat_err, "#font[22]{(#frac{#delta#sigma}{#sigma})_{mes}^{Stat.}}", "f");
    // lg1->AddEntry(hrdiff_c, "#font[22]{ #times 0.5}", "lp");
    lg1->AddEntry(hrdiff_c, "#font[22]{Ratio(= #frac{d#sigma_{pyth} - d#sigma_{mes}}{d#sigma_{mes}})}", "lp");
    lg1->SetTextSize(0.07);
    lg1->Draw();
    gPad->Update();

    TLatex tex1;
    tex1.SetTextSize(0.08);
    // tex1.DrawLatex(0.45, hsys_err->GetMinimum() + 0.25, "#font[12]{10\% luminosity measurement uncertainty not shown.}");
    tex1.DrawLatex(0.35, -0.48, "#font[22]{10\% luminosity uncertainty not shown.}");

    can->Update();
}

#endif