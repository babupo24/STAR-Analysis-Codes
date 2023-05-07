#include <iostream>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iomanip>

using namespace std;
void drawXsec(TH1D *hcom, TH1D *hpyth, TH1D *hjp2, TH1D *hjp1, TH1D *hjp0, TCanvas *can);
void setBinError(TH1D *hRatio, TH1D *hData, TH1D *hEmbed);
double getBinError(double dataBinCt, double dataBinErr, double embedBinCt, double embedBinErr);

void drawXsecMoneyPlot(TH1D *hcom, TH1D *hpyth, TH1D *hsys_err, TH1D *hstat_err, TH1D *htotal_err, TCanvas *can);

void drawErrors(TH1D *htotal, TH1D *hsys, TH1D *hpid, TH1D *heff, TH1D *htrgspr, TH1D *hstat, TH1D *hbias);

void xsecAna()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    // gStyle->SetLabelFont(13, "xyz");
    gStyle->SetLegendBorderSize(0);

    const int nBins = 13;
    Double_t nBinsEdges[nBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};
    // binning for cross-section money plot. Add bins at underflow and overflow
    const int nBinsU = 15;
    Double_t nBinsEdgesU[nBinsU + 1] = {0.05, 0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0, 4.25};

    // trigger prescales
    Double_t presJP0 = 0.00707; //(1. / 141.276);
    Double_t presJP1 = 0.3892;  //(1. / 2.56938);
    Double_t presJP2 = 1.0;

    TFile *funfolding = new TFile("../BinResolutionAndUnfolding/UnfoldingResults/unfoldingOutput.root", "r");
    // unfolded in |eta|<1
    TH1D *hun_jp0 = (TH1D *)funfolding->Get("UnfoldedJP0");
    TH1D *hun_jp1 = (TH1D *)funfolding->Get("UnfoldedJP1");
    TH1D *hun_jp2 = (TH1D *)funfolding->Get("UnfoldedJP2");

    // get particle level cross section (for comparison)
    TFile *fparticle = new TFile("hist_XsecParticle.root", "r");
    TH1D *hXsecParticle = (TH1D *)fparticle->Get("hUnbiasedMinv");
    Double_t scalePart = hXsecParticle->Integral();

    // luminosity
    Double_t L_jp0 = 0.16;
    Double_t L_jp1 = 7.68;
    Double_t L_jp2 = 18.80;
    Double_t L_jptotal = 26.64;
    // tracking efficiency:
    Double_t eff_trk_jp0 = 0.92;
    Double_t eff_trk_jp1 = 0.92;
    Double_t eff_trk_jp2 = 0.92;

    // tracking efficiency combined
    Double_t trkPairEff = 0.7632;
    Double_t trkPairEffErr = 0.00356;
    // tracking efficiency per bin
    double trkeff_jp0[13] = {0.888706, 0.890627, 0.869874, 0.880854, 0.889219, 0.889851, 0.895486, 0.904407, 0.90956, 0.931546, 0.916641, 0.911637, 0.891178};
    double trkefferr_jp0[13] = {1.06954e-07, 7.77831e-08, 7.15327e-08, 1.0936e-07, 1.48747e-07, 2.8056e-07, 4.68323e-07, 7.20706e-07, 1.2522e-06, 2.07082e-06, 4.45592e-06, 1.12821e-05, 3.39753e-05};
    double trkeff_jp1[13] = {0.884145, 0.884121, 0.86267, 0.883308, 0.890884, 0.893681, 0.900803, 0.900202, 0.920559, 0.924294, 0.92829, 0.908445, 0.875614};
    double trkefferr_jp1[13] = {3.15813e-07, 2.27766e-07, 2.03141e-07, 2.79039e-07, 3.42844e-07, 5.74394e-07, 9.00969e-07, 1.29421e-06, 1.88018e-06, 3.26229e-06, 5.97537e-06, 1.4349e-05, 4.19778e-05};
    double trkeff_jp2[13] = {0.876888, 0.882371, 0.854922, 0.880348, 0.884055, 0.892542, 0.899577, 0.908812, 0.914405, 0.913439, 0.932468, 0.916973, 0.865696};
    double trkefferr_jp2[13] = {7.20508e-07, 5.08041e-07, 4.53362e-07, 5.80138e-07, 6.77582e-07, 1.05714e-06, 1.54391e-06, 2.06328e-06, 3.09668e-06, 5.15335e-06, 8.37119e-06, 1.777e-05, 5.24983e-05};
    // dipion triggering efficiency
    // jp0 efficiency and error
    Double_t eff_trg_jp0[13] = {0.0647609, 0.0663649, 0.0734078, 0.0840551, 0.109116, 0.150331, 0.197242, 0.261445, 0.408169, 0.447309, 0.501666, 0.668775, 0.759649};
    Double_t efferr_trg_jp0[13] = {0.000446209, 0.000330274, 0.000310635, 0.000505293, 0.000819134, 0.00188765, 0.00383486, 0.00661998, 0.00748233, 0.0149041, 0.0275889, 0.0519418, 0.0985095};
    // jp1 efficiency and error
    Double_t eff_trg_jp1[13] = {0.00762761, 0.00806872, 0.00943897, 0.0127224, 0.0203221, 0.0349688, 0.0512922, 0.0836607, 0.162998, 0.195088, 0.246034, 0.425724, 0.533614};
    Double_t efferr_trg_jp1[13] = {8.31865e-05, 6.43366e-05, 6.41374e-05, 0.000102156, 0.000170519, 0.00041934, 0.00085, 0.00192461, 0.00316734, 0.00682422, 0.0121178, 0.0320969, 0.0574833};
    // jp2 efficiency and error
    Double_t eff_trg_jp2[13] = {0.00151809, 0.00163691, 0.0019597, 0.0029936, 0.00543582, 0.0104283, 0.0176766, 0.0307442, 0.0640216, 0.0869714, 0.119685, 0.258225, 0.361955};
    Double_t efferr_trg_jp2[13] = {2.11924e-05, 1.42757e-05, 1.59223e-05, 1.75488e-05, 3.76751e-05, 0.000108269, 0.000299922, 0.00065822, 0.00114225, 0.00334942, 0.00557087, 0.0197869, 0.0413919};
    // trigger  bias
    double sys_trigBias[13] = {0.0217904, 0.020958, 0.0230904, 0.0138005, 0.00169444, 0.00280884, 0.00213932, 0.00525629, 0.0149015, 0.0155368, 0.0102302, 0.00550931, 0.033424};
    double sys_trigBiasErr[13] = {0.00694856, 0.0051036, 0.00451617, 0.00715329, 0.00999482, 0.0193596, 0.0318765, 0.0450885, 0.0296183, 0.0527495, 0.0867533, 0.151848, 0.207589};
    // cross-section factor (1/(eff_trk x eff_trg x L))
    Double_t binWidth[13] = {0};
    Double_t f_jp0[13] = {0};
    Double_t f_jp1[13] = {0};
    Double_t f_jp2[13] = {0};

    for (int i = 0; i < 13; i++)
    {
        binWidth[i] = hun_jp0->GetBinWidth(i + 1);
        f_jp0[i] = (Double_t)(binWidth[i] * trkeff_jp0[i] * eff_trg_jp0[i] * L_jp0);
        f_jp1[i] = (Double_t)(binWidth[i] * trkeff_jp1[i] * eff_trg_jp1[i] * L_jp1);
        f_jp2[i] = (Double_t)(binWidth[i] * trkeff_jp2[i] * eff_trg_jp2[i] * L_jp2);
    }
    // yields per bin (dN)
    Double_t dN_jp0[nBins] = {0};            // unfolded dipion raw yields
    Double_t stat_err_jp0[nBins] = {0};      // stat error from unfolding
    Double_t xsec_stat_err_jp0[nBins] = {0}; // stat error from unfolding
    Double_t norm_stat_err_jp0[nBins] = {0}; // total error normalized by cross section yields
    Double_t xsec_jp0[nBins] = {0};          // cross section
    Double_t dN_jp1[nBins] = {0};
    Double_t stat_err_jp1[nBins] = {0};
    Double_t xsec_stat_err_jp1[nBins] = {0};
    Double_t norm_stat_err_jp1[nBins] = {0};
    Double_t xsec_jp1[nBins] = {0};
    Double_t dN_jp2[nBins] = {0};
    Double_t stat_err_jp2[nBins] = {0};
    Double_t xsec_stat_err_jp2[nBins] = {0};
    Double_t norm_stat_err_jp2[nBins] = {0};
    Double_t xsec_jp2[nBins] = {0};

    // xsec for combined trigger
    Double_t dN_jp[nBins] = {0};
    Double_t stat_err_jp[nBins] = {0};
    Double_t xsec_stat_err_jp[nBins] = {0};
    Double_t total_stat_err_jp[nBins] = {0};
    Double_t norm_stat_err_jp[nBins] = {0};
    Double_t xsec_jp[nBins] = {0};

    Double_t sys_trg_err[13] = {0};
    Double_t sys_eff_err[13] = {0};
    double sys_pid_err[13] = {0.161342, 0.164592, 0.185448, 0.257799, 0.221166, 0.215697, 0.187357, 0.14168, 0.122346, 0.0988125, 0.111084, 0.136663, 0.153807};
    Double_t sys_err[13] = {0};
    Double_t norm_stat_err[13] = {0};

    // pythia cross ection
    Double_t xsec_pyth[13] = {0};
    Double_t stat_err_pyth[13] = {0};

    int ww = 20;
    // cout << "total" << setw(ww) << "jp2" << setw(ww) << "jp1" << setw(ww) << "jp0" << endl;

    for (int i = 0; i < nBins; i++)
    {
        if (hun_jp0->GetBinContent(i + 1) == 0)
            continue;
        // jp0
        dN_jp0[i] = hun_jp0->GetBinContent(i + 1);
        stat_err_jp0[i] = hun_jp0->GetBinError(i + 1);
        xsec_jp0[i] = (Double_t)dN_jp0[i] / f_jp0[i];
        xsec_stat_err_jp0[i] = sqrt(xsec_jp0[i]);
        norm_stat_err_jp0[i] = (Double_t)(xsec_stat_err_jp0[i] / xsec_jp0[i]);
        // jp1
        dN_jp1[i] = hun_jp1->GetBinContent(i + 1);
        stat_err_jp1[i] = hun_jp1->GetBinError(i + 1);
        xsec_jp1[i] = (Double_t)dN_jp1[i] / f_jp1[i];
        xsec_stat_err_jp1[i] = sqrt(xsec_jp1[i]);
        norm_stat_err_jp1[i] = (Double_t)(xsec_stat_err_jp1[i] / xsec_jp1[i]);
        // jp2
        dN_jp2[i] = hun_jp2->GetBinContent(i + 1);
        stat_err_jp2[i] = hun_jp2->GetBinError(i + 1);
        xsec_jp2[i] = (Double_t)dN_jp2[i] / f_jp2[i];
        xsec_stat_err_jp2[i] = sqrt(xsec_jp2[i]);
        norm_stat_err_jp2[i] = (Double_t)(xsec_stat_err_jp2[i] / xsec_jp2[i]);

        // weights for the combining triggers
        double w_jp0 = 1.0 / xsec_stat_err_jp0[i];
        double w2_jp0 = w_jp0 * w_jp0;
        double w_jp1 = 1.0 / xsec_stat_err_jp1[i];
        double w2_jp1 = w_jp1 * w_jp1;
        double w_jp2 = 1.0 / xsec_stat_err_jp2[i];
        double w2_jp2 = w_jp2 * w_jp2;
        // cout << i << "  " << xsec_jp2[i] << "  " << stat_err_jp2[i] << "  " << sqrt(xsec_jp2[i]) << endl;
        //  weighted average of all triggers
        xsec_jp[i] = (w2_jp0 * xsec_jp0[i] + w2_jp1 * xsec_jp1[i] + w2_jp2 * xsec_jp2[i]) / (w2_jp0 + w2_jp1 + w2_jp2);
        stat_err_jp[i] = 1. / sqrt(w2_jp0 + w2_jp1 + w2_jp2);
        norm_stat_err[i] = stat_err_jp[i] / xsec_jp[i];

        // xsec pythia
        xsec_pyth[i] = (hXsecParticle->GetBinContent(i + 1) / binWidth[i]);
        stat_err_pyth[i] = hXsecParticle->GetBinError(i + 1);

        // find a trigger spread relative to the combined cross section and apply as systematic error
        double spr_jp0 = (xsec_jp[i] - xsec_jp0[i]) / xsec_jp[i];
        double spr_jp1 = (xsec_jp[i] - xsec_jp1[i]) / xsec_jp[i];
        double spr_jp2 = (xsec_jp[i] - xsec_jp2[i]) / xsec_jp[i];
        double spr_avg = (spr_jp0 + spr_jp1 + spr_jp2) / 3.0;
        sys_trg_err[i] = spr_avg;
        // cout << spr_avg << ", ";
        // weights for combined systematics
        double sysw_pEffjp0 = 1.0 / trkefferr_jp0[i]; // tracking efficiency
        double sysw2_pEffjp0 = sysw_pEffjp0 * sysw_pEffjp0;
        double sysw_trgjp0 = 1.0 / efferr_trg_jp0[i]; // triggering efficiency
        double sysw2_trgjp0 = sysw_trgjp0 * sysw_trgjp0;

        double sysw_pEffjp1 = 1.0 / trkefferr_jp1[i];
        double sysw2_pEffjp1 = sysw_pEffjp1 * sysw_pEffjp1;
        double sysw_trgjp1 = 1.0 / efferr_trg_jp1[i];
        double sysw2_trgjp1 = sysw_trgjp1 * sysw_trgjp1;

        double sysw_pEffjp2 = 1.0 / trkefferr_jp2[i];
        double sysw2_pEffjp2 = sysw_pEffjp2 * sysw_pEffjp2;
        double sysw_trgjp2 = 1.0 / efferr_trg_jp2[i];
        double sysw2_trgjp2 = sysw_trgjp2 * sysw_trgjp2;

        double sysw_err_trk = 1. / (sqrt(sysw2_pEffjp0 + sysw2_pEffjp1 + sysw2_pEffjp2));
        double sysw_err_trg = 1. / (sqrt(sysw2_trgjp0 + sysw2_trgjp1 + sysw2_trgjp2));
        sys_eff_err[i] = sqrt(sysw_err_trg * sysw_err_trg + sysw_err_trk * sysw_err_trk);
        // systematic error
        // sys_err[i] = sqrt(pow(sys_trigBias[i], 2) + pow(sysw_err_trg, 2) + pow(sysw_err_trk, 2) + pow(spr_avg, 2) + pow(sys_pid_err[i], 2));
        sys_err[i] = sqrt(pow(sys_trigBias[i], 2) + pow(spr_avg, 2) + pow(sys_pid_err[i], 2));
    }
    // cout << endl;
    // book histograms for money plot
    TH1D *hxsec_jp0 = new TH1D("hxsec_jp0", "", nBinsU, nBinsEdgesU);
    TH1D *hxsec_jp1 = new TH1D("hxsec_jp1", "", nBinsU, nBinsEdgesU);
    TH1D *hxsec_jp2 = new TH1D("hxsec_jp2", "", nBinsU, nBinsEdgesU);
    TH1D *hxsec_pyth = new TH1D("hxsec_pyth", "", nBinsU, nBinsEdgesU);
    TH1D *hxsec_comb = new TH1D("hxsec_comb", "", nBinsU, nBinsEdgesU);

    // for systematic error
    TH1D *htrg_spr_err = new TH1D("htrg_spr_err", "", nBinsU, nBinsEdgesU); // trig spread and efficiencies
    TH1D *hpid_err = new TH1D("hpid_err", "", nBinsU, nBinsEdgesU);         // efficiencies error
    TH1D *hbias_err = new TH1D("hbias_err", "", nBinsU, nBinsEdgesU);       // efficiencies error
    TH1D *heff_err = new TH1D("heff_err", "", nBinsU, nBinsEdgesU);         // efficiencies error
    TH1D *hsys_err = new TH1D("hsys_err", "", nBinsU, nBinsEdgesU);         // tota systematic error
    TH1D *hstat_err = new TH1D("hstat_err", "", nBinsU, nBinsEdgesU);       // xsec statistical error
    TH1D *htotal_err = new TH1D("htotal_err", "", nBinsU, nBinsEdgesU);     // stat + sys combined

    for (int i = 0; i < nBinsU; i++)
    {
        if (i < 1 || i > nBins)
        {
            continue;
        }
        hxsec_jp0->SetBinContent(i + 1, xsec_jp0[i - 1]);
        hxsec_jp0->SetBinError(i + 1, stat_err_jp0[i - 1]);

        hxsec_jp1->SetBinContent(i + 1, xsec_jp1[i - 1]);
        hxsec_jp1->SetBinError(i + 1, stat_err_jp1[i - 1]);

        hxsec_jp2->SetBinContent(i + 1, xsec_jp2[i - 1]);
        hxsec_jp2->SetBinError(i + 1, stat_err_jp2[i - 1]);

        hxsec_comb->SetBinContent(i + 1, xsec_jp[i - 1]);
        hxsec_comb->SetBinError(i + 1, stat_err_jp[i - 1]);

        hxsec_pyth->SetBinContent(i + 1, xsec_pyth[i - 1]);
        hxsec_pyth->SetBinError(i + 1, stat_err_pyth[i - 1]);

        // fill error histograms
        hpid_err->SetBinContent(i + 1, 0.0);
        hpid_err->SetBinError(i + 1, sys_pid_err[i - 1]);

        hbias_err->SetBinContent(i + 1, 0.0);
        hbias_err->SetBinError(i + 1, sys_trigBias[i - 1]);

        htrg_spr_err->SetBinContent(i + 1, 0.0);
        htrg_spr_err->SetBinError(i + 1, sys_trg_err[i - 1]);

        heff_err->SetBinContent(i + 1, 0.0);
        heff_err->SetBinError(i + 1, sys_eff_err[i - 1]);

        hsys_err->SetBinContent(i + 1, 0.0);
        hsys_err->SetBinError(i + 1, sys_err[i - 1]);

        hstat_err->SetBinContent(i + 1, 0.0);
        hstat_err->SetBinError(i + 1, norm_stat_err[i - 1]);

        htotal_err->SetBinContent(i + 1, 0.0);
        htotal_err->SetBinError(i + 1, sqrt(pow(norm_stat_err[i - 1], 2) + pow(sys_err[i - 1], 2)));
    }

    //----
    // draw money plot
    TCanvas *cmpt = new TCanvas("cmpt", "", 0, 0, 700, 800);
    drawXsecMoneyPlot(hxsec_comb, hxsec_pyth, hsys_err, hstat_err, htotal_err, cmpt);
    cmpt->SaveAs("Results/xsec_MoneyPlot.pdf");
    // delete cmpt;

    // TCanvas *cxsec = new TCanvas("cxsec", "", 0, 0, 700, 800);
    // drawXsec(hxsec_comb, hxsec_pyth, hxsec_jp2, hxsec_jp1, hxsec_jp0, cxsec);
    // cxsec->SaveAs("Results/xsec_all.pdf");
    // delete cxsec;

    // drawErrors(htotal_err, hsys_err, hpid_err, heff_err, htrg_spr_err, hstat_err, hbias_err);

} // main

void drawXsec(TH1D *hcom, TH1D *hpyth, TH1D *hjp2, TH1D *hjp1, TH1D *hjp0, TCanvas *can)
{
    can->cd();
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogy();

    hpyth->GetXaxis()->SetLabelSize(0);
    hpyth->GetYaxis()->SetTitle("#font[22]{#frac{d#sigma^{#pi^{+}#pi^{-}}}{dM^{#pi^{+}#pi^{-}}}(= #frac{1}{L #upoint  #epsilon_{trk}^{#pi^{+}} #upoint #epsilon_{trk}^{#pi^{-}} #upoint #epsilon_{trg}^{#pi^{+}#pi^{-}}} #times #frac{dN^{#pi^{+}#pi^{-}}_{true}}{dM_{inv}^{#pi^{+}#pi^{-}}})}");
    hpyth->GetYaxis()->CenterTitle();
    hpyth->GetYaxis()->SetTitleOffset(1.5);
    hpyth->GetXaxis()->CenterTitle();
    hpyth->Scale(hcom->Integral() / hpyth->Integral());
    hpyth->SetMaximum(5e10);
    hpyth->SetMarkerStyle(20);
    hpyth->SetMarkerColor(1);
    hpyth->SetLineColor(1);
    hpyth->SetLineWidth(2);
    hpyth->SetLineStyle(1);
    hpyth->Draw("e1 ");

    hcom->SetMarkerStyle(20);
    hcom->SetMarkerColor(2);
    hcom->SetLineColor(2);
    hcom->SetLineWidth(2);
    hcom->Draw("e1 SAME");

    hjp0->SetMarkerStyle(20);
    hjp0->SetMarkerColor(3);
    hjp0->SetLineColor(3);
    hjp0->SetLineWidth(2);
    hjp0->SetLineStyle(2);
    hjp0->Draw("e1 SAME");

    hjp1->SetMarkerStyle(20);
    hjp1->SetMarkerColor(4);
    hjp1->SetLineColor(4);
    hjp1->SetLineWidth(2);
    hjp1->SetLineStyle(2);
    hjp1->Draw(" e1 SAME");

    hjp2->SetMarkerStyle(22);
    hjp2->SetMarkerColor(6);
    hjp2->SetLineColor(6);
    hjp2->SetLineWidth(2);
    hjp2->SetLineStyle(2);
    hjp2->Draw("e1 SAME");

    TLegend *lx = new TLegend(0.7, 0.5, 0.85, 0.8);
    // lx->SetNColumns(4);
    lx->SetTextSize(0.04);
    lx->AddEntry(hcom, " Comb.", "lp");
    lx->AddEntry(hpyth, " Pythia", "lp");
    lx->AddEntry(hjp2, " jp2", "lp");
    lx->AddEntry(hjp1, " jp1", "lp");
    lx->AddEntry(hjp0, " jp0", "lp");
    lx->Draw();

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.04);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.3, "#font[22]{p + p #rightarrow #pi^{+}#pi^{-} + X  at #sqrt{s} = 200 GeV }");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.1, "#font[22]{|#eta^{#pi^{+}#pi^{-}}| < 1}");

    gPad->Update();

    can->cd();
    TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    gPad->SetGrid(1, 1);

    TH1D *h_R2 = (TH1D *)hpyth->Clone();
    h_R2->Add(hjp2, -1);
    h_R2->Divide(hjp2);
    h_R2->SetBit(TH1::kNoTitle);
    setBinError(h_R2, hpyth, hjp2);

    h_R2->GetYaxis()->SetTitle("#font[22]{#frac{Pythia - X}{X}}");
    h_R2->GetXaxis()->SetTitle("#font[22]{M^{#pi^{+}#pi^{-}}_{inv} (GeV/c^{2})}");
    h_R2->GetYaxis()->CenterTitle();
    h_R2->GetYaxis()->SetRangeUser(-1.1, 1.1);
    h_R2->GetYaxis()->SetLabelSize(0.11);
    h_R2->GetYaxis()->SetTitleOffset(0.55);
    h_R2->GetYaxis()->SetTitleSize(0.11);
    h_R2->GetXaxis()->SetLabelOffset(0.02);
    h_R2->GetXaxis()->SetTitleOffset(1.1);
    h_R2->GetXaxis()->SetTitleSize(0.11);
    h_R2->GetXaxis()->SetLabelSize(0.11);
    h_R2->GetXaxis()->SetTickLength(0.08);
    h_R2->GetYaxis()->SetNdivisions(505);
    h_R2->SetLineColor(6);
    h_R2->SetMarkerStyle(0);
    h_R2->SetLineStyle(2);
    h_R2->Draw(" e1");

    TH1D *h_R1 = (TH1D *)hpyth->Clone();
    h_R1->Add(hjp1, -1);
    h_R1->Divide(hjp1);
    h_R1->SetBit(TH1::kNoTitle);
    setBinError(h_R1, hpyth, hjp1);
    h_R1->SetLineColor(4);
    h_R1->SetLineStyle(2);
    h_R1->SetMarkerStyle(0);
    h_R1->Draw("e1 same ");

    TH1D *h_R0 = (TH1D *)hpyth->Clone();
    h_R0->Add(hjp0, -1);
    h_R0->Divide(hjp0);
    h_R0->SetBit(TH1::kNoTitle);
    setBinError(h_R0, hpyth, hjp0);
    h_R0->SetLineColor(3);
    h_R0->SetLineStyle(2);
    h_R0->SetMarkerStyle(0);
    // h_R0->SetMarkerColor(1);
    h_R0->Draw("same ");

    TH1D *h_RC = (TH1D *)hpyth->Clone();
    h_RC->Add(hcom, -1);
    h_RC->Divide(hcom);
    setBinError(h_RC, hpyth, hcom);
    h_RC->SetLineColor(2);
    h_RC->SetLineStyle(1);
    h_RC->SetMarkerStyle(0);
    h_RC->Draw("e1 same ");

    gPad->Update();

    can->Update();
}
void drawXsecMoneyPlot(TH1D *hcom, TH1D *hpyth, TH1D *hsys_err, TH1D *hstat_err, TH1D *htotal_err, TCanvas *can)
{
    can->cd();
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogy();

    hpyth->GetXaxis()->SetLabelSize(0);
    hpyth->GetYaxis()->SetTitle("#font[22]{#frac{d#sigma^{#pi^{+}#pi^{-}}}{dM^{#pi^{+}#pi^{-}}}  (pb)}");
    hpyth->GetYaxis()->CenterTitle();
    hpyth->GetYaxis()->SetTitleOffset(1.5);
    hpyth->GetXaxis()->CenterTitle();
    hpyth->Scale(hcom->Integral() / hpyth->Integral());
    hpyth->SetMaximum(5e10);
    // hpyth->SetMarkerStyle(20);
    // hpyth->SetMarkerColor(4);
    hpyth->SetLineColor(1);
    hpyth->SetLineWidth(2);
    hpyth->SetLineStyle(1);
    hpyth->Draw("e1");

    hcom->SetMarkerStyle(20);
    hcom->SetMarkerColor(2);
    hcom->SetLineColor(2);
    hcom->SetLineWidth(2);
    hcom->Draw("E1 SAME");

    TLegend *lx = new TLegend(0.56, 0.7, 0.85, 0.74);
    lx->SetNColumns(2);
    lx->SetTextSize(0.033);
    lx->AddEntry(hcom, "#font[22]{ Measured}", "lp");
    lx->AddEntry(hpyth, "#font[22]{ Pythia}", "lp");
    lx->Draw();

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.045);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.4, "#font[22]{#color[2]{STAR Preliminary 2012}}");
    tex.SetTextSize(0.033);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.15, "#font[22]{p + p #rightarrow #pi^{+}#pi^{-} + X  at #sqrt{s} = 200 GeV }");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.05, "#font[22]{|#eta^{#pi^{+}#pi^{-}}| < 1, 1 < p_{T}^{#pi^{+}#pi^{-}} < 15 (GeV/c), 0.02 < cone < 0.7}");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.02, "#font[22]{L_{int} = 26 (pb)^{-1}, }");

    gPad->Update();

    can->cd();
    TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    pad2->SetGrid(0, 0);

    pad2->Update();

    htotal_err->GetYaxis()->SetRangeUser(-0.58, 0.58);
    htotal_err->GetYaxis()->SetLabelSize(0.11);
    htotal_err->GetYaxis()->SetTitleOffset(0.55);
    htotal_err->GetYaxis()->SetTitle("#font[22]{#frac{#delta#sigma}{#sigma}}");
    htotal_err->GetYaxis()->SetTitleSize(0.11);
    htotal_err->GetXaxis()->SetLabelOffset(0.02);
    htotal_err->GetXaxis()->SetTitle("#font[22]{M_{inv}^{#pi^{+}#pi^{-}} (GeV/c^{2})}");
    htotal_err->GetXaxis()->SetTitleOffset(1.1);
    htotal_err->GetXaxis()->SetTitleSize(0.11);
    htotal_err->GetXaxis()->CenterTitle();
    htotal_err->GetYaxis()->CenterTitle();
    htotal_err->GetXaxis()->SetLabelSize(0.11);
    htotal_err->GetXaxis()->SetTickLength(0.08);
    htotal_err->GetYaxis()->SetNdivisions(505);
    htotal_err->SetFillColorAlpha(28, 1.0);
    htotal_err->SetFillStyle(1000);
    htotal_err->Draw("E2");

    pad2->Update();
    TLine *l1 = new TLine(pad2->GetUxmin(), 0.0, pad2->GetUxmax(), 0.0);
    l1->SetLineWidth(2);
    l1->SetLineStyle(2);
    l1->Draw();
    pad2->Update();
    /*
 hsys_err->SetFillColorAlpha(8, 1);
 hsys_err->SetFillStyle(3144);
 hsys_err->Draw("E2 SAME");
 */
    // hstat_err->SetFillColorAlpha(2, 1);
    hstat_err->SetLineColor(2);
    hstat_err->Draw("E1 SAME");

    // draw relative difference between pythia and measured
    TH1D *hrdiff = (TH1D *)hpyth->Clone();
    hrdiff->Add(hcom, -1);
    hrdiff->Divide(hcom);
    hrdiff->SetLineColor(1);
    hrdiff->SetLineWidth(2);
    hrdiff->SetLineStyle(1);
    hrdiff->Draw("same E1");

    TLegend *lg1 = new TLegend(0.3, 0.88, 0.8, 0.95);
    lg1->SetNColumns(4);
    lg1->AddEntry(htotal_err, "#font[22]{ #delta#sigma_{stat #oplus sys}}", "f");
    // lg1->AddEntry(hsys_err, "#font[22]{ #sigma_{sys}}", "f");
    lg1->AddEntry(hstat_err, "#font[22]{ #delta#sigma_{stat}}", "lep");
    lg1->AddEntry(hrdiff, "#font[22]{ Rel. diff. w/ Pythia}", "lep");
    lg1->SetTextSize(0.09);
    lg1->Draw();
    gPad->Update();

    can->Update();
}
void setBinError(TH1D *hRatio, TH1D *hData, TH1D *hEmbed)
{
    for (int i = 1; i < hData->GetNbinsX(); i++)
    {
        double binError = 0;
        if (hData->GetBinContent(i) == 0 && hEmbed->GetBinContent(i) == 0)
            continue;
        binError = getBinError(hData->GetBinContent(i), hData->GetBinError(i), hEmbed->GetBinContent(i), hEmbed->GetBinError(i));
        // cout << "bin " << i << "Content = " << hData->GetBinContent(i) << " Error = " << hData->GetBinError(i) << " Embed Content = " << hEmbed->GetBinContent(i) << " Err = " << hEmbed->GetBinError(i) << endl;
        hRatio->SetBinError(i, binError);
    }
}

double getBinError(double dataBinCt, double dataBinErr, double embedBinCt, double embedBinErr)
{

    return (dataBinCt / embedBinCt) * sqrt(pow(dataBinErr / dataBinCt, 2) + pow(embedBinErr / embedBinCt, 2));
}
// drawErrors(htotal_err, hsys_err, hpid_err, heff_err, htrg_spr_err, hstat_err);
void drawErrors(TH1D *htotal, TH1D *hsys, TH1D *hpid, TH1D *heff, TH1D *htrgspr, TH1D *hstat, TH1D *hbias)
{
    TCanvas *cerr = new TCanvas("cerr", "", 500, 350);
    cerr->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    htotal->GetXaxis()->SetTitle("M_{inv}");
    htotal->GetXaxis()->SetTitleSize(0.06);
    htotal->GetYaxis()->SetTitleSize(0.06);
    htotal->GetXaxis()->SetLabelSize(0.06);
    htotal->GetYaxis()->SetLabelSize(0.06);
    htotal->GetYaxis()->SetTitle("Uncertainty");
    htotal->GetYaxis()->SetRangeUser(-0.3, 0.3);
    htotal->SetFillColorAlpha(42, 1);
    htotal->Draw("E2");

    TLegend *leg = new TLegend(0.2, 0.85, 0.9, 0.9);
    leg->SetNColumns(7);
    leg->AddEntry(htotal, " total", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total.pdf");

    hsys->SetLineColor(2);
    hsys->Draw("E1 SAME");
    leg->AddEntry(hsys, " comb.sys.", "lep");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total_sys.pdf");

    hpid->SetFillStyle(3144);
    hpid->SetFillColorAlpha(6, 0.5);
    hpid->Draw("E2 SAME");
    leg->AddEntry(hpid, " pid", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total_sys_pid.pdf");
    /*
    heff->SetFillColorAlpha(4, 1);
    heff->Draw("E2 SAME");
    leg->AddEntry(heff, " eff.", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total_sys_pid_eff.pdf");
    */
    htrgspr->SetFillColorAlpha(3, 1);
    htrgspr->Draw("E2 SAME");
    leg->AddEntry(htrgspr, " trg.spr", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total_sys_pid_trgspr.pdf");

    hstat->SetFillColorAlpha(2, 1);
    hstat->Draw("E2 SAME");
    leg->AddEntry(hstat, " stat", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total_sys_pid_trgspr_stat.pdf");

    hbias->SetLineColor(8);
    hbias->Draw("E SAME");
    leg->AddEntry(hbias, " bias", "lep");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_total_sys_pid_trgspr_stat_bias.pdf");

    leg->Draw();
    cerr->Update();
    cerr->SaveAs("Results/sys_err_all.pdf");
    delete cerr;
}

void printArrayValues(int size, double val[], const char *valName)
{
    cout << "double " << valName << "[" << size << "] = { ";
    for (int i = 0; i < size; i++)
    {
        if (i == size - 1)
            cout << val[i] << "};" << endl;
        else
            cout << val[i] << ",";
    }
}
