#include <iostream>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;
void drawXsec(TH1D *hjp2, TH1D *hjp1, TH1D *hjp0, TCanvas *can, const char *cname);
void setBinError(TH1D *hRatio, TH1D *hData, TH1D *hEmbed);
double getBinError(double dataBinCt, double dataBinErr, double embedBinCt, double embedBinErr);
void xsecAna()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    // gStyle->SetLegendBorderSize(0);

    const int nBins = 13;
    Double_t nBinsEdges[nBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};
    // binning for cross-section money plot. Add bins at underflow and overflow
    const int nBinsU = 15;
    Double_t nBinsEdgesU[nBinsU + 1] = {0.05, 0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0, 4.25};

    // trigger prescales
    Double_t presJP0 = 0.00707; //(1. / 141.276);
    Double_t presJP1 = 0.3892;  //(1. / 2.56938);
    Double_t presJP2 = 1.0;

    TFile *fin = new TFile("../BinResolutionAndUnfolding/UnfoldingResults/unfoldingOutput.root", "r");
    // unfolded in |eta|<1
    TH1D *hun_jp0 = (TH1D *)fin->Get("UnfoldedJP0");
    TH1D *hun_jp1 = (TH1D *)fin->Get("UnfoldedJP1");
    TH1D *hun_jp2 = (TH1D *)fin->Get("UnfoldedJP2");

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
    // new values from embed trees that has correct no. of events and good data-embed match includeing more end runs
    Double_t pairEffJP0[13] = {0.755596, 0.745754, 0.757203, 0.759054, 0.775456, 0.765619, 0.794246, 0.788862, 0.787895, 0.828457, 0.782415, 0.691203, 0.635486};

    Double_t pairEffErrJP0[13] = {0.018504, 0.012761, 0.0109636, 0.0132696, 0.0130157, 0.0190241, 0.0226193, 0.0251595, 0.0389486, 0.0520254, 0.0587974, 0.0663692, 0.0928496};

    Double_t pairEffJP1[13] = {0.754459, 0.740175, 0.748992, 0.758559, 0.775749, 0.781586, 0.771335, 0.794878, 0.833484, 0.834591, 0.805139, 0.785028, 0.645855};

    Double_t pairEffErrJP1[13] = {0.018504, 0.012761, 0.0109636, 0.0132696, 0.0130157, 0.0190241, 0.0226193, 0.0251595, 0.0389486, 0.0520254, 0.0587974, 0.0663692, 0.0928496};

    Double_t pairEffJP2[13] = {0.743954, 0.728791, 0.747529, 0.75281, 0.770241, 0.784688, 0.743925, 0.807385, 0.822877, 0.813964, 0.812033, 0.81501, 0.705472};

    Double_t pairEffErrJP2[13] = {0.0256901, 0.0131463, 0.0115316, 0.00946569, 0.0107372, 0.0134516, 0.0398716, 0.0207971, 0.024248, 0.0524755, 0.0496478, 0.0758663, 0.115594};

    // dipion triggering efficiency
    // new values with good embedding runs after less event issue resolved (435 embedding runs)
    /*// jp0 efficiency and error
    Double_t eff_trg_jp0[13] = {0.0931162, 0.0943369, 0.103936, 0.131643, 0.158628, 0.2184, 0.272021, 0.316554, 0.351829, 0.34123, 0.591596, 0.984124, 1.13258};
    Double_t efferr_trg_jp0[13] = {0.000550693, 0.000405033, 0.000385203, 0.000687935, 0.00109193, 0.0025018, 0.00490356, 0.0074884, 0.0144634, 0.0211906, 0.0355467, 0.0717775, 0.13699};
    // jp1 efficiency and error
    Double_t eff_trg_jp1[13] = {0.0109936, 0.0114659, 0.0129404, 0.0193034, 0.0284457, 0.0499181, 0.0714804, 0.10004, 0.128901, 0.157633, 0.30484, 0.652834, 0.813984};
    Double_t efferr_trg_jp1[13] = {0.00011011, 8.25671e-05, 7.62766e-05, 0.000140926, 0.00022971, 0.000593137, 0.00117423, 0.00211926, 0.0043121, 0.00974095, 0.0159981, 0.046846, 0.080441};
    // jp2 efficiency and error
    Double_t eff_trg_jp2[13] = {0.00208857, 0.00224209, 0.00263264, 0.00433279, 0.00750271, 0.0145892, 0.0244818, 0.037786, 0.0523582, 0.0735169, 0.156924, 0.380886, 0.582972};
    Double_t efferr_trg_jp2[13] = {2.59113e-05, 1.90692e-05, 1.77113e-05, 2.9782e-05, 7.69039e-05, 0.000139034, 0.000376674, 0.000759807, 0.00163107, 0.00460581, 0.00808878, 0.0273548, 0.0585053};
   */
    // jp0 efficiency and error
    Double_t eff_trg_jp0[13] = {0.0923402, 0.0942138, 0.103686, 0.131355, 0.15839, 0.218144, 0.272181, 0.316211, 0.35184, 0.513864, 0.588915, 0.985694, 1.12947};
    Double_t efferr_trg_jp0[13] = {0.000547047, 0.000403385, 0.000382652, 0.000683076, 0.00108522, 0.0024899, 0.0048747, 0.00742807, 0.0143286, 0.0157066, 0.0351342, 0.0713016, 0.135711};
    // jp1 efficiency and error
    Double_t eff_trg_jp1[13] = {0.0108946, 0.0114603, 0.0129201, 0.0192628, 0.0284066, 0.0499675, 0.0715268, 0.0999069, 0.128872, 0.237605, 0.304006, 0.657327, 0.815849};
    Double_t efferr_trg_jp1[13] = {0.000109331, 8.20771e-05, 7.58507e-05, 0.000140109, 0.00022887, 0.00060527, 0.00116741, 0.00210084, 0.00427554, 0.00711265, 0.0158261, 0.0467808, 0.0801746};
    // jp2 efficiency and error
    Double_t eff_trg_jp2[13] = {0.00206183, 0.00224183, 0.00262727, 0.00431819, 0.00748315, 0.0145654, 0.0244884, 0.0377554, 0.0522194, 0.110814, 0.156778, 0.384113, 0.58051};
    Double_t efferr_trg_jp2[13] = {2.56258e-05, 1.9039e-05, 1.75471e-05, 2.94599e-05, 7.61407e-05, 0.000137931, 0.000374207, 0.000754832, 0.00161248, 0.00350601, 0.00802174, 0.0273927, 0.0578907};

    // cross-section factor (1/(eff_trk x eff_trg x L))
    Double_t f_jp0[13] = {0};
    Double_t f_jp1[13] = {0};
    Double_t f_jp2[13] = {0};
    for (int i = 0; i < 13; i++)
    {
        // f_jp0[i] = (Double_t)(trkPairEff * eff_trg_jp0[i] * L_jp0);
        // f_jp1[i] = (Double_t)(trkPairEff * eff_trg_jp1[i] * L_jp1);
        // f_jp2[i] = (Double_t)(trkPairEff * eff_trg_jp2[i] * L_jp2);
        f_jp0[i] = (Double_t)(pairEffJP0[i] * eff_trg_jp0[i] * L_jp0);
        f_jp1[i] = (Double_t)(pairEffJP1[i] * eff_trg_jp1[i] * L_jp1);
        f_jp2[i] = (Double_t)(pairEffJP2[i] * eff_trg_jp2[i] * L_jp2);
        cout << "Bin = " << i << " fjp0 = " << f_jp0[i] << " fjp1 = " << f_jp1[i] << " fjp2 = " << f_jp2[i] << endl;
    }
    // yields per bin (dN)
    Double_t dN_jp0[nBins] = {0};
    Double_t stat_err_jp0[nBins] = {0};
    Double_t total_stat_err_jp0[nBins] = {0}; // errors added in quadrature
    Double_t norm_stat_err_jp0[nBins] = {0};  // total error normalized by cross section yields
    Double_t xsec_jp0[nBins] = {0};
    Double_t dN_jp1[nBins] = {0};
    Double_t stat_err_jp1[nBins] = {0};
    Double_t total_stat_err_jp1[nBins] = {0};
    Double_t norm_stat_err_jp1[nBins] = {0};
    Double_t xsec_jp1[nBins] = {0};
    Double_t dN_jp2[nBins] = {0};
    Double_t stat_err_jp2[nBins] = {0};
    Double_t total_stat_err_jp2[nBins] = {0};
    Double_t norm_stat_err_jp2[nBins] = {0};
    Double_t xsec_jp2[nBins] = {0};

    for (int i = 0; i < nBins; i++)
    {
        if (hun_jp0->GetBinContent(i + 1) == 0)
            continue;
        dN_jp0[i] = hun_jp0->GetBinContent(i + 1);
        stat_err_jp0[i] = hun_jp0->GetBinError(i + 1);
        xsec_jp0[i] = (Double_t)dN_jp0[i] / f_jp0[i];
        total_stat_err_jp0[i] = (Double_t)sqrt(pow(stat_err_jp0[i], 2) + pow(efferr_trg_jp0[i], 2) + pow(pairEffErrJP0[i], 2));
        norm_stat_err_jp0[i] = (Double_t)(total_stat_err_jp0[i] / xsec_jp0[i]);

        dN_jp1[i] = hun_jp1->GetBinContent(i + 1);
        stat_err_jp1[i] = hun_jp1->GetBinError(i + 1);
        xsec_jp1[i] = (Double_t)dN_jp1[i] / f_jp1[i];
        total_stat_err_jp1[i] = (Double_t)sqrt(pow(stat_err_jp1[i], 2) + pow(efferr_trg_jp1[i], 2) + pow(pairEffErrJP1[i], 2));
        norm_stat_err_jp1[i] = (Double_t)(total_stat_err_jp1[i] / xsec_jp1[i]);

        dN_jp2[i] = hun_jp2->GetBinContent(i + 1);
        stat_err_jp2[i] = hun_jp2->GetBinError(i + 1);
        xsec_jp2[i] = (Double_t)dN_jp2[i] / f_jp2[i];
        total_stat_err_jp2[i] = (Double_t)sqrt(pow(stat_err_jp2[i], 2) + pow(efferr_trg_jp2[i], 2) + pow(pairEffErrJP2[i], 2));
        norm_stat_err_jp2[i] = (Double_t)(total_stat_err_jp2[i] / xsec_jp2[i]);
    }

    // book histograms for money plot
    TH1D *hxsec_jp0 = new TH1D("hxsec_jp0", "", nBinsU, nBinsEdgesU);
    TH1D *hxsec_jp1 = new TH1D("hxsec_jp1", "", nBinsU, nBinsEdgesU);
    TH1D *hxsec_jp2 = new TH1D("hxsec_jp2", "", nBinsU, nBinsEdgesU);
    for (int i = 0; i < nBinsU; i++)
    {
        if (i < 1 || i > nBins)
        {
            continue;
        }
        hxsec_jp0->SetBinContent(i + 1, xsec_jp0[i - 1]);
        hxsec_jp0->SetBinError(i + 1, total_stat_err_jp0[i - 1]);
        hxsec_jp1->SetBinContent(i + 1, xsec_jp1[i - 1]);
        hxsec_jp1->SetBinError(i + 1, total_stat_err_jp1[i - 1]);
        hxsec_jp2->SetBinContent(i + 1, xsec_jp2[i - 1]);
        hxsec_jp2->SetBinError(i + 1, total_stat_err_jp2[i - 1]);

        // cout << "Bin = " << i << " Yields = " << dN_jp0[i - 1] << endl;
    }

    // for combined triggers
    TH1D *hxsecCom = (TH1D *)hxsec_jp2->Clone();
    hxsecCom->Add(hxsec_jp1, presJP1);
    hxsecCom->Add(hxsec_jp0, presJP0);
    //----
    const char *cname = "can_Xsec";
    TCanvas *cxsec;
    drawXsec(hxsec_jp2, hxsec_jp1, hxsec_jp0, cxsec, cname);

} // main

void drawXsec(TH1D *hjp2, TH1D *hjp1, TH1D *hjp0, TCanvas *can, const char *cname)
{
    can = new TCanvas(Form("%s", cname), "", 0, 0, 700, 800);
    can->cd();
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogy();

    hjp2->SetMaximum(3e9);
    hjp2->GetXaxis()->SetLabelSize(0);
    hjp2->GetYaxis()->SetTitle("#frac{d#sigma^{#pi^{+}#pi^{-}}}{dM^{#pi^{+}#pi^{-}}}(= #frac{1}{L #upoint #epsilon_{trk}^{#pi^{+}} #upoint #epsilon_{trk}^{#pi^{-}} #upoint #epsilon_{trg}^{#pi^{+}#pi^{-}}} #times #frac{dN^{#pi^{+}#pi^{-}}_{true}}{dM_{inv}^{#pi^{+}#pi^{-}}})");
    hjp2->GetYaxis()->CenterTitle();
    hjp2->GetYaxis()->SetTitleOffset(1.5);
    hjp2->GetXaxis()->CenterTitle();
    hjp2->SetMarkerStyle(22);
    hjp2->SetMarkerColor(2);
    hjp2->SetLineColor(2);
    hjp2->SetLineWidth(2);
    hjp2->Draw("e1");

    hjp1->SetMarkerStyle(20);
    hjp1->SetMarkerColor(4);
    hjp1->SetLineColor(4);
    hjp1->SetLineWidth(2);
    hjp1->SetLineStyle(2);
    hjp1->Draw(" e1 SAME");

    hjp0->SetMarkerStyle(20);
    hjp0->SetMarkerColor(1);
    hjp0->SetLineColor(1);
    hjp0->SetLineWidth(2);
    hjp0->SetLineStyle(2);
    hjp0->Draw("e1 SAME");

    TLegend *lx = new TLegend(0.55, 0.75, 0.8, 0.8);
    lx->SetNColumns(3);
    lx->SetTextSize(0.04);
    lx->AddEntry(hjp2, " jp2", "lp");
    lx->AddEntry(hjp1, " jp1", "lp");
    lx->AddEntry(hjp0, " jp0", "lp");
    lx->Draw();

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.04);
    tex.DrawLatex(1.3, 1e9, "p + p #rightarrow #pi^{+}#pi^{-} + X  at #sqrt{s} = 200 GeV ");
    tex.DrawLatex(1.3, 3e8, "|#eta^{#pi^{+}#pi^{-}}| < 1");

    gPad->Update();

    can->cd();
    TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    gPad->SetGrid(1, 1);

    TH1D *h_R21 = (TH1D *)hjp2->Clone();
    h_R21->Add(hjp1, -1);
    h_R21->Divide(hjp1);
    h_R21->SetBit(TH1::kNoTitle);
    setBinError(h_R21, hjp2, hjp1);
    h_R21->SetLineColor(4);
    // h_R21->GetYaxis()->SetTitle("X_{jp2} - X_{jp1(0)}/X_{jp1(0)}");
    h_R21->GetYaxis()->SetTitle("Rel. diff w/ jp2");
    h_R21->GetXaxis()->SetTitle("M^{#pi^{+}#pi^{-}}_{inv} (GeV/c^{2})");
    h_R21->GetYaxis()->CenterTitle();
    h_R21->GetYaxis()->SetRangeUser(-.8, 0.8);
    h_R21->GetYaxis()->SetLabelSize(0.11);
    h_R21->GetYaxis()->SetTitleOffset(0.55);
    h_R21->GetYaxis()->SetTitleSize(0.11);
    h_R21->GetXaxis()->SetLabelOffset(0.02);
    h_R21->GetXaxis()->SetTitleOffset(1.1);
    h_R21->GetXaxis()->SetTitleSize(0.11);
    h_R21->GetXaxis()->SetLabelSize(0.11);
    h_R21->GetXaxis()->SetTickLength(0.08);
    h_R21->GetYaxis()->SetNdivisions(505);
    h_R21->SetMarkerStyle(8);
    h_R21->SetMarkerColor(4);
    h_R21->SetLineColor(4);
    h_R21->Draw(" e1");

    TH1D *h_R20 = (TH1D *)hjp2->Clone();
    h_R20->Add(hjp0, -1);
    h_R20->Divide(hjp0);
    h_R20->SetBit(TH1::kNoTitle);
    setBinError(h_R20, hjp2, hjp0);
    h_R20->SetLineColor(1);
    h_R20->SetLineStyle(2);
    h_R20->SetMarkerStyle(8);
    h_R20->SetMarkerColor(1);
    h_R20->Draw("same ");

    gPad->Update();

    can->Update();
    can->SaveAs(Form("./Results/%s.pdf", cname));
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
