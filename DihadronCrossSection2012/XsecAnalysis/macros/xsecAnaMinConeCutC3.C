#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iomanip>
#include <cmath>
#include "helper_functions.h"

using namespace std;

void drawXsec(TH1D *hcom, TH1D *hpyth, TH1D *hjp2, TH1D *hjp1, TH1D *hjp0);
// void drawXsec(TH1D *hcom, TH1D *hpyth, TH1D *hjp2, TH1D *hjp1, TH1D *hjp0, TCanvas *can);
void setBinError(TH1D *hRatio, TH1D *hData, TH1D *hEmbed);
double getBinError(double dataBinCt, double dataBinErr, double embedBinCt, double embedBinErr);

void drawXsecMoneyPlot(TH1D *hcom, TH1D *hpyth, TH1D *hsys_err, TH1D *hstat_err, TH1D *htotal_err, TCanvas *can);

void drawErrors(TH1D *htotal, TH1D *hsys, TH1D *hpid, TH1D *hloss, TH1D *htrgspr, TH1D *hstat, TH1D *hbias);

void xsecAnaMinConeCutC3()
{
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptDate(0);
    // gStyle->SetLabelFont(13, "xyz");
    gStyle->SetLegendBorderSize(0);

    TFile *fxout = new TFile("ResultsC3/fxsecOutputC3.root", "Recreate");
    // cross-section binning
    const int nBins = 13;
    Double_t nBinsEdges[nBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

    const int nBinsC3 = 8;
    Double_t nBinsEdgesC3[nBinsC3 + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60};

    // binning for cross-section money plot. Add bins at underflow and overflow
    const int nBinsUC3 = 10;
    Double_t nBinsEdgesUC3[nBinsUC3 + 1] = {0.05, 0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.65};

    // const int nBinsU = 15;
    // Double_t nBinsEdgesU[nBinsU + 1] = {0.05, 0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0, 4.25};

    // trigger prescales
    Double_t presjp0 = 0.00707; //(1. / 141.276);
    Double_t presjp1 = 0.3892;  //(1. / 2.56938);
    Double_t presjp2 = 1.0;
    // pythia cross section

    // partonic pt parameters
    double partPtRange[nBins + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., -1};

    Double_t crossSection[nBins] = {9.0012e+00, 1.46253e+00, 3.54466e-01, 1.51622e-01, 2.49062e-02, 5.84527e-03, 2.30158e-03, 3.42755e-04, 4.57002e-05, 9.72535e-06, 4.69889e-07, 2.69202e-08, 1.43453e-09}; // in mb

    // TFile *frawX = new TFile("PythiaXsec/pythiaXsec.root", "R"); // good embed runs
    TFile *frawX = new TFile("PythiaXsec/pythiaXsecAll.root", "R"); // all embed runs

    // get raw cross section in wach pt bins
    TH1D *hrawX[13];
    TH1D *hevProcess[13];
    double lum[13] = {0};
    double evProcessed[13] = {0};
    for (int i = 0; i < 13; i++)
    {
        hrawX[i] = (TH1D *)frawX->Get(Form("hpythiaMinv_%g_%g", partPtRange[i], partPtRange[i + 1]));
        hevProcess[i] = (TH1D *)frawX->Get(Form("hevProcessed_%g_%g", partPtRange[i], partPtRange[i + 1]));
        evProcessed[i] = hevProcess[i]->Integral();
        lum[i] = evProcessed[i] / (crossSection[i] * 1e9);
        hrawX[i]->Scale(1. / lum[i]);
    }

    TH1D *hpythiaXsec = (TH1D *)hrawX[0]->Clone();
    for (int i = 1; i < 13; i++)
    {
        hpythiaXsec->Add(hrawX[i]);
    }
    TH1D *hpythiaXsecW = (TH1D *)hpythiaXsec->Clone();
    hpythiaXsecW->Reset();

    for (int i = 1; i <= 13; i++)
    {
        double nn = hpythiaXsec->GetBinContent(i);
        double bw = hpythiaXsec->GetBinWidth(i);
        double ne = hpythiaXsec->GetBinError(i);
        hpythiaXsecW->SetBinContent(i, nn / bw);
        hpythiaXsecW->SetBinError(i, ne / bw);
    }
    //----

    TFile *funfolding = new TFile("../BinResolutionAndUnfolding/UnfoldingResultsMinConeCut/unfoldingOutputC7C3.root", "r");
    //  unfolded in |eta|<1
    TH1D *hun_jp0 = (TH1D *)funfolding->Get("UnfoldedJP0C3");
    TH1D *hun_jp1 = (TH1D *)funfolding->Get("UnfoldedJP1C3");
    TH1D *hun_jp2 = (TH1D *)funfolding->Get("UnfoldedJP2C3");

    // get particle level cross section (for comparison)
    // TFile *fparticle = new TFile("hist_XsecParticle.root", "r");
    // TH1D *hXsecParticle = (TH1D *)fparticle->Get("hUnbiasedMinv");
    TH1D *hXsecParticle = (TH1D *)hpythiaXsecW->Clone();

    // luminosity
    Double_t L_jp0 = 0.16;
    Double_t L_jp1 = 7.68;
    Double_t L_jp2 = 18.80;
    Double_t L_jptotal = 26.64;

    // tracking efficiency per bin
    // pair efficiency with the minimum cone cut
    double trkeff_jp0[13] = {0.888295, 0.890521, 0.868916, 0.880734, 0.888951, 0.889116, 0.895198, 0.905515, 0.909525, 0.934664, 0.914252, 0.920243, 0.888233};
    double trkefferr_jp0[13] = {1.10635e-07, 8.02667e-08, 7.39238e-08, 1.12719e-07, 1.53353e-07, 2.89457e-07, 4.82991e-07, 7.38947e-07, 1.28912e-06, 2.0917e-06, 4.61979e-06, 1.12867e-05, 3.49742e-05};
    double trkeff_jp1[13] = {0.884958, 0.884586, 0.862004, 0.882875, 0.89147, 0.89712, 0.900171, 0.900808, 0.919683, 0.925122, 0.92678, 0.921612, 0.871139};
    double trkefferr_jp1[13] = {3.25733e-07, 2.3539e-07, 2.09923e-07, 2.88629e-07, 3.53001e-07, 5.84872e-07, 9.34296e-07, 1.33235e-06, 1.95069e-06, 3.35083e-06, 6.22082e-06, 1.40152e-05, 4.35365e-05};
    double trkeff_jp2[13] = {0.876053, 0.881871, 0.854031, 0.879692, 0.883777, 0.892923, 0.900998, 0.909497, 0.912793, 0.911193, 0.931605, 0.918144, 0.86194};
    double trkefferr_jp2[13] = {7.47627e-07, 5.25193e-07, 4.69177e-07, 5.99939e-07, 6.99555e-07, 1.08637e-06, 1.58333e-06, 2.12297e-06, 3.23232e-06, 5.37244e-06, 8.66935e-06, 1.81999e-05, 5.42626e-05};

    // dipion triggering efficiency
    /*
    // jp0 efficiency and error
    Double_t eff_trg_jp0[13] = {0.0647609, 0.0663649, 0.0734078, 0.0840551, 0.109116, 0.150331, 0.197242, 0.261445, 0.408169, 0.447309, 0.501666, 0.668775, 0.759649};
    Double_t efferr_trg_jp0[13] = {0.000446209, 0.000330274, 0.000310635, 0.000505293, 0.000819134, 0.00188765, 0.00383486, 0.00661998, 0.00748233, 0.0149041, 0.0275889, 0.0519418, 0.0985095};
    // jp1 efficiency and error
    Double_t eff_trg_jp1[13] = {0.00762761, 0.00806872, 0.00943897, 0.0127224, 0.0203221, 0.0349688, 0.0512922, 0.0836607, 0.162998, 0.195088, 0.246034, 0.425724, 0.533614};
    Double_t efferr_trg_jp1[13] = {8.31865e-05, 6.43366e-05, 6.41374e-05, 0.000102156, 0.000170519, 0.00041934, 0.00085, 0.00192461, 0.00316734, 0.00682422, 0.0121178, 0.0320969, 0.0574833};
    // jp2 efficiency and error
    Double_t eff_trg_jp2[13] = {0.00151809, 0.00163691, 0.0019597, 0.0029936, 0.00543582, 0.0104283, 0.0176766, 0.0307442, 0.0640216, 0.0869714, 0.119685, 0.258225, 0.361955};
    Double_t efferr_trg_jp2[13] = {2.11924e-05, 1.42757e-05, 1.59223e-05, 1.75488e-05, 3.76751e-05, 0.000108269, 0.000299922, 0.00065822, 0.00114225, 0.00334942, 0.00557087, 0.0197869, 0.0413919};
    */
    // jp0 efficiency and error
    Double_t eff_trg_jp0[13] = {0.06218, 0.0635835, 0.0666106, 0.0842533, 0.115018, 0.162886, 0.227979, 0.318971, 0.500256, 0.5584, 0.638504, 0.760584, 0.901282};
    Double_t efferr_trg_jp0[13] = {0.000358476, 0.000264799, 0.000237964, 0.000419787, 0.000761594, 0.00177198, 0.00397223, 0.00726696, 0.00820762, 0.0168248, 0.0306366, 0.052663, 0.103569};
    // jp1 efficiency and error
    Double_t eff_trg_jp1[13] = {0.00731371, 0.00773927, 0.00830041, 0.0123315, 0.0206927, 0.0372235, 0.0600657, 0.103605, 0.211771, 0.257325, 0.347275, 0.507202, 0.641553};
    Double_t efferr_trg_jp1[13] = {7.22544e-05, 5.52032e-05, 4.79204e-05, 8.71347e-05, 0.000165822, 0.000433534, 0.000934156, 0.00208306, 0.00359666, 0.00750071, 0.0155309, 0.0346532, 0.0605443};
    // jp2 efficiency and error
    Double_t eff_trg_jp2[13] = {0.00138597, 0.00151376, 0.00168918, 0.00276861, 0.0054136, 0.010864, 0.0205764, 0.0396581, 0.0873528, 0.121622, 0.174098, 0.300651, 0.438564};
    Double_t efferr_trg_jp2[13] = {1.70826e-05, 1.3147e-05, 1.11464e-05, 1.6876e-05, 5.42097e-05, 9.61059e-05, 0.000298852, 0.000770725, 0.00137986, 0.00379984, 0.00705385, 0.0203718, 0.0420877};

    // trigger  bias
    double sys_trigBias[13] = {0.0217904, 0.020958, 0.0230904, 0.0138005, 0.00169444, 0.00280884, 0.00213932, 0.00525629, 0.0149015, 0.0155368, 0.0102302, 0.00550931, 0.033424};

    //  signal pid fraction (apply as a correction)
    // true fractions -----------------
    double fraction_truejp2[13] = {0.84499, 0.842992, 0.838184, 0.805859, 0.833143, 0.837585, 0.847277, 0.859741, 0.860903, 0.876257, 0.874977, 0.878856, 0.845706};
    double fraction_truejp1[13] = {0.832808, 0.831881, 0.821748, 0.765399, 0.808291, 0.810806, 0.825833, 0.857797, 0.867719, 0.886949, 0.88463, 0.867285, 0.812604};
    double fraction_truejp0[13] = {0.838803, 0.835594, 0.814567, 0.742295, 0.778612, 0.784243, 0.812714, 0.858089, 0.877497, 0.901201, 0.888658, 0.863375, 0.846732};
    // Missed fractions-----------------
    double fraction_missedjp2[13] = {1.38673, 1.39851, 1.3703, 1.35543, 1.33698, 1.32011, 1.30906, 1.29129, 1.29228, 1.29154, 1.31678, 1.2389, 1.25844};
    double fraction_missedjp1[13] = {1.38834, 1.3936, 1.38832, 1.3708, 1.33332, 1.33016, 1.32992, 1.31654, 1.29919, 1.29851, 1.26968, 1.2463, 1.28287};
    double fraction_missedjp0[13] = {1.40607, 1.40268, 1.39037, 1.37083, 1.33549, 1.31754, 1.30873, 1.29688, 1.29907, 1.27603, 1.2869, 1.25255, 1.27453};

    // cross-section factor
    Double_t binWidth[13] = {0};
    Double_t f_jp0[13] = {0};
    Double_t f_jp1[13] = {0};
    Double_t f_jp2[13] = {0};

    for (int i = 0; i < nBinsC3; i++)
    { // factor that should be multiplied with the true yields
        binWidth[i] = hun_jp0->GetBinWidth(i + 1);
        f_jp0[i] = (Double_t)((fraction_truejp0[i] * fraction_missedjp0[i]) / (binWidth[i] * trkeff_jp0[i] * eff_trg_jp0[i] * L_jp0));
        f_jp1[i] = (Double_t)((fraction_truejp1[i] * fraction_missedjp1[i]) / (binWidth[i] * trkeff_jp1[i] * eff_trg_jp1[i] * L_jp1));
        f_jp2[i] = (Double_t)((fraction_truejp2[i] * fraction_missedjp2[i]) / (binWidth[i] * trkeff_jp2[i] * eff_trg_jp2[i] * L_jp2));
    }
    // yields per bin (dN)
    Double_t dN_jp0[nBinsC3] = {0};            // unfolded dipion raw yields
    Double_t stat_err_jp0[nBinsC3] = {0};      // stat error from unfolding
    Double_t xsec_stat_err_jp0[nBinsC3] = {0}; // stat error from unfolding
    Double_t norm_stat_err_jp0[nBinsC3] = {0}; // total error normalized by cross section yields
    Double_t xsec_jp0[nBins] = {0};            // cross section
    Double_t dN_jp1[nBins] = {0};
    Double_t stat_err_jp1[nBinsC3] = {0};
    Double_t xsec_stat_err_jp1[nBinsC3] = {0};
    Double_t norm_stat_err_jp1[nBinsC3] = {0};
    Double_t xsec_jp1[nBinsC3] = {0};
    Double_t dN_jp2[nBinsC3] = {0};
    Double_t stat_err_jp2[nBinsC3] = {0};
    Double_t xsec_stat_err_jp2[nBinsC3] = {0};
    Double_t norm_stat_err_jp2[nBinsC3] = {0};
    Double_t xsec_jp2[nBinsC3] = {0};

    // xsec for combined trigger
    Double_t dN_jp[nBinsC3] = {0};
    Double_t stat_err_jp[nBinsC3] = {0};
    Double_t xsec_stat_err_jp[nBinsC3] = {0};
    Double_t total_stat_err_jp[nBinsC3] = {0};
    Double_t norm_stat_err_jp[nBinsC3] = {0};
    Double_t xsec_jp[nBinsC3] = {0};
    // systematic error from trigger spread
    Double_t sys_trg_err[nBinsC3] = {0};
    double sys_loss_jp[nBinsC3] = {0};
    // systematic errors for combinatorial background correction
    // f_fake  and f_loss systematic uncertainties

    Double_t sys_pid_err[nBinsC3] = {0};

    double sys_pid_jp0[13] = {0.0319655, 0.0380019, 0.0381703, 0.0625011, 0.0790663, 0.103529, 0.106695, 0.102389, 0.087327, 0.0863025, 0.0847956, 0.0899803, 0.0622231};
    double sys_pid_jp1[13] = {0.0431741, 0.0416838, 0.0420531, 0.063713, 0.0772958, 0.0827711, 0.110822, 0.0956831, 0.0950896, 0.0811664, 0.106994, 0.108509, 0.0468354};
    double sys_pid_jp2[13] = {0.0355776, 0.0318222, 0.0465825, 0.0663322, 0.0641217, 0.0763546, 0.0851183, 0.0913476, 0.0994355, 0.0767705, 0.105245, 0.098161, 0.0938268};

    // double sys_pidloss[13] = {0.00994507, 0.00698479, 0.00580664, 0.00561845, 0.00350876, 0.00667169, 0.0078051, 0.00315009, 0.00317296, 0.00157208, 0.0114905, 0.0199317, 0.0665472};

    double sys_loss_jp0[13] = {0.0685003, 0.0630261, 0.0639657, 0.0616306, 0.060329, 0.060322, 0.0567175, 0.0546299, 0.0613261, 0.0690677, 0.067528, 0.0683821, 0.0870841};
    double sys_loss_jp1[13] = {0.0610002, 0.0657847, 0.0642313, 0.0592798, 0.0583507, 0.0643571, 0.0648272, 0.0670913, 0.0580218, 0.0627106, 0.0542626, 0.0644079, 0.0817957};
    double sys_loss_jp2[13] = {0.0656017, 0.0649353, 0.0598546, 0.0623264, 0.05982, 0.061688, 0.0620466, 0.0585203, 0.0659548, 0.0589355, 0.0658444, 0.0651514, 0.0632301};
    // total errors
    Double_t sys_err[nBinsC3] = {0};
    Double_t norm_stat_err[nBinsC3] = {0};

    // pythia cross ection
    Double_t xsec_pyth[13] = {0};
    Double_t stat_err_pyth[13] = {0};

    int ww = 15;
    // cout << "total" << "\t" << "jp2" << "\t" << "jp1" << "\t" << "jp0" << endl;

    ofstream fxsec0;
    fxsec0.open("ResultsC3/fxsec_jp0_numbers.tex");
    fxsec0 << "\%Cross section results for jp0 " << endl;

    ofstream fxsec1;
    fxsec1.open("ResultsC3/fxsec_jp1_numbers.tex");
    fxsec1 << "\%Cross section results for jp1" << endl;

    ofstream fxsec2;
    fxsec2.open("ResultsC3/fxsec_jp2_numbers.tex");
    fxsec2 << "\%Cross section results for jp2" << endl;

    ofstream fxsec;
    fxsec.open("ResultsC3/fxsec_values.tex");
    fxsec << "% 1. Bin Range 2. Xsec 3. StatErr 4. Xsec_Pyth 5. Pyth_Err 6. TrgEffErr  7. SysPID 8. Sys_TrgBias 9. Sys_Total " << endl;

    ofstream fsyserr;
    fsyserr.open("ResultsC3/syserr_numbers.tex");
    fsyserr << "%\ total systematic errors by numbers " << endl;

    for (int i = 0; i < nBinsC3; i++)
    {
        if (hun_jp0->GetBinContent(i + 1) == 0)
            continue;

        // jp0
        dN_jp0[i] = hun_jp0->GetBinContent(i + 1);
        stat_err_jp0[i] = hun_jp0->GetBinError(i + 1);

        xsec_jp0[i] = (Double_t)(dN_jp0[i] * f_jp0[i]);
        xsec_stat_err_jp0[i] = stat_err_jp0[i] * f_jp0[i];
        norm_stat_err_jp0[i] = (Double_t)(xsec_stat_err_jp0[i] / xsec_jp0[i]);
        // jp1
        dN_jp1[i] = hun_jp1->GetBinContent(i + 1);
        stat_err_jp1[i] = hun_jp1->GetBinError(i + 1);
        xsec_jp1[i] = (Double_t)(dN_jp1[i] * f_jp1[i]);
        xsec_stat_err_jp1[i] = stat_err_jp1[i] * f_jp1[i];
        norm_stat_err_jp1[i] = (Double_t)(xsec_stat_err_jp1[i] / xsec_jp1[i]);
        // jp2
        dN_jp2[i] = hun_jp2->GetBinContent(i + 1);
        stat_err_jp2[i] = hun_jp2->GetBinError(i + 1);
        xsec_jp2[i] = (Double_t)(dN_jp2[i] * f_jp2[i]);
        xsec_stat_err_jp2[i] = stat_err_jp2[i] * f_jp2[i];
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
        xsec_stat_err_jp[i] = 1. / sqrt(w2_jp0 + w2_jp1 + w2_jp2);
        norm_stat_err[i] = xsec_stat_err_jp[i] / xsec_jp[i];

        // xsec pythia
        xsec_pyth[i] = (hXsecParticle->GetBinContent(i + 1));
        stat_err_pyth[i] = hXsecParticle->GetBinError(i + 1);

        // find a trigger spread relative to the combined cross section and apply as systematic error
        double spr_jp0 = fabs(xsec_jp[i] - xsec_jp0[i]) / xsec_jp[i];
        double spr_jp1 = fabs(xsec_jp[i] - xsec_jp1[i]) / xsec_jp[i];
        double spr_jp2 = fabs(xsec_jp[i] - xsec_jp2[i]) / xsec_jp[i];
        if (i >= 10)
        {
            spr_jp0 = 0.03;
            spr_jp1 = 0.02;
            spr_jp2 = 0.01;
        }
        // sys_trg_err[i] = (w2_jp0 * spr_jp0 + w2_jp1 * spr_jp1 + w2_jp2 * spr_jp2) / (w2_jp0 + w2_jp1 + w2_jp2); // weighted average

        sys_trg_err[i] = TMath::Max(spr_jp0, TMath::Max(spr_jp1, spr_jp2)); // max offsed of all three triggers is assigned as systematic error

        sys_pid_err[i] = (w2_jp0 * sys_pid_jp0[i] + w2_jp1 * sys_pid_jp1[i] + w2_jp2 * sys_pid_jp0[i]) / ((w2_jp0 + w2_jp1 + w2_jp2));

        sys_loss_jp[i] = (w2_jp0 * sys_loss_jp0[i] + w2_jp1 * sys_loss_jp1[i] + w2_jp2 * sys_loss_jp0[i]) / ((w2_jp0 + w2_jp1 + w2_jp2));

        // total systematic error
        sys_err[i] = sqrt(pow(sys_trigBias[i], 2) + pow(sys_trg_err[i], 2) + pow(sys_pid_err[i], 2) + pow(sys_loss_jp[i], 2));

        // systematic errors
        if (i == 0)
            fsyserr << "% 1. Bin No. 2. trigeff 3. f_fake 4. f_loss 5. trigger bias  6. sys_total" << endl;

        fsyserr << i << "&" << setprecision(3) << sys_trg_err[i] << "&" << setprecision(3) << sys_pid_err[i] << "&" << setprecision(3) << sys_loss_jp[i] << "&" << setprecision(3) << sys_trigBias[i] << "&" << setprecision(3) << sys_err[i] << "\\\\" << endl;

        // jp0 numbers
        if (i == 0)
            fxsec0 << "% 1. Bin No.2. BW 3. L 4. Etrg 5. Etrk 6.ff 7. fl 8. N_{un} 9. #delta_{un} 9. #sigma_{jp0} 10. #delta_{stat,jp0} " << endl;

        fxsec0 << i << "&" << binWidth[i] << "&" << L_jp0 << "&" << setprecision(3) << eff_trg_jp0[i] << "&" << setprecision(3) << trkeff_jp0[i] << "&" << setprecision(3) << fraction_truejp0[i] << "&" << setprecision(3) << fraction_missedjp0[i] << "&" << setprecision(3) << dN_jp0[i] << "&" << setprecision(3) << stat_err_jp0[i] << "&" << setprecision(3) << xsec_jp0[i] << "&" << setprecision(3) << xsec_stat_err_jp0[i] << "\\\\" << endl;
        // jp1 numbers
        if (i == 0)
            fxsec1 << "% 1. Bin No.2. BW 3. L 4. Etrg 5. Etrk 6.ff 7. fl 8. N_{un} 9. #delta_{un} 9. #sigma_{jp0} 10. #delta_{stat,jp0} " << endl;

        fxsec1 << i << "&" << binWidth[i] << "&" << L_jp1 << "&" << setprecision(3) << eff_trg_jp1[i] << "&" << setprecision(3) << trkeff_jp1[i] << "&" << setprecision(3) << fraction_truejp1[i] << "&" << setprecision(3) << fraction_missedjp1[i] << "&" << setprecision(3) << dN_jp1[i] << "&" << setprecision(3) << stat_err_jp1[i] << "&" << setprecision(3) << xsec_jp1[i] << "&" << setprecision(3) << xsec_stat_err_jp1[i] << "\\\\" << endl;
        // jp2 numbers
        if (i == 0)
            fxsec2 << "% 1. Bin No.2. BW 3. L 4. Etrg 5. Etrk 6.ff 7. fl 8. N_{un} 9. #delta_{un} 9. #sigma_{jp0} 10. #delta_{stat,jp0} " << endl;

        fxsec2 << i << "&" << binWidth[i] << "&" << L_jp2 << "&" << setprecision(3) << eff_trg_jp2[i] << "&" << setprecision(3) << trkeff_jp2[i] << "&" << setprecision(3) << fraction_truejp2[i] << "&" << setprecision(3) << fraction_missedjp2[i] << "&" << setprecision(3) << dN_jp2[i] << "&" << setprecision(3) << stat_err_jp2[i] << "&" << setprecision(3) << xsec_jp2[i] << "&" << setprecision(3) << xsec_stat_err_jp2[i] << "\\\\" << endl;

        // total cross section
        fxsec << nBinsEdges[i] << "-" << nBinsEdges[i + 1] << "&" << setprecision(3) << xsec_jp[i] << "&" << setprecision(3) << xsec_stat_err_jp[i] << "&" << setprecision(3) << xsec_pyth[i] << "&" << setprecision(3) << stat_err_pyth[i] << "&" << setprecision(3) << sys_err[i] << "\\\\" << endl;
    }

    // book histograms for money plot
    TH1D *hxsec_jp0 = new TH1D("hxsec_jp0", "", nBinsUC3, nBinsEdgesUC3);
    TH1D *hxsec_jp1 = new TH1D("hxsec_jp1", "", nBinsUC3, nBinsEdgesUC3);
    TH1D *hxsec_jp2 = new TH1D("hxsec_jp2", "", nBinsUC3, nBinsEdgesUC3);
    TH1D *hxsec_pyth = new TH1D("hxsec_pyth", "", nBinsUC3, nBinsEdgesUC3);
    TH1D *hxsec_comb = new TH1D("hxsec_comb", "", nBinsUC3, nBinsEdgesUC3);

    // for systematic error
    TH1D *htrg_spr_err = new TH1D("htrg_spr_err", "", nBinsUC3, nBinsEdgesUC3); // trig spread and efficiencies
    TH1D *hpid_err = new TH1D("hpid_err", "", nBinsUC3, nBinsEdgesUC3);         // efficiencies error
    TH1D *hloss_err = new TH1D("hloss_err", "", nBinsUC3, nBinsEdgesUC3);       // efficiencies error
    TH1D *hbias_err = new TH1D("hbias_err", "", nBinsUC3, nBinsEdgesUC3);       // efficiencies error
    // TH1D *hdatemb_err = new TH1D("hdatemb_err", "", nBinsUC3, nBinsEdgesUC3);   // efficiencies error
    TH1D *hsys_err = new TH1D("hsys_err", "", nBinsUC3, nBinsEdgesUC3);     // tota systematic error
    TH1D *hstat_err = new TH1D("hstat_err", "", nBinsUC3, nBinsEdgesUC3);   // xsec statistical error
    TH1D *htotal_err = new TH1D("htotal_err", "", nBinsUC3, nBinsEdgesUC3); // stat + sys combined

    for (int i = 0; i < nBinsUC3; i++)
    {
        if (i < 1 || i > nBins)
        {
            continue;
        }
        hxsec_jp0->SetBinContent(i + 1, xsec_jp0[i - 1]);
        hxsec_jp0->SetBinError(i + 1, xsec_stat_err_jp0[i - 1]);

        hxsec_jp1->SetBinContent(i + 1, xsec_jp1[i - 1]);
        hxsec_jp1->SetBinError(i + 1, xsec_stat_err_jp1[i - 1]);

        hxsec_jp2->SetBinContent(i + 1, xsec_jp2[i - 1]);
        hxsec_jp2->SetBinError(i + 1, xsec_stat_err_jp2[i - 1]);

        hxsec_comb->SetBinContent(i + 1, xsec_jp[i - 1]);
        hxsec_comb->SetBinError(i + 1, xsec_stat_err_jp[i - 1]);

        // cout << xsec_stat_err_jp[i - 1] << endl;

        hxsec_pyth->SetBinContent(i + 1, xsec_pyth[i - 1]);
        hxsec_pyth->SetBinError(i + 1, stat_err_pyth[i - 1]);

        // fill error histograms
        hpid_err->SetBinContent(i + 1, 0.0);
        hpid_err->SetBinError(i + 1, sys_pid_err[i - 1]);

        hloss_err->SetBinContent(i + 1, 0.0);
        hloss_err->SetBinError(i + 1, sys_loss_jp[i - 1]);

        hbias_err->SetBinContent(i + 1, 0.0);
        hbias_err->SetBinError(i + 1, sys_trigBias[i - 1]);

        htrg_spr_err->SetBinContent(i + 1, 0.0);
        htrg_spr_err->SetBinError(i + 1, sys_trg_err[i - 1]);

        // hdatemb_err->SetBinContent(i + 1, 0.0);
        // hdatemb_err->SetBinError(i + 1, sys_datemb[i - 1]);

        hsys_err->SetBinContent(i + 1, 0.0);
        hsys_err->SetBinError(i + 1, sys_err[i - 1]);

        hstat_err->SetBinContent(i + 1, 0.0);
        hstat_err->SetBinError(i + 1, norm_stat_err[i - 1]);

        htotal_err->SetBinContent(i + 1, 0.0);
        htotal_err->SetBinError(i + 1, sqrt(pow(norm_stat_err[i - 1], 2) + pow(sys_err[i - 1], 2)));
    }

    // trigger efficiency systematic studies
    // trigEff(hxsec_comb, hxsec_jp0, hxsec_jp1, hxsec_jp2);

    // draw money plot

    TCanvas *cmpt = new TCanvas("cmpt", "", 0, 0, 700, 800);
    drawXsecMoneyPlot(hxsec_comb, hxsec_pyth, hsys_err, hstat_err, htotal_err, cmpt);

    cmpt->SaveAs("ResultsC3/xsec_MoneyPlot.pdf");
    cmpt->SaveAs("ResultsC3/xsec_MoneyPlot.png");
    cmpt->SaveAs("ResultsC3/xsec_MoneyPlot.eps");
    delete cmpt;

    drawXsec(hxsec_comb, hxsec_pyth, hxsec_jp2, hxsec_jp1, hxsec_jp0);

    drawErrors(htotal_err, hsys_err, hpid_err, hloss_err, htrg_spr_err, hstat_err, hbias_err);

    fxout->cd();
    hxsec_comb->Write();
    hxsec_jp0->Write();
    hxsec_jp1->Write();
    hxsec_jp2->Write();
    hxsec_pyth->Write();
    // hdat->Write();
    // hemb->Write();

    //----

} // main

// void drawXsec(TH1D *hcom, TH1D *hpyth, TH1D *hjp2, TH1D *hjp1, TH1D *hjp0, TCanvas *can)
void drawXsec(TH1D *hcom, TH1D *hpyth, TH1D *hjp2, TH1D *hjp1, TH1D *hjp0)
{
    TCanvas *cxsec = new TCanvas("cxsec", "", 0, 0, 700, 800);
    cxsec->cd();
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLogy();

    hjp0->GetXaxis()->SetLabelSize(0);
    hjp0->GetYaxis()->SetTitle("#font[22]{#frac{d#sigma^{#pi^{+}#pi^{-}}}{dM^{#pi^{+}#pi^{-}}}(= #frac{f_{fake}*f_{loss}}{L #upoint  #epsilon_{trk}^{#pi^{+}} #upoint #epsilon_{trk}^{#pi^{-}} #upoint #epsilon_{trg}^{#pi^{+}#pi^{-}}} #times #frac{dN^{#pi^{+}#pi^{-}}_{true}}{dM_{inv}^{#pi^{+}#pi^{-}}})}");
    hjp0->GetYaxis()->CenterTitle();
    hjp0->GetYaxis()->SetTitleOffset(1.5);
    hjp0->GetXaxis()->CenterTitle();
    hjp0->SetMaximum(5e10);
    hjp0->SetMinimum(10);
    // hjp0->SetMarkerStyle(20);
    // hjp0->SetMarkerColor(1);
    hjp0->SetLineColor(1);
    hjp0->SetLineWidth(2);
    hjp0->SetLineStyle(1);
    hjp0->Draw("E1");

    /*
    hjp0->SetMarkerStyle(20);
    hjp0->SetMarkerColor(3);
    hjp0->SetLineColor(3);
    hjp0->SetLineWidth(2);
    hjp0->SetLineStyle(2);
    hjp0->Draw("E1 SAME");
    */
    hcom->SetMarkerStyle(24);
    hcom->SetMarkerColor(4);
    hcom->SetLineColor(4);
    hcom->SetLineStyle(4);
    hcom->SetLineWidth(2);
    hcom->Draw("E1 SAME");

    // hjp1->SetMarkerStyle(20);
    // hjp1->SetMarkerColor(4);
    hjp1->SetLineColor(2);
    hjp1->SetLineWidth(2);
    hjp1->SetLineStyle(2);
    hjp1->Draw("E1 SAME");

    // hjp2->SetMarkerStyle(22);
    // hjp2->SetMarkerColor(6);
    hjp2->SetLineColor(3);
    hjp2->SetLineWidth(2);
    hjp2->SetLineStyle(3);
    hjp2->Draw("E1 SAME");

    TLegend *lx = new TLegend(0.7, 0.6, 0.85, 0.8);
    // lx->SetNColumns(4);
    lx->SetTextSize(0.04);
    lx->AddEntry(hcom, " Comb.", "lep");
    // lx->AddEntry(hpyth, " Pythia", "lp");
    lx->AddEntry(hjp2, " jp2", "lp");
    lx->AddEntry(hjp1, " jp1", "lp");
    lx->AddEntry(hjp0, " jp0", "lp");
    lx->Draw();

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.04);
    tex.DrawLatex(1.3, hjp0->GetMaximum() * 0.3, "#font[22]{p + p #rightarrow #pi^{+}#pi^{-} + X  at #sqrt{s} = 200 GeV }");
    tex.DrawLatex(1.3, hjp0->GetMaximum() * 0.1, "#font[22]{|#eta^{#pi^{+}#pi^{-}}| < 1}");

    gPad->Update();

    cxsec->cd();
    TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    gPad->SetGrid(1, 1);

    TH1D *h_R2 = (TH1D *)hjp2->Clone();
    h_R2->Add(hcom, -1);
    h_R2->Divide(hcom);
    h_R2->SetBit(TH1::kNoTitle);
    // setBinError(h_R2, hcom, hjp2);

    h_R2->GetYaxis()->SetTitle("#font[22]{#frac{#sigma_{jp0,1,2} - #sigma_{comb}}{#sigma_{comb}}}");
    h_R2->GetXaxis()->SetTitle("#font[22]{M^{#pi^{+}#pi^{-}}_{inv} (GeV/c^{2})}");
    h_R2->GetYaxis()->CenterTitle();
    h_R2->GetYaxis()->SetRangeUser(-0.4, 0.4);
    h_R2->GetYaxis()->SetLabelSize(0.11);
    h_R2->GetYaxis()->SetTitleOffset(0.55);
    h_R2->GetYaxis()->SetTitleSize(0.11);
    h_R2->GetXaxis()->SetLabelOffset(0.02);
    h_R2->GetXaxis()->SetTitleOffset(1.1);
    h_R2->GetXaxis()->SetTitleSize(0.11);
    h_R2->GetXaxis()->SetLabelSize(0.11);
    h_R2->GetXaxis()->SetTickLength(0.08);
    h_R2->GetYaxis()->SetNdivisions(505);
    h_R2->SetLineColor(3);
    h_R2->SetMarkerStyle(0);
    h_R2->SetLineStyle(1);
    h_R2->Draw("e1");

    TH1D *h_R1 = (TH1D *)hjp1->Clone();
    h_R1->Add(hcom, -1);
    h_R1->Divide(hcom);
    h_R1->SetBit(TH1::kNoTitle);
    // setBinError(h_R1, hcom, hjp1);
    h_R1->SetLineColor(2);
    h_R1->SetLineStyle(1);
    h_R1->SetMarkerStyle(0);
    h_R1->Draw("e1 same ");

    TH1D *h_R0 = (TH1D *)hjp0->Clone();
    h_R0->Add(hcom, -1);
    h_R0->Divide(hcom);
    h_R0->SetBit(TH1::kNoTitle);
    // setBinError(h_R0, hcom, hjp0);
    h_R0->SetLineColor(1);
    h_R0->SetLineStyle(1);
    h_R0->SetMarkerStyle(0);
    // h_R0->SetMarkerColor(1);
    h_R0->Draw("e1 same ");

    TLegend *leg2 = new TLegend(0.3, 0.85, 0.6, 0.95);
    leg2->SetNColumns(3);
    leg2->AddEntry(h_R2, " jp2", "le");
    leg2->AddEntry(h_R1, " jp1", "le");
    leg2->AddEntry(h_R0, " jp0", "le");
    leg2->Draw();

    gPad->Update();

    cxsec->Update();
    cxsec->SaveAs("ResultsC3/xsec_all.pdf");
    cxsec->SaveAs("ResultsC3/xsec_all.png");
    // cxsec->SaveAs("ResultsC3/xsec_all.eps");
    delete cxsec;
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
    hpyth->SetMaximum(5e10);
    hpyth->SetMinimum(10);
    hpyth->SetLineColor(1);
    hpyth->SetLineWidth(2);
    hpyth->SetLineStyle(1);
    hpyth->Draw("E1");

    hcom->SetLineColor(2);
    hcom->SetLineWidth(2);
    // hcom->SetMarkerStyle(20);
    hcom->SetMarkerColor(2);
    hcom->Draw("SAME E1");

    TLegend *lx = new TLegend(0.45, 0.61, 0.75, 0.65);
    lx->SetNColumns(2);
    lx->SetTextSize(0.033);
    lx->AddEntry(hcom, "#font[22]{ Measured}", "lep");
    lx->AddEntry(hpyth, "#font[22]{ Pythia 6.4.28 @ Perugia 2012, PARP(90)=0.213, CTEQ6 PDFs}", "lep");
    lx->Draw();

    TLatex tex;
    tex.SetTextAlign(12);
    tex.SetTextSize(0.045);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.42, "#font[22]{#color[2]{STAR Preliminary 2012}}");
    tex.SetTextSize(0.033);
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.07, "#font[22]{p + p #rightarrow #pi^{+}#pi^{-} + X  at #sqrt{s} = 200 GeV }");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.02, "#font[22]{|#eta^{#pi^{+}#pi^{-}}| < 1, 1 < p_{T}^{#pi^{+}#pi^{-}} < 15 (GeV/c), 0.02 < cone < 0.7}");
    tex.DrawLatex(1.3, hpyth->GetMaximum() * 0.006, "#font[22]{0.27 < M_{inv}^{#pi^{+}#pi^{-}} < 4.0 (GeV/c^{2}),  L_{int} = 26 (pb)^{-1} }");

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

    hsys_err->GetYaxis()->SetRangeUser(-0.6, 0.6);
    hsys_err->GetYaxis()->SetLabelSize(0.11);
    hsys_err->GetYaxis()->SetTitleOffset(0.55);
    // hsys_err->GetYaxis()->SetTitle("#font[22]{#frac{#sigma_{pyth} - #sigma_{mes}}{#sigma_{mes}}}");
    hsys_err->GetYaxis()->SetTitle("#font[22]{#delta#sigma}");
    hsys_err->GetYaxis()->SetTitleSize(0.11);
    hsys_err->GetXaxis()->SetLabelOffset(0.02);
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
    hrdiff_c->SetLineColor(2);
    hrdiff_c->SetLineWidth(2);
    hrdiff_c->Draw("e1 SAME");

    pad2->Update();
    TGaxis *axis = new TGaxis(pad2->GetUxmax(), pad2->GetUymin(), pad2->GetUxmax(), pad2->GetUymax(), -0.6, 0.6, 505, "+L");
    axis->SetTitle("Ratio");
    axis->SetTitleSize(0.11);
    axis->SetLabelSize(0.11);
    axis->SetLabelOffset(0.01);
    axis->SetTitleOffset(0.08);
    axis->CenterTitle();
    axis->Draw();

    TLegend *lg1 = new TLegend(0.18, 0.85, 0.38, 0.95);
    lg1->SetNColumns(3);
    // lg1->AddEntry(hsys_err, "#font[22]{ #delta#sigma_{sys}}", "f");
    lg1->AddEntry(hsys_err, "#font[22]{ #delta#sigma_{sys}}", "f");
    // lg1->AddEntry(hstat_err, "#font[22]{ #delta#sigma_{stat}}", "f");
    lg1->AddEntry(hstat_err, "#font[22]{ #delta#sigma_{stat}}", "f");
    // lg1->AddEntry(hrdiff_c, "#font[22]{ #times 0.5}", "lp");
    lg1->SetTextSize(0.1);
    lg1->Draw();
    gPad->Update();

    TLatex tex1;
    tex1.SetTextSize(0.08);
    // tex1.DrawLatex(0.45, hsys_err->GetMinimum() + 0.25, "#font[12]{10\% luminosity measurement uncertainty not shown.}");
    tex1.DrawLatex(1.35, 0.4, "#font[22]{10\% luminosity measurement uncertainty not shown.}");

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
void drawErrors(TH1D *htotal, TH1D *hsys, TH1D *hpid, TH1D *hloss, TH1D *htrgspr, TH1D *hstat, TH1D *hbias)
{
    TCanvas *cerr = new TCanvas("cerr", "", 500, 350);
    cerr->cd();
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    hsys->GetXaxis()->SetTitle("M_{inv}");
    hsys->GetXaxis()->SetTitleSize(0.06);
    hsys->GetYaxis()->SetTitleSize(0.06);
    hsys->GetXaxis()->SetLabelSize(0.06);
    hsys->GetYaxis()->SetLabelSize(0.06);
    hsys->GetYaxis()->SetTitle("Uncertainty");
    hsys->GetYaxis()->SetTitleOffset(1.);
    hsys->GetYaxis()->SetRangeUser(-0.29, 0.29);
    hsys->SetFillColorAlpha(42, 1);
    hsys->Draw("E2");

    TLegend *leg = new TLegend(0.2, 0.82, 0.85, 0.87);
    leg->SetNColumns(6);
    /*
        leg->AddEntry(hsys, "#delta_{sys}", "f");
        leg->Draw();
        cerr->Update();
        cerr->SaveAs("ResultsC3/sys_err_total.pdf");
        cerr->SaveAs("ResultsC3/sys_err_total.png");

    hsys->SetLineColor(2);
    hsys->Draw("E1 SAME");
    */
    leg->AddEntry(hsys, "#delta_{sys}", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("ResultsC3/sys_err_total_sys.pdf");
    cerr->SaveAs("ResultsC3/sys_err_total_sys.png");

    hpid->SetFillStyle(3144);
    hpid->SetFillColorAlpha(6, 0.5);
    hpid->Draw("E2 SAME");
    leg->AddEntry(hpid, "#delta_{fake}", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid.pdf");
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid.png");

    hloss->SetFillStyle(3144);
    hloss->SetFillColorAlpha(1, 0.5);
    hloss->Draw("E2 SAME");
    leg->AddEntry(hloss, "#delta_{loss}", "f");
    leg->Draw();
    cerr->Update();

    hstat->SetFillColorAlpha(2, 1);
    hstat->Draw("E2 SAME");
    leg->AddEntry(hstat, " #delta_{stat}", "f");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid_stat.pdf");
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid_stat.png");

    htrgspr->SetLineColorAlpha(3, 1);
    htrgspr->SetLineWidth(2);
    htrgspr->Draw("E1 SAME");
    leg->AddEntry(htrgspr, " #delta_{trg}", "lep");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid_stat_trgspr.pdf");
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid_stat_trgspr.png");

    hbias->SetLineColor(4);
    hbias->SetLineWidth(2);
    hbias->SetLineStyle(2);
    hbias->Draw("E1 SAME");
    leg->AddEntry(hbias, "#delta_{bias}", "lep");
    leg->Draw();
    cerr->Update();
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid_trgspr_stat_bias.pdf");
    cerr->SaveAs("ResultsC3/sys_err_total_sys_pid_trgspr_stat_bias.png");

    leg->Draw();
    cerr->Update();
    cerr->SaveAs("ResultsC3/sys_err_all.pdf");
    cerr->SaveAs("ResultsC3/sys_err_all.png");
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
