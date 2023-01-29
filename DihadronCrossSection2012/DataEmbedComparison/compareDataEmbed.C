#include <iostream>
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TStyle.h"

//--->START MAIN PROGRAM
//________________________________________________________________________________
// function to draw tracks quantities
void drawTracks(TH1D *hdat, TH1D *hemb, TCanvas *can, const char *trig);
// function to draw vertex and pair event quantities per trigger
void drawEvents(TH1D *hdat, TH1D *hemb, TCanvas *can, const char *trig);
void drawEventsTest(TH1D *hdat, TH1D *hemb, TH1D *hemball, TCanvas *can, const char *trig);
// function to draw track quantities for combined trigger
double getBinError(double dataBinCt, double dataBinErr, double embedBinCt, double embedBinErr);
void setBinError(TH1D *hRatio, TH1D *hData, TH1D *hEmbed);

void compareDataEmbed()
{

        gROOT->cd();

        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
        gStyle->SetStatX(0.43);
        gStyle->SetStatY(0.99);
        gStyle->SetStatW(0.15);
        gStyle->SetStatH(0.33);
        gStyle->SetLegendBorderSize(0);

        // TFile *embed = new TFile("ReadEmbedTrees_Test/histEmbed341.root ");
        TFile *embed = new TFile("./Embed/histEmbed.root");
        TFile *embedall = new TFile("../ReadEmbedTrees/histEmbed_all.root ");
        TFile *data = new TFile("./Data/histData.root ");
        // TFile *data = new TFile("ReadDataTrees_Test/histData.root ");
        if (!embed && !data)
        {
                cout << "At least one input file doesn't exist. "
                     << "\n";
                cout << "Exiting....." << endl;
                // break;
        }

        // Double_t presJP0 = 0.00845; //(1. / 118.335);
        Double_t presJP0 = 0.00707; //(1. / 141.35);
        // Double_t presJP0 = 1; //(1. / 141.35);
        Double_t presJP1 = 0.3978; //(1. / 2.514);
        // Double_t presJP1 = 1; //(1. / 2.514);
        Double_t presJP2 = 1.0;

        const int nTrig = 3;
        const char *trigType[3] = {"JP0", "JP1", "JP2"};
        const int nCh = 2;
        const char *chType[2] = {"Pos", "Neg"};

        TH1D *hpartonicPt;
        TCanvas *cpartonicPt;

        // histograms for combined trigger
        // for data
        TH1D *hVZDatJP;
        TCanvas *cVZJP;

        TH1D *hxsecMinvDatJP0;
        TH1D *hxsecMinvDatJP1;
        TH1D *hxsecMinvDatJP2;
        TH1D *hxsecMinvEmbJP0;
        TH1D *hxsecMinvEmbJP1;
        TH1D *hxsecMinvEmbJP2;
        TCanvas *cxsecMinvJP;

        TH1D *hMinvDatJP;
        TCanvas *cMinvJP;

        TH1D *hConeDatJP;
        TCanvas *cConeJP;

        TH1D *hEtaPairDatJP;
        TCanvas *cEtaPairJP;

        TH1D *hpTPairDatJP;
        TCanvas *cpTPairJP;

        TH1D *htrackEtaDatJP;
        TCanvas *ctrackEtaJP;

        TH1D *htrackPhiDatJP;
        TCanvas *ctrackPhiJP;

        TH1D *htrackpTDatJP;
        TCanvas *ctrackpTJP;

        TH1D *htracknSigmaPionDatJP;
        TCanvas *ctracknSigmaPionJP;

        TH1D *hpionEtaDatJP;
        TCanvas *cpionEtaJP;

        TH1D *hpionPhiDatJP;
        TCanvas *cpionPhiJP;

        TH1D *hpionpTDatJP;
        TCanvas *cpionpTJP;

        // for embedding
        TH1D *hVZEmbJP;
        TH1D *hMinvEmbJP;
        TH1D *hConeEmbJP;
        TH1D *hEtaPairEmbJP;
        TH1D *hpTPairEmbJP;

        TH1D *htrackEtaEmbJP;
        TH1D *htrackPhiEmbJP;
        TH1D *htrackpTEmbJP;
        TH1D *htracknSigmaPionEmbJP;

        TH1D *hpionEtaEmbJP;
        TH1D *hpionPhiEmbJP;
        TH1D *hpionpTEmbJP;

        // trigger wise histograms
        TH1D *hVZDat[nTrig];
        TH1D *hVZEmb[nTrig];
        TCanvas *cVZ[nTrig];

        TH1D *hMinvDat[nTrig];
        TH1D *hMinvEmb[nTrig];
        TCanvas *cMinv[nTrig];

        TH1D *hConeDat[nTrig];
        TH1D *hConeEmb[nTrig];
        TH1D *hConeEmbU[nTrig];
        TCanvas *cCone[nTrig];

        TH1D *hpTPairDat[nTrig];
        TH1D *hpTPairEmb[nTrig];
        TH1D *hpTPairEmbU[nTrig];
        TCanvas *cpTPair[nTrig];

        TH1D *hEtaPairDat[nTrig];
        TH1D *hEtaPairDatU[nTrig];
        TH1D *hEtaPairEmb[nTrig];
        TH1D *hEtaPairEmbU[nTrig];
        TCanvas *cEtaPair[nTrig];

        TH1D *htrackpTDat[nTrig];
        TH1D *htrackpTEmb[nTrig];
        TCanvas *ctrackpT[nTrig];

        TH1D *htrackpTEmbAll[nTrig];
        TCanvas *ctrackpTAll[nTrig];

        TH1D *htrackEtaDat[nTrig];
        TH1D *htrackEtaEmb[nTrig];
        TCanvas *ctrackEta[nTrig];

        TH1D *htrackPhiDat[nTrig];
        TH1D *htrackPhiEmb[nTrig];
        TCanvas *ctrackPhi[nTrig];

        TH1D *htracknSigmaPionDat[nTrig];
        TH1D *htracknSigmaPionEmb[nTrig];
        TCanvas *ctracknSigmaPion[nTrig];

        TH1D *hpionpTDat[nTrig];
        TH1D *hpionpTEmb[nTrig];
        TCanvas *cpionpT[nTrig];

        TH1D *hpionEtaDat[nTrig];
        TH1D *hpionEtaNegDatU[nTrig];
        TH1D *hpionEtaPosDatU[nTrig];
        TH1D *hpionEtaEmb[nTrig];
        TCanvas *cpionEta[nTrig];

        TH1D *hpionPhiDat[nTrig];
        TH1D *hpionPhiEmb[nTrig];
        TH1D *hpionPhiPosEmbU[nTrig];
        TH1D *hpionPhiNegEmbU[nTrig];
        TCanvas *cpionPhi[nTrig];
        TCanvas *cpionPhiPosU[nTrig];
        TCanvas *cpionPhiNegU[nTrig];

        // get partonic pt histogram and draw
        hpartonicPt = (TH1D *)embed->Get("hpartonicPt");
        TH1D *hpartonicPtAll = (TH1D *)embedall->Get("hpartonicPt");
        cpartonicPt = new TCanvas("cpartonicPt", "", 500, 500);
        cpartonicPt->cd();
        hpartonicPt->SetLineColor(2);
        hpartonicPtAll->SetMarkerColor(4);
        hpartonicPtAll->SetTitle("");
        hpartonicPtAll->GetXaxis()->SetTitle("partonic p_{T} GeV/c");
        hpartonicPtAll->Draw(" hist E ");
        hpartonicPt->Draw("hist E same");
        TLegend *lpart = new TLegend(0.75, 0.75, 0.9, 0.9);
        lpart->AddEntry(hpartonicPt, "Matched sample", "l");
        lpart->AddEntry(hpartonicPtAll, "Full sample", "lp");
        lpart->Draw();
        cpartonicPt->SetLogy();
        cpartonicPt->Update();
        cpartonicPt->SaveAs("./Plots/partonicPt.pdf");

        // get hist in xsec bins
        hxsecMinvDatJP0 = (TH1D *)data->Get("hxsecMJP0");
        hxsecMinvDatJP1 = (TH1D *)data->Get("hxsecMJP1");
        hxsecMinvDatJP2 = (TH1D *)data->Get("hxsecMJP2");
        hxsecMinvEmbJP0 = (TH1D *)embed->Get("hxsecMJP0");
        // hxsecMinvEmbJP0->Scale(presJP0);
        hxsecMinvEmbJP1 = (TH1D *)embed->Get("hxsecMJP1");
        // hxsecMinvEmbJP1->Scale(presJP1);
        hxsecMinvEmbJP2 = (TH1D *)embed->Get("hxsecMJP2");
        // hxsecMinvEmbJP2->Scale(presJP2);
        // plot data and embedding mass distribution in xsec bins

        TCanvas *cm = new TCanvas("cm", "", 700, 900);
        cm->cd();
        TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
        pad1->SetLeftMargin(0.15);
        pad1->SetBottomMargin(0);
        pad1->Draw();
        pad1->cd();
        pad1->SetLogy();

        hxsecMinvEmbJP0->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
        hxsecMinvEmbJP0->GetYaxis()->SetTitle("#pi^{+#pi^{-}} Yields");
        hxsecMinvEmbJP0->GetYaxis()->SetTitleOffset(1.3);
        hxsecMinvEmbJP0->SetLineColor(2);
        hxsecMinvEmbJP0->SetLineStyle(2);
        hxsecMinvEmbJP0->SetMarkerStyle(8);
        hxsecMinvEmbJP0->SetMarkerColor(2);
        hxsecMinvEmbJP0->DrawNormalized("hist P E");
        hxsecMinvEmbJP1->SetLineColor(3);
        hxsecMinvEmbJP1->SetLineStyle(2);
        hxsecMinvEmbJP1->SetMarkerStyle(8);
        hxsecMinvEmbJP1->SetMarkerColor(3);
        hxsecMinvEmbJP1->DrawNormalized("hist P E SAME");
        hxsecMinvEmbJP2->SetLineColor(4);
        hxsecMinvEmbJP2->SetLineStyle(2);
        hxsecMinvEmbJP2->SetMarkerStyle(8);
        hxsecMinvEmbJP2->SetMarkerColor(4);
        hxsecMinvEmbJP2->DrawNormalized("hist P E SAME");

        hxsecMinvDatJP0->SetLineStyle(1);
        hxsecMinvDatJP0->SetLineColor(2);
        hxsecMinvDatJP0->DrawNormalized("hist E SAME");
        hxsecMinvDatJP1->SetLineStyle(1);
        hxsecMinvDatJP1->SetLineColor(3);
        hxsecMinvDatJP1->DrawNormalized("hist E SAME");
        hxsecMinvDatJP2->SetLineStyle(1);
        hxsecMinvDatJP2->SetLineColor(4);
        hxsecMinvDatJP2->DrawNormalized("hist E SAME");
        TLegend *lx = new TLegend(0.7, 0.7, 0.9, 0.9);
        lx->AddEntry(hxsecMinvDatJP0, "Data-JP0", "l");
        lx->AddEntry(hxsecMinvDatJP1, "Data-JP1", "l");
        lx->AddEntry(hxsecMinvDatJP2, "Data-JP2", "l");
        lx->AddEntry(hxsecMinvEmbJP0, "Embed-JP0", "lp");
        lx->AddEntry(hxsecMinvEmbJP1, "Embed-JP1", "lp");
        lx->AddEntry(hxsecMinvEmbJP2, "Embed-JP2", "lp");
        lx->Draw();
        cm->cd();
        TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->SetLeftMargin(0.15);
        pad2->Draw();
        pad2->cd();
        gPad->SetGrid(1, 1);
        TH1D *rJP0 = (TH1D *)hxsecMinvDatJP0->Clone();
        rJP0->Scale(1 / rJP0->Integral());
        TH1D *emJP0 = (TH1D *)hxsecMinvEmbJP0->Clone();
        emJP0->Scale(1 / emJP0->Integral());
        rJP0->Add(emJP0, -1.0);
        rJP0->Divide(emJP0);
        rJP0->GetYaxis()->SetRangeUser(-1.5, 1.5);
        rJP0->GetYaxis()->SetTitle("#frac{Dat - Emb}{Emb}");
        rJP0->GetXaxis()->SetTitle("M_{inv} GeV/c^{2}");
        rJP0->GetXaxis()->SetTitleSize(0.08);
        rJP0->GetXaxis()->SetLabelSize(0.08);
        rJP0->GetYaxis()->SetTitleSize(0.08);
        rJP0->GetYaxis()->SetLabelSize(0.08);
        rJP0->GetYaxis()->SetTitleOffset(0.5);
        rJP0->Draw();
        TH1D *rJP1 = (TH1D *)hxsecMinvDatJP1->Clone();
        rJP1->Scale(1 / rJP1->Integral());
        TH1D *emJP1 = (TH1D *)hxsecMinvEmbJP1->Clone();
        emJP1->Scale(1 / emJP1->Integral());
        rJP1->Add(emJP1, -1.0);
        rJP1->Divide(emJP1);
        rJP1->Draw("SAME");

        TH1D *rJP2 = (TH1D *)hxsecMinvDatJP2->Clone();
        rJP2->Scale(1 / rJP2->Integral());
        TH1D *emJP2 = (TH1D *)hxsecMinvEmbJP2->Clone();
        emJP2->Scale(1 / emJP2->Integral());
        rJP2->Add(emJP2, -1.0);
        rJP2->Divide(emJP2);
        rJP2->Draw("SAME");
        cm->Update();
        cm->SaveAs("./Plots/Minv_PerTrig.pdf");

        // compare Minv for combined trigger in xsec bins
        TH1D *hxsecMJPDat = (TH1D *)hxsecMinvDatJP2->Clone();
        hxsecMJPDat->Add(hxsecMinvDatJP1);
        hxsecMJPDat->Add(hxsecMinvDatJP0);
        hxsecMJPDat->SetName("InvMassJP");
        TH1D *hxsecMJPEmb = (TH1D *)hxsecMinvEmbJP2->Clone();
        hxsecMJPEmb->Add(hxsecMinvEmbJP1, presJP1);
        hxsecMJPEmb->Add(hxsecMinvEmbJP0, presJP0);
        drawEvents(hxsecMJPDat, hxsecMJPEmb, cxsecMinvJP, "JP");
        ///--------

        // get histograms per trigger per trigger type
        for (int i = 0; i < nTrig; i++)
        {

                // vertex
                hVZDat[i] = (TH1D *)data->Get(Form("hVZ_%s", trigType[i]));
                hVZDat[i]->Rebin(2);
                // hVZDat[i]->Sumw2();
                hVZEmb[i] = (TH1D *)embed->Get(Form("hVZ_%s", trigType[i]));
                hVZEmb[i]->Rebin(2);
                // opening angle between pion pair
                hConeDat[i] = (TH1D *)data->Get(Form("hCone_%s", trigType[i]));
                hConeDat[i]->Rebin(2);
                // hConeDat[i]->Sumw2();
                hConeEmb[i] = (TH1D *)embed->Get(Form("hCone_%s", trigType[i]));
                hConeEmb[i]->Rebin(2);
                // invariant mass
                hMinvDat[i] = (TH1D *)data->Get(Form("hMinv_%s", trigType[i]));
                hMinvDat[i]->Rebin(2);
                // hMinvDat[i]->Sumw2();
                hMinvEmb[i] = (TH1D *)embed->Get(Form("hMinv_%s", trigType[i]));
                hMinvEmb[i]->Rebin(2);
                // pTPair
                hpTPairDat[i] = (TH1D *)data->Get(Form("hpTPair_%s", trigType[i]));
                hpTPairDat[i]->Rebin(2);
                // hpTPairDat[i]->Sumw2();
                hpTPairEmb[i] = (TH1D *)embed->Get(Form("hpTPair_%s", trigType[i]));
                hpTPairEmb[i]->Rebin(2);
                // EtaPair
                hEtaPairDat[i] = (TH1D *)data->Get(Form("hEtaPair_%s", trigType[i]));
                hEtaPairDat[i]->Rebin(2);
                // hEtaPairDat[i]->Sumw2();
                hEtaPairEmb[i] = (TH1D *)embed->Get(Form("hEtaPair_%s", trigType[i]));
                hEtaPairEmb[i]->Rebin(2);

                // track eta
                htrackEtaDat[i] = (TH1D *)data->Get(Form("htrackEta_%s", trigType[i]));
                htrackEtaDat[i]->Rebin(2);
                // htrackEtaDat[i]->Sumw2();
                htrackEtaEmb[i] = (TH1D *)embed->Get(Form("htrackEta_%s", trigType[i]));
                htrackEtaEmb[i]->Rebin(2);
                // track pt
                htrackpTDat[i] = (TH1D *)data->Get(Form("htrackpT_%s", trigType[i]));
                htrackpTDat[i]->Rebin(2);
                // htrackpTDat[i]->Sumw2();
                htrackpTEmb[i] = (TH1D *)embed->Get(Form("htrackpT_%s", trigType[i]));
                htrackpTEmb[i]->Rebin(2);

                htrackpTEmbAll[i] = (TH1D *)embedall->Get(Form("htrackpT_%s", trigType[i]));
                htrackpTEmbAll[i]->Rebin(2);
                // track phi
                htrackPhiDat[i] = (TH1D *)data->Get(Form("htrackPhi_%s", trigType[i]));
                htrackPhiDat[i]->Rebin(2);
                // htrackPhiDat[i]->Sumw2();
                htrackPhiEmb[i] = (TH1D *)embed->Get(Form("htrackPhi_%s", trigType[i]));
                htrackPhiEmb[i]->Rebin(2);

                // track nSigmaPion
                htracknSigmaPionDat[i] = (TH1D *)data->Get(Form("htracknSigmaPion_%s", trigType[i]));
                htracknSigmaPionDat[i]->Rebin(2);
                // htracknSigmaPionDat[i]->Sumw2();
                htracknSigmaPionEmb[i] = (TH1D *)embed->Get(Form("htracknSigmaPion_%s", trigType[i]));
                htracknSigmaPionEmb[i]->Rebin(2);

                // pion eta
                hpionEtaDat[i] = (TH1D *)data->Get(Form("hpionEta_%s", trigType[i]));
                hpionEtaDat[i]->Rebin(2);
                // hpionEtaDat[i]->Sumw2();
                hpionEtaEmb[i] = (TH1D *)embed->Get(Form("hpionEta_%s", trigType[i]));
                hpionEtaEmb[i]->Rebin(2);
                // pion pt
                hpionpTDat[i] = (TH1D *)data->Get(Form("hpionpT_%s", trigType[i]));
                hpionpTDat[i]->Rebin(2);
                // hpionpTDat[i]->Sumw2();
                hpionpTEmb[i] = (TH1D *)embed->Get(Form("hpionpT_%s", trigType[i]));
                hpionpTEmb[i]->Rebin(2);
                // pion phi
                hpionPhiDat[i] = (TH1D *)data->Get(Form("hpionPhi_%s", trigType[i]));
                hpionPhiDat[i]->Rebin(2);
                // hpionPhiDat[i]->Sumw2();
                hpionPhiEmb[i] = (TH1D *)embed->Get(Form("hpionPhi_%s", trigType[i]));
                hpionPhiEmb[i]->Rebin(2);

                // for combined trigger
                if (i == 0)
                {
                        hVZDatJP = (TH1D *)hVZDat[i]->Clone();
                        hVZDatJP->SetName("VZ_JP");

                        hMinvDatJP = (TH1D *)hMinvDat[i]->Clone();
                        hMinvDatJP->SetName("Minv_JP");
                        hConeDatJP = (TH1D *)hConeDat[i]->Clone();
                        hConeDatJP->SetName("Cone_JP");
                        hEtaPairDatJP = (TH1D *)hEtaPairDat[i]->Clone();
                        hEtaPairDatJP->SetName("Pair_Eta_JP");
                        hpTPairDatJP = (TH1D *)hpTPairDat[i]->Clone();
                        hpTPairDatJP->SetName("Pair_pT_JP");

                        htrackEtaDatJP = (TH1D *)htrackEtaDat[i]->Clone();
                        htrackEtaDatJP->SetName("Track_Eta_JP");
                        htrackPhiDatJP = (TH1D *)htrackPhiDat[i]->Clone();
                        htrackPhiDatJP->SetName("Track_Phi_JP");
                        htrackpTDatJP = (TH1D *)htrackpTDat[i]->Clone();
                        htrackpTDatJP->SetName("Track_pT_JP");
                        htracknSigmaPionDatJP = (TH1D *)htracknSigmaPionDat[i]->Clone();
                        htracknSigmaPionDatJP->SetName("Track_nSigmaPion_JP");

                        hpionEtaDatJP = (TH1D *)hpionEtaDat[i]->Clone();
                        hpionEtaDatJP->SetName("Pion_Eta_JP");
                        hpionPhiDatJP = (TH1D *)hpionPhiDat[i]->Clone();
                        hpionPhiDatJP->SetName("Pion_Phi_JP");
                        hpionpTDatJP = (TH1D *)hpionpTDat[i]->Clone();
                        hpionpTDatJP->SetName("Pion_pT_JP");
                }
                else
                {
                        hVZDatJP->Add(hVZDat[i]);
                        hMinvDatJP->Add(hMinvDat[i]);
                        hConeDatJP->Add(hConeDat[i]);
                        hEtaPairDatJP->Add(hEtaPairDat[i]);
                        hpTPairDatJP->Add(hpTPairDat[i]);

                        htrackEtaDatJP->Add(htrackEtaDat[i]);
                        htrackPhiDatJP->Add(htrackPhiDat[i]);
                        htrackpTDatJP->Add(htrackpTDat[i]);
                        htracknSigmaPionDatJP->Add(htracknSigmaPionDat[i]);

                        hpionEtaDatJP->Add(hpionEtaDat[i]);
                        hpionPhiDatJP->Add(hpionPhiDat[i]);
                        hpionpTDatJP->Add(hpionpTDat[i]);
                }
        }
        // combine triggers for embedding
        // get JP2, which had prescale = 1
        hVZEmbJP = (TH1D *)hVZEmb[2]->Clone();
        hVZEmbJP->SetName("VZ (JP)");
        hMinvEmbJP = (TH1D *)hMinvEmb[2]->Clone();
        hMinvEmbJP->SetName("Minv (JP)");
        hConeEmbJP = (TH1D *)hConeEmb[2]->Clone();
        hConeEmbJP->SetName("Cone (JP)");
        hpTPairEmbJP = (TH1D *)hpTPairEmb[2]->Clone();
        hpTPairEmbJP->SetName("Pair p_{T} (JP)");
        hEtaPairEmbJP = (TH1D *)hEtaPairEmb[2]->Clone();
        hEtaPairEmbJP->SetName("Pair #eta (JP)");

        htrackpTEmbJP = (TH1D *)htrackpTEmb[2]->Clone();
        htrackpTEmbJP->SetName("track p_{T} (JP)");
        htrackEtaEmbJP = (TH1D *)htrackEtaEmb[2]->Clone();
        htrackEtaEmbJP->SetName("track #eta (JP)");
        htrackPhiEmbJP = (TH1D *)htrackPhiEmb[2]->Clone();
        htrackPhiEmbJP->SetName("track #phi (JP)");
        htracknSigmaPionEmbJP = (TH1D *)htracknSigmaPionEmb[2]->Clone();
        htracknSigmaPionEmbJP->SetName("track n#sigma#pi (JP)");

        hpionpTEmbJP = (TH1D *)hpionpTEmb[2]->Clone();
        hpionpTEmbJP->SetName("pion p_{T} (JP)");
        hpionEtaEmbJP = (TH1D *)hpionEtaEmb[2]->Clone();
        hpionEtaEmbJP->SetName("pion #eta  (JP)");
        hpionPhiEmbJP = (TH1D *)hpionPhiEmb[2]->Clone();
        hpionPhiEmbJP->SetName("pion #phi  (JP)");
        // combined all triggers distributions for embedding using prescale
        for (int i = 0; i < 2; i++)
        {
                if (i == 0)
                {
                        // add prescaled JP0 to JP2
                        hVZEmbJP->Add(hVZEmb[i], presJP0);
                        hMinvEmbJP->Add(hMinvEmb[i], presJP0);
                        hConeEmbJP->Add(hConeEmb[i], presJP0);
                        hpTPairEmbJP->Add(hpTPairEmb[i], presJP0);
                        hEtaPairEmbJP->Add(hEtaPairEmb[i], presJP0);

                        htrackEtaEmbJP->Add(htrackEtaEmb[i], presJP0);
                        htrackPhiEmbJP->Add(htrackPhiEmb[i], presJP0);
                        htrackpTEmbJP->Add(htrackpTEmb[i], presJP0);
                        htracknSigmaPionEmbJP->Add(htracknSigmaPionEmb[i], presJP0);

                        hpionEtaEmbJP->Add(hpionEtaEmb[i], presJP0);
                        hpionPhiEmbJP->Add(hpionPhiEmb[i], presJP0);
                        hpionpTEmbJP->Add(hpionpTEmb[i], presJP0);
                }
                else if (i == 1)
                {
                        // add prescaled JP1 to JP2
                        hVZEmbJP->Add(hVZEmb[i], presJP1);
                        hMinvEmbJP->Add(hMinvEmb[i], presJP1);
                        hConeEmbJP->Add(hConeEmb[i], presJP1);
                        hpTPairEmbJP->Add(hpTPairEmb[i], presJP1);
                        hEtaPairEmbJP->Add(hEtaPairEmb[i], presJP1);

                        htrackEtaEmbJP->Add(htrackEtaEmb[i], presJP1);
                        htrackPhiEmbJP->Add(htrackPhiEmb[i], presJP1);
                        htrackpTEmbJP->Add(htrackpTEmb[i], presJP1);
                        htracknSigmaPionEmbJP->Add(htracknSigmaPionEmb[i], presJP1);

                        hpionEtaEmbJP->Add(hpionEtaEmb[i], presJP1);
                        hpionPhiEmbJP->Add(hpionPhiEmb[i], presJP1);
                        hpionpTEmbJP->Add(hpionpTEmb[i], presJP1);
                }
        }
        // draw events and tracks quantities per trigger
        for (int i = 0; i < nTrig; i++)
        {
                drawEvents(hVZDat[i], hVZEmb[i], cVZ[i], trigType[i]);
                drawEvents(hConeDat[i], hConeEmb[i], cCone[i], trigType[i]);
                drawEvents(hpTPairDat[i], hpTPairEmb[i], cpTPair[i], trigType[i]);
                drawEvents(hEtaPairDat[i], hEtaPairEmb[i], cEtaPair[i], trigType[i]);
                drawEvents(hMinvDat[i], hMinvEmb[i], cMinv[i], trigType[i]);

                drawEvents(htrackEtaDat[i], htrackEtaEmb[i], ctrackEta[i], trigType[i]);
                drawEvents(htrackPhiDat[i], htrackPhiEmb[i], ctrackPhi[i], trigType[i]);
                drawEvents(htrackpTDat[i], htrackpTEmb[i], ctrackpT[i], trigType[i]);
                // drawEventsTest(htrackpTDat[i], htrackpTEmb[i], htrackpTEmbAll[i], ctrackpTAll[i], trigType[i]);
                drawEvents(htracknSigmaPionDat[i], htracknSigmaPionEmb[i], ctracknSigmaPion[i], trigType[i]);

                drawEvents(hpionEtaDat[i], hpionEtaEmb[i], cpionEta[i], trigType[i]);
                drawEvents(hpionPhiDat[i], hpionPhiEmb[i], cpionPhi[i], trigType[i]);
                drawEvents(hpionpTDat[i], hpionpTEmb[i], cpionpT[i], trigType[i]);
        }

        // draw combined trigger
        drawEvents(hVZDatJP, hVZEmbJP, cVZJP, "JP");
        drawEvents(hMinvDatJP, hMinvEmbJP, cMinvJP, "JP");

        drawEvents(hpTPairDatJP, hpTPairEmbJP, cpTPairJP, "JP");
        drawEvents(hEtaPairDatJP, hEtaPairEmbJP, cEtaPairJP, "JP");
        drawEvents(hConeDatJP, hConeEmbJP, cConeJP, "JP");

        drawEvents(htrackEtaDatJP, htrackEtaEmbJP, ctrackEtaJP, "JP");
        drawEvents(htrackPhiDatJP, htrackPhiEmbJP, ctrackPhiJP, "JP");
        drawEvents(htrackpTDatJP, htrackpTEmbJP, ctrackpTJP, "JP");
        drawEvents(htracknSigmaPionDatJP, htracknSigmaPionEmbJP, ctracknSigmaPionJP, "JP");

        drawEvents(hpionEtaDatJP, hpionEtaEmbJP, cpionEtaJP, "JP");
        drawEvents(hpionPhiDatJP, hpionPhiEmbJP, cpionPhiJP, "JP");
        drawEvents(hpionpTDatJP, hpionpTEmbJP, cpionpTJP, "JP");
}

void drawEvents(TH1D *hdat, TH1D *hemb, TCanvas *can, const char *trig)
{
        const char *hnames[37] = {"htrackpT_JP0", "htrackpT_JP1", "htrackpT_JP2", "htracknSigmaPion_JP0", "htracknSigmaPion_JP1", "htracknSigmaPion_JP2", "hpionpT_JP0", "hpionpT_JP1", "hpionpT_JP2", "hMinv_JP0", "hMinv_JP1", "hMinv_JP2", "hpTPair_JP0", "hpTPair_JP1", "hpTPair_JP2", "Track_pT_JP", "Track_nSigmaPion_JP", "Pion_pT_JP", "Minv_JP", "Pair_pT_JP", "UniquePionNeg_pT_JP", "UniquePionPos_pT_JP", "UniquePair_pT_JP", "UniquePair_Minv_JP", "hpionpTNegU_JP2", "hpionpTNegU_JP1", "hpionpTNegU_JP0", "hpionpTPosU_JP2", "hpionpTPosU_JP1", "hpionpTPosU_JP0", "hMinvU_JP2", "hMinvU_JP1", "hMinvU_JP0", "hpTPairU_JP2", "hpTPairU_JP1", "hpTPairU_JP0", "InvMassJP"};

        can = new TCanvas(Form("c%s_%s", hdat->GetName(), trig), "", 0, 0, 600, 800);
        can->cd();
        TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
        pad1->SetLeftMargin(0.15);
        pad1->SetBottomMargin(0);
        pad1->Draw();
        pad1->cd();
        gPad->SetGrid(0, 0);

        for (int i = 0; i < 37; i++)
        {
                if (strcmp(hnames[i], hdat->GetName()) == 0)
                {
                        gPad->SetLogy();
                }
        }

        hdat->SetLineWidth(2);
        hdat->SetLineColor(4);
        hdat->Scale(1 / hdat->Integral());
        // hdat->GetXaxis()->SetTitleOffset(0.5);
        hemb->GetXaxis()->SetTitleOffset(0.5);
        // hdat->GetXaxis()->SetTitleSize(0.6);
        hemb->GetXaxis()->SetTitleSize(0.6);
        // hdat->GetYaxis()->SetLabelSize(0.04);
        hemb->GetYaxis()->SetLabelSize(0.06);
        // hdat->SetTitle(Form("%s", hdat->GetName()));
        hemb->SetTitle("");
        hemb->SetMarkerStyle(24);
        hemb->SetMarkerSize(.5);
        hemb->SetMarkerColor(2);
        hemb->SetLineColor(2);
        hemb->Scale(1 / hemb->Integral());

        hemb->Draw("hist  E");
        hdat->Draw("hist E SAME");
        gPad->Update();
        // TLegend *leg = new TLegend(0.78, 0.77, 0.88, 0.88);
        TLegend *leg = new TLegend(0.45, 0.05, 0.65, 0.2);
        leg->AddEntry(hdat, "Data", "l");
        leg->AddEntry(hemb, "Embed.", "lp");
        leg->SetTextSize(0.05);
        leg->Draw();
        gPad->Update();

        can->cd();
        TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->SetLeftMargin(0.15);
        pad2->Draw();
        pad2->cd();
        gPad->SetGrid(1, 1);

        TH1D *h_R = (TH1D *)hdat->Clone();
        h_R->Add(hemb, -1);
        h_R->Divide(hemb);
        h_R->SetBit(TH1::kNoTitle);
        setBinError(h_R, hdat, hemb);
        h_R->SetLineColor(2);
        // h_R->SetMarkerStyle(1);
        // h_R->SetMarkerColor(1);
        h_R->GetYaxis()->SetTitle("#frac{Data-Embed}{Embed}");
        h_R->GetYaxis()->SetTitleSize(0.08);
        h_R->GetYaxis()->CenterTitle();
        h_R->GetYaxis()->SetRangeUser(-1.0, 1.0);
        h_R->GetYaxis()->SetLabelSize(0.15);
        h_R->GetYaxis()->SetTitleOffset(0.42);
        h_R->GetYaxis()->SetTitleSize(0.15);
        h_R->GetXaxis()->SetTitle(Form("%s", hdat->GetName()));
        h_R->GetXaxis()->SetTitleOffset(1.0);
        h_R->GetXaxis()->SetTitleSize(0.15);
        h_R->GetXaxis()->SetLabelSize(0.15);
        h_R->GetYaxis()->SetNdivisions(505);
        h_R->SetMarkerStyle(2);
        h_R->Draw();
        //     h_R->Fit("pol0");

        hdat->SetDirectory(gROOT);
        hemb->SetDirectory(gROOT);
        can->Draw();
        can->SaveAs(Form("./Plots/%s.pdf", hdat->GetName()));
}
void drawEventsTest(TH1D *hdat, TH1D *hemb, TH1D *hemball, TCanvas *can, const char *trig)
{
        const char *hnames[37] = {"htrackpT_JP0", "htrackpT_JP1", "htrackpT_JP2", "htracknSigmaPion_JP0", "htracknSigmaPion_JP1", "htracknSigmaPion_JP2", "hpionpT_JP0", "hpionpT_JP1", "hpionpT_JP2", "hMinv_JP0", "hMinv_JP1", "hMinv_JP2", "hpTPair_JP0", "hpTPair_JP1", "hpTPair_JP2", "Track_pT_JP", "Track_nSigmaPion_JP", "Pion_pT_JP", "Minv_JP", "Pair_pT_JP", "UniquePionNeg_pT_JP", "UniquePionPos_pT_JP", "UniquePair_pT_JP", "UniquePair_Minv_JP", "hpionpTNegU_JP2", "hpionpTNegU_JP1", "hpionpTNegU_JP0", "hpionpTPosU_JP2", "hpionpTPosU_JP1", "hpionpTPosU_JP0", "hMinvU_JP2", "hMinvU_JP1", "hMinvU_JP0", "hpTPairU_JP2", "hpTPairU_JP1", "hpTPairU_JP0", "InvMassJP"};

        can = new TCanvas(Form("c%s_%s", hdat->GetName(), trig), "", 0, 0, 600, 800);
        can->cd();
        TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
        pad1->SetLeftMargin(0.15);
        pad1->SetBottomMargin(0);
        pad1->Draw();
        pad1->cd();
        gPad->SetGrid(0, 0);

        for (int i = 0; i < 37; i++)
        {
                if (strcmp(hnames[i], hdat->GetName()) == 0)
                {
                        gPad->SetLogy();
                }
        }

        hdat->SetLineWidth(2);
        hdat->SetLineColor(4);
        hdat->Scale(1 / hdat->Integral());
        // hdat->GetXaxis()->SetTitleOffset(0.5);
        hemball->GetXaxis()->SetTitleOffset(0.5);
        // hdat->GetXaxis()->SetTitleSize(0.6);
        hemball->GetXaxis()->SetTitleSize(0.6);
        // hdat->GetYaxis()->SetLabelSize(0.04);
        hemball->GetYaxis()->SetLabelSize(0.06);
        // hdat->SetTitle(Form("%s", hdat->GetName()));
        hemball->SetTitle("");
        hemb->SetMarkerStyle(24);
        hemb->SetMarkerSize(.5);
        hemb->SetMarkerColor(2);
        hemb->SetLineColor(2);
        hemb->Scale(1 / hemb->Integral());

        hemball->SetLineColor(2);
        hemball->SetLineWidth(2);
        hemball->Scale(1 / hemball->Integral());

        hemball->Draw("hist  E");
        hemb->Draw("hist  E SAME");
        hdat->Draw("hist E SAME");
        gPad->Update();
        // TLegend *leg = new TLegend(0.78, 0.77, 0.88, 0.88);
        TLegend *leg = new TLegend(0.45, 0.05, 0.65, 0.25);
        leg->AddEntry(hdat, "Data", "l");
        leg->AddEntry(hemb, "Embed. Good", "lp");
        leg->AddEntry(hemball, "Embed. All", "lp");
        leg->SetTextSize(0.05);
        leg->Draw();
        gPad->Update();

        can->cd();
        TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->SetLeftMargin(0.15);
        pad2->Draw();
        pad2->cd();
        gPad->SetGrid(1, 1);

        TH1D *h_R = (TH1D *)hdat->Clone();
        h_R->Add(hemb, -1);
        h_R->Divide(hemb);
        h_R->SetBit(TH1::kNoTitle);
        setBinError(h_R, hdat, hemb);
        h_R->SetLineColor(2);
        // h_R->SetMarkerStyle(1);
        // h_R->SetMarkerColor(1);
        h_R->GetYaxis()->SetTitle("#frac{Data-Embed}{Embed}");
        h_R->GetYaxis()->SetTitleSize(0.08);
        h_R->GetYaxis()->CenterTitle();
        h_R->GetYaxis()->SetRangeUser(-1.0, 1.0);
        h_R->GetYaxis()->SetLabelSize(0.15);
        h_R->GetYaxis()->SetTitleOffset(0.42);
        h_R->GetYaxis()->SetTitleSize(0.15);
        h_R->GetXaxis()->SetTitle(Form("%s", hdat->GetName()));
        h_R->GetXaxis()->SetTitleOffset(1.0);
        h_R->GetXaxis()->SetTitleSize(0.15);
        h_R->GetXaxis()->SetLabelSize(0.15);
        h_R->GetYaxis()->SetNdivisions(505);
        h_R->SetMarkerStyle(2);
        h_R->Draw();
        //     h_R->Fit("pol0");

        hdat->SetDirectory(gROOT);
        hemb->SetDirectory(gROOT);
        can->Draw();
        can->SaveAs(Form("./Plots/%s.pdf", hdat->GetName()));
}
double getBinError(double dataBinCt, double dataBinErr, double embedBinCt, double embedBinErr)
{

        return (dataBinCt / embedBinCt) * sqrt(pow(dataBinErr / dataBinCt, 2) + pow(embedBinErr / embedBinCt, 2));
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

void drawTracks(TH1D *hdat, TH1D *hemb, TCanvas *can, const char *trig)
{
        const char *hnames[24] = {"htrackpT_JP0", "htrackpT_JP1", "htrackpT_JP2", "htracknSigmaPion_JP0", "htracknSigmaPion_JP1", "htracknSigmaPion_JP2", "hpionpT_JP0", "hpionpT_JP1", "hpionpT_JP2", "hMinv_JP0", "hMinv_JP1", "hMinv_JP2", "hpTPair_JP0", "hpTPair_JP1", "hpTPair_JP2", "Track_pT_JP", "Track_nSigmaPion_JP", "Pion_pT_JP", "Minv_JP", "Pair_pT_JP", "UniquePionNeg_pT_JP", "UniquePionPos_pT_JP", "UniquePair_pT_JP", "UniquePair_Minv_JP"};

        can = new TCanvas(Form("c%s_%s", hdat->GetName(), trig), "", 0, 0, 600, 700);
        can->cd();

        TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0);
        pad1->SetLeftMargin(0.15);
        pad1->Draw();
        pad1->cd();
        gPad->SetGrid(0, 0);
        for (int i = 0; i < 24; i++)
        {
                if (strcmp(hnames[i], hdat->GetName()) == 0)
                {
                        gPad->SetLogy();
                        continue;
                }
        }
        hdat->SetLineWidth(2);
        hdat->SetLineColor(4);
        hdat->GetXaxis()->SetTitleOffset(0.8);
        hdat->GetXaxis()->SetTitleSize(0.6);
        hdat->Scale(1 / hdat->Integral());
        // hdat->Sumw2(kFALSE);
        hdat->SetTitle(Form("%s", hdat->GetName()));
        // hemb->SetLineWidth(2);
        hemb->SetLineColor(2);
        hemb->Scale(1 / hemb->Integral());

        hdat->Draw("hist E");
        hemb->Draw("hist E same");

        TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg->AddEntry(hdat, "Data", "l");
        leg->AddEntry(hemb, "Embed.", "l");
        leg->SetTextSize(0.03);
        leg->Draw();

        can->cd();
        TPad *pad2 = new TPad("pad2", "", 0, 0.05, 1.0, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetLeftMargin(0.2);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();
        gPad->SetGrid(1, 1);

        TH1D *h_R = (TH1D *)hdat->Clone();
        h_R->Add(hemb, -1);
        h_R->Divide(hemb);
        h_R->SetBit(TH1::kNoTitle);
        h_R->SetLineColor(1);
        h_R->GetYaxis()->SetTitle("#frac{Data-Embed}{Embed}");
        h_R->GetYaxis()->SetTitleSize(0.08);
        h_R->GetYaxis()->CenterTitle();
        h_R->GetYaxis()->SetRangeUser(-.5, .5);
        h_R->GetYaxis()->SetLabelSize(0.15);
        h_R->GetYaxis()->SetTitleOffset(0.42);
        h_R->GetYaxis()->SetTitleSize(0.15);
        h_R->GetXaxis()->SetTitle(Form("%s", hdat->GetName()));
        h_R->GetXaxis()->SetTitleOffset(1.0);
        h_R->GetXaxis()->SetTitleSize(0.15);
        h_R->GetXaxis()->SetLabelSize(0.15);
        h_R->GetYaxis()->SetNdivisions(505);
        h_R->SetMarkerStyle(2);
        h_R->Draw("E");

        hdat->SetDirectory(gROOT);
        hemb->SetDirectory(gROOT);
        can->Draw();
        can->SaveAs(Form("./Plots/%s.pdf", hdat->GetName()));
}
