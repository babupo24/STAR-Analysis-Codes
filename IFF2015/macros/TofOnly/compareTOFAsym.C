#include <iostream>
#include "TGraphErrors.h"
#include <cmath>
#include "TCanvas.h"

using namespace std;


void compareTOFAsym(){
	gStyle->SetLegendBorderSize(0);
	double Mbkg[5]={0.386641, 0.582974, 0.756721, 0.969758, 1.32812};
	double Msig[5]={0.390506, 0.582592, 0.756392, 0.962861, 1.42839};
	double errx[5]={0};  

	double pTsig[5]    = {3.42273, 4.00489, 4.67533, 5.67678, 8.30727};
	double Agt_pTsig[5]={0.00465894, 0.00824099,0.00992505,0.0130413, 0.0302641};
        double Egt_pTsig[5]={0.00227235, 0.00208303,0.00196724,0.00190284, 0.00192238};
	double Alt_pTsig[5]={-0.00107789, 0.0012272,0.00254389,-0.000104231, 0.00496254};
        double Elt_pTsig[5]={0.00226486, 0.0020729,0.0019545,0.00187831, 0.00182148};

	double pTbkg[5]    ={3.37934, 3.97461, 4.63947, 5.62811, 8.14033}; 
        double Agt_pTbkg[5]={0.00592812, -0.00984245,-0.00994387,0.0162177, 0.0245844};
        double Egt_pTbkg[5]={0.00832525, 0.0108387,0.0138235,0.0171862, 0.0199603};
        double Alt_pTbkg[5]={-0.0120163, -0.00231737,-0.0080582,0.00331296, 0.0198817};
        double Elt_pTbkg[5]={0.0082977, 0.0107857,0.0137826,0.0171486, 0.0199587};
	
	double Agt_Msig[5]={0.00518938, 0.014124,0.0237483,0.0175566, 0.0106613};
	double Egt_Msig[5]={0.00196768, 0.00198494,0.00199649,0.00204516, 0.0020565};
	double Alt_Msig[5]={2.34039e-06, 0.000233546,0.00297367,0.00379058, 0.00164068};
	double Elt_Msig[5]={0.00196853, 0.00196012,0.00193228,0.00201406, 0.0020383};

	double Agt_Mbkg[5]={-0.0156053, 0.0072731,-0.00203504,0.00586201, 0.0138005};
	double Egt_Mbkg[5]={0.0119505, 0.0120655,0.0122357,0.0113117, 0.0131677};
	double Alt_Mbkg[5]={-0.0116171, -0.000518376,-0.014758,0.00101275, -0.00166904};
	double Elt_Mbkg[5]={0.0119087, 0.0120466,0.0122094,0.0112728, 0.0130882};



	TGraphErrors * grMsiggt = new TGraphErrors(5, Msig, Agt_Msig, errx, Egt_Msig);
	grMsiggt->SetTitle(0);
	grMsiggt->GetXaxis()->SetTitle("M_{inv}(p_{T} integrated)");
	grMsiggt->GetYaxis()->SetTitle("A_{UT}");
	grMsiggt->GetYaxis()->SetTitleOffset(1.2);
	grMsiggt->GetYaxis()->SetRangeUser(-0.03, 0.04);
	grMsiggt->GetXaxis()->SetLimits(0.25, 2.0);
	grMsiggt->GetXaxis()->SetTitleOffset(1.2);
	grMsiggt->SetMarkerStyle(20);
	grMsiggt->SetMarkerColor(2);
	grMsiggt->SetLineColor(2);

	TGraphErrors * grMsiglt = new TGraphErrors(5, Msig, Alt_Msig, errx, Elt_Msig);
	grMsiglt->SetMarkerStyle(24);
	grMsiglt->SetMarkerColor(2);
	grMsiglt->SetLineColor(2);


	TGraphErrors * grMbkggt = new TGraphErrors(5, Mbkg, Agt_Mbkg, errx, Egt_Mbkg);
	grMbkggt->SetMarkerStyle(20);
	grMbkggt->SetMarkerColor(4);
	grMbkggt->SetLineColor(4);

	TGraphErrors * grMbkglt = new TGraphErrors(5, Mbkg, Alt_Mbkg, errx, Elt_Mbkg);
	grMbkglt->SetMarkerStyle(24);
	grMbkglt->SetMarkerColor(4);
	grMbkglt->SetLineColor(4);

	TCanvas *canM = new TCanvas("canM","canM", 700,500);
	canM->cd();
	canM->SetLeftMargin(0.12);
	canM->SetBottomMargin(0.12);
	canM->SetGrid(0,0);
	grMsiggt->Draw("AP");
	grMsiglt->Draw("SAME P");
	grMbkggt->Draw("SAME P");
	grMbkglt->Draw("SAME P");
	
	canM->Update();
	TLine *lineM = new TLine(canM->GetUxmin(),0,canM->GetUxmax(),0);
	lineM->SetLineStyle(2);
	lineM->Draw();
	
	TLegend *legM = new TLegend(0.65,0.65,0.85,0.85);
	legM->AddEntry(grMsiggt, "Signal, #eta>0","lp");
	legM->AddEntry(grMsiglt, "Signal, #eta<0","lp");
	legM->AddEntry(grMbkggt, "Bkg, #eta>0","lp");
	legM->AddEntry(grMbkglt, "Bkg, #eta<0","lp");
	legM->SetTextSize(0.04);
	legM->Draw();
	
	canM->Update();

	TGraphErrors * grpTsiggt = new TGraphErrors(5, pTsig, Agt_pTsig, errx, Egt_pTsig);
	grpTsiggt->SetTitle(0);
	grpTsiggt->GetXaxis()->SetTitle("p_{T}(M_{inv} integrated)");
	grpTsiggt->GetYaxis()->SetTitle("A_{UT}");
	grpTsiggt->GetYaxis()->SetTitleOffset(1.2);
	grpTsiggt->GetYaxis()->SetRangeUser(-0.03, 0.04);
	grpTsiggt->GetXaxis()->SetLimits(3.0, 10.0);
	grpTsiggt->GetXaxis()->SetTitleOffset(1.2);
	grpTsiggt->SetMarkerStyle(20);
	grpTsiggt->SetMarkerColor(2);
	grpTsiggt->SetLineColor(2);

	TGraphErrors * grpTsiglt = new TGraphErrors(5, pTsig, Alt_pTsig, errx, Elt_pTsig);
	grpTsiglt->SetMarkerStyle(24);
	grpTsiglt->SetMarkerColor(2);
	grpTsiglt->SetLineColor(2);


	TGraphErrors * grpTbkggt = new TGraphErrors(5, pTbkg, Agt_pTbkg, errx, Egt_pTbkg);
	grpTbkggt->SetMarkerStyle(20);
	grpTbkggt->SetMarkerColor(4);
	grpTbkggt->SetLineColor(4);

	TGraphErrors * grpTbkglt = new TGraphErrors(5, pTbkg, Alt_pTbkg, errx, Elt_pTbkg);
	grpTbkglt->SetMarkerStyle(24);
	grpTbkglt->SetMarkerColor(4);
	grpTbkglt->SetLineColor(4);

	TCanvas *canpT = new TCanvas("canpT","canpT", 700,500);
	canpT->cd();
	canpT->SetLeftMargin(0.12);
	canpT->SetBottomMargin(0.12);
	canpT->SetGrid(0,0);
	grpTsiggt->Draw("AP");
	grpTsiglt->Draw("SAME P");
	grpTbkggt->Draw("SAME P");
	grpTbkglt->Draw("SAME P");
	
	canpT->Update();
	TLine *linepT = new TLine(canpT->GetUxmin(),0,canpT->GetUxmax(),0);
	linepT->SetLineStyle(2);
	linepT->Draw();
	
	TLegend *legpT = new TLegend(0.65,0.15,0.85,0.35);
	legpT->AddEntry(grpTsiggt, "Signal, #eta>0","lp");
	legpT->AddEntry(grpTsiglt, "Signal, #eta<0","lp");
	legpT->AddEntry(grpTbkggt, "Bkg, #eta>0","lp");
	legpT->AddEntry(grpTbkglt, "Bkg, #eta<0","lp");
	legpT->SetTextSize(0.04);
	legpT->Draw();
	
	canpT->Update();

	canpT->SaveAs("TOFAsymmetrySignalBkgMbin.pdf");
	canM->SaveAs("TOFAsymmetrySignalBkgpTbin.pdf");

	TCanvas *canpT1 = new TCanvas("canpT1","canpT1", 700,500);
	canpT1->cd();
	canpT1->SetLeftMargin(0.12);
	canpT1->SetBottomMargin(0.12);
	canpT1->SetGrid(0,0);

	TGraphErrors * grpTbkggt = new TGraphErrors(5, pTbkg, Agt_pTbkg, errx, Egt_pTbkg);
	grpTbkggt->SetTitle(0);
	grpTbkggt->GetXaxis()->SetTitle("p_{T}(M_{inv} integrated)");
	grpTbkggt->GetYaxis()->SetTitle("A_{UT}");
	grpTbkggt->GetYaxis()->SetTitleOffset(1.2);
	grpTbkggt->GetYaxis()->SetRangeUser(-0.03, 0.04);
	grpTbkggt->GetXaxis()->SetLimits(3.0, 10.0);
	grpTbkggt->GetXaxis()->SetTitleOffset(1.2);
	grpTbkggt->SetMarkerStyle(20);
	grpTbkggt->SetMarkerColor(2);
	grpTbkggt->SetLineColor(2);
}	
