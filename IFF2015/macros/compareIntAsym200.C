#include <iostream>
#include <iomanip>
#include "TGraphErrors.h"
#include "TMath.h"
#include "TCanvas.h"

using namespace std;

void compareIntAsym200()
{

	// integrated pT result
	double M7[9] = {0.3487, 0.453749, 0.564329, 0.662024, 0.756702, 0.874893, 1.03266, 1.2362, 1.68868};
	double asym7[9] = {0.0065842, 0.00506351, 0.0126739, 0.0143018, 0.022683, 0.0237341, 0.0158502, 0.0102208, 0.00869371};
	double err7[9] = {0.001514, 0.00150577, 0.00151905, 0.00152489, 0.00155356, 0.00142348, 0.00156034, 0.0015549, 0.00179993};
	// pion pair fraction
	double pairPurityMGt[9] = {0.768775, 0.815674, 0.79902, 0.796715, 0.825616, 0.805815, 0.785091, 0.758806, 0.718494};
	double pairpurityMLt[9] = {0.781753, 0.826413, 0.81047, 0.809187, 0.834432, 0.81532, 0.795783, 0.772011, 0.732798};
	// bias ratio
	double trigBiasGt[9] = {1.07026, 1.01152, 1.03012, 1.00929, 0.987516, 1.14202, 1.05634, 1.24891, 1.08074};
	double trigBiasGtErr[9] = {0.0522475, 0.0436093, 0.0512523, 0.0500996, 0.0469354, 0.0518734, 0.0542341, 0.0689986, 0.0443317};
	double trigBiasLt[9] = {1.15514, 1.10721, 0.886907, 0.908908, 1.0177, 1.00483, 1.15098, 1.0318, 1.16321};
	double trigBiasLtErr[9] = {0.0561346, 0.051315, 0.0472823, 0.0488303, 0.0450848, 0.052088, 0.0680179, 0.0571762, 0.0591917};

	// calculate systematic error
	double pidSys[9] = {0};
	double biasSys[9] = {0};
	double combSys[9] = {0};
	for (int i = 0; i < 9; i++)
	{
		pidSys[i] = (1 - pairPurityMGt[i]) * TMath::Max(asym7[i], err7[i]);
		biasSys[i] = (1 - trigBiasGt[i]) * TMath::Max(asym7[i], err7[i]);
		combSys[i] = sqrt(pow(pidSys[i], 2) + pow(biasSys[i], 2));
	}

	// run 11 results
	double run11M[7] = {0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 2.2};
	double run11asym[7] = {0.0023, 0.0269, 0.0207, 0.0249, 0.0298, 0.0196, 0.0104};
	double run11err[7] = {0.0090, 0.0075, 0.0072, 0.0060, 0.0050, 0.0060, 0.0061};
	double run11pid[7] = {0.0003, 0.0040, 0.0031, 0.0037, 0.0045, 0.0029, 0.0016};
	double run11bias[7] = {0.0002, 0.0028, 0.0022, 0.0026, 0.0031, 0.0021, 0.0011};
	double run11os[7] = {0, 0, 0, 0, 0, 0, 0};
	double run11sys[7] = {0};
	for (int i = 0; i < 7; i++)
	{
		run11sys[i] = sqrt(pow(run11pid[i], 2) + pow(run11bias[i], 2));
	}
	//-------------------
	double avg_pT7 = 5.25;

	double run06M[5] = {0.36, 0.50, 0.69, 0.88, 1.19};
	double run06asym[5] = {0.0054, 0.018, 0.023, 0.070, 0.039};
	double run06err[5] = {0.010, 0.0068, 0.0081, 0.013, 0.020};

	// theory from paper (500 GeV) extracted using xyscan
	double thM500[15] = {0.371, 0.433, 0.493, 0.523, 0.571, 0.633, 0.704, 0.751, 0.811, 0.886, 0.944, 0.996, 1.06, 1.13, 1.2};
	double thA500[15] = {0.00788, 0.0101, 0.013, 0.0166, 0.0225, 0.029, 0.0368, 0.0413, 0.0454, 0.0503, 0.0547, 0.051, 0.0478, 0.0441, 0.0402};
	double thE500[15] = {0.00365, 0.00406, 0.00522, 0.00591, 0.00771, 0.00986, 0.0126, 0.0141, 0.0157, 0.0173, 0.0186, 0.0179, 0.0167, 0.0155, 0.0142};

	// theory from Radici et. al. @√s = 200 GeV

	/*Mh     lower    upper
	  0.36   0.0027   0.0076
	  0.5    0.0077   0.0132
	  0.69   0.0173   0.0286
	  0.88   0.0262   0.0438
	  1.19   0.0200   0.0360
	 */
	double thM200[5] = {0.36, 0.5, 0.69, 0.88, 1.19};
	double thEl200[5] = {0.0027, 0.0077, 0.0173, 0.0262, 0.02};
	double thEh200[5] = {0.0076, 0.0132, 0.0286, 0.0438, 0.0360};

	double thA200[5] = {0};
	double thE200[5] = {0};
	double thEtest[5] = {0};
	for (int i = 0; i < 5; i++)
	{
		thA200[i] = 0.5 * (thEl200[i] + thEh200[i]);
		thE200[i] = thEh200[i] - thA200[i];
		thEtest[i] = (-1) * thEl200[i] + thA200[i];
		// cout<<thE200[i]<<" == "<<thEtest[i]<<endl;
	}

	ofstream outfile;
	outfile.open("intAutVsM200GeV.txt");
	for (int i = 0; i < 9; i++)
	{
		if (i < 5)
		{
			outfile << setprecision(4) << thM200[i] << " & " << thEl200[i] << " & " << thEh200[i] << " & " << M7[i] << " & " << asym7[i] << " & " << err7[i] << " & " << combSys[i] << " & " << run06M[i] << " & " << run06asym[i] << " & " << run06err[i] << endl;
		}
		else
		{
			outfile << setprecision(4) << " & "
					<< " & "
					<< " & " << M7[i] << " & " << asym7[i] << " & " << err7[i] << " & " << combSys[i] << " & "
					<< " & "
					<< " & " << endl;
		}
	}

	TCanvas *myCan = new TCanvas("myCan", "myCan", 700, 500);

	gStyle->SetOptStat(0);
	gStyle->SetOptDate(0);
	gStyle->SetTitle("");
	gStyle->SetLegendBorderSize(0);
	// myCan -> Divide(2,1);
	// myCan -> cd(1);
	// gPad->SetPad(0.0,.35,1.,1.);
	gPad->SetGrid(0, 0);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.12);

	TGraphErrors *th200 = new TGraphErrors(5, thM200, thA200, 0, thE200);
	// TGraphErrors *th200 = new TGraphErrors(15,thM500,thA500,0,thE500);
	th200->SetTitle("");
	th200->SetFillStyle(1001);
	th200->SetFillColorAlpha(39, 0.8);
	// th200->SetFillColor(kViolet);

	th200->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})  }}");
	th200->SetTitle("");
	th200->GetYaxis()->SetLabelSize(0.04);
	th200->GetYaxis()->SetLabelFont(22);
	th200->GetYaxis()->SetTitleSize(0.04);
	th200->GetYaxis()->SetTitleOffset(1.2);
	th200->GetYaxis()->SetRangeUser(-0.01, 0.099);
	th200->GetXaxis()->SetLimits(0.21, 2.4);
	// th200->GetXaxis()->SetNdivisions(505);
	// th200->GetYaxis()->SetNdivisions(500);
	th200->GetXaxis()->SetTitle("#font[22]{M^{#pi^{+}#pi^{-}}_{inv}(GeV/c^{2})}");
	th200->GetXaxis()->SetTitleSize(0.04);
	th200->GetXaxis()->SetLabelSize(0.04);
	th200->GetXaxis()->SetLabelFont(22);
	th200->GetXaxis()->SetTitleOffset(1.1);
	th200->Draw("AE3");

	double err[9] = {0};
	double errxps[9] = {0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
	double errxps11[7] = {0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};

	TGraphAsymmErrors *Run11pids = new TGraphAsymmErrors(7, run11M, run11asym, errxps, errxps, run11os, run11pid);
	Run11pids->SetFillStyle(3144);
	Run11pids->SetFillColor(6);
	Run11pids->Draw("same 2");
	TGraphAsymmErrors *Run11trigs = new TGraphAsymmErrors(7, run11M, run11asym, errxps, errxps, run11bias, run11os);
	Run11trigs->SetFillStyle(3144);
	Run11trigs->SetFillColor(9);
	Run11trigs->Draw("same 2");

	TGraphErrors *Run11 = new TGraphErrors(7, run11M, run11asym, err, run11err);
	Run11->SetMarkerStyle(20);
	Run11->SetMarkerSize(1.);
	Run11->SetMarkerColor(1);
	Run11->SetLineColor(1);
	Run11->SetLineStyle(1);
	// Run11->GetXaxis()->SetLimits(0.2, 2.5);
	Run11->Draw("P same");
	// run11----done

	TGraphErrors *Run06 = new TGraphErrors(5, run06M, run06asym, 0, run06err);
	Run06->SetMarkerStyle(20);
	Run06->SetMarkerColor(4);
	Run06->SetLineColor(4);
	Run06->SetMarkerSize(1.);
	Run06->GetXaxis()->SetLimits(0.2, 2.5);
	Run06->Draw("same P");

	TGraphErrors *Run15s = new TGraphErrors(9, M7, asym7, errxps, combSys);
	Run15s->SetFillStyle(0);
	Run15s->SetLineColor(2);
	Run15s->Draw("same 2");

	TGraphErrors *Run15 = new TGraphErrors(9, M7, asym7, err, err7);
	Run15->SetMarkerStyle(20);
	Run15->SetMarkerSize(1.);
	Run15->SetMarkerColor(kRed);
	Run15->SetLineColor(kRed);
	// Run15-> SetLineWidth(2);
	Run15->GetXaxis()->SetLimits(0.2, 2.5);
	// Run15-> GetYaxis()->SetRangeUser(-0.01, 0.09);
	Run15->Draw("P same");
	myCan->Update();
	TLine *line1 = new TLine(gPad->GetUxmin(), 0., gPad->GetUxmax(), 0.);
	line1->SetLineStyle(2);
	line1->Draw();

	TLegend *legend = new TLegend(0.58, 0.5, 0.85, 0.94);
	legend->AddEntry(th200, "Radici et. al., #sqrt{s} = 200 GeV", "f");
	legend->AddEntry(Run15, "Run 15, Cone < 0.7", "lpf");
	legend->AddEntry("", "#sqrt{s} = 200 GeV, #LT p_{T} #GT = 5.25 GeV/c", "");
	legend->AddEntry(Run11, "Run 11, Cone < 0.7", "lp");
	legend->AddEntry("", "#sqrt{s} = 510 GeV, #LT p_{T} #GT = 13 GeV/c", "");
	legend->AddEntry(Run06, "Run 06, Cone < 0.3", "lpf");
	legend->AddEntry("", "#sqrt{s} = 200 GeV, #LT p_{T} #GT = 6 GeV/c", "");
	legend->AddEntry(Run15s, " Run 15 Combined Syst. Error", "f");
	legend->AddEntry(Run11pids, " Run 11 PID Syst. Error", "f");
	legend->AddEntry(Run11trigs, " Run 11 Trig. Bias Syst. Error", "f");
	legend->AddEntry("", "#eta^{#pi^{+}#pi^{-}} > 0 ", " ");
	legend->SetTextSize(0.032);
	legend->SetTextFont(22);
	legend->Draw();

	TLatex latex;
	latex.SetTextFont(22);
	latex.SetTextSize(0.028);
	latex.SetTextAlign(13); // align at top
	latex.SetTextSize(0.055);
	latex.DrawLatex(0.27, 0.096, "#color[2]{#font[22]{STAR Preliminary 2015}}");
	latex.SetTextSize(0.04);
	latex.DrawLatex(0.27, 0.0895, "#color[1]{#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X} at #sqrt{s} = 200 GeV}");
	// latex.DrawLatex(0.27, 0.0895, "#color[1]{#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X}}");
	latex.SetTextSize(0.03);
	latex.DrawLatex(0.6, -0.003, "#color[1]{#font[22]{#pm 3% scale uncertainty from beam polarization (not shown)}}");

	myCan->SaveAs("STARAut200GeV_Wsys.pdf");
	myCan->SaveAs("STARAut200GeV_Wsys.png");
	// myCan->SaveAs("STARAut200GeV_WOsys.pdf");
}
