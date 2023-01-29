#include <iostream>
#include <iomanip>
#include "TGraphErrors.h"
#include "TMath.h"
#include "TCanvas.h"


using namespace std;

void asymIntMinv_AfterPIDRelativeDiff()
{
	
//integrated pT result ----Preliminary
double M7[9] = {0.3487, 0.453749, 0.564329, 0.662024, 0.756702, 0.874893, 1.03266, 1.2362, 1.68868};
double asym7[9]={0.0065842, 0.00506351,0.0126739,0.0143018, 0.022683, 0.0237341, 0.0158502, 0.0102208, 0.00869371};
double err7[9]={0.001514, 0.00150577,0.00151905,0.00152489, 0.00155356, 0.00142348, 0.00156034, 0.0015549, 0.00179993};
//pion pair fraction --Preliminary
double pairPurityMGt[9]={0.768775,0.815674,0.79902,0.796715,0.825616,0.805815,0.785091,0.758806,0.718494};
double pairpurityMLt[9]={0.781753,0.826413,0.81047,0.809187,0.834432,0.81532,0.795783,0.772011,0.732798};

//tpc only
double M[9] = {0.348699, 0.453745, 0.564318, 0.662012, 0.756711, 0.874902, 1.03267, 1.23621, 1.68853};
double A_tpc[9]={0.00634007, 0.00555659,0.0129098,0.0159221, 0.0223458, 0.0238164, 0.0146982, 0.00995199, 0.00804356};
double Aerr_tpc[9]={0.00153057, 0.00152348,0.00153763,0.00154682, 0.00156767, 0.00143773, 0.00157216, 0.00156993, 0.00181772};
//tpc or tof
double A_F[9]={0.00633879, 0.00530626,0.0139225,0.0166477, 0.0233493, 0.0251667, 0.0152965, 0.0109216, 0.00791919};
double Aerr_F[9]={0.00157303, 0.00156351,0.00158448,0.00159235, 0.00161233, 0.0014911, 0.00163803, 0.00163575, 0.00186892};
//tpc && tof
double A_tof[9]={0.00384879, 0.00416929,0.0129391,0.0164233, 0.0221693, 0.021228, 0.0119228, 0.00971422, 0.00941017};
double Aerr_tof[9]={0.00241906, 0.00242958,0.00242384,0.00242292, 0.00242979, 0.00219865, 0.00247791, 0.0024986, 0.00289007};
double xerr[9]={0};
double rdiff[9]={0};
double rdifferr[9]={0};
for(int i=0; i<9; i++){
rdiff[i]=(A_tpc[i]-A_tof[i])/A_tpc[i];
rdifferr[i]=sqrt(pow((A_tof[i]/A_tpc[i])*Aerr_tpc[i],2)+pow(Aerr_tof[i]/A_tpc[i],2));
}
TGraphErrors *gr_diff=new TGraphErrors(9,M,rdiff,xerr, rdifferr);
gr_diff->SetTitle("");
gr_diff->GetYaxis()->SetLabelSize(0.04);
gr_diff->GetYaxis()->SetLabelFont(22);
gr_diff->GetYaxis()->SetTitle("A_{UT}^{Tpc} - A_{UT}^{Tpc&&Tof}");
gr_diff->GetYaxis()->SetTitleSize(0.04);
gr_diff->GetYaxis()->SetTitleOffset(1.2);
gr_diff->GetYaxis()->SetRangeUser(-2,2);
gr_diff->GetXaxis()->SetLimits(0.21, 1.9);
gr_diff->GetXaxis()->SetTitle("#font[22]{M^{#pi^{+}#pi^{-}}_{inv}(GeV/c^{2})}");
gr_diff->GetXaxis()->SetTitleSize(0.04);
gr_diff-> GetXaxis()->SetLabelSize(0.04);
gr_diff-> GetXaxis()->SetLabelFont(22);
gr_diff->GetXaxis()->SetTitleOffset(1.1);
gr_diff->SetMarkerStyle(20);
gr_diff->SetMarkerColor(2);
gr_diff->SetLineColor(2);
TCanvas *cdiff=new TCanvas("cdiff","",700,500);
cdiff->cd();
gPad->SetGrid(0,0);
gPad->SetTopMargin(0.05);
gPad->SetLeftMargin(0.11);
gPad->SetRightMargin(0.05);
gPad->SetBottomMargin(0.12);
gr_diff->Draw("AP");
gr_diff->Fit("pol0");
gPad->Update();
 TLine *line=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
line->SetLineStyle(2);
line->Draw();
TLatex tex;
tex.SetTextSize(0.04);
tex.SetTextAlign(13);
tex.DrawLatex(0.4,1.8,Form("#font[22]{#color[2]{chi2/ndf = %g / %i}}",gr_diff->GetFunction("pol0")->GetChisquare(), gr_diff->GetFunction("pol0")->GetNDF()));
tex.DrawLatex(0.4,1.5,Form("#font[22]{#color[2]{p0 = %g #pm %g}}",gr_diff->GetFunction("pol0")->GetParameter(0), gr_diff->GetFunction("pol0")->GetParError(0)));
double fp0=gr_diff->GetFunction("pol0")->GetParameter(0);
cdiff->Update();
cdiff->SaveAs("./Plots/intMinv_Difference.pdf");

//pionpair fraction ----After PID only TOF tracks
 double pionPairMGtT[9]={0.975651,0.987897,0.985057,0.981433,0.984766,0.983844,0.98232,0.979519,0.967469};
//pion pair fraction ----After PID TPC or TOF
 double pionPairMGt[9]={0.866901,0.896771,0.888716,0.886221,0.903343,0.892411,0.879951,0.86373,0.837623};
//----TPC or TOF -----<<<
//bias ratio
double trigBiasGt[9]={1.07026, 1.01152, 1.03012, 1.00929, 0.987516, 1.14202, 1.05634, 1.24891, 1.08074};
double trigBiasGtErr[9]={0.0522475, 0.0436093, 0.0512523, 0.0500996, 0.0469354, 0.0518734, 0.0542341, 0.0689986, 0.0443317};
double trigBiasLt[9]={1.15514, 1.10721, 0.886907, 0.908908, 1.0177, 1.00483, 1.15098, 1.0318, 1.16321};
double trigBiasLtErr[9]={0.0561346, 0.051315, 0.0472823, 0.0488303, 0.0450848, 0.052088, 0.0680179, 0.0571762, 0.0591917};

//calculate systematic error
double pidSys[9]={0};
double biasSys[9]={0};
double combSys[9]={0};
for(int i=0; i<9; i++){
	//pidSys[i]=sqrt(pow(fp0,2)+pow(((1-pionPairMGtT[i])*max(A_F[i],Aerr_F[i])),2));
	pidSys[i]=fp0*A_tpc[i];
	biasSys[i]=(1-trigBiasGt[i])*TMath::Max(A_tpc[i],Aerr_tpc[i]);
	combSys[i]=sqrt(pow(pidSys[i],2)+pow(biasSys[i],2));
}


	
	double avg_pT7=5.25;

	double run06M[5]={0.36, 0.50, 0.69, 0.88, 1.19};
	double run06asym[5]={0.0054, 0.018, 0.023, 0.070, 0.039};
	double run06err[5]={0.010, 0.0068, 0.0081, 0.013, 0.020};

	//theory from paper (500 GeV) extracted using xyscan
	double thM500[15]={0.371, 0.433,  0.493, 0.523, 0.571, 0.633, 0.704, 0.751, 0.811, 0.886, 0.944, 0.996, 1.06, 1.13, 1.2};
	double thA500[15]={0.00788, 0.0101, 0.013, 0.0166, 0.0225, 0.029, 0.0368, 0.0413, 0.0454, 0.0503, 0.0547, 0.051, 0.0478, 0.0441, 0.0402};
	double thE500[15]={0.00365, 0.00406, 0.00522, 0.00591, 0.00771, 0.00986, 0.0126, 0.0141, 0.0157, 0.0173, 0.0186, 0.0179, 0.0167, 0.0155, 0.0142};

	//theory from Radici et. al. @√s = 200 GeV

	/*Mh     lower    upper
	  0.36   0.0027   0.0076
	  0.5    0.0077   0.0132
	  0.69   0.0173   0.0286
	  0.88   0.0262   0.0438
	  1.19   0.0200   0.0360
	 */
	double thM200[5]={0.36, 0.5,  0.69, 0.88, 1.19};
	double thEl200[5]={0.0027, 0.0077, 0.0173, 0.0262, 0.02};
	double thEh200[5]={0.0076, 0.0132, 0.0286, 0.0438, 0.0360};

	double thA200[5]={0};
	double thE200[5]={0};
	double thEtest[5]={0};
	for(int i=0; i<5; i++){
		thA200[i]=0.5*(thEl200[i]+thEh200[i]);
		thE200[i]=thEh200[i]-thA200[i];
		thEtest[i]=(-1)*thEl200[i]+thA200[i];
		//cout<<thE200[i]<<" == "<<thEtest[i]<<endl;
	}

ofstream outfile;
outfile.open("intAutVsM200GeV.txt");
for(int i=0; i<9; i++){
if(i<5){
outfile<<setprecision(4)<<thM200[i]<<" & "<<thEl200[i]<<" & "<<thEh200[i]<<" & "<<M7[i]<<" & "<<asym7[i]<<" & "<<err7[i]<<" & "<<combSys[i]<<" & "<< run06M[i]<<" & "<<run06asym[i]<<" & "<<run06err[i]<<endl;
} else {
outfile<<setprecision(4)<<" & "<<" & "<<" & "<<M7[i]<<" & "<<asym7[i]<<" & "<<err7[i]<<" & "<<combSys[i]<<" & "<<" & "<<" & "<<endl;
}
} 




	TCanvas *myCan = new TCanvas("myCan","myCan",700,500);


	gStyle -> SetOptStat(0);
	gStyle -> SetOptDate(0);
	gStyle -> SetTitle("");
	gStyle->SetLegendBorderSize(0);
	//myCan -> Divide(2,1);
	//myCan -> cd(1);
	//gPad->SetPad(0.0,.35,1.,1.);
	gPad->SetGrid(0,0);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.12);

	TGraphErrors *th200 = new TGraphErrors(5,thM200,thA200,0,thE200);
	//TGraphErrors *th200 = new TGraphErrors(15,thM500,thA500,0,thE500);
	th200->SetTitle("");
	th200->SetFillStyle(1001);
	th200->SetFillColorAlpha(39,0.8);
	//th200->SetFillColor(kViolet);

	th200->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})  }}");
	th200->SetTitle("");
	th200->GetYaxis()->SetLabelSize(0.04);
	th200->GetYaxis()->SetLabelFont(22);
	th200->GetYaxis()->SetTitleSize(0.04);
	th200->GetYaxis()->SetTitleOffset(1.2);
	th200->GetYaxis()->SetRangeUser(-0.01,0.099);
	th200->GetXaxis()->SetLimits(0.21, 1.9);
	//th200->GetXaxis()->SetNdivisions(501);
	//th200->GetYaxis()->SetNdivisions(500);
	th200->GetXaxis()->SetTitle("#font[22]{M^{#pi^{+}#pi^{-}}_{inv}(GeV/c^{2})}");
	th200->GetXaxis()->SetTitleSize(0.04);
	th200-> GetXaxis()->SetLabelSize(0.04);
	th200-> GetXaxis()->SetLabelFont(22);
	th200->GetXaxis()->SetTitleOffset(1.1);
	th200->Draw("AE3");




	double errx[9]={0};
	double errxps[9]={0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015};
	double errxps11[7]={0.015,0.015,0.015,0.015,0.015,0.015,0.015};

	TGraphErrors *Run06= new TGraphErrors(5,run06M,run06asym,0,run06err);
	Run06-> SetMarkerStyle(20);
	Run06-> SetMarkerColor(4);
	Run06-> SetLineColor(4);
	Run06-> SetMarkerSize(1.2);
	Run06-> GetXaxis()->SetLimits(0.2,2.5);
	Run06->Draw("same P");
	
	TGraphErrors *Run15s = new TGraphErrors(9,M,A_tpc,errxps,combSys);
	Run15s->SetFillStyle(0);
	//Run15s->SetFillColorAlpha(kRed-5,0.8);
	Run15s->SetLineColor(2);
	Run15s->Draw("same 2");
	
	TGraphErrors *Run15 = new TGraphErrors(9,M,A_tpc,errx,Aerr_tpc);
	Run15-> SetMarkerStyle(20);
	Run15-> SetMarkerSize(1.2);
	Run15-> SetMarkerColor(kRed);
	Run15-> SetLineColor(kRed);
	//Run15-> SetLineWidth(2);
	Run15-> GetXaxis()->SetLimits(0.2,2.5);
	//Run15-> GetYaxis()->SetRangeUser(-0.01, 0.09);
	Run15->Draw("P same");
	myCan->Update();
	TLine *line1=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
	line1->SetLineStyle(2);
	line1->Draw();

	TLegend *legend=new TLegend(0.68,0.56,0.93, 0.94);
	//legend->AddEntry(th200, "Radici et. al., #sqrt{s} = 200 GeV","f");
	legend->AddEntry(th200, "Radici et. al.","f");
	legend->AddEntry(Run15, "Run 15, Cone < 0.7","lP");
	legend->AddEntry("", "#LT p_{T} #GT = 5.25 GeV/c","");
	legend->AddEntry(Run06, "Run 06, Cone < 0.3","lP");
	legend->AddEntry("", "#LT p_{T} #GT = 6 GeV/c","");
	legend->AddEntry(Run15s, " Syst. Error","f");
	legend->AddEntry("", "#eta^{#pi^{+}#pi^{-}} > 0 "," ");
	//legend->AddEntry(Run15s, "Run 15 Syst. Error","f");
	//legend->AddEntry(Run11s, "Run 11 Syst. Error","f");
	legend->SetTextSize(0.036);
	legend->SetTextFont(22);
	legend->Draw();
	/*
	TLegend *legend2=new TLegend(0.67,0.55,0.9, 0.75);
	legend2->AddEntry("","#eta^{#pi^{+}#pi^{-}} > 0","");
	legend2->AddEntry(Run15,"#LT p_{T} #GT = 5.25 GeV/c","P");
	legend2->AddEntry(Run06,"#LT p_{T} #GT = 6 GeV/c","P");
	legend2->SetTextSize(0.036);
	legend2->SetTextFont(22);
	//legend2->SetNColumns(2);
	legend2->Draw();
	*/
	TLatex latex;
	latex.SetTextFont(22);
	latex.SetTextSize(0.028);
	latex.SetTextAlign(13);  //align at top
	//latex.DrawLatex(1.42,0.06,"#color[2]{#font[22]{Run15, <p_{T}^{#pi^{+}#pi^{-}}> = 8.3 GeV/c, Cone < 0.7}}");
	//latex.DrawLatex(1.42,0.053,"#color[1]{#font[22]{Run11, <p_{T}^{#pi^{+}#pi^{-}}> = 13 GeV/c}, Cone < 0.7}");
	//latex.DrawLatex(1.42,0.046,"#color[4]{#font[22]{Run06, <p_{T}^{#pi^{+}#pi^{-}}> = 6 GeV/c}, Cone < 0.3}");
	latex.SetTextSize(0.055);
	latex.DrawLatex(0.27,0.096,"#color[2]{#font[22]{STAR Preliminary 2015}}");
	latex.SetTextSize(0.04);
	latex.DrawLatex(0.27,0.0895,"#color[1]{#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X} at #sqrt{s} = 200 GeV}");
	//latex.DrawLatex(0.27,0.08,"#color[1]{#font[12]{#eta^{#pi^{+}#pi^{-}} > 0 }}");
	latex.SetTextSize(0.03);
	latex.DrawLatex(0.6,-0.003,"#color[1]{#font[22]{#pm 3% scale uncertainty from beam polarization (not shown)}}");
	myCan->Update();
	myCan->Update();



	myCan->SaveAs("./Plots/AutvsMIntPt_AfterPID.pdf");
	//myCan->SaveAs("STARAut200GeV_Wsys.eps");
	//myCan->SaveAs("STARAut200GeV_Wsys.png");



}
