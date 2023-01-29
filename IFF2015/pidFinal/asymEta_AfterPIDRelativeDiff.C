#include <iostream>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include <iomanip> 
#include <cmath> 

using namespace std;

void asymEta_AfterPIDRelativeDiff(){
//for X and Z as a function of eta
 TFile *geant=new TFile("GEANT_combinedTrees_Sigmam1to2.root");

TProfile *hx1geant = (TProfile*)geant->Get("hx1vsEta");
TProfile *hx2geant = (TProfile*)geant->Get("hx2vsEta");

TProfile *hz1geant = (TProfile*)geant->Get("hz1vsEta");
TProfile *hz2geant = (TProfile*)geant->Get("hz2vsEta");

TProfile *avgXgeant = (TProfile*)hx1geant->Clone();
avgXgeant->Add(hx2geant);

TProfile *avgZgeant = (TProfile*)hz1geant->Clone();
avgZgeant->Add(hz2geant);

double avg_x[9]={0};
double avg_xerr[9]={0};
double errx[9]={0};
double avg_z[9]={0};
double avg_zerr[9]={0};
double binCenter[9]={0};
ofstream fmc;
fmc.open("avgXandZ.txt");
fmc<<"//average x and z values in 9 eta bins from MC simulation."<<endl;
fmc<<"double avg_x[9]={";
//draw average x and z in the same pad 
TH1F* hx=(TH1F*)avgXgeant->Clone();
TH1F* hz=(TH1F*)avgZgeant->Clone();
for(int i=0; i<9; i++){
 binCenter[i]=hx->GetBinCenter(i+1);
 avg_x[i]=hx->GetBinContent(i+1);
 if(i==8){
  fmc<<avg_x[i]<<"};"<<endl;
 }else{
  fmc<<avg_x[i]<<",";
 }
 avg_xerr[i]=hx->GetBinError(i+1);
 avg_z[i]=hz->GetBinContent(i+1);
 avg_zerr[i]=hz->GetBinError(i+1);
 }
fmc<<"double avg_xerr[9]={";
for(int i=0; i<9; i++){
 if(i==8){
  fmc<<avg_xerr[i]<<"};"<<endl;
 }else{
  fmc<<avg_xerr[i]<<",";
 }
}
fmc<<"double avg_z[9]={";
for(int i=0; i<9; i++){
 if(i==8){
  fmc<<avg_z[i]<<"};"<<endl;
 }else{
  fmc<<avg_z[i]<<",";
 }
}
fmc<<"double avg_zerr[9]={";
for(int i=0; i<9; i++){
 if(i==8){
  fmc<<avg_zerr[i]<<"};"<<endl;
 }else{
  fmc<<avg_zerr[i]<<",";
 }
}

//A_UT vs Eta----Tpc only
double Eta[9]={-0.78988, -0.567138,-0.374526,-0.188095,-0.00349934,0.181359, 0.371768,0.571623, 0.797184};
double A_tpc[9]={0.00237124, 0.00178508,0.00332655,0.00282821,0.00457845, 0.00754294, 0.0121409, 0.0163978, 0.0241117};
double Aerr_tpc[9]={0.00109241, 0.00109075,0.00109383,0.00109019,0.00108671, 0.00109823, 0.00111833, 0.00113991, 0.00119084};
//---------------
//A_UT vs Eta Tpc Or Tof
double A_F[9]={0.00270787, 0.0019989,0.00356342,0.00266537,0.00542049, 0.0082484, 0.0130979, 0.0165652, 0.0242299};
double Aerr_F[9]={0.00111228, 0.00112171,0.00113351,0.00113446,0.00113317, 0.00114367, 0.00116181, 0.00117042, 0.00121023};
//--------------
//A_UT vs Eta Tof only
double A_tof[9]={-0.000459428, -0.000506288,0.00182886,0.00229372,0.00339382, 0.00698642, 0.00957079, 0.0168978, 0.0248798};
double Aerr_tof[9]={0.00185453, 0.00166434,0.0016759,0.00173602,0.00176355, 0.00174082, 0.00168958, 0.00170398, 0.00191897};
//-------------

//flat rate
double ratioEta[9]={1.184,1.184,1.184,1.184,1.184,1.184,1.184,1.184,1.184};
double ratioEtaErr[9]={0.0626708,0.0531546,0.0418376,0.0523667,0.0605381,0.0520935,0.0694973,0.0488803,0.041424};

//priliminary fractions-----
double pairpurityPtGt[5]={0.67307,0.817963,0.794931,0.765595,0.793264};
double pairpurityPtLt[5]={0.660449,0.8359,0.808052,0.739741,0.80209};
double pairpurityMGt[5]={0.524506,0.767845,0.81264,0.827867,0.842939};
double pairpurityMLt[5]={0.528519,0.779716,0.83881,0.83958,0.848677};
double pairpurityEta[9]={0.777056,0.798529,0.800554,0.803539,0.802013,0.798089,0.79153,0.784865,0.760883};
//--------------------------
//PID fractions TOF only
double pairpurityEtaT[9]={0.97593,0.98328,0.982627,0.978905,0.975671,0.978282,0.983113,0.983608,0.976228};
//---------------------
//calculate difference in asymmetry 
double rdiff[9]={0};
 double rdifferr[9]={0};
 double xerr[9]={0};
 for(int i=0; i<9; i++){
 rdiff[i]=(A_tpc[i]-A_tof[i])/A_tpc[i];
 rdifferr[i]=sqrt(pow((A_tof[i]/A_tpc[i])*Aerr_tpc[i],2)+pow(Aerr_tof[i]/A_tpc[i],2));
 }
 TGraphErrors *gr_diff=new TGraphErrors(9,Eta,rdiff,xerr,rdifferr);
 gr_diff->SetTitle("");
 gr_diff->GetYaxis()->SetLabelSize(0.04);
 gr_diff->GetYaxis()->SetLabelFont(22);
 gr_diff->GetYaxis()->SetTitle("A_{UT}^{Tpc} - A_{UT}^{Tpc&&Tof}");
 gr_diff->GetYaxis()->SetTitleSize(0.04);
 gr_diff->GetYaxis()->SetTitleOffset(1.2);
 gr_diff->GetYaxis()->SetRangeUser(-2,2);
 gr_diff->GetXaxis()->SetLimits(-1, 1.);
 gr_diff->GetXaxis()->SetTitle("#font[22]{#eta^{#pi^{+}#pi^{-}}}");
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
 tex.DrawLatex(-0.5,1.8,Form("#font[22]{#color[2]{chi2/ndf = %g / %i}}",gr_diff->GetFunction("pol0")->GetChisquare(), gr_diff->GetFunction("pol0")->GetNDF()));
 tex.DrawLatex(-0.5,1.5,Form("#font[22]{#color[2]{p0 = %g #pm %g}}",gr_diff->GetFunction("pol0")->GetParameter(0), gr_diff->GetFunction("pol0")->GetParError(0)));
 double fp0=gr_diff->GetFunction("pol0")->GetParameter(0);
 cdiff->Update();
 cdiff->SaveAs("./Plots/intEta_Difference.pdf");
 




//theory from skoby 
double Mth[6]={0.298412698412699,     0.50188679245283,       0.703144654088051,      0.89433962264151,       1.0994708994709,        1.3308176100629};
  double Ath[6]={0.00925619834710749,   0.0148760330578513,     0.0401652892561984,     0.047603305785124,      0.0285950413223141,     0.0204958677685951};
  double dAth[6]={0.00297520661157025,  0.0028099173553719,     0.00826446280991735,    0.0135537190082645,     0.0112396694214876,     0.00892561983471074};


//bias calculation
double biasM1[9]={0};
double biasM2[9]={0};
double biasM3[9]={0};
double biasM4[9]={0};
double biasM5[9]={0};

double biaspT1[9]={0};
double biaspT2[9]={0};
double biaspT3[9]={0};
double biaspT4[9]={0};
double biaspT5[9]={0};

double biasEta[9]={0};

double ex[9]={0};
double exs[9]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
double eyl[9]={0};
double eyh[9]={0};
double pideyl[9]={0};
double pideyh[9]={0};
double pidEta[9]={0};
double combSys[9]={0};
for (int pbin=0; pbin<9; pbin++){
 //bias for A_UT vs Eta
  biasEta[pbin]=(1-ratioEta[pbin])*max(A_tpc[pbin],Aerr_tpc[pbin]);
  pidEta[pbin]=fp0*A_tpc[pbin];
  combSys[pbin]=sqrt(pow(biasEta[pbin],2) + pow(pidEta[pbin],2));
}

for (int pbin=0; pbin<9; pbin++){
if(A_F[pbin]>0){eyl[pbin]=abs(biasEta[pbin]); pideyh[pbin]=(1-pairpurityEta[pbin])*abs(A_F[pbin]);}
if(A_F[pbin]<0){eyh[pbin]=abs(biasEta[pbin]);pideyl[pbin]=(1-pairpurityEta[pbin])*abs(A_F[pbin]);}
}
ofstream table;
table.open("table_AutvsEta.text");
table<<"<#eta>"<<"\t "<<"A_{ut}"<<"\t "<<"#sigma_{stat}\t"<<"#sigma_{pid}\t "<<"#sigma_{bias}\t"<<endl;
for(int i=0; i<9; i++){
table<<setprecision(4)<<Eta[i]<<" & "<<avg_x[i]<<" & "<<avg_z[i]<<" & "<< A_F[i]<<" & "<<Aerr_F[i]<<" & "<<combSys[i]<<endl;
}
/*
for (int pbin=0; pbin<9; pbin++){
if(ratioEta[pbin]>1){eyl[pbin]=abs(biasEta[pbin]);}
if(ratioEta[pbin]<1){eyh[pbin]=abs(biasEta[pbin]);}
}
*/
gStyle->SetLegendBorderSize(0);
const int n=9;
 //TCanvas *myCanA = new TCanvas("myCanA","myCanA",150,10,990,660);
 TCanvas *myCanA = new TCanvas("myCanA","myCanA",600,500);
gStyle -> SetOptStat(0);
//gStyle->SetLegendBorderSize(0);
myCanA -> Divide(2,1);
myCanA -> cd(2);
gPad->SetPad(0.0,0.0,1,0.35);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.3);
gPad->SetTopMargin(0.0);
gPad->SetGrid(0,0);
TGraphErrors *gZ = new TGraphErrors(n,Eta,avg_z,errx,avg_zerr);
gZ->SetTitle("");
gZ->GetYaxis()->SetRangeUser(0,0.65);
gZ->GetXaxis()->SetLimits(-1.,1.);
gZ->SetMarkerStyle(21);
gZ->SetMarkerColor(2);

gZ->GetXaxis()->SetTitle("#eta^{#pi^{+}#pi^{-}}");
gZ->GetXaxis()->SetTitleOffset(0.8);
gZ->GetXaxis()->SetTickLength(0.1);
gZ->GetYaxis()->SetTickLength(0.04);
gZ->GetYaxis()->SetNdivisions(505);
//gZ->GetXaxis()->SetNdivisions(505);
gZ->GetXaxis()->SetTitleFont(22);
gZ->GetXaxis()->SetLabelFont(22);
gZ->GetXaxis()->SetLabelSize(0.12);
gZ->GetXaxis()->SetTitleSize(0.15);
gZ->GetYaxis()->SetLabelSize(0.13);
gZ->GetYaxis()->SetLabelFont(22);
//gZ->GetYaxis()->SetNdivisions(5, 5,0,kTRUE);
gZ->Draw("AP");

gPad->Update();
TGraphErrors *gX = new TGraphErrors(n,Eta,avg_x,errx,avg_xerr);
gX->SetMarkerStyle(20);
gX->SetMarkerColor(4);
gX->Draw("P same");

TLegend *leg1=new TLegend(0.3,0.53,0.5,0.65);
leg1->SetNColumns(2);
leg1->AddEntry(gZ,"#font[22]{#color[2]{< z >}}", "P");
leg1->AddEntry(gX,"#font[22]{#color[4]{< x >}}", "P");
leg1->SetTextSize(0.12);
leg1->Draw();
gPad->Update();

myCanA -> cd(1);
gPad->SetPad(0.0,0.35,1.,1.);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.00);
gPad->SetGrid(0,0);

TGraphErrors *syst = new TGraphErrors(n, Eta, A_tpc, exs,combSys);//Forward average
  syst->SetTitle("");
syst->GetYaxis()-> SetTitle("A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}  ");
syst->GetYaxis()-> SetTitleFont(22);
syst->SetTitle("");
syst->GetYaxis()->SetTitleOffset(1.);
syst->GetYaxis()->SetRangeUser(-0.005,0.035);
syst->GetYaxis()->SetNdivisions(5, 5,0,kTRUE);
syst->GetXaxis()->SetLabelSize(0);

syst->GetYaxis()->SetTitleSize(0.07);
syst->GetXaxis()->SetTickLength(0.06);
syst->GetYaxis()->SetLabelSize(0.07);
syst->GetYaxis()->SetLabelFont(22);
syst->SetFillStyle(0); 
syst->SetLineColor(1); 
syst->SetLineWidth(1); 
syst->Draw("AE2");


TGraphErrors *data = new TGraphErrors(n, Eta, A_tpc, ex, Aerr_tpc);//Forward average
data-> SetMarkerStyle(20);
data-> SetMarkerColor(2);
data-> SetLineColor(2);
data->Draw("P same ");
myCanA->Update();
TLine *lineA=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
lineA->SetLineStyle(2);
lineA->Draw();

TLatex text;
text.SetTextSize(0.07);
text.SetTextFont(22);
text.DrawLatex(-0.905, 0.032,"#font[22]{#color[2]{STAR Preliminary 2015}}");
text.SetTextSize(0.06);
text.DrawLatex(-0.905, 0.0285, "#font[22]{ p^{#uparrow}+ p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 200 GeV}");
text.DrawLatex(-0.77, -0.0022, "#scale[0.8]{#pm 3% scale uncertainty from beam polarization (not shown)}");
//gStyle->SetLegendBorderSize(0);
TLegend *leg = new TLegend(0.18,0.64, 0.38, 0.71);
leg->AddEntry(syst," #font[22]{Syst. Error}","f");
//leg->AddEntry("","cone > 0.7","");
leg->SetTextSize(.06);
leg->Draw();
myCanA->Update();
myCanA->SaveAs("./Plots/AutVsEta_AfterPID.pdf");
//myCanA->SaveAs("moneyplot_AutVsEtaXZ_CombSys.png");
//myCanA->SaveAs("moneyplot_AutVsEtaXZ_CombSys.eps");
}
