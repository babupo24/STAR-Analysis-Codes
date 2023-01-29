//This macro estimates the particle fractions in respective asymmetry bins. 
//The structure of this macro is as follows:
//1. Reads histograms file for PID, which contains all histograms required 
//   for particle fraction calculation in all asymmetry bins. 
//2. In 5 pT bins, A_UT as a function of Minv in 9 bins is calculated. 
//   I have calculated particle fractions in 5 those pT bins, which will be applied in all 9 minv bins within each 
//   pT bins.     
//3. In 5 Minv bins, A_UT as a function of pT in 9 bins is calculated. 
//   I have calculated particle fractions in those 5 Minv bins, which will be applied in all 9 pT bins within each 
//   Minv bins.
//4. A_UT is extracted in 9 eta bins, integrated over pT and Minv. Particle fractions is calculated in all 9 eta 
//   bins. 
//5. A_UT is extracted in 9 Minv bins, integrated over pT. Particle fractions is calculated in all 9 Minv bins. 
//
//General idea of PID calculation in each asymmetry bins: 
//1. Read 2D hist- nSigmaPion vs nSigmaProton/Kaon/Electron to estimate realtive position of proton/kaon/electron
//   from pion peak position. 
//2. Pion mean is estimated using TOF independently. 
//3. Multigaussian function composite of 3 gaussians: for pion, kaon+proton,  and electron.
//   Since proton and kaon are inseparable, both are accounted by a single gaussian. Koan peak is used for gaussian
//   mean for kaon+ proton.  
//4. nSigmaPion distribution is fitted with multigaussian fit. In the multigaus fit, the means are fixed and the 
//   width are varied for the best fit convergence.        
//5. Particle fraction is calculated as: 
//         f = particle_fit->Integral(-1,2)/total_fit->Integral(-1,2); (-1,2)
//6. Fractions are calculated in two eta bins: eta>0, eta<0 and for positive and negative charge separately. 
//7. Since this analysis uses π+π- pair, the fraction of which is needed. 
//8. The π+π- pair fraction is calculated as a product of π+ fraction and π- fraction in respective A_UT and 
//   eta bins. 
//9. All the numerical values are printed in text files. 
//10.All plots are saved in "Plots", text files are saved in "textFiles" folders. 
// 
//How to run? 
// 1. have your histogram file in Histos folder.
// 2. mkdir Plots textFiles
// 3. root4star piddEdxFinal.C
//////----------------------------------------------------------------------------

#include <iostream>
#include "TH1.h" 
#include "TH2.h" 
#include "TCanvas.h" 
#include "TTree.h" 
#include "TFile.h" 
#include <iomanip> 

using namespace std;
void piddEdxFinal(){
 TFile *infile = new TFile("Histos/pidHist_StTOFWithTOFCut_02.root","R");

 const int nBins=9;
 const char *etaDir[2]={"Gt","Lt"};
 const char *etas[2]={"> 0","< 0"};
 const char *charge[2]={"Pos","Neg"};
 const char *charges[2]={"+ve","-ve"};
 const char *ch[2]={"+","-"};
 const char *kbins[2]={"pTbin","mBin"};
 const char *spec[4]={"Pion","Proton","Kaon","Electron"};
const char *part[3]={"P","K","E"};
double avgPt[5]={3.4, 4.0, 4.7, 5.7, 8.3};
double avgM[5]={0.39, 0.58, 0.76, 0.96, 1.42};
 double pT[6]={2.80,   3.71,  4.303 ,   5.084,   6.404,  15.00}; //pT bin range for A_UT vs Minv
 double M[6]={0.20, 0.4795, 0.6714, 0.8452, 1.1100, 4.00};// Minv bin range for A_UT vs pT
 double eta_range[10]={-1.00, -0.668, -0.469, -0.281, -0.096, 0.089, 0.275, 0.470, 0.675, 1.00};//eta bins for A_UT vs eta
double avgEta[9]={-0.789853, -0.567055,-0.374604,-0.188188,-0.00350474,0.181505, 0.37192,0.571621, 0.797284};

//integrated binning
double pT_in[10]= {2.800, 3.470,  3.790, 4.120, 4.490, 4.938,  5.505, 6.300, 7.660, 15.00 };
double M_in[10]= {0.250, 0.403, 0.516,   0.612,  0.711,  0.803,  0.921,  1.070,  1.286,  4.000};
double eta_in[10]={-1.00, -0.668, -0.469, -0.281, -0.096, 0.089, 0.275, 0.470, 0.675, 1.00};
double avgIntM[9] = {0.3487, 0.453749, 0.564329, 0.662024, 0.756702, 0.874893, 1.03266, 1.2362, 1.68868};

//find proton, kaon, electron peak positions in the pion signal region.
 double peakPt[3][2][2][5];  double sigmaPt[3][2][2][5];
 double peakM[3][2][2][5];   double sigmaM[3][2][2][5];
 double peakEta[3][2][9];    double sigmaEta[3][2][9];
 double peakIntM[3][2][2][9];    double sigmaIntM[3][2][9];
//pion mean and width from K^0_s-> π^+π^-
 //Pion Mean
Double_t piPosGtM[5]={0.00129139,-0.00850172,-0.0241273,-0.0353401,-0.0445746};
Double_t piPosLtM[5]={0.00746207,0.0174659,-0.000166883,-0.00153431,-0.0151362};
Double_t piNegGtM[5]={-0.0288659,-0.0114434,-0.0171198,-0.0148943,-0.0170946};
Double_t piNegLtM[5]={0.0288481,0.0175343,0.0163743,0.0259707,0.00641877};
//Pion Width
Double_t piPosGtW[5]={0.977472,0.986781,0.982175,0.978301,0.995949};
Double_t piPosLtW[5]={0.954947,0.951614,0.932312,0.951007,0.951007};
Double_t piNegGtW[5]={0.979396,0.977307,0.989981,0.974688,0.971426};
Double_t piNegLtW[5]={0.952709,0.949421,0.938136,0.937794,0.938797};
///----------------------------------------------------------
//pion maen and width pt bins  using tof
double piMeanPt[5]={-0.0362011,-0.0612811,-0.0686669,-0.0676157,-0.0878517};
double piWidthPt[5]={1.0176,1.01859,1.00694,0.992572,0.976291};
//pion mean and width eta bins
double piMeanEta[9]={-0.0423608,-0.0460702,-0.0272924,0.00520594,0.0237695,0.00194142,-0.0429511,-0.0790397,-0.0900327};
double piWidthEta[9]={0.94073,0.957204,0.96034,0.956218,0.959724,0.977369,1.00071,1.00447,0.994622};
//pion mean and width Minv bins
double piMeanM[5]={-0.0455023,-0.033907,-0.0274415,-0.034777,-0.0732631};
double piWidthM[5]={0.991953,0.989967,0.983732,0.985867,0.992272};
//--------------------------------------------------------------
  //for pT-bins
  TH2D *hist2DpT[3][2][2][5];
   TH1D *hprofilepT[3][2][2][5];
  //for M bins 
  TH2D *hist2DM[3][2][2][5];
  TH1D *hprofileM[3][2][2][5];

for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
    for(int netaDir=0; netaDir<2; netaDir++){
      for(int nbin=0; nbin<5; nbin++){
      hprofilepT[npart][ncharge][netaDir][nbin]=new TH1D(Form("%s%s%s%i_pTbin",part[npart],charge[ncharge],etaDir[netaDir],nbin),"",100,-1,1);
      hprofileM[npart][ncharge][netaDir][nbin]=new TH1D(Form("%s%s%s%i_Mbin",part[npart],charge[ncharge],etaDir[netaDir],nbin),"",100,-1,1);
      }
    }
  }
 }
//Get histograms from the root file 

 for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
    for(int netaDir=0; netaDir<2; netaDir++){
      for(int nbin=0; nbin<5; nbin++){

         //pTbins 
       hist2DpT[npart][ncharge][netaDir][nbin]=(TH2D*)infile->Get(Form("hPivs%s%s%s_ptBin%i",part[npart], charge[ncharge],etaDir[netaDir],nbin));
       hprofilepT[npart][ncharge][netaDir][nbin]=hist2DpT[npart][ncharge][netaDir][nbin]->ProfileX(Form("%s%s%s%ipT",part[npart],charge[ncharge],etaDir[netaDir],nbin),1,-1,"s");
       hprofilepT[npart][ncharge][netaDir][nbin]->Fit("pol3");
       peakPt[npart][ncharge][netaDir][nbin]=hprofilepT[npart][ncharge][netaDir][nbin]->GetFunction("pol3")->GetParameter(0);
       sigmaPt[npart][ncharge][netaDir][nbin]=hprofilepT[npart][ncharge][netaDir][nbin]->GetBinError((int)hprofilepT[npart][ncharge][netaDir][nbin]->GetXaxis()->FindBin(0.0));	
cout<<"sigma pT: particle: "<<part[npart]<<", charge: "<<charge[ncharge]<<" eta: "<< etaDir[netaDir]<<", pad: "<<nbin<<" "<<sigmaPt[npart][ncharge][netaDir][nbin]<<endl;
       //Mbins
       hist2DM[npart][ncharge][netaDir][nbin]=(TH2D*)infile->Get(Form("hPivs%s%s%s_mBin%i",part[npart], charge[ncharge],etaDir[netaDir],nbin));
       hprofileM[npart][ncharge][netaDir][nbin]=hist2DM[npart][ncharge][netaDir][nbin]->ProfileX();
       hprofileM[npart][ncharge][netaDir][nbin]->Fit("pol3");
       peakM[npart][ncharge][netaDir][nbin]=hprofileM[npart][ncharge][netaDir][nbin]->GetFunction("pol3")->GetParameter(0);	

      }
    }
  }
}

//for integrated eta and integrated pT bins
 TH2D* histEta[3][2][9];
 TH1D* hprofileEta[3][2][9];
 for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
      for(int nbin=0; nbin<9; nbin++){
       histEta[npart][ncharge][nbin]=(TH2D*)infile->Get(Form("hPivs%s%s_etaBin%i",part[npart],charge[ncharge],nbin));
       hprofileEta[npart][ncharge][nbin]=histEta[npart][ncharge][nbin]->ProfileX();
       hprofileEta[npart][ncharge][nbin]->Fit("pol3");
       peakEta[npart][ncharge][nbin]=hprofileEta[npart][ncharge][nbin]->GetFunction("pol3")->GetParameter(0);
      }
    }
 }

 //for integrated p_T bins  added 10/29/2021
 TH2D* histIntM[3][2][2][9];
 TH1D* hprofileIntM[3][2][2][9];
 TH1D *hnsigmaIntM[3][2][2][9];//species(4 - Pi, P, K, E),charge(2  +/-), eta direction (2 ><0), number of bins(5)
 for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
    for(int netaDir=0; netaDir<2;netaDir++){
      for(int nbin=0; nbin<9; nbin++){
       histIntM[npart][ncharge][netaDir][nbin]=(TH2D*)infile->Get(Form("inthPivs%s%s%s_mBin%i",part[npart],charge[ncharge],etaDir[netaDir],nbin));
       hprofileIntM[npart][ncharge][netaDir][nbin]=histIntM[npart][ncharge][netaDir][nbin]->ProfileX();
       hprofileIntM[npart][ncharge][netaDir][nbin]->Fit("pol3");
       peakIntM[npart][ncharge][netaDir][nbin]=hprofileIntM[npart][ncharge][netaDir][nbin]->GetFunction("pol3")->GetParameter(0);
     hnsigmaIntM[npart][ncharge][netaDir][nbin]=(TH1D*)infile->Get(Form("inthsigma%s%s_%s_mBin%i",spec[npart],charge[ncharge],etaDir[netaDir],nbin));

      }
    }
  }
 }
//---------------------
ofstream fpeak;
fpeak.open("textFiles/peakPosition.txt");
fpeak<<"npart: 0 - proton      1 - Kaon     2 - Electron"<<endl;
fpeak<<"ncharge: 0 - +ve      1 - -Ve"<<endl;
fpeak<<"netaDir: 0 - >0      1 - <0"<<endl;
fpeak<<"###pT Bins: "<<endl;

 for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
    for(int netaDir=0; netaDir<2; netaDir++){
      for(int nbin=0; nbin<5; nbin++){
	fpeak<<npart<<" "<<ncharge<<" "<<netaDir<<" "<<nbin<<" "<<peakPt[npart][ncharge][netaDir][nbin]<<endl;	
      }
    }
  }
}

	fpeak<<"###M Bins: "<<endl;
 for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
    for(int netaDir=0; netaDir<2; netaDir++){
      for(int nbin=0; nbin<5; nbin++){
	fpeak<<npart<<" "<<ncharge<<" "<<netaDir<<" "<<nbin<<" "<<peakM[npart][ncharge][netaDir][nbin]<<endl;	
      }
    }
  }
}

	fpeak<<"###Eta Bins: "<<endl;
 for(int npart=0; npart<3; npart++){
  for(int ncharge=0; ncharge<2;ncharge++){
      for(int nbin=0; nbin<9; nbin++){
	fpeak<<npart<<" "<<ncharge<<" "<<nbin<<" "<<peakEta[npart][ncharge][nbin]<<endl;	
      }
    }
  }
//end peak finder
//////////////////////////////////////////////////////


//multi-gaussian fit on nSigmaPion distribution for the particle fractions
 TH1D *hnsigmaM[4][2][2][5];//species(4 - Pi, P, K, E),charge(2  +/-), eta direction (2 ><0), number of bins(5)
  TH1D *hnsigmapT[4][2][2][5];//species(4 - Pi, P, K, E),charge(2  +/-), eta direction (2 ><0), number of bins(5)
  TH1D *hnsigmaEta[4][2][9];//species(4 - Pi, P, K, E),charge(2  +/-),  number of bins(9)
 
for(int nspec=0; nspec<4;nspec++){
  for(int ncharge=0; ncharge<2;ncharge++){
   for(int netaDir=0; netaDir<2; netaDir++){
    for(int nbin=0; nbin<5; nbin++){
     hnsigmaM[nspec][ncharge][netaDir][nbin]=(TH1D*)infile->Get(Form("hsigma%s%s_%s_mBin%i",spec[nspec],charge[ncharge],etaDir[netaDir],nbin));
     hnsigmaM[nspec][ncharge][netaDir][nbin]->Rebin(1);
     //hnsigmaM[nspec][ncharge][netaDir][nbin]->Sumw2();
     hnsigmapT[nspec][ncharge][netaDir][nbin]=(TH1D*)infile->Get(Form("hsigma%s%s_%s_pTbin%i",spec[nspec],charge[ncharge],etaDir[netaDir],nbin));
     hnsigmapT[nspec][ncharge][netaDir][nbin]->Rebin(1);
     //hnsigmapT[nspec][ncharge][netaDir][nbin]->Sumw2();
    }
   }
  }
 }

for(int nspec=0; nspec<4;nspec++){
  for(int ncharge=0; ncharge<2;ncharge++){
    for(int nbin=0; nbin<9; nbin++){
     //hnsigmaEta[nspec][ncharge][nbin]=(TH1D*)infile->Get(Form("hsigma%s%s_Etabin%i",spec[nspec],charge[ncharge],nbin));
     hnsigmaEta[nspec][ncharge][nbin]=(TH1D*)infile->Get(Form("hSigma%s%s_etaBin%i",spec[nspec],charge[ncharge],nbin));
     hnsigmaEta[nspec][ncharge][nbin]->Rebin(1);
     //hnsigmaEta[nspec][ncharge][nbin]->Sumw2();

   }
  }
 }//integrated binning
  gStyle->SetStatW(0.21); gStyle->SetStatH(0.17); 
  //gStyle->SetStatX(1.05); //gStyle->SetStatH(0.16); 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptDate(0);
  gStyle->SetLegendBorderSize(0);
 TF1 *fitFuncPt;//fit each species
 TF1 *fitFuncM;
 TF1 *fitFuncIntM;
 TF1 *fitFuncEta;
const char *percent="%";
//for signal fractions
 double pionIntMbin[2][2][9]={0};
 double protonIntMbin[2][2][9]={0};
 double electronIntMbin[2][2][9]={0};
 
 double pionPtbin[2][2][5]={0};       double pionMbin[2][2][5]={0};       double pionEtabin[2][9]={0};    
 double kaonPtbin[2][2][5]={0};       double kaonMbin[2][2][5]={0};       double kaonEtabin[2][9]={0};    
 double protonPtbin[2][2][5]={0};     double protonMbin[2][2][5]={0};     double protonEtabin[2][9]={0};  
 double electronPtbin[2][2][5]={0};   double electronMbin[2][2][5]={0};   double electronEtabin[2][9]={0};

 double pionPtbinErr[2][2][5]={0};       double pionMbinErr[2][2][5]={0};       double pionEtabinErr[2][9]={0};    
 double kaonPtbinErr[2][2][5]={0};       double kaonMbinErr[2][2][5]={0};       double kaonEtabinErr[2][9]={0};    
 double protonPtbinErr[2][2][5]={0};     double protonMbinErr[2][2][5]={0};     double protonEtabinErr[2][9]={0};  
 double electronPtbinErr[2][2][5]={0};   double electronMbinErr[2][2][5]={0};   double electronEtabinErr[2][9]={0};


 double fparsPt[11];
 double fparsM[11];
 char *funcName;
 TLegend *legend; 
 TCanvas *cvs[2][2][2];//2-M/pT bin, 2-charge +/-, 2- eta direction ><0
 
 for(int nkbins=0; nkbins<2; nkbins++){
   for(int ncharge=0; ncharge<2; ncharge++){
     for(int netaDir=0; netaDir<2; netaDir++){
      //if(nkbins!=0 && ncharge!=0 && netaDir!=0) continue; //for fit test, only for pion +ve eta>0
      cvs[nkbins][ncharge][netaDir]=new TCanvas(Form("cvs%s%s%s",kbins[nkbins],charge[ncharge],etaDir[netaDir]),"",1200,900);
   cvs[nkbins][ncharge][netaDir]->SetFrameFillStyle(0);
   cvs[nkbins][ncharge][netaDir]->Divide(3,2);
  cvs[nkbins][ncharge][netaDir]->Print("Plots/pionfitgaus_pTMBin.pdf[");  
  if(nkbins==0){ //pT_bins
   ofstream outputPtbin;
   outputPtbin.open(Form("textFiles/purityPtbin_%s%s.txt",charge[ncharge],etaDir[netaDir]));
   for(int pad=0; pad<6; pad++){
    if( pad==5) continue; 

    cvs[nkbins][ncharge][netaDir]->cd(pad+1);
    cout<<"pT Bin, charge: "<<charge[ncharge]<<" , eta: "<<etas[netaDir]<<", Pad:  "<<pad<<endl;

    gPad->SetGrid(0,0);
    gPad->SetLogy();
    hnsigmapT[0][ncharge][netaDir][pad]->SetLineColor(1);
    hnsigmapT[0][ncharge][netaDir][pad]->GetXaxis()->SetTitle("nSigmaPion"); 
    hnsigmapT[0][ncharge][netaDir][pad]->SetTitle(Form("%g < p_{T} < %g, #eta %s, charge = %s ", pT[pad], pT[pad+1], etas[netaDir], charges[ncharge]));
    //hnSigmaTotalpT[ncharge][netaDir][pad]->Sumw2(kFALSE);
    hnsigmapT[0][ncharge][netaDir][pad]->Draw();
     //fitFuncPt=new TF1("fitFuncPt","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])+[9]*TMath::Gaus(x,[10],[11])",-10,7);         		
     fitFuncPt=new TF1("fitFuncPt","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])",-10,7);         		
     //fitFuncPt=new TF1("fitFuncPt","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",-10,7);         		
     fitFuncPt->SetLineColor(1);
     fitFuncPt->SetLineWidth(2);

     //proton+kaon
     //fitFuncPt->SetParameter(0,aproton*0.1);
     //fitFuncPt->SetParLimits(0,aproton*0.01, aproton*10);
     fitFuncPt->SetParameter(1,peakPt[1][ncharge][netaDir][pad]);
     fitFuncPt->SetParLimits(1,peakPt[1][ncharge][netaDir][pad],peakPt[1][ncharge][netaDir][pad]);
     fitFuncPt->SetParameter(2,1);
     fitFuncPt->SetParLimits(2,0.8, 1.3);
     fitFuncPt->SetParName(0, "A_{K+p}");
     fitFuncPt->SetParName(1, "#mu_{K+p}");
     fitFuncPt->SetParName(2, "#sigma_{K+p}");


     //pion
     fitFuncPt->SetParameter(4,piMeanPt[pad]);
     fitFuncPt->SetParLimits(4,piMeanPt[pad],piMeanPt[pad]);
     //fitFuncPt->FixParameter(5,piWidthPt[pad]);
     fitFuncPt->SetParameter(5,1);
     fitFuncPt->SetParLimits(5,0.8,1.2);

     fitFuncPt->SetParName(3, "A_{#pi}");
     fitFuncPt->SetParName(4, "#mu_{#pi}");
     fitFuncPt->SetParName(5, "#sigma_{#pi}");

     //electron
     //fitFuncPt->SetParameter(9,aelectron*0.1);
     //fitFuncPt->SetParLimits(9,aelectron*0.01,aelectron*10);
     fitFuncPt->SetParameter(7,peakPt[2][ncharge][netaDir][pad]);
     fitFuncPt->SetParLimits(7,peakPt[2][ncharge][netaDir][pad], peakPt[2][ncharge][netaDir][pad]);
     fitFuncPt->SetParameter(8,1.5);
     fitFuncPt->SetParLimits(8,1,2.2);

     fitFuncPt->SetParName(6, "A_{e}");
     fitFuncPt->SetParName(7, "#mu_{e}");
     fitFuncPt->SetParName(8, "#sigma_{e}");

    /* //pileup
     //fitFuncPt->SetParameter(9,apileup);
     //fitFuncPt->SetParLimits(9,apileup*0.1,apileup*10);
     fitFuncPt->SetParameter(10,6.5);
     fitFuncPt->SetParLimits(10,6,10);
     fitFuncPt->SetParameter(11,2);
     fitFuncPt->SetParLimits(11,1.5,3);
    */
     //fitFuncPt->SetParNames("A_{#pi}","#mu_{#pi}","#sigma_{#pi^{+}}","A_{P}","#mu_{p}","#sigma_{P}","A_{K}","#mu_{K}","#sigma_{K}","A_{e}","#mu_{e}");
    /* //fitFuncPt->SetParName(11, "#sigma_{e}");
     fitFuncPt->SetParName(9, "A_{pileup}");
     fitFuncPt->SetParName(10, "#mu_{pileup}");
     fitFuncPt->SetParName(11, "#sigma_{pileup}");
    */

     TFitResultPtr fitPt= hnsigmapT[0][ncharge][netaDir][pad]->Fit(fitFuncPt,"SR");
     gPad->Update();


     TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
     ps->SetName("fitstats");
     ps->SetTextSize(0.031);

     TList *listOfLines = ps->GetListOfLines();

     // hnsigmapT[0][ncharge][netaDir][pad]->SetStats(0); 
     //gPad->Modified();

     double fitPars[15]={0};
     fitFuncPt->GetParameters(fitPars);

     TF1 *fpion= new TF1("fpion",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncPt->GetParameter(3),fitFuncPt->GetParameter(4),fitFuncPt->GetParameter(5)),-10,7);
     fpion->SetLineColor(2);
     fpion->SetLineWidth(2);
     fpion->Draw("same ");

     TF1 *fpk= new TF1("fpk",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncPt->GetParameter(0),fitFuncPt->GetParameter(1),fitFuncPt->GetParameter(2)),-10,7);
     fpk->SetLineColor(6);
     fpk->SetLineWidth(2);
     fpk->Draw("same ");
     TF1 *felectron= new TF1("felectron",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncPt->GetParameter(6),fitFuncPt->GetParameter(7),fitFuncPt->GetParameter(8)),-10,7);
     felectron->SetLineColor(7);
     felectron->SetLineWidth(2);
     felectron->Draw("same");
    /* //TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncPt->GetParameter(9),fitFuncPt->GetParameter(10),fitFuncPt->GetParameter(11)),-10,7);
     TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncPt->GetParameter(9),fitFuncPt->GetParameter(10),fitFuncPt->GetParameter(11)),-10,7);
     fpileup->SetLineColor(8);
     fpileup->SetLineWidth(2);
     fpileup->Draw("same");
    */
     gPad->Update();
     double pion=(double)fpion->Integral(-1,2)/(double)fitFuncPt->Integral(-1,2);
     pionPtbin[ncharge][netaDir][pad]=pion; 
     double pk=(double)fpk->Integral(-1,2)/(double)fitFuncPt->Integral(-1,2);
     protonPtbin[ncharge][netaDir][pad]=pk; 
     double electron=(double)felectron->Integral(-1,2)/(double)fitFuncPt->Integral(-1,2);
     electronPtbin[ncharge][netaDir][pad]=electron; 
     const char *per="%";
     TLatex purity;
     purity.SetTextAlign(13);
     purity.SetTextSize(0.038);
     purity.DrawLatex(-14,2e5,"#font[12]{#color[1]{Total fit}}");
     purity.DrawLatex(-14,1e5,Form("#font[12]{#color[2]{f(#pi) = %f}}",pion));
     purity.DrawLatex(-14,5e4,Form("#font[12]{#color[6]{f(p+k) = %f}}",pk));
     purity.DrawLatex(-14,2e4,Form("#font[12]{#color[7]{f(e) = %f}}",electron));
     purity.DrawLatex(-14,8e3,Form("#font[12]{#color[2]{#mu_{#pi}^{fixed}= %f}}", piMeanPt[pad]));
     purity.DrawLatex(-14,2e3,Form("#font[12]{#color[6]{#mu_{p+k}^{fixed} = %f}}",peakPt[1][ncharge][netaDir][pad]));
     purity.DrawLatex(-14,5e2,Form("#font[12]{#color[7]{#mu_{e}^{fixed} = %f}}",peakPt[2][ncharge][netaDir][pad]));
     //purity.DrawLatex(-14,1e3,Form("#font[12]{#color[6]{K_{peak} = %f}}",peakPt[1][ncharge][netaDir][pad]));

   }//pad loop
   cvs[nkbins][ncharge][netaDir]->Print("Plots/pionfitgaus_pTMBin.pdf");  
  }//pT-bins

	//purity in Minv bins 
  if(nkbins==1){ 
   ofstream outputMbin;
   outputMbin.open(Form("textFiles/purityMbin_%s%s.txt",charge[ncharge],etaDir[netaDir]));

   for(int pad=0; pad<6; pad++){
    if(pad==5) continue; 
    cvs[nkbins][ncharge][netaDir]->cd(pad+1);
    gPad->SetGrid(0,0);
    gPad->SetLogy();
    cout<<"Mbin,  "<<"charge: "<<charge[ncharge]<<" , eta: "<<etas[netaDir]<<", Pad:  "<<pad<<endl;
    hnsigmaM[0][ncharge][netaDir][pad]->SetLineColor(1);
    hnsigmaM[0][ncharge][netaDir][pad]->GetXaxis()->SetTitle("nSigmaPion"); 
    hnsigmaM[0][ncharge][netaDir][pad]->SetTitle(Form("%g < M_{inv} < %g, #eta %s, charge = %s ", M[pad], M[pad+1], etas[netaDir], charges[ncharge]));
    hnsigmaM[0][ncharge][netaDir][pad]->Draw();


    //proton + kaon combined fit
    //fitFuncM=new TF1("fitFuncM","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])+[9]*TMath::Gaus(x,[10],[11])",-10,7);         		
    fitFuncM=new TF1("fitFuncM","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])",-10,7);         		
    //fitFuncM=new TF1("fitFuncM","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",-10,7);         		
    fitFuncM->SetLineColor(1);
    fitFuncM->SetLineWidth(2);

    //proton
    //fitFuncM->SetParameter(0,aproton*0.1);
    //fitFuncM->SetParLimits(0,aproton*0.01, aproton*10);
    fitFuncM->SetParameter(1,peakM[1][ncharge][netaDir][pad]);
    fitFuncM->SetParLimits(1,peakM[1][ncharge][netaDir][pad],peakM[1][ncharge][netaDir][pad]);
    fitFuncM->SetParameter(2,1);
    fitFuncM->SetParLimits(2,0.5, 1.5);
    fitFuncM->SetParName(0, "A_{K}");
    fitFuncM->SetParName(1, "#mu_{K}");
    fitFuncM->SetParName(2, "#sigma_{K}");


    //pion
    fitFuncM->FixParameter(4,piMeanM[pad]);
    //fitFuncM->FixParameter(5,piWidthM[pad]);
    fitFuncM->SetParameter(5,1);
    fitFuncM->SetParLimits(5,0.8,1.2);

    fitFuncM->SetParName(3, "A_{#pi}");
    fitFuncM->SetParName(4, "#mu_{#pi}");
    fitFuncM->SetParName(5, "#sigma_{#pi}");

    //electron
    //fitFuncM->SetParameter(9,aelectron*0.1);
    //fitFuncM->SetParLimits(9,aelectron*0.01,aelectron*10);
    fitFuncM->FixParameter(7,peakM[2][ncharge][netaDir][pad]);
    fitFuncM->SetParameter(8,1.5);
    fitFuncM->SetParLimits(8,1,2.2);

    fitFuncM->SetParName(6, "A_{e}");
    fitFuncM->SetParName(7, "#mu_{e}");
    fitFuncM->SetParName(8, "#sigma_{e}");

 /*   //pileup
    //fitFuncM->SetParameter(9,apileup);
    //fitFuncM->SetParLimits(9,apileup*0.1,apileup*10);
    fitFuncM->SetParameter(10,6.5);
    fitFuncM->SetParLimits(10,6,10);
    fitFuncM->SetParameter(11,2);
    fitFuncM->SetParLimits(11,1.5,3);

    //fitFuncM->SetParNames("A_{#pi}","#mu_{#pi}","#sigma_{#pi^{+}}","A_{P}","#mu_{p}","#sigma_{P}","A_{K}","#mu_{K}","#sigma_{K}","A_{e}","#mu_{e}");
    //fitFuncM->SetParName(11, "#sigma_{e}");
    fitFuncM->SetParName(9, "A_{pileup}");
    fitFuncM->SetParName(10, "#mu_{pileup}");
    fitFuncM->SetParName(11, "#sigma_{pileup}");
*/

    TFitResultPtr fitM= hnsigmaM[0][ncharge][netaDir][pad]->Fit(fitFuncM,"SR");
    gPad->Update();


    TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
    ps->SetName("fitstats");
    ps->SetTextSize(0.031);

    TList *listOfLines = ps->GetListOfLines();

    // hnsigmaM[0][ncharge][netaDir][pad]->SetStats(0); 
    //gPad->Modified();

    double fitPars[15]={0};
    fitFuncM->GetParameters(fitPars);

    TF1 *fpion= new TF1("fpion",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(3),fitFuncM->GetParameter(4),fitFuncM->GetParameter(5)),-10,7);
    fpion->SetLineColor(2);
    fpion->SetLineWidth(2);
    fpion->Draw("same ");

    //TF1 *fproton= new TF1("fproton",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(3),peakM[0][ncharge][netaDir][pad],fitFuncM->GetParameter(4)),-10,7);
    TF1 *fpk= new TF1("fpk",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(0),fitFuncM->GetParameter(1),fitFuncM->GetParameter(2)),-10,7);
    fpk->SetLineColor(6);
    fpk->SetLineWidth(2);
    fpk->Draw("same ");
    //TF1 *fkaon= new TF1("fkaon",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(0),fitFuncM->GetParameter(1),fitFuncM->GetParameter(2)),-10,7);
    //fkaon->SetLineColor(6);
    //fkaon->SetLineWidth(2);
    //fkaon->Draw("same ");
    TF1 *felectron= new TF1("felectron",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(6),fitFuncM->GetParameter(7),fitFuncM->GetParameter(8)),-10,7);
    felectron->SetLineColor(7);
    felectron->SetLineWidth(2);
    felectron->Draw("same");
    /*//TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(9),fitFuncM->GetParameter(10),fitFuncM->GetParameter(11)),-10,7);
    TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncM->GetParameter(9),fitFuncM->GetParameter(10),fitFuncM->GetParameter(11)),-10,7);
    fpileup->SetLineColor(8);
    fpileup->SetLineWidth(2);
    fpileup->Draw("same");
    */
    gPad->Update();
    double pion=(double)fpion->Integral(-1,2)/(double)fitFuncM->Integral(-1,2);
    pionMbin[ncharge][netaDir][pad]=pion; 
    double pk=(double)fpk->Integral(-1,2)/(double)fitFuncM->Integral(-1,2);
    protonMbin[ncharge][netaDir][pad]=pk; 
    double electron=(double)felectron->Integral(-1,2)/(double)fitFuncM->Integral(-1,2);
    electronMbin[ncharge][netaDir][pad]=electron; 
    const char *per="%";
    TLatex purity;
    purity.SetTextAlign(13);
    purity.SetTextSize(0.038);

     purity.DrawLatex(-14,2e5,"#font[12]{#color[1]{Total fit}}");
     purity.DrawLatex(-14,1e5,Form("#font[12]{#color[2]{f(#pi) = %f}}",pion));
     purity.DrawLatex(-14,5e4,Form("#font[12]{#color[6]{f(p+k) = %f}}",pk));
     purity.DrawLatex(-14,2e4,Form("#font[12]{#color[7]{f(e) = %f}}",electron));
     purity.DrawLatex(-14,8e3,Form("#font[12]{#color[2]{#mu_{#pi}^{fixed}= %f}}", piMeanM[pad]));
     purity.DrawLatex(-14,2e3,Form("#font[12]{#color[6]{#mu_{p+k}^{fixed} = %f}}",peakM[1][ncharge][netaDir][pad]));
     purity.DrawLatex(-14,5e2,Form("#font[12]{#color[7]{#mu_{e}^{fixed} = %f}}",peakM[2][ncharge][netaDir][pad]));

   }//pad selection
  }//mbin loop
	cvs[nkbins][ncharge][netaDir]->Print("Plots/pionfitgaus_pTMBin.pdf");  
   cvs[nkbins][ncharge][netaDir]->SaveAs(Form("Plots/pionFit_%s%s%s.pdf",kbins[nkbins],charge[ncharge],etaDir[netaDir]));
    cvs[nkbins][ncharge][netaDir]->Print("Plots/pionfitgaus_pTMBin.pdf");
//if(nkbins==1 && ncharge==1 && netaDir==1)cvs[nkbins][ncharge][netaDir]->Print("Plots/pionfitgaus_pTMBin.pdf]");	 

      }//eta loop
    }//charge loop

   }//nkbins loop 

//pion fit in eta bins started 
 TCanvas *ceta[2];//for + and - charge
   for(int ncharge=0; ncharge<2; ncharge++){
      ceta[ncharge]=new TCanvas(Form("ceta%s",charge[ncharge]),"",1200,900);
ceta[ncharge]->Print("Plots/pionfitgaus_EtaBin.pdf[");
   ceta[ncharge]->Divide(3,3);
   ofstream outputEtabin;
   outputEtabin.open(Form("textFiles/purityEtabin_%s.txt",charge[ncharge]));
     for(int pad=0; pad<9; pad++){
       ceta[ncharge]->cd(pad+1);
       cout<<"charge: "<<charge[ncharge]<<", Pad:  "<<pad<<endl;
       gPad->SetGrid(0,0);
       gPad->SetLogy();
	  hnsigmaEta[0][ncharge][pad]->SetLineColor(1);
 	  hnsigmaEta[0][ncharge][pad]->GetXaxis()->SetTitle("nSigmaPion"); 
 	  hnsigmaEta[0][ncharge][pad]->GetXaxis()->SetTitleSize(0.04); 
 	  hnsigmaEta[0][ncharge][pad]->GetXaxis()->SetLabelSize(0.06); 
 	  hnsigmaEta[0][ncharge][pad]->GetYaxis()->SetLabelSize(0.06); 
 	  hnsigmaEta[0][ncharge][pad]->SetTitle(Form("%g < #eta^{#pi^{+}#pi^{-}} < %g, charge = %s ", eta_range[pad], eta_range[pad+1],charges[ncharge]));
	  hnsigmaEta[0][ncharge][pad]->Draw();

       	  double nproton= hnsigmaEta[0][ncharge][pad]->GetXaxis()->FindBin(peakEta[0][ncharge][pad]);//find bin where proton mean falls in nSigmaPion region
          double aproton= hnsigmaEta[0][ncharge][pad]->GetBinContent(nproton);//get pion amplitude at proton mean position
          double nkaon= hnsigmaEta[0][ncharge][pad]->GetXaxis()->FindBin(peakEta[1][ncharge][pad]);//bin for kaon mean
          double akaon= hnsigmaEta[0][ncharge][pad]->GetBinContent(nkaon);//pion amp at kaon mean
          double nelectron= hnsigmaEta[0][ncharge][pad]->GetXaxis()->FindBin(peakEta[2][ncharge][pad]);//bin at electron mean
          double aelectron= hnsigmaEta[0][ncharge][pad]->GetBinContent(nelectron);//pion amp at electron mean
          double npileup= hnsigmaEta[0][ncharge][pad]->GetXaxis()->FindBin(9.5);//pileup mean
          double apileup= hnsigmaEta[0][ncharge][pad]->GetBinContent(npileup);//pileup amp
        
     //proton + kaon combined fit

     //fitFuncEta=new TF1("fitFuncEta","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])+[9]*TMath::Gaus(x,[10],[11])",-10,7);         		
     fitFuncEta=new TF1("fitFuncEta","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])",-10,7);         		
     //fitFuncEta=new TF1("fitFuncEta","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",-10,7);         		
     fitFuncEta->SetLineColor(1);
     fitFuncEta->SetLineWidth(2);

     //proton+kaon
     fitFuncEta->SetParameter(1,peakEta[1][ncharge][pad]);
     fitFuncEta->SetParLimits(1,peakEta[1][ncharge][pad],peakEta[1][ncharge][pad]);
     fitFuncEta->SetParameter(2,1);
     fitFuncEta->SetParLimits(2,0.8, 1.2);
     fitFuncEta->SetParLimits(2,0.8, 1.3);
     if(ncharge==1&&pad==4)fitFuncEta->SetParLimits(2,0.8, 1.3);
     fitFuncEta->SetParName(0, "A_{K}");
     fitFuncEta->SetParName(1, "#mu_{K}");
     fitFuncEta->SetParName(2, "#sigma_{K}");

     //pion
     fitFuncEta->FixParameter(4,piMeanEta[pad]);
     //fitFuncEta->FixParameter(5,piWidthEta[pad]);
     fitFuncEta->SetParameter(5,1);
     fitFuncEta->SetParLimits(5,0.8,1.2);

     fitFuncEta->SetParName(3, "A_{#pi}");
     fitFuncEta->SetParName(4, "#mu_{#pi}");
     fitFuncEta->SetParName(5, "#sigma_{#pi}");

     //electron
     fitFuncEta->FixParameter(7,peakEta[2][ncharge][pad]);
     fitFuncEta->SetParameter(8,1.5);
     fitFuncEta->SetParLimits(8,1,2.2);

     fitFuncEta->SetParName(6, "A_{e}");
     fitFuncEta->SetParName(7, "#mu_{e}");
     fitFuncEta->SetParName(8, "#sigma_{e}");

     /*//pileup
     fitFuncEta->SetParameter(10,6.5);
     fitFuncEta->SetParLimits(10,6,10);
     fitFuncEta->SetParameter(11,2);
     fitFuncEta->SetParLimits(11,1.5,3);

     //fitFuncEta->SetParNames("A_{#pi}","#mu_{#pi}","#sigma_{#pi^{+}}","A_{P}","#mu_{p}","#sigma_{P}","A_{K}","#mu_{K}","#sigma_{K}","A_{e}","#mu_{e}");
     //fitFuncEta->SetParName(11, "#sigma_{e}");
     fitFuncEta->SetParName(9, "A_{pileup}");
     fitFuncEta->SetParName(10, "#mu_{pileup}");
     fitFuncEta->SetParName(11, "#sigma_{pileup}");
     */

     TFitResultPtr fitEta= hnsigmaEta[0][ncharge][pad]->Fit(fitFuncEta,"SR");
     gPad->Update();


     TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
     ps->SetName("fitstats");
     ps->SetTextSize(0.031);

     TList *listOfLines = ps->GetListOfLines();

     // hnsigmaEta[0][ncharge][pad]->SetStats(0); 
     //gPad->Modified();

     double fitPars[15]={0};
     fitFuncEta->GetParameters(fitPars);

     TF1 *fpion= new TF1("fpion",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(3),fitFuncEta->GetParameter(4),fitFuncEta->GetParameter(5)),-10,7);
     fpion->SetLineColor(2);
     fpion->SetLineWidth(2);
     fpion->Draw("same ");

    /* //TF1 *fproton= new TF1("fproton",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(3),peakEta[0][ncharge][pad],fitFuncEta->GetParameter(4)),-10,7);
     TF1 *fproton= new TF1("fproton",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(0),fitFuncEta->GetParameter(1),fitFuncEta->GetParameter(2)),-10,7);
     fproton->SetLineColor(4);
     fproton->SetLineWidth(2);
     fproton->Draw("same ");
    */
     TF1 *fpk= new TF1("fpk",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(0),fitFuncEta->GetParameter(1),fitFuncEta->GetParameter(2)),-10,7);
     fpk->SetLineColor(6);
     fpk->SetLineWidth(2);
     fpk->Draw("same ");
     TF1 *felectron= new TF1("felectron",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(6),fitFuncEta->GetParameter(7),fitFuncEta->GetParameter(8)),-10,7);
     felectron->SetLineColor(7);
     felectron->SetLineWidth(2);
     felectron->Draw("same");
   /*  //TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(9),fitFuncEta->GetParameter(10),fitFuncEta->GetParameter(11)),-10,7);
     TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncEta->GetParameter(9),fitFuncEta->GetParameter(10),fitFuncEta->GetParameter(11)),-10,7);
     fpileup->SetLineColor(8);
     fpileup->SetLineWidth(2);
     fpileup->Draw("same");
   */
     gPad->Update();
     double pion=(double)fpion->Integral(-1,2)/(double)fitFuncEta->Integral(-1,2);
     pionEtabin[ncharge][pad]=pion; 
     double pk=(double)fpk->Integral(-1,2)/(double)fitFuncEta->Integral(-1,2);
     protonEtabin[ncharge][pad]=pk; 
     double electron=(double)felectron->Integral(-1,2)/(double)fitFuncEta->Integral(-1,2);
     electronEtabin[ncharge][pad]=electron; 
     const char *per="%";
       TLatex purity;
       purity.SetTextAlign(13);
       purity.SetTextSize(0.05);

     purity.DrawLatex(-14,2e5,"#font[12]{#color[1]{Total fit}}");
     purity.DrawLatex(-14,8e4,Form("#font[12]{#color[2]{f(#pi) = %f}}",pion));
     purity.DrawLatex(-14,2e4,Form("#font[12]{#color[6]{f(p+k) = %f}}",pk));
     purity.DrawLatex(-14,8e3,Form("#font[12]{#color[7]{f(e) = %f}}",electron));
     purity.DrawLatex(-14,3e3,Form("#font[12]{#color[2]{#mu_{#pi}^{fixed}= %f}}", piMeanEta[pad]));
     purity.DrawLatex(-14,8e2,Form("#font[12]{#color[6]{#mu_{p+k}^{fixed} = %f}}",peakEta[1][ncharge][pad]));
     purity.DrawLatex(-14,2e2,Form("#font[12]{#color[7]{#mu_{e}^{fixed} = %f}}",peakEta[2][ncharge][pad]));

   }//pad
   ceta[ncharge]->SaveAs(Form("Plots/pionFit_EtaBin_%s.pdf",charge[ncharge]));
ceta[ncharge]->Print("Plots/pionfitgaus_EtaBin.pdf");
if(ncharge==1)ceta[ncharge]->Print("Plots/pionfitgaus_EtaBin.pdf]");
 }//charge
//--------------------------------------PID in eta ends


//PID for AUT vs Minv, integrated p_T bins
TCanvas *icvs[2][2];
   for(int ncharge=0; ncharge<2; ncharge++){
     for(int netaDir=0; netaDir<2; netaDir++){
      icvs[ncharge][netaDir]=new TCanvas(Form("icvs%s%s",charge[ncharge],etaDir[netaDir]),"",1200,900);
      icvs[ncharge][netaDir]->SetFrameFillStyle(0);
      icvs[ncharge][netaDir]->Divide(3,3);
      icvs[ncharge][netaDir]->Print("Plots/pionfitgaus_IntMBin.pdf[");  
   ofstream outputIntMbin;
   outputIntMbin.open(Form("textFiles/purityIntMbin_%s%s.txt",charge[ncharge],etaDir[netaDir]));
   for(int pad=0; pad<9; pad++){

    icvs[ncharge][netaDir]->cd(pad+1);
    cout<<"int M  Bin, charge: "<<charge[ncharge]<<" , eta: "<<etas[netaDir]<<", Pad:  "<<pad<<endl;

    gPad->SetGrid(0,0);
    gPad->SetLogy();
    hnsigmaIntM[0][ncharge][netaDir][pad]->SetLineColor(1);
    hnsigmaIntM[0][ncharge][netaDir][pad]->GetXaxis()->SetTitle("nSigmaPion"); 
    hnsigmaIntM[0][ncharge][netaDir][pad]->SetTitle(Form("%g < M_{inv} < %g, #eta %s, charge = %s ", M_in[pad], M_in[pad+1], etas[netaDir], charges[ncharge]));
    //hnSigmaTotalIntM[ncharge][netaDir][pad]->Sumw2(kFALSE);
    hnsigmaIntM[0][ncharge][netaDir][pad]->Draw();
     //fitFuncIntM=new TF1("fitFuncIntM","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])+[9]*TMath::Gaus(x,[10],[11])",-10,7);         		
     fitFuncIntM=new TF1("fitFuncIntM","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])",-10,7);         		
     //fitFuncIntM=new TF1("fitFuncIntM","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",-10,7);         		
     fitFuncIntM->SetLineColor(1);
     fitFuncIntM->SetLineWidth(2);

     //proton+kaon
     //fitFuncIntM->SetParameter(0,aproton*0.1);
     //fitFuncIntM->SetParLimits(0,aproton*0.01, aproton*10);
     fitFuncIntM->SetParameter(1,peakIntM[1][ncharge][netaDir][pad]);
     fitFuncIntM->SetParLimits(1,peakIntM[1][ncharge][netaDir][pad],peakIntM[1][ncharge][netaDir][pad]);
     fitFuncIntM->SetParameter(2,1);
     fitFuncIntM->SetParLimits(2,0.5, 1.3);
     fitFuncIntM->SetParName(0, "A_{K+p}");
     fitFuncIntM->SetParName(1, "#mu_{K+p}");
     fitFuncIntM->SetParName(2, "#sigma_{K+p}");


     //pion
     //fitFuncIntM->SetParameter(4,piMeanIntM[pad]);
     fitFuncIntM->SetParameter(4,0.);
     //fitFuncIntM->SetParLimits(4,piMeanIntM[pad],piMeanIntM[pad]);
     fitFuncIntM->SetParLimits(4,-0.1,0.1);
     //fitFuncIntM->FixParameter(5,piWidthIntM[pad]);
     fitFuncIntM->SetParameter(5,1);
     fitFuncIntM->SetParLimits(5,0.8,1.2);

     fitFuncIntM->SetParName(3, "A_{#pi}");
     fitFuncIntM->SetParName(4, "#mu_{#pi}");
     fitFuncIntM->SetParName(5, "#sigma_{#pi}");

     //electron
     //fitFuncIntM->SetParameter(9,aelectron*0.1);
     //fitFuncIntM->SetParLimits(9,aelectron*0.01,aelectron*10);
     fitFuncIntM->SetParameter(7,peakIntM[2][ncharge][netaDir][pad]);
     fitFuncIntM->SetParLimits(7,peakIntM[2][ncharge][netaDir][pad], peakIntM[2][ncharge][netaDir][pad]);
     fitFuncIntM->SetParameter(8,1.5);
     fitFuncIntM->SetParLimits(8,1,2.2);

     fitFuncIntM->SetParName(6, "A_{e}");
     fitFuncIntM->SetParName(7, "#mu_{e}");
     fitFuncIntM->SetParName(8, "#sigma_{e}");

    /* //pileup
     //fitFuncIntM->SetParameter(9,apileup);
     //fitFuncIntM->SetParLimits(9,apileup*0.1,apileup*10);
     fitFuncIntM->SetParameter(10,6.5);
     fitFuncIntM->SetParLimits(10,6,10);
     fitFuncIntM->SetParameter(11,2);
     fitFuncIntM->SetParLimits(11,1.5,3);
    */
     //fitFuncIntM->SetParNames("A_{#pi}","#mu_{#pi}","#sigma_{#pi^{+}}","A_{P}","#mu_{p}","#sigma_{P}","A_{K}","#mu_{K}","#sigma_{K}","A_{e}","#mu_{e}");
    /* //fitFuncIntM->SetParName(11, "#sigma_{e}");
     fitFuncIntM->SetParName(9, "A_{pileup}");
     fitFuncIntM->SetParName(10, "#mu_{pileup}");
     fitFuncIntM->SetParName(11, "#sigma_{pileup}");
    */

     TFitResultPtr fitIntM= hnsigmaIntM[0][ncharge][netaDir][pad]->Fit(fitFuncIntM,"SR");
     gPad->Update();


     TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
     ps->SetName("fitstats");
     ps->SetTextSize(0.031);

     TList *listOfLines = ps->GetListOfLines();

     // hnsigmaIntM[0][ncharge][netaDir][pad]->SetStats(0); 
     //gPad->Modified();

     double fitPars[15]={0};
     fitFuncIntM->GetParameters(fitPars);

     TF1 *fpion= new TF1("fpion",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncIntM->GetParameter(3),fitFuncIntM->GetParameter(4),fitFuncIntM->GetParameter(5)),-10,7);
     fpion->SetLineColor(2);
     fpion->SetLineWidth(2);
     fpion->Draw("same ");

     TF1 *fpk= new TF1("fpk",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncIntM->GetParameter(0),fitFuncIntM->GetParameter(1),fitFuncIntM->GetParameter(2)),-10,7);
     fpk->SetLineColor(6);
     fpk->SetLineWidth(2);
     fpk->Draw("same ");
     TF1 *felectron= new TF1("felectron",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncIntM->GetParameter(6),fitFuncIntM->GetParameter(7),fitFuncIntM->GetParameter(8)),-10,7);
     felectron->SetLineColor(7);
     felectron->SetLineWidth(2);
     felectron->Draw("same");
    /* //TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncIntM->GetParameter(9),fitFuncIntM->GetParameter(10),fitFuncIntM->GetParameter(11)),-10,7);
     TF1 *fpileup= new TF1("fpileup",Form("%f*exp(-0.5*pow(x-%f,2)/pow(%f,2))",fitFuncIntM->GetParameter(9),fitFuncIntM->GetParameter(10),fitFuncIntM->GetParameter(11)),-10,7);
     fpileup->SetLineColor(8);
     fpileup->SetLineWidth(2);
     fpileup->Draw("same");
    */
     gPad->Update();
     double pion=(double)fpion->Integral(-1,2)/(double)fitFuncIntM->Integral(-1,2);
     pionIntMbin[ncharge][netaDir][pad]=pion; 
     double pk=(double)fpk->Integral(-1,2)/(double)fitFuncIntM->Integral(-1,2);
     protonIntMbin[ncharge][netaDir][pad]=pk; 
     double electron=(double)felectron->Integral(-1,2)/(double)fitFuncIntM->Integral(-1,2);
     electronIntMbin[ncharge][netaDir][pad]=electron; 
     const char *per="%";
     TLatex purity;
     purity.SetTextAlign(13);
     purity.SetTextSize(0.038);
     purity.DrawLatex(-14,2e5,"#font[12]{#color[1]{Total fit}}");
     purity.DrawLatex(-14,1e5,Form("#font[12]{#color[2]{f(#pi) = %f}}",pion));
     purity.DrawLatex(-14,5e4,Form("#font[12]{#color[6]{f(p+k) = %f}}",pk));
     purity.DrawLatex(-14,2e4,Form("#font[12]{#color[7]{f(e) = %f}}",electron));
     //purity.DrawLatex(-14,8e3,Form("#font[12]{#color[2]{#mu_{#pi}^{fixed}= %f}}", piMeanIntM[pad]));
     purity.DrawLatex(-14,2e3,Form("#font[12]{#color[6]{#mu_{p+k}^{fixed} = %f}}",peakIntM[1][ncharge][netaDir][pad]));
     purity.DrawLatex(-14,5e2,Form("#font[12]{#color[7]{#mu_{e}^{fixed} = %f}}",peakIntM[2][ncharge][netaDir][pad]));
     //purity.DrawLatex(-14,1e3,Form("#font[12]{#color[6]{K_{peak} = %f}}",peakIntM[1][ncharge][netaDir][pad]));

   }//pad loop
   icvs[ncharge][netaDir]->Print("Plots/pionfitgaus_IntMBin.pdf");  
  }
}
//--------------------------------------
//----Calculate and plot pair-purity
TGraphErrors *grPtbinPosGt[4];
TGraphErrors *grPtbinPosLt[4];
TGraphErrors *grPtbinNegGt[4];
TGraphErrors *grPtbinNegLt[4];
TGraphErrors *grMbinPosGt[4];
TGraphErrors *grMbinPosLt[4];
TGraphErrors *grMbinNegGt[4];
TGraphErrors *grMbinNegLt[4];

TGraphErrors *grEtaBinN;  TGraphErrors *grEtaBinKaonN;  TGraphErrors *grEtaBinProtonN;  TGraphErrors *grEtaBinElectronN; 
TGraphErrors *grEtaBinP;  TGraphErrors *grEtaBinKaonP;  TGraphErrors *grEtaBinProtonP;  TGraphErrors *grEtaBinElectronP;
double pairpurityEtabin[9]={0};
double pionPairEta[9]={0};
double kpPairEta[9]={0}; 
double ePairEta[9]={0}; 

double pionPairPtGt[5]={0};                                        
double pionPairPtLt[5]={0};                                        
double pionPairMGt[5]={0};                                         
double pionPairMLt[5]={0};                                         
double pionPairIntMGt[9]={0};                                         
double pionPairIntMLt[9]={0};                                         
                                                                   
double kpPairPtGt[5]={0};                                        
double kpPairPtLt[5]={0};                                        
double kpPairMGt[5]={0};                                         
double kpPairMLt[5]={0};                                         
double kpPairIntMGt[9]={0};                                         
double kpPairIntMLt[9]={0};                                         
                                                                   
double ePairPtGt[5]={0};                                        
double ePairPtLt[5]={0};                                        
double ePairMGt[5]={0};                                         
double ePairMLt[5]={0};                                         
double ePairIntMGt[9]={0};                                         
double ePairIntMLt[9]={0};                                         
//pair fractions from TOF only PID--------------------- 
//pion
double pionPairPtGtTO[5]={0.995685,0.99413,0.992819,0.981137,0.951615};
double pionPairPtLtTO[5]={0.997217,0.995865,0.994613,0.982959,0.954391};
double pionPairMGtTO[5]={0.978901,0.98193,0.981563,0.980261,0.971465};
double pionPairMLtTO[5]={0.981487,0.984854,0.98391,0.982889,0.973322};
//kaon+proton
double kpPairPtGtTO[5]={0,0,0,2.08895e-05,0.000131916};
double kpPairPtLtTO[5]={0,0,0,1.90933e-05,0.000122819};
double kpPairMGtTO[5]={1.88864e-05,3.60047e-05,3.40103e-05,4.00047e-05,7.8545e-05};
double kpPairMLtTO[5]={1.49303e-05,2.6513e-05,2.64161e-05,3.2977e-05,7.32954e-05};
// pair purity in eta bins
double pionPairEtaTO[9]={0.97593,0.98328,0.982627,0.978905,0.975671,0.978282,0.983113,0.983608,0.976228};
double kpPairEtaTO[9]={3.55615e-05,2.98475e-05,3.96205e-05,5.86822e-05,7.06504e-05,5.12718e-05,2.85773e-05,2.04411e-05,2.54495e-05};
double electronPairEtaTO[9]={3.69162e-05,8.40101e-06,5.82411e-06,8.56687e-06,1.40317e-05,1.3517e-05,9.81749e-06,1.35601e-05,4.67372e-05};
// pair purity in integrated pT bins
double pionPairIntMGtT[9]={0.975651,0.987897,0.985057,0.981433,0.984766,0.983844,0.98232,0.979519,0.967469};
double pionPairIntMLtT[9]={0,0,0,0,0,0,0,0,0};
double kpPairIntMGtT[9]={1.03398e-05,1.02371e-05,1.9807e-05,3.59803e-05,1.7999e-05,2.05885e-05,2.80107e-05,4.46095e-05,9.39217e-05};
double kpPairIntMLtT[9]={0,0,0,0,0,0,0,0,0};
double electronPairIntMGtT[9]={8.00411e-05,8.2071e-06,9.2425e-06,1.09599e-05,1.15676e-05,1.26905e-05,1.28323e-05,1.30567e-05,4.49852e-05};
double electronPairIntMLtT[9]={0,0,0,0,0,0,0,0,0};


//-------------------------------------------------------
                                                                   
                                                                    
for(int n=0; n<5; n++){
//signal pair fractions
pionPairPtGt[n]=pionPtbin[0][0][n]*pionPtbin[1][0][n];
pionPairPtLt[n]=pionPtbin[0][1][n]*pionPtbin[1][1][n];
pionPairMGt[n]=pionMbin[0][0][n]*pionMbin[1][0][n];
pionPairMLt[n]=pionMbin[0][1][n]*pionMbin[1][1][n];

kpPairPtGt[n]=protonPtbin[0][0][n]*protonPtbin[1][0][n];
kpPairPtLt[n]=protonPtbin[0][1][n]*protonPtbin[1][1][n];
kpPairMGt[n]=protonMbin[0][0][n]*protonMbin[1][0][n];
kpPairMLt[n]=protonMbin[0][1][n]*protonMbin[1][1][n];

ePairPtGt[n]=electronPtbin[0][0][n]*electronPtbin[1][0][n];
ePairPtLt[n]=electronPtbin[0][1][n]*electronPtbin[1][1][n];
ePairMGt[n]=electronMbin[0][0][n]*electronMbin[1][0][n];
ePairMLt[n]=electronMbin[0][1][n]*electronMbin[1][1][n];
}


// signal  pair fractions
TCanvas *cpt = new TCanvas("cpt","",500,500);
cpt->cd();
gPad->SetGrid(1,1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);
//pion in signal
TGraphErrors *pionSigGtPt=new TGraphErrors(5,avgPt,pionPairPtGt,0,0);
pionSigGtPt->SetTitle("");
pionSigGtPt->GetYaxis()->SetTitle("pair fractions");
pionSigGtPt->GetYaxis()->SetRangeUser(-0.05,1.05);
pionSigGtPt->GetXaxis()->SetLimits(2,10);
pionSigGtPt->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}} bin center");
pionSigGtPt->SetMarkerStyle(20);
pionSigGtPt->SetMarkerSize(1);
pionSigGtPt->SetMarkerColor(2);
pionSigGtPt->GetXaxis()->SetNdivisions(5, 5,0,kTRUE);
pionSigGtPt->GetYaxis()->SetTitleOffset(0.8);
pionSigGtPt->GetXaxis()->SetTitleOffset(1.1);
pionSigGtPt->GetYaxis()->SetLabelSize(0.05);
pionSigGtPt->GetXaxis()->SetLabelSize(0.05);
pionSigGtPt->GetYaxis()->SetTitleSize(0.05);
pionSigGtPt->GetXaxis()->SetTitleSize(0.05);
pionSigGtPt->Draw("AP");

TGraphErrors *pionSigLtPt=new TGraphErrors(5,avgPt,pionPairPtLt,0,0);
pionSigLtPt->SetMarkerStyle(24);
pionSigLtPt->SetMarkerColor(2);
pionSigLtPt->Draw("same P");
//TOF only pion
TGraphErrors *pionSigGtPtTO=new TGraphErrors(5,avgPt,pionPairPtGtTO,0,0);
pionSigGtPtTO->SetMarkerStyle(23);
pionSigGtPtTO->SetMarkerColor(4);
pionSigGtPtTO->Draw("same P");
TGraphErrors *pionSigLtPtTO=new TGraphErrors(5,avgPt,pionPairPtLtTO,0,0);
pionSigLtPtTO->SetMarkerStyle(26);
pionSigLtPtTO->SetMarkerColor(4);
pionSigLtPtTO->Draw("same P");
//------------

TGraphErrors *kpSigGtPt=new TGraphErrors(5,avgPt,kpPairPtGt,0,0);
kpSigGtPt->SetMarkerStyle(22);
kpSigGtPt->SetMarkerColor(4);
kpSigGtPt->Draw("same P");

TGraphErrors *kpSigLtPt=new TGraphErrors(5,avgPt,kpPairPtLt,0,0);
kpSigLtPt->SetMarkerStyle(23);
kpSigLtPt->SetMarkerColor(7);
kpSigLtPt->Draw("same P");

TGraphErrors *eSigGtPt=new TGraphErrors(5,avgPt,ePairPtGt,0,0);
eSigGtPt->SetMarkerStyle(25);
eSigGtPt->SetMarkerColor(6);
eSigGtPt->Draw("same P");

TGraphErrors *eSigLtPt=new TGraphErrors(5,avgPt,ePairPtLt,0,0);
eSigLtPt->SetMarkerStyle(26);
eSigLtPt->SetMarkerColor(8);
eSigLtPt->Draw("same P");

TLegend *l1=new TLegend(0.4,0.4,0.7,0.7);
l1->SetNColumns(2);
l1->AddEntry(""," #eta > 0 "," ");
l1->AddEntry(""," #eta < 0 "," ");
l1->AddEntry(pionSigGtPt," Pion TOF cut", " lp ");
l1->AddEntry(pionSigGtPt," Pion TOF cut"," lp ");
l1->AddEntry(kpSigGtPt," K+p TOF cut"," lp ");
l1->AddEntry(kpSigLtPt," K+P TOF cut"," lp ");
l1->AddEntry(eSigGtPt," e TOF cut "," lp ");
l1->AddEntry(eSigLtPt," e TOF cut "," lp ");
l1->AddEntry(pionSigGtPtTO," Pion TOF only ", " lp ");
l1->AddEntry(pionSigGtPtTO," Pion TOF only"," lp ");
l1->Draw();
cpt->Update();
cpt->SaveAs("Plots/newpionpurityStTOFPtWithTOFonly.pdf");

TCanvas *cM = new TCanvas("cM","",500,500);
cM->cd();
gPad->SetGrid(1,1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);
//pion in signal
TGraphErrors *pionSigGtM=new TGraphErrors(5,avgM,pionPairMGt,0,0);
pionSigGtM->SetTitle("");
pionSigGtM->GetYaxis()->SetTitle("pair fractions");
pionSigGtM->GetYaxis()->SetRangeUser(-0.05,1.05);
pionSigGtM->GetXaxis()->SetTitle("M^{#pi^{+}#pi^{-}} bin center");
pionSigGtM->SetMarkerStyle(20);
pionSigGtM->SetMarkerSize(1);
pionSigGtM->SetMarkerColor(2);
pionSigGtM->GetXaxis()->SetNdivisions(5, 5,0,kTRUE);
pionSigGtM->GetYaxis()->SetTitleOffset(0.8);
pionSigGtM->GetXaxis()->SetTitleOffset(1.1);
pionSigGtM->GetYaxis()->SetLabelSize(0.05);
pionSigGtM->GetXaxis()->SetLabelSize(0.05);
pionSigGtM->GetYaxis()->SetTitleSize(0.05);
pionSigGtM->GetXaxis()->SetTitleSize(0.05);
pionSigGtM->Draw("AP");

TGraphErrors *pionSigLtM=new TGraphErrors(5,avgM,pionPairMLt,0,0);
pionSigLtM->SetMarkerStyle(24);
pionSigLtM->SetMarkerColor(3);
pionSigLtM->Draw("same P");
//TOF only pion
TGraphErrors *pionSigGtMTO=new TGraphErrors(5,avgM,pionPairMGtTO,0,0);
pionSigGtMTO->SetMarkerStyle(23);
pionSigGtMTO->SetMarkerColor(4);
pionSigGtMTO->Draw("same P");
TGraphErrors *pionSigLtMTO=new TGraphErrors(5,avgM,pionPairMLtTO,0,0);
pionSigLtMTO->SetMarkerStyle(26);
pionSigLtMTO->SetMarkerColor(4);
pionSigLtMTO->Draw("same P");
//------------

TGraphErrors *kpSigLtM=new TGraphErrors(5,avgM,kpPairMLt,0,0);
kpSigLtM->SetMarkerStyle(22);
kpSigLtM->SetMarkerColor(4);
kpSigLtM->Draw("same P");

TGraphErrors *kpSigGtM=new TGraphErrors(5,avgM,kpPairMGt,0,0);
kpSigGtM->SetMarkerStyle(23);
kpSigGtM->SetMarkerColor(5);
kpSigGtM->Draw("same P");

TGraphErrors *eSigLtM=new TGraphErrors(5,avgM,ePairMLt,0,0);
eSigLtM->SetMarkerStyle(26);
eSigLtM->SetMarkerColor(6);
eSigLtM->Draw("same P");

TGraphErrors *eSigGtM=new TGraphErrors(5,avgM,ePairMGt,0,0);
eSigGtM->SetMarkerStyle(27);
eSigGtM->SetMarkerColor(7);
eSigGtM->Draw("same P");

TLegend *l2=new TLegend(0.4,0.4,0.7,0.7);
l2->SetNColumns(2);
l2->AddEntry(""," #eta > 0 "," ");
l2->AddEntry(""," #eta < 0 "," ");
l2->AddEntry(pionSigGtM," Pion TOF cut"," lp");
l2->AddEntry(pionSigLtM," Pion TOF cut","lp");
l2->AddEntry(kpSigGtM," K+p TOF cut","lp");
l2->AddEntry(kpSigLtM," K+p TOF cut ","lp");
l2->AddEntry(eSigGtM," e TOF cut ","lp");
l2->AddEntry(eSigLtM," e TOF cut ","lp");
l2->AddEntry(pionSigGtMTO," Pion TOF only"," lp");
l2->AddEntry(pionSigLtMTO," Pion TOF only","lp");
l2->Draw();
cM->Update();
cM->SaveAs("Plots/pionpurityStTOF_MWithTOFonly.pdf");

//M-Bins,  integated pT bins------------------- 
//pair purity in integrated pT bins
for(int n=0; n<9; n++){
pionPairIntMGt[n]=pionIntMbin[0][0][n]*pionIntMbin[1][0][n];
pionPairIntMLt[n]=pionIntMbin[0][1][n]*pionIntMbin[1][1][n];
kpPairIntMGt[n]=protonIntMbin[0][0][n]*protonIntMbin[1][0][n];
kpPairIntMLt[n]=protonIntMbin[0][1][n]*protonIntMbin[1][1][n];
ePairIntMGt[n]=electronIntMbin[0][0][n]*electronIntMbin[1][0][n];
ePairIntMLt[n]=electronIntMbin[0][1][n]*electronIntMbin[1][1][n];
}
TCanvas *cIntM = new TCanvas("cIntM","",500,500);
cIntM->cd();
gPad->SetGrid(1,1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);
//pion in signal
TGraphErrors *pionSigGtIntM=new TGraphErrors(9,avgIntM,pionPairIntMGt,0,0);
pionSigGtIntM->SetTitle("");
pionSigGtIntM->GetYaxis()->SetTitle("pair fractions");
pionSigGtIntM->GetYaxis()->SetRangeUser(-0.05,1.05);
pionSigGtIntM->GetXaxis()->SetTitle("IntM^{#pi^{+}#pi^{-}} bin center");
pionSigGtIntM->SetMarkerStyle(20);
pionSigGtIntM->SetMarkerSize(1);
pionSigGtIntM->SetMarkerColor(2);
pionSigGtIntM->GetXaxis()->SetNdivisions(9, 5,0,kTRUE);
pionSigGtIntM->GetYaxis()->SetTitleOffset(0.8);
pionSigGtIntM->GetXaxis()->SetTitleOffset(1.1);
pionSigGtIntM->GetYaxis()->SetLabelSize(0.05);
pionSigGtIntM->GetXaxis()->SetLabelSize(0.05);
pionSigGtIntM->GetYaxis()->SetTitleSize(0.05);
pionSigGtIntM->GetXaxis()->SetTitleSize(0.05);
pionSigGtIntM->Draw("AP");

TGraphErrors *pionSigGtIntMT=new TGraphErrors(9,avgIntM,pionPairIntMGtT,0,0);//TOF only PID
pionSigGtIntMT->SetMarkerStyle(24);
pionSigGtIntMT->SetMarkerSize(1);
pionSigGtIntMT->SetMarkerColor(2);
pionSigGtIntMT->Draw("same P");

TGraphErrors *pionSigLtIntM=new TGraphErrors(9,avgIntM,pionPairIntMLt,0,0);
pionSigLtIntM->SetMarkerStyle(24);
pionSigLtIntM->SetMarkerColor(3);
pionSigLtIntM->Draw("same P");

TGraphErrors *kpSigLtIntM=new TGraphErrors(9,avgIntM,kpPairIntMLt,0,0);
kpSigLtIntM->SetMarkerStyle(22);
kpSigLtIntM->SetMarkerColor(4);
kpSigLtIntM->Draw("same P");

TGraphErrors *kpSigGtIntM=new TGraphErrors(9,avgIntM,kpPairIntMGt,0,0);
kpSigGtIntM->SetMarkerStyle(23);
kpSigGtIntM->SetMarkerColor(5);
kpSigGtIntM->Draw("same P");

TGraphErrors *eSigLtIntM=new TGraphErrors(9,avgIntM,ePairIntMLt,0,0);
eSigLtIntM->SetMarkerStyle(26);
eSigLtIntM->SetMarkerColor(6);
eSigLtIntM->Draw("same P");

TGraphErrors *eSigGtIntM=new TGraphErrors(9,avgIntM,ePairIntMGt,0,0);
eSigGtIntM->SetMarkerStyle(27);
eSigGtIntM->SetMarkerColor(7);
eSigGtIntM->Draw("same P");

TLegend *l2m=new TLegend(0.4,0.4,0.7,0.7);
l2m->SetNColumns(2);
l2m->AddEntry(""," #eta > 0 "," ");
l2m->AddEntry(""," #eta < 0 "," ");
l2m->AddEntry(pionSigGtIntM," Pion TOF cut"," lp");
l2m->AddEntry(pionSigLtIntM," Pion TOF cut","lp");
l2m->AddEntry(kpSigGtIntM," K+p TOF cut","lp");
l2m->AddEntry(kpSigLtIntM," K+p TOF cut ","lp");
l2m->AddEntry(eSigGtIntM," e TOF cut ","lp");
l2m->AddEntry(eSigLtIntM," e TOF cut ","lp");
l2m->AddEntry(pionSigGtIntMT," Pion TOF only"," lp");
l2m->Draw();
cIntM->Update();
cIntM->SaveAs("Plots/pionpurityStTOF_IntMWithTOFCut.pdf");


////////////////////////
//-----Eta Bin-------------//
for(int n=0; n<9; n++){
//signal pair fractions
pionPairEta[n]=pionEtabin[0][n]*pionEtabin[1][n];
kpPairEta[n]=protonEtabin[0][n]*protonEtabin[1][n];
ePairEta[n]=electronEtabin[0][n]*electronEtabin[1][n];
}

TCanvas *cpairEta = new TCanvas("cpairEta","",500,500);
cpairEta->cd();
gPad->SetGrid(1,1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);
//pion in signal
TGraphErrors *grEta = new TGraphErrors(9,avgEta,pionPairEta, 0, 0);
grEta->SetTitle("");
grEta->GetYaxis()->SetTitle("pair fractions");
grEta->GetYaxis()->SetRangeUser(-0.05,1.05);
grEta->GetXaxis()->SetLimits(-1,1);
grEta->GetXaxis()->SetTitle("#eta^{#pi^{+}#pi^{-}} bin center");
grEta->SetMarkerStyle(20);
grEta->SetMarkerSize(1);
grEta->SetMarkerColor(2);
grEta->GetXaxis()->SetNdivisions(5, 5,0,kTRUE);
grEta->GetYaxis()->SetTitleOffset(0.8);
grEta->GetXaxis()->SetTitleOffset(1.1);
grEta->GetYaxis()->SetLabelSize(0.05);
grEta->GetXaxis()->SetLabelSize(0.05);
grEta->GetYaxis()->SetTitleSize(0.05);
grEta->GetXaxis()->SetTitleSize(0.05);
grEta->Draw("ap");

//TOF only pion
TGraphErrors *pionSigEtaTO=new TGraphErrors(9,avgEta,pionPairEtaTO,0,0);
pionSigEtaTO->SetMarkerStyle(23);
pionSigEtaTO->SetMarkerColor(4);
pionSigEtaTO->Draw("same P");
//------------

//kaon+proton in signal
TGraphErrors *grkpEta = new TGraphErrors(9,avgEta,kpPairEta,0,0);
grkpEta->SetMarkerStyle(21);
grkpEta->SetMarkerSize(1);
grkpEta->SetMarkerColor(2);
grkpEta->Draw("p same");
//electron 
TGraphErrors *grEEta = new TGraphErrors(9,avgEta,ePairEta,0, 0);
grEEta->SetMarkerStyle(22);
grEEta->SetMarkerSize(1);
grEEta->SetMarkerColor(2);
grEEta->Draw("p same");

TLegend *legEta=new TLegend(0.4,0.45,0.7,0.7);
legEta->AddEntry(grEta," #pi TOF cut","lp");
legEta->AddEntry(pionSigEtaTO," #pi TOF only","lp");
legEta->AddEntry(grkpEta," k+p TOF cut","lp");
legEta->AddEntry(grEEta," e TOF cut","lp");
legEta->SetTextSize(0.035);
legEta->Draw();
cpairEta->Update();
cpairEta->SaveAs("Plots/pairpurityStTOF_EtaWithTOFonly.pdf");




//-------------------------------


ofstream fpairpurity;
fpairpurity.open("textFiles/pairpurity.txt");
fpairpurity<<"/////////// Signal fractions ///////////////// "<<endl;
fpairpurity<<"//pion"<<endl;
fpairpurity<<"double pionPairPtGt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<pionPairPtGt[i]<<"};"<<endl;
  else {fpairpurity<<pionPairPtGt[i]<<",";}
}
fpairpurity<<"double pionPairPtLt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<pionPairPtLt[i]<<"};"<<endl;
  else {fpairpurity<<pionPairPtLt[i]<<",";}
}
fpairpurity<<"double pionPairMGt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<pionPairMGt[i]<<"};"<<endl;
  else {fpairpurity<<pionPairMGt[i]<<",";}
}

fpairpurity<<"double pionPairMLt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<pionPairMLt[i]<<"};"<<endl;
  else {fpairpurity<<pionPairMLt[i]<<",";}
}

fpairpurity<<"//kaon+proton"<<endl;
fpairpurity<<"double kpPairPtGt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<kpPairPtGt[i]<<"};"<<endl;
  else {fpairpurity<<kpPairPtGt[i]<<",";}
}
fpairpurity<<"double kpPairPtLt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<kpPairPtLt[i]<<"};"<<endl;
  else {fpairpurity<<kpPairPtLt[i]<<",";}
}
fpairpurity<<"double kpPairMGt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<kpPairMGt[i]<<"};"<<endl;
  else {fpairpurity<<kpPairMGt[i]<<",";}
}
fpairpurity<<"double kpPairMLt[5]={";
for(int i=0; i<5; i++){
  if(i==4)fpairpurity<<kpPairMLt[i]<<"};"<<endl;
  else {fpairpurity<<kpPairMLt[i]<<",";}
}

fpairpurity<<"// pair purity in eta bins"<<endl;
fpairpurity<<"double pionPairEta[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<pionPairEta[i]<<"};"<<endl;
  else {fpairpurity<<pionPairEta[i]<<",";}
}

fpairpurity<<"double kpPairEta[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<kpPairEta[i]<<"};"<<endl;
  else {fpairpurity<<kpPairEta[i]<<",";}
}
fpairpurity<<"double electronPairEta[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<ePairEta[i]<<"};"<<endl;
  else {fpairpurity<<ePairEta[i]<<",";}
}
fpairpurity<<"// pair purity in integrated pT bins"<<endl;
fpairpurity<<"double pionPairIntMGt[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<pionPairIntMGt[i]<<"};"<<endl;
  else {fpairpurity<<pionPairIntMGt[i]<<",";}
}
fpairpurity<<"double pionPairIntMLt[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<pionPairIntMLt[i]<<"};"<<endl;
  else {fpairpurity<<pionPairIntMLt[i]<<",";}
}

fpairpurity<<"double kpPairIntMGt[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<kpPairIntMGt[i]<<"};"<<endl;
  else {fpairpurity<<kpPairIntMGt[i]<<",";}
}
fpairpurity<<"double kpPairIntMLt[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<kpPairIntMLt[i]<<"};"<<endl;
  else {fpairpurity<<kpPairIntMLt[i]<<",";}
}
fpairpurity<<"double electronPairIntMGt[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<ePairIntMGt[i]<<"};"<<endl;
  else {fpairpurity<<ePairIntMGt[i]<<",";}
}
fpairpurity<<"double electronPairIntMLt[9]={";
for(int i=0; i<9; i++){
  if(i==8)fpairpurity<<ePairIntMLt[i]<<"};"<<endl;
  else {fpairpurity<<ePairIntMLt[i]<<",";}
}
//PID table 
ofstream ftable;
ftable.open("textFiles/pidtable.txt");

ftable<<"//Columns: 1. range 2. particle 3. +ve fraction 4. -ve fraction 5. pair fraction 6. +ve fraction 7. -ve fraction 8. pair fraction"<<endl;
ftable<<"//Columns 3,4,5 for eta>0 and columns 6,7,8 for eta<0 "<<endl;
ftable<<"// pT bins-------- "<<endl;
for(int i=0; i<5; i++){
  ftable<<pT[i]<<"-"<<pT[i+1]<<" & "<<" pi "<<" & "<< Form("%4.4f",pionPtbin[0][0][i])<<" & "<<  Form("%4.4f",pionPtbin[1][0][i])<<" & "<< Form("%4.4f",pionPairPtGt[i])<< " & "<<Form("%4.4f",pionPtbin[0][1][i])<<" & "<<  Form("%4.4f",pionPtbin[1][1][i])<<" & "<< Form("%4.4f",pionPairPtLt[i])<<" \\\\"<<endl;
  ftable<<"\t & "<<" K+p "<<" & "<< Form("%4.4f",protonPtbin[0][0][i])<<" & "<<  Form("%4.4f",protonPtbin[1][0][i])<<" & "<< Form("%4.4f",kpPairPtGt[i])<< " & "<<Form("%4.4f",protonPtbin[0][1][i])<<" & "<<  Form("%4.4f",protonPtbin[1][1][i])<<" & "<< Form("%4.4f",kpPairPtLt[i])<<" \\\\"<<endl;
  ftable<<" \t & "<<" e "<<" & "<< Form("%4.4f",electronPtbin[0][0][i])<<" & "<<  Form("%4.4f",electronPtbin[1][0][i])<<" & "<< Form("%4.4f",ePairPtGt[i])<< " & "<<Form("%4.4f",electronPtbin[0][1][i])<<" & "<<  Form("%4.4f",electronPtbin[1][1][i])<<" & "<< Form("%4.4f",ePairPtLt[i])<<" \\\\"<<endl;
}
ftable<<"// M bins-------- "<<endl;
for(int i=0; i<5; i++){
  ftable<<M[i]<<"-"<<M[i+1]<<" & "<<" pi "<<" & "<< Form("%4.4f",pionMbin[0][0][i])<<" & "<<  Form("%4.4f",pionMbin[1][0][i])<<" & "<< Form("%4.4f",pionPairMGt[i])<< " & "<<Form("%4.4f",pionMbin[0][1][i])<<" & "<<  Form("%4.4f",pionMbin[1][1][i])<<" & "<< Form("%4.4f",pionPairMLt[i])<<" \\\\"<<endl;
  ftable<<"\t & "<<" K+p "<<" & "<< Form("%4.4f",protonMbin[0][0][i])<<" & "<<  Form("%4.4f",protonMbin[1][0][i])<<" & "<< Form("%4.4f",kpPairMGt[i])<< " & "<<Form("%4.4f",protonMbin[0][1][i])<<" & "<<  Form("%4.4f",protonMbin[1][1][i])<<" & "<< Form("%4.4f",kpPairMLt[i])<<" \\\\"<<endl;
  ftable<<" \t & "<<" e "<<" & "<< Form("%4.4f",electronMbin[0][0][i])<<" & "<<  Form("%4.4f",electronMbin[1][0][i])<<" & "<< Form("%4.4f",ePairMGt[i])<< " & "<<Form("%4.4f",electronMbin[0][1][i])<<" & "<<  Form("%4.4f",electronMbin[1][1][i])<<" & "<< Form("%4.4f",ePairMLt[i])<<" \\\\"<<endl;
}
ftable<<"// Eta bins-------- "<<endl;
for(int i=0; i<9; i++){
  ftable<<eta_range[i]<<"-"<<eta_range[i+1]<<" & "<<" pi "<<" & "<< Form("%4.4f",pionEtabin[0][i])<<" & "<<  Form("%4.4f",pionEtabin[1][i])<<" & "<< Form("%4.4f",pionPairEta[i])<<" \\\\"<<endl;
  ftable<<"\t & "<<" K+p "<<" & "<< Form("%4.4f",protonEtabin[0][i])<<" & "<<  Form("%4.4f",protonEtabin[1][i])<<" & "<< Form("%4.4f",kpPairEta[i])<<" \\\\"<<endl;
  ftable<<" \t & "<<" e "<<" & "<< Form("%4.4f",electronEtabin[0][i])<<" & "<<  Form("%4.4f",electronEtabin[1][i])<<" & "<< Form("%4.4f",ePairEta[i])<<" \\\\"<<endl;
}

ftable<<"// Int pT  bins-------- "<<endl;
for(int i=0; i<9; i++){
  ftable<<M_in[i]<<"-"<<M_in[i+1]<<" & "<<" pi "<<" & "<< Form("%4.4f",pionIntMbin[0][0][i])<<" & "<<  Form("%4.4f",pionIntMbin[1][0][i])<<" & "<< Form("%4.4f",pionPairIntMGt[i])<< " & "<<Form("%4.4f",pionIntMbin[0][1][i])<<" & "<<  Form("%4.4f",pionIntMbin[1][1][i])<<" & "<< Form("%4.4f",pionPairIntMLt[i])<<" \\\\"<<endl;
  ftable<<"\t & "<<" K+p "<<" & "<< Form("%4.4f",protonIntMbin[0][0][i])<<" & "<<  Form("%4.4f",protonIntMbin[1][0][i])<<" & "<< Form("%4.4f",kpPairIntMGt[i])<< " & "<<Form("%4.4f",protonIntMbin[0][1][i])<<" & "<<  Form("%4.4f",protonIntMbin[1][1][i])<<" & "<< Form("%4.4f",kpPairIntMLt[i])<<" \\\\"<<endl;
  ftable<<" \t & "<<" e "<<" & "<< Form("%4.4f",electronIntMbin[0][0][i])<<" & "<<  Form("%4.4f",electronIntMbin[1][0][i])<<" & "<< Form("%4.4f",ePairIntMGt[i])<< " & "<<Form("%4.4f",electronIntMbin[0][1][i])<<" & "<<  Form("%4.4f",electronIntMbin[1][1][i])<<" & "<< Form("%4.4f",ePairIntMLt[i])<<" \\\\"<<endl;
}


}//main
