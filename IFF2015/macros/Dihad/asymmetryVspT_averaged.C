//This code is modefied to calculate asymmetry for individual beams and averaged for total asymmetry
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace std;

ofstream Output;


void asymmetryVspT_averaged(const char *ifile="/star/u/pokhrel/GPFS/IFF_TREES/StartLessTOF/dihadNtuples.root"){
	//Histograms to get mean polarization for yellow and blue beam 
	TH1D *hpolB = new TH1D("hpolB", "", 100, 0, 1);
	TH1D *hpolY = new TH1D("hpolY", "", 100, 0, 1);
	//pT distribution for each pT bin for average pT
	TFile *chi2f=new TFile("chi2.root","RECREATE");
	TH1D *hMinv[5];
	TH1D *hMinvGtB[5];
        TH1D *hMinvLtB[5];
        TH1D *hMinvGtY[5];
        TH1D *hMinvLtY[5];
	for(Int_t n=0; n<5; n++){hMinv[n]=new TH1D(Form("Minv_mbin_%i",n), "",100, 0, 5);}
        for(Int_t n=0; n<5; n++){hMinvGtB[n]=new TH1D(Form("MinvGtB_mbin_%i",n), "",100, 0, 5);}
        for(Int_t n=0; n<5; n++){hMinvLtB[n]=new TH1D(Form("MinvLtB_mbin_%i",n), "",100, 0, 5);}
        for(Int_t n=0; n<5; n++){hMinvGtY[n]=new TH1D(Form("MinvGtY_mbin_%i",n), "",100, 0, 5);}
        for(Int_t n=0; n<5; n++){hMinvLtY[n]=new TH1D(Form("MinvLtY_mbin_%i",n), "",100, 0, 5);}
	TH1F *chi2NdfBGt=new TH1F("chi2NdfBGt","Blue eta > 0", 100, 0, 50);
	TH1F *chi2NdfBLt=new TH1F("chi2NdfBLt","Blue eta < 0", 100, 0, 50);
	TH1F *chi2NdfYGt=new TH1F("chi2NdfYGt","Yellow eta > 0", 100, 0, 50);
	TH1F *chi2NdfYLt=new TH1F("chi2NdfYLt","Yellow eta < 0", 100, 0, 50);
		
	TFile *f = new TFile(ifile);
	TTree *ntuple1tpc = (TTree*)f->Get("ntuple1tpc"); //get trees 
	TTree *ntuple2tpc = (TTree*)f->Get("ntuple2tpc");
	TTree *ntuple4tpc = (TTree*)f->Get("ntuple4tpc");
	TTree *ntuple6 = (TTree*)f->Get("ntuple6");
	//define variables to hold trees variables 
Output.open("AsymVspTTPConly.txt");
	float eta_pair;
	float fspinconfig;
	float cone;
	float pT_pair;
	float Minv;
	float PhiRSB;
	float PhiRSY;
	float fitPts_min_pair;
	float polB_corr, polY_corr;
	double pi = 3.14159265359;

	double NpairsBGt[5][9]={0};     double NpairsBLt[5][9] ={0};
        double pTpairsBGt[5][9]={0};    double pTpairsBLt[5][9]={0};
        double MpairsBGt[5][9]={0};     double MpairsBLt[5][9]={0};
        double NpairsYGt[5][9]={0};     double NpairsYLt[5][9]={0};
        double pTpairsYGt[5][9]={0};    double pTpairsYLt[5][9]={0};
        double MpairsYGt[5][9]={0};     double MpairsYLt[5][9]={0};


	double Npairs[9];
	double pTpairs[9];
	double Mpairs[9];
	double pT1[9];//for pT  bin average 
	double pT2[9];
	double pT3[9];
	double pT4[9];
	double pT5[9];
	double avg_Minv[5]={0};
	for(int k=0; k<9; k++) 
	{	
		Npairs[k]=0.;
		pTpairs[k]=0.;
		Mpairs[k]=0.;
		pT1[k]=0;
		pT2[k]=0;
		pT3[k]=0;
		pT4[k]=0;
		pT5[k]=0;
	}

	//get the values of variables from tree 


	ntuple1tpc->SetBranchAddress("fspinconfig",&fspinconfig);
	ntuple1tpc->SetBranchAddress("cone",&cone);
	ntuple2tpc->SetBranchAddress("Minv",&Minv);
	ntuple2tpc->SetBranchAddress("pT_pair",&pT_pair);
	ntuple2tpc->SetBranchAddress("eta_pair",&eta_pair);
	ntuple4tpc->SetBranchAddress("PhiRSB",&PhiRSB);
	ntuple4tpc->SetBranchAddress("PhiRSY",&PhiRSY);
	ntuple6->SetBranchAddress("polB_corr", &polB_corr);
	ntuple6->SetBranchAddress("polY_corr", &polY_corr);


	//To store average polarization values from histograms
	double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;

	Int_t nentries = (Int_t)ntuple1tpc->GetEntries(); //entries of trees 
	cout<<"Entries: "<<nentries<<endl;
	//add friend to the tree. no root file to add friend means the tree is on the same root file 
	ntuple1tpc->AddFriend("ntuple2tpc"); 
	ntuple1tpc->AddFriend("ntuple4tpc");
	//double pT[10]={2.80, 3.47,  3.79, 4.12, 4.49, 4.938,  5.505, 6.300, 7.66, 15 };	
	//double eta_range[10]={-1.200,   -0.668,  -0.469, -0.281, -0.096, 0.089, 0.275, 0.47,   0.675,  1.2 };
	Double_t avgMinv[5]={0};
	//double pT[6]={2.80, 3.727, 4.343, 5.157, 6.529, 15.00};
	//double M[10]={0.250, 0.403, 0.516, 2, 0.711, 0.803, 0.921, 1.070, 1.286, 4.000};
	double M[6]={0.20, 0.4795, 0.6714, 0.8452, 1.1100, 4.00}; 
	//pT binning for each Minv-Bin
	double p1[10]={2.80, 3.432, 3.709, 3.9884, 4.3098, 4.6781, 5.147, 5.809, 6.966, 15};
	double p2[10]={2.80, 3.432, 3.7192, 4.0147, 4.3474, 4.745, 5.255, 5.973, 7.217, 15};    
	double p3[10]={2.80, 3.4294, 3.734, 4.0515, 4.4108, 4.8396, 5.386, 6.153, 7.468, 15}; 
	double p4[10]={2.80, 3.3631, 3.6523, 3.9558, 4.3046, 4.729, 5.2784, 6.056, 7.392, 15}; 
	double p5[10]={2.80, 3.7464, 4.144, 4.538, 4.9705, 5.4766, 6.111, 6.9860, 8.459, 15}; 

/*	
	int ne=(int)ntuple6->GetEntries();
	for (int pol=0; pol<ne;pol++)
	{
		ntuple6->GetEntry(pol);
		hpolB->Fill(polB_corr);			
		hpolY->Fill(polY_corr);			
	}
	
*/
        TCanvas *fitCanvas[5];
	for(Int_t mbin=0;mbin<5;mbin++)
	{
		//if (mbin>0) continue;//for test


		int NpTGtUpB[165]={0};
		int NpTLtUpB[165]={0};
		int NpTGtDnB[165]={0};
		int NpTLtDnB[165]={0};
		int NMinvGtUpB[165]={0};
		int NMinvLtUpB[165]={0};
		int NMinvGtDnB[165]={0};
		int NMinvLtDnB[165]={0};
		int NetaUpB[165]={0};
		int NetaDnB[165]={0};
		int NpTGtUpY[165]={0};
		int NpTLtUpY[165]={0};
		int NpTGtDnY[165]={0};
		int NpTLtDnY[165]={0};
		int NMinvGtUpY[165]={0};
		int NMinvLtUpY[165]={0};
		int NMinvGtDnY[165]={0};
		int NMinvLtDnY[165]={0};
		int NetaUpY[165]={0};
		int NetaDnY[165]={0};
		int NpTGtUpT[165]={0};
		int NpTLtUpT[165]={0};
		int NpTGtDnT[165]={0};
		int NpTLtDnT[165]={0};
		int NMinvGtUpT[165]={0};
		int NMinvLtUpT[165]={0};
		int NMinvGtDnT[165]={0};
		int NMinvLtDnT[165]={0};
		int NetaUpT[165]={0};
		int NetaDnT[165]={0};
	
		double pT[10]={0};
		//choose correspondig invariant mass bin for different  pT-Bin
		if(mbin==0) for(int k=0; k<10; k++){pT[k]=p1[k];}
		if(mbin==1) for(int k=0; k<10; k++){pT[k]=p2[k];}
		if(mbin==2) for(int k=0; k<10; k++){pT[k]=p3[k];}
		if(mbin==3) for(int k=0; k<10; k++){pT[k]=p4[k];}
		if(mbin==4) for(int k=0; k<10; k++){pT[k]=p5[k];}
		
			
			//cout << "mbin: "<< mbin<<", pT[10]={"<<pT[0]<<", "<< pT[1]<<", "<< pT[2]<<", "<< pT[3]<<", "<< pT[4]<<", "<< pT[5]<<", "<< pT[6]<<", "<< pT[7]<<", "<< pT[8]<<", "<< pT[9]<<"}"<<endl;;
		cout << "Entries: "<<nentries<<endl;
		for (Int_t i=0; i < nentries; i++ )
		//for (Int_t i=0; i < 10000; i++ )
		{

			ntuple1tpc -> GetEntry(i);
			if(Minv<M[mbin] || Minv>=M[mbin+1])continue;
			if(cone>0.7)continue;	
			if(Minv>4.)continue;
			if(Minv>0.4876 && Minv<0.5076)  continue;//dodge K0 mass range, doesn't cause asymmetry. Implemention was not working before.This is fine now.
			//if(fitPts_min_pair<15)continue;
			
			hMinv[mbin]->Fill(Minv);
			//cout << "polarizztion values   "<< polB_corr << "  "<<polY_corr <<endl;
			//cout << "selection cuts are good ...working on phi loop.... " << endl; 



			//-------------BLUE-------------------------------------
			//Phi
			for(int phi=0;phi<16;phi++)
			{
				if(PhiRSB>=(phi-8.)/8.*pi && PhiRSB<=(phi-7.)/8.*pi)
				{
					//Invariant mass loop 
					for(int m=0;m<9;m++)
					{

						if(pT_pair>=pT[m] && pT_pair<pT[m+1])
						{
							Npairs[m]=Npairs[m]+1;
							pTpairs[m]=pTpairs[m]+pT_pair;
							Mpairs[m]=Mpairs[m]+Minv;
					       if(eta_pair>0){
                                                        hMinvGtB[mbin]->Fill(Minv);
                                                        NpairsBGt[mbin][m]=NpairsBGt[mbin][m]+1;
                                                        pTpairsBGt[mbin][m]=pTpairsBGt[mbin][m]+pT_pair;
                                                        MpairsBGt[mbin][m]=MpairsBGt[mbin][m]+Minv;
                                                }
                                                if(eta_pair<0){
                                                        hMinvLtB[mbin]->Fill(Minv);
                                                        NpairsBLt[mbin][m]=NpairsBLt[mbin][m]+1;
                                                        pTpairsBLt[mbin][m]=pTpairsBLt[mbin][m]+pT_pair;
                                                        MpairsBLt[mbin][m]=MpairsBLt[mbin][m]+Minv;
                                                }
					
							

							if(fspinconfig==51 || fspinconfig==53)
							{
								if(eta_pair>0)
								{
									NpTGtUpB[m*16+phi]++;
								}
								if(eta_pair<0)
								{
									NpTLtUpB[m*16+phi]++;
								}
							}
							if(fspinconfig==83 || fspinconfig==85)
							{
								if(eta_pair>0)
								{
									NpTGtDnB[m*16+phi]++;
								}
								if(eta_pair<0)
								{
									NpTLtDnB[m*16+phi]++;
								}
							}
						}
					}//pT loop

				}//phi-range loop 
			}//phi loop
			/**************** BLUE BEAM ENDS *****************************////////

			////***************** YELLOW BEAM *****************************////////////
			//Phi loop 
			for(int phi=0;phi<16;phi++)
			{
				if(PhiRSY>=(phi-8.)/8.*pi && PhiRSY<=(phi-7.)/8.*pi)
				{	
					//pT loop
					for(int m=0;m<9;m++)
					{
						if(pT_pair>=pT[m] && pT_pair<pT[m+1])
						{
					       if(eta_pair>0){
                                                        hMinvGtY[mbin]->Fill(Minv);
                                                        NpairsYGt[mbin][m]=NpairsYGt[mbin][m]+1;
                                                        pTpairsYGt[mbin][m]=pTpairsYGt[mbin][m]+pT_pair;
                                                        MpairsYGt[mbin][m]=MpairsYGt[mbin][m]+Minv;
                                                }
                                                if(eta_pair<0){
                                                        hMinvLtY[mbin]->Fill(Minv);
                                                        NpairsYLt[mbin][m]=NpairsYLt[mbin][m]+1;
                                                        pTpairsYLt[mbin][m]=pTpairsYLt[mbin][m]+pT_pair;
                                                        MpairsYLt[mbin][m]=MpairsYLt[mbin][m]+Minv;
                                                }
							if(fspinconfig==51 || fspinconfig==83)
							{
								//if(eta_pair>0)
								if(eta_pair<0)
								{
									NpTGtUpY[m*16+phi]++;
								}
								//if(eta_pair<0)
								if(eta_pair>0)
								{
									NpTLtUpY[m*16+phi]++;
								}
							}
							if(fspinconfig==53 || fspinconfig==85)
							{
								//if(eta_pair>0)
								if(eta_pair<0)
								{
									NpTGtDnY[m*16+phi]++;
								}
								//if(eta_pair<0)
								if(eta_pair>0)
								{
									NpTLtDnY[m*16+phi]++;
								}
							}
						}
					}//pT loop 
				}//PhiRS range loop
			}//phi loop
		}//entries  for loop 
		/******************** YELLOW BEAM ENDS ************************//////
			
		//Get average polarizations and rms from histograms 
		/*avgPolB = hpolB->GetMean();
		avgPolY = hpolY->GetMean();
		rmsB    = hpolB->GetRMS();
		rmsY    = hpolY->GetRMS();
		*/	
		avgPolB=0.57534; rmsB=0.0370551;    //avgPolB = hpolB->GetMean();
                avgPolY=0.58560; rmsY=0.0386434;     //avgPolY = hpolY->GetMean();	
		
	
		//Store pT bin averages for each pT bin in an array for plotting purpose 
		 
                if(mbin==0){for(int k0=0;k0<9;k0++){pT1[k0]=0.5*((pTpairsBGt[0][k0]/(double)NpairsBGt[0][k0])+(pTpairsYGt[0][k0]/(double)NpairsYGt[0][k0]));}}
                if(mbin==1){for(int k0=0;k0<9;k0++){pT2[k0]=0.5*((pTpairsBGt[1][k0]/(double)NpairsBGt[1][k0])+(pTpairsYGt[1][k0]/(double)NpairsYGt[1][k0]));}}
                if(mbin==2){for(int k0=0;k0<9;k0++){pT3[k0]=0.5*((pTpairsBGt[2][k0]/(double)NpairsBGt[2][k0])+(pTpairsYGt[2][k0]/(double)NpairsYGt[2][k0]));}}
                if(mbin==3){for(int k0=0;k0<9;k0++){pT4[k0]=0.5*((pTpairsBGt[3][k0]/(double)NpairsBGt[3][k0])+(pTpairsYGt[3][k0]/(double)NpairsYGt[3][k0]));}}
                if(mbin==4){for(int k0=0;k0<9;k0++){pT5[k0]=0.5*((pTpairsBGt[4][k0]/(double)NpairsBGt[4][k0])+(pTpairsYGt[4][k0]/(double)NpairsYGt[4][k0]));}}
		
		/*if(mbin==0){for(int k0=0;k0<9;k0++){pT1[k0]=pTpairs[k0]/(double)Npairs[k0];}}	
		if(mbin==1){for(int k1=0;k1<9;k1++){pT2[k1]=pTpairs[k1]/(double)Npairs[k1];}}	
		if(mbin==2){for(int k2=0;k2<9;k2++){pT3[k2]=pTpairs[k2]/(double)Npairs[k2];}}	
		if(mbin==3){for(int k3=0;k3<9;k3++){pT4[k3]=pTpairs[k3]/(double)Npairs[k3];}}	
		if(mbin==4){for(int k4=0;k4<9;k4++){pT5[k4]=pTpairs[k4]/(double)Npairs[k4];}}	
		*/
		double dAdB,dAdC,dAdD,dAdE,dAdP;//variables for error calculation
		//Use rms for polarization error from distribution!
		double dP_B = rmsB; // Polarization errors blue
		double dP_Y = rmsY; // Polarization errors yellow 
		double dA[165];// Asymmetry amplitude from sine fit 
		double B,C,D,E; 
		double a,b;
		double Asym[165];
		double pi = 3.14159265359;		
 		fitCanvas[mbin]=new TCanvas(Form("fitCanvas_MBin%i",mbin),"",650,550);
                fitCanvas[mbin]->Print(Form("FitPlots_Mbin_%i_AsymVspTTPConly.pdf(", mbin));
		
		// Asymmetry for BLUE eta >0
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					a = sqrt(NpTGtUpB[m*16+ang]*NpTGtDnB[m*16+ang+8]);
					b = sqrt(NpTGtDnB[m*16+ang]*NpTGtUpB[m*16+ang+8]);
					B = NpTGtUpB[m*16+ang];
					C = NpTGtUpB[m*16+ang+8];
					D = NpTGtDnB[m*16+ang];
					E = NpTGtDnB[m*16+ang+8];
				}
				if(ang>7)
				{
					a = sqrt(NpTGtUpB[m*16+ang]*NpTGtDnB[m*16+ang-8]);
					b = sqrt(NpTGtDnB[m*16+ang]*NpTGtUpB[m*16+ang-8]);
					B = NpTGtUpB[m*16+ang];
					C = NpTGtUpB[m*16+ang-8];
					D = NpTGtDnB[m*16+ang];
					E = NpTGtDnB[m*16+ang-8];		
				}
				Asym[ang]=(1./avgPolB)*((a-b)/(a+b));
				dAdB = (1./avgPolB)*E*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdE = (1./avgPolB)*B*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdD = (-1./avgPolB)*C*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdC = (-1./avgPolB)*D*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdP = (-1./(avgPolB*avgPolB))*(sqrt(B*E)-sqrt(D*C))/(sqrt(B*E)+sqrt(D*C));
				dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_B)**2);
		//cout <<"polB:"<<avgPolB<<", polY: "<<avgPolY<< ", a="<<a<<", b="<<b<<", B="<<B<<", C="<<C<<", D="<<D<<", E="<<E<<", Asym["<<ang<<"]="<<Asym[ang]<<", dA["<<ang<<"]="<<dA[ang]<<endl;	
		}//angle loop

			double Abg[9];// Here Abg refers to Asymmetry for Blue in eta >0
			double deltaAbg[9];
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			gStyle->SetOptDate(0);
			grt = new TGraphErrors(8,angle,Asym,ex,dA);
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
			fit->SetParameter(0,0.0001);
			grt->Fit(fit,"R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			sprintf(title,"Minv bin %i, pT bin %i, BLUE, eta > 0",mbin,m);
			grt->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			cout<<"ChiSquare = "<<chi2Ndf[m]<<endl;
			chi2NdfBGt->Fill(chi2Ndf[m]);
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->CenterTitle();
			grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Abg[m]=fit->GetParameter(0);
			deltaAbg[m]=fit->GetParError(0);
			//sprintf(name,"Minv%i_pTbin%i_etaGt_BLUETPConly.pdf",mbin, m);
			//c1->SaveAs(name);
                fitCanvas[mbin]->Print(Form("FitPlots_Mbin_%i_AsymVspTTPConly.pdf", mbin));

		}//invariant mass loop
		
		// YELLOW eta > 0
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					a = sqrt(NpTGtUpY[m*16+ang]*NpTGtDnY[m*16+ang+8]);
					b = sqrt(NpTGtDnY[m*16+ang]*NpTGtUpY[m*16+ang+8]);
					B = NpTGtUpY[m*16+ang];
					C = NpTGtUpY[m*16+ang+8];
					D = NpTGtDnY[m*16+ang];
					E = NpTGtDnY[m*16+ang+8];
				}
				if(ang>7)
				{
					a = sqrt(NpTGtUpY[m*16+ang]*NpTGtDnY[m*16+ang-8]);
					b = sqrt(NpTGtDnY[m*16+ang]*NpTGtUpY[m*16+ang-8]);
					B = NpTGtUpY[m*16+ang];
					C = NpTGtUpY[m*16+ang-8];
					D = NpTGtDnY[m*16+ang];
					E = NpTGtDnY[m*16+ang-8];
				}
				Asym[ang]=(1./avgPolY)*((a-b)/(a+b));

				dAdB = (1./avgPolY)*E*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdE = (1./avgPolY)*B*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdD = (-1./avgPolY)*C*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdC = (-1./avgPolY)*D*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdP = (-1./(avgPolY*avgPolY))*(sqrt(B*E)-sqrt(D*C))/(sqrt(B*E)+sqrt(D*C));
				dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_Y)**2);
			}
			double Ayg[9];
			double deltaAyg[9];
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			gStyle->SetOptDate(0);
			gr = new TGraphErrors(8,angle,Asym,ex,dA);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
			fit->SetParameter(0,0.0001);
			gr->Fit(fit,"R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			sprintf(title,"Minv bin %i, pT bin %i, YELLOW, eta > 0",mbin, m);
			gr->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			cout<<"ChiSquare= "<<chi2Ndf[m]<<endl;
			chi2NdfYGt->Fill(chi2Ndf[m]);
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle(kTRUE);
			gr->GetYaxis()->SetRangeUser(-0.08,0.08);
			Ayg[m]=fit->GetParameter(0);
			deltaAyg[m]=fit->GetParError(0);
                        fitCanvas[mbin]->Print(Form("FitPlots_Mbin_%i_AsymVspTTPConly.pdf", mbin));

		}//pT loop
			//BLUE beam in eta < 0 
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					a = sqrt(NpTLtUpB[m*16+ang]*NpTLtDnB[m*16+ang+8]);
					b = sqrt(NpTLtDnB[m*16+ang]*NpTLtUpB[m*16+ang+8]);
					B = NpTLtUpB[m*16+ang];
					C = NpTLtUpB[m*16+ang+8];
					D = NpTLtDnB[m*16+ang];
					E = NpTLtDnB[m*16+ang+8];
				}
				if(ang>7)
				{
					a = sqrt(NpTLtUpB[m*16+ang]*NpTLtDnB[m*16+ang-8]);
					b = sqrt(NpTLtDnB[m*16+ang]*NpTLtUpB[m*16+ang-8]);
					B = NpTLtUpB[m*16+ang];
					C = NpTLtUpB[m*16+ang-8];
					D = NpTLtDnB[m*16+ang];
					E = NpTLtDnB[m*16+ang-8];		
				}
				Asym[ang]=(1./avgPolB)*((a-b)/(a+b));
				dAdB = (1./avgPolB)*E*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdE = (1./avgPolB)*B*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdD = (-1./avgPolB)*C*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdC = (-1./avgPolB)*D*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdP = (-1./(avgPolB*avgPolB))*(sqrt(B*E)-sqrt(D*C))/(sqrt(B*E)+sqrt(D*C));
				dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_B)**2);
			}//angle loop

			double Abl[9];
			double deltaAbl[9];  /// bl refers to blue beam eta less than 0
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			gStyle->SetOptDate(0);
			grt = new TGraphErrors(8,angle,Asym,ex,dA);
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
			fit->SetParameter(0,0.0001);
			grt->Fit(fit,"R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			sprintf(title,"Minv bin %i, pT bin %i, BLUE, #eta < 0",mbin, m);
			grt->SetTitle(title);
			chi2Ndf[m] = (double)fit->GetChisquare()/(double)fit->GetNDF();
			chi2NdfBLt->Fill(chi2Ndf[m]);
			cout<<"ChiSquare = "<<chi2Ndf[m]<<endl;
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->SetTitleOffset(1);
			grt->GetYaxis()->CenterTitle(kTRUE);
			grt->GetYaxis()->SetRangeUser(-0.08,0.08);
			Abl[m]=fit->GetParameter(0);
			deltaAbl[m]=fit->GetParError(0);
                	fitCanvas[mbin]->Print(Form("FitPlots_Mbin_%i_AsymVspTTPConly.pdf", mbin));

		}//invariant mass loop

	
		//YELLOW Beam Eta<0 
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					a = sqrt(NpTLtUpY[m*16+ang]*NpTLtDnY[m*16+ang+8]);
					b = sqrt(NpTLtDnY[m*16+ang]*NpTLtUpY[m*16+ang+8]);
					B = NpTLtUpY[m*16+ang];
					C = NpTLtUpY[m*16+ang+8];
					D = NpTLtDnY[m*16+ang];
					E = NpTLtDnY[m*16+ang+8];
				}
				if(ang>7)
				{
					a = sqrt(NpTLtUpY[m*16+ang]*NpTLtDnY[m*16+ang-8]);
					b = sqrt(NpTLtDnY[m*16+ang]*NpTLtUpY[m*16+ang-8]);
					B = NpTLtUpY[m*16+ang];
					C = NpTLtUpY[m*16+ang-8];
					D = NpTLtDnY[m*16+ang];
					E = NpTLtDnY[m*16+ang-8];
				}
				Asym[ang]=(1./avgPolY)*((a-b)/(a+b));

				dAdB = (1./avgPolY)*E*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdE = (1./avgPolY)*B*sqrt(D*C)/(sqrt(B*E)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdD = (-1./avgPolY)*C*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdC = (-1./avgPolY)*D*sqrt(B*E)/(sqrt(D*C)*((sqrt(B*E)+sqrt(D*C))**2));
				dAdP = (-1./(avgPolY*avgPolY))*(sqrt(B*E)-sqrt(D*C))/(sqrt(B*E)+sqrt(D*C));
				dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_Y)**2);
			}
			double Ayl[9];
			double deltaAyl[9];
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			gStyle->SetOptDate(0);
			gr = new TGraphErrors(8,angle,Asym,ex,dA);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
			fit->SetParameter(0,0.0001);
			gr->Fit(fit,"R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			sprintf(title,"Minv bin %i, pT bin %i, YELLOW, #eta < 0",mbin, m);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2NdfYLt->Fill(chi2Ndf[m]);
			cout<<"ChiSquare= "<<chi2Ndf[m]<<endl;
			gr->SetTitle(title);
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle(kTRUE);
			gr->GetYaxis()->SetRangeUser(-0.08,0.08);
			Ayl[m]=fit->GetParameter(0);
			deltaAyl[m]=fit->GetParError(0);
                	fitCanvas[mbin]->Print(Form("FitPlots_Mbin_%i_AsymVspTTPConly.pdf", mbin));
	
		}//pT loop
                fitCanvas[mbin]->Print(Form("FitPlots_Mbin_%i_AsymVspTTPConly.pdf)", mbin));

		//Yellow beam asymmetry ends
	
	
		if(mbin==0)
		{	//for Blue beam eta > 0
			double A_pT1bg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4],Abg[5],Abg[6], Abg[7], Abg[8]};
			double deltaA_pT1bg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4],deltaAbg[5],deltaAbg[6], deltaAbg[7], deltaAbg[8]};
			//for Blue beam eta < 0
			double A_pT1bl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4],Abl[5],Abl[6], Abl[7], Abl[8]};
			double deltaA_pT1bl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4],deltaAbl[5],deltaAbl[6], deltaAbl[7], deltaAbl[8]};
			// for Yellow beam eta >0
			double A_pT1yg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4],Ayg[5],Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT1yg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4],deltaAyg[5],deltaAyg[6], deltaAyg[7], deltaAyg[8]};
			
			double A_pT1yl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4],Ayl[5],Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT1yl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4],deltaAyl[5],deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		
		if(mbin==1)
		{	 //for Blue beam eta > 0
			double A_pT2bg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4],Abg[5],Abg[6], Abg[7], Abg[8]};
			double deltaA_pT2bg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4],deltaAbg[5],deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT2bl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4],Abl[5],Abl[6], Abl[7], Abl[8]};
			double deltaA_pT2bl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4],deltaAbl[5],deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT2yg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4],Ayg[5],Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT2yg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4],deltaAyg[5],deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT2yl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4],Ayl[5],Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT2yl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4],deltaAyl[5],deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		if(mbin==2)
		{
			double A_pT3bg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4],Abg[5],Abg[6], Abg[7], Abg[8]};
			double deltaA_pT3bg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4],deltaAbg[5],deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT3bl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4],Abl[5],Abl[6], Abl[7], Abl[8]};
			double deltaA_pT3bl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4],deltaAbl[5],deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT3yg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4],Ayg[5],Ayg[6],  Ayg[7], Ayg[8]};
			double deltaA_pT3yg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4],deltaAyg[5],deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT3yl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4],Ayl[5],Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT3yl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4],deltaAyl[5],deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		if(mbin==3)
		{	
			double A_pT4bg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4],Abg[5],Abg[6], Abg[7], Abg[8]};
			double deltaA_pT4bg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4],deltaAbg[5],deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT4bl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4],Abl[5],Abl[6], Abl[7], Abl[8]};
			double deltaA_pT4bl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4],deltaAbl[5],deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT4yg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4],Ayg[5],Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT4yg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4],deltaAyg[5],deltaAyg[6], deltaAyg[7], deltaAyg[8]};   

			double A_pT4yl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4],Ayl[5],Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT4yl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4],deltaAyl[5],deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		if(mbin==4)
		{	
			double A_pT5bg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4],Abg[5],Abg[6], Abg[7], Abg[8]};
			double deltaA_pT5bg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4],deltaAbg[5],deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT5bl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4],Abl[5],Abl[6], Abl[7], Abl[8]};
			double deltaA_pT5bl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4],deltaAbl[5],deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT5yg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4],Ayg[5],Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT5yg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4],deltaAyg[5],deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT5yl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4],Ayl[5],Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT5yl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4],deltaAyl[5],deltaAyl[6], deltaAyg[7], deltaAyg[8]};
		}

	}//mbin bin loop
		double MinvBGt[5][9]={0};            double MinvBLt[5][9]={0};
                double MinvYGt[5][9]={0};            double MinvYLt[5][9]={0};
                double MinvAvgGt[5][9]={0};          double MinvAvgLt[5][9]={0};
                double pTBGt[5][9]={0};         double pTBLt[5][9]={0};
                double pTYGt[5][9]={0};         double pTYLt[5][9]={0};
                double pTAvgGt[5][9]={0};       double pTAvgLt[5][9]={0};
                double avgMinvGt[5]={0};
                double avgMinvLt[5]={0};
                double avg_MinvGtB[5]={0};
                double avg_MinvLtB[5]={0};
                double avg_MinvGtY[5]={0};
                double avg_MinvLtY[5]={0};
		for(Int_t n=0; n<5; n++){avg_Minv[n]=hMinv[n]->GetMean();}
		for(Int_t n=0; n<5; n++){
                  avg_Minv[n]=hMinv[n]->GetMean();
                  avg_MinvGtB[n]=hMinvGtB[n]->GetMean();
                  avg_MinvLtB[n]=hMinvLtB[n]->GetMean();
                  avg_MinvGtY[n]=hMinvGtY[n]->GetMean();
                  avg_MinvLtY[n]=hMinvLtY[n]->GetMean();
                } //cout << "avg_Minv["<<n<<"]="<<avg_Minv[n]<< endl;}
                for(int pbin=0; pbin<5; pbin++){
                  avgMinvGt[pbin]=0.5*(avg_MinvGtB[pbin]+avg_MinvGtY[pbin]);
                  avgMinvLt[pbin]=0.5*(avg_MinvLtB[pbin]+avg_MinvLtY[pbin]);

                  for(int mbin=0; mbin<9; mbin++){
                    MinvBGt[pbin][mbin]=MpairsBGt[pbin][mbin]/(double)NpairsBGt[pbin][mbin];pTBGt[pbin][mbin]=pTpairsBGt[pbin][mbin]/(double)NpairsBGt[pbin][mbin];
                    MinvYGt[pbin][mbin]=MpairsYGt[pbin][mbin]/(double)NpairsYGt[pbin][mbin];pTYGt[pbin][mbin]=pTpairsYGt[pbin][mbin]/(double)NpairsYGt[pbin][mbin];
                    MinvAvgGt[pbin][mbin]=0.5*(MinvBGt[pbin][mbin]+MinvYGt[pbin][mbin]);
                    pTAvgGt[pbin][mbin]=0.5*(pTBGt[pbin][mbin]+pTYGt[pbin][mbin]);

                    MinvBLt[pbin][mbin]=MpairsBLt[pbin][mbin]/(double)NpairsBLt[pbin][mbin];pTBLt[pbin][mbin]=pTpairsBLt[pbin][mbin]/(double)NpairsBLt[pbin][mbin];
                    MinvYLt[pbin][mbin]=MpairsYLt[pbin][mbin]/(double)NpairsYLt[pbin][mbin];pTYLt[pbin][mbin]=pTpairsYLt[pbin][mbin]/(double)NpairsYLt[pbin][mbin];
                    MinvAvgLt[pbin][mbin]=0.5*(MinvBLt[pbin][mbin]+MinvYLt[pbin][mbin]);
                    pTAvgLt[pbin][mbin]=0.5*(pTBLt[pbin][mbin]+pTYLt[pbin][mbin]);
                   }
                  }
        
		Output<<"avgMinvGtB[5]={";
                for(int pbin=0; pbin<5; pbin++){
                        Output<<avg_MinvGtB[pbin]<<",";
                        if(pbin==4)Output<<"};"<<endl;
                  }
                Output<<"avgMinvGtY[5]={";
                for(int pbin=0; pbin<5; pbin++){
                        Output<<avg_MinvGtY[pbin]<<",";
                        if(pbin==4)Output<<"};"<<endl;
                  }
                Output<<"avgMinvLtB[5]={";
                for(int pbin=0; pbin<5; pbin++){
                        Output<<avg_MinvLtB[pbin]<<",";
                        if(pbin==4)Output<<"};"<<endl;
                  }
                Output<<"avgMinvLtY[5]={";
                for(int pbin=0; pbin<5; pbin++){
                        Output<<avg_MinvLtY[pbin]<<",";
                        if(pbin==4)Output<<"};"<<endl;
                  }

                Output<<"avgMinvGt[5]={";
                for(int pbin=0; pbin<5; pbin++){
                        Output<<avgMinvGt[pbin]<<",";
                        if(pbin==4)Output<<"};"<<endl;
                  }
                Output<<"avgMinvLt[5]={";
                for(int pbin=0; pbin<5; pbin++){
                        Output<<avgMinvLt[pbin]<<",";
                        if(pbin==4)Output<<"};"<<endl;
                  }
		Output<<"AvgMinvGt[5][9]={";
                for(int pbin=0; pbin<5; pbin++){
                   Output<<"{";
                  for(int mbin=0; mbin<9; mbin++){
                        Output<<MinvAvgGt[pbin][mbin]<<",";
                        if(mbin==8)Output<<"};"<<endl;
                  }
                }
                Output<<"};"<<endl;
                Output<<"AvgpTGt[5][9]={";
                for(int pbin=0; pbin<5; pbin++){
                   Output<<"{";
                  for(int mbin=0; mbin<9; mbin++){
                        Output<<pTAvgGt[pbin][mbin]<<",";
                        if(mbin==8)Output<<"};"<<endl;
                  }
                }
                Output<<"};"<<endl;
                Output<<"AvgMinvLt[5][9]={";
                for(int pbin=0; pbin<5; pbin++){
                   Output<<"{";
                  for(int mbin=0; mbin<9; mbin++){
                        Output<<MinvAvgLt[pbin][mbin]<<",";
                        if(mbin==8)Output<<"};"<<endl;
                  }
                }
                Output<<"};"<<endl;
                Output<<"AvgpTLt[5][9]={";
                for(int pbin=0; pbin<5; pbin++){
                   Output<<"{";
                  for(int mbin=0; mbin<9; mbin++){
                        Output<<pTAvgLt[pbin][mbin]<<",";
                        if(mbin==8)Output<<"};"<<endl;
                  }
                }
                Output<<"};"<<endl;




			//calculate average asymmetry  
			double avgA_pT1g[9] ={0};
			double avgA_pT1l[9] ={0};
			double avgA_pT2g[9] ={0};
			double avgA_pT2l[9] ={0};
			double avgA_pT3g[9] ={0};
			double avgA_pT3l[9] ={0};
			double avgA_pT4g[9] ={0};
			double avgA_pT4l[9] ={0};
			double avgA_pT5g[9] ={0};
			double avgA_pT5l[9] ={0};

			double WavgA_pT1g[9] ={0};   double WerrA_pT1g[9] ={0};
                        double WavgA_pT1l[9] ={0};   double WerrA_pT1l[9] ={0};
                        double WavgA_pT2g[9] ={0};   double WerrA_pT2g[9] ={0};
                        double WavgA_pT2l[9] ={0};   double WerrA_pT2l[9] ={0};
                        double WavgA_pT3g[9] ={0};   double WerrA_pT3g[9] ={0};
                        double WavgA_pT3l[9] ={0};   double WerrA_pT3l[9] ={0};
                        double WavgA_pT4g[9] ={0};   double WerrA_pT4g[9] ={0};
                        double WavgA_pT4l[9] ={0};   double WerrA_pT4l[9] ={0};
                        double WavgA_pT5g[9] ={0};   double WerrA_pT5g[9] ={0};
                        double WavgA_pT5l[9] ={0};   double WerrA_pT5l[9] ={0};


				for (Int_t ii=0; ii<9; ii++)
				{
					avgA_pT1g[ii]= (A_pT1bg[ii]+A_pT1yg[ii])/2.;
					avgA_pT1l[ii]= (A_pT1bl[ii]+A_pT1yl[ii])/2.;
					avgA_pT2g[ii]= (A_pT2bg[ii]+A_pT2yg[ii])/2.;
					avgA_pT2l[ii]= (A_pT2bl[ii]+A_pT2yl[ii])/2.;
					avgA_pT3g[ii]= (A_pT3bg[ii]+A_pT3yg[ii])/2.;
					avgA_pT3l[ii]= (A_pT3bl[ii]+A_pT3yl[ii])/2.;
					avgA_pT4g[ii]= (A_pT4bg[ii]+A_pT4yg[ii])/2.;
					avgA_pT4l[ii]= (A_pT4bl[ii]+A_pT4yl[ii])/2.;
					avgA_pT5g[ii]= (A_pT5bg[ii]+A_pT5yg[ii])/2.;
					avgA_pT5l[ii]= (A_pT5bl[ii]+A_pT5yl[ii])/2.;

					 WavgA_pT1g[ii]= (A_pT1bg[ii]*(1/pow(deltaA_pT1bg[ii],2))+A_pT1yg[ii]*(1/pow(deltaA_pT1yg[ii],2)))/((1/pow(deltaA_pT1bg[ii],2))+(1/pow(deltaA_pT1yg[ii],2)));
                                        WavgA_pT2g[ii]= (A_pT2bg[ii]*(1/pow(deltaA_pT2bg[ii],2))+A_pT2yg[ii]*(1/pow(deltaA_pT2yg[ii],2)))/((1/pow(deltaA_pT2bg[ii],2))+(1/pow(deltaA_pT2yg[ii],2)));
                                        WavgA_pT3g[ii]= (A_pT3bg[ii]*(1/pow(deltaA_pT3bg[ii],2))+A_pT3yg[ii]*(1/pow(deltaA_pT3yg[ii],2)))/((1/pow(deltaA_pT3bg[ii],2))+(1/pow(deltaA_pT3yg[ii],2)));
                                        WavgA_pT4g[ii]= (A_pT4bg[ii]*(1/pow(deltaA_pT4bg[ii],2))+A_pT4yg[ii]*(1/pow(deltaA_pT4yg[ii],2)))/((1/pow(deltaA_pT4bg[ii],2))+(1/pow(deltaA_pT4yg[ii],2)));
                                        WavgA_pT5g[ii]= (A_pT5bg[ii]*(1/pow(deltaA_pT5bg[ii],2))+A_pT5yg[ii]*(1/pow(deltaA_pT5yg[ii],2)))/((1/pow(deltaA_pT5bg[ii],2))+(1/pow(deltaA_pT5yg[ii],2)));
					
					WerrA_pT1g[ii]= 1/sqrt((1/pow(deltaA_pT1bg[ii],2))+(1/pow(deltaA_pT1yg[ii],2)));//total asym err, pT bin 1, eta > 0
                                        WerrA_pT2g[ii]= 1/sqrt((1/pow(deltaA_pT2bg[ii],2))+(1/pow(deltaA_pT2yg[ii],2)));
                                        WerrA_pT3g[ii]= 1/sqrt((1/pow(deltaA_pT3bg[ii],2))+(1/pow(deltaA_pT3yg[ii],2)));
                                        WerrA_pT4g[ii]= 1/sqrt((1/pow(deltaA_pT4bg[ii],2))+(1/pow(deltaA_pT4yg[ii],2)));
                                        WerrA_pT5g[ii]= 1/sqrt((1/pow(deltaA_pT5bg[ii],2))+(1/pow(deltaA_pT5yg[ii],2)));

                                        WavgA_pT1l[ii]= (A_pT1bl[ii]*(1/pow(deltaA_pT1bl[ii],2))+A_pT1yl[ii]*(1/pow(deltaA_pT1yl[ii],2)))/((1/pow(deltaA_pT1bl[ii],2))+(1/pow(deltaA_pT1yl[ii],2)));
                                        WavgA_pT2l[ii]= (A_pT2bl[ii]*(1/pow(deltaA_pT2bl[ii],2))+A_pT2yl[ii]*(1/pow(deltaA_pT2yl[ii],2)))/((1/pow(deltaA_pT2bl[ii],2))+(1/pow(deltaA_pT2yl[ii],2)));
                                        WavgA_pT3l[ii]= (A_pT3bl[ii]*(1/pow(deltaA_pT3bl[ii],2))+A_pT3yl[ii]*(1/pow(deltaA_pT3yl[ii],2)))/((1/pow(deltaA_pT3bl[ii],2))+(1/pow(deltaA_pT3yl[ii],2)));
                                        WavgA_pT4l[ii]= (A_pT4bl[ii]*(1/pow(deltaA_pT4bl[ii],2))+A_pT4yl[ii]*(1/pow(deltaA_pT4yl[ii],2)))/((1/pow(deltaA_pT4bl[ii],2))+(1/pow(deltaA_pT4yl[ii],2)));
                                        WavgA_pT5l[ii]= (A_pT5bl[ii]*(1/pow(deltaA_pT5bl[ii],2))+A_pT5yl[ii]*(1/pow(deltaA_pT5yl[ii],2)))/((1/pow(deltaA_pT5bl[ii],2))+(1/pow(deltaA_pT5yl[ii],2)));


                                        WerrA_pT1l[ii]= 1/sqrt((1/pow(deltaA_pT1bl[ii],2))+(1/pow(deltaA_pT1yl[ii],2)));//total asym err, pT bin 1, eta < 0
                                        WerrA_pT2l[ii]= 1/sqrt((1/pow(deltaA_pT2bl[ii],2))+(1/pow(deltaA_pT2yl[ii],2)));
                                        WerrA_pT3l[ii]= 1/sqrt((1/pow(deltaA_pT3bl[ii],2))+(1/pow(deltaA_pT3yl[ii],2)));
                                        WerrA_pT4l[ii]= 1/sqrt((1/pow(deltaA_pT4bl[ii],2))+(1/pow(deltaA_pT4yl[ii],2)));
                                        WerrA_pT5l[ii]= 1/sqrt((1/pow(deltaA_pT5bl[ii],2))+(1/pow(deltaA_pT5yl[ii],2)));



				}	
			//total asymmetry error calculation
			//since two asymmetry(BLUE+YELLOW) are averaged, error propagation formaula is used for total asymmetry error 
			//if function, f = (a+b)/2, Then Error, df =1/2 √{(∂f/∂a)^2*(∆a)^2+(∂f/∂b)^2*(∆b)^2}
			double errA_pT1g[9] ={0};
			double errA_pT1l[9] ={0};
			double errA_pT2g[9] ={0};
			double errA_pT2l[9] ={0};
			double errA_pT3g[9] ={0};
			double errA_pT3l[9] ={0};
			double errA_pT4g[9] ={0};
			double errA_pT4l[9] ={0};
			double errA_pT5g[9] ={0};
			double errA_pT5l[9] ={0};
				for (Int_t iii=0; iii<9; iii++)
				{
					errA_pT1g[iii]= .5*sqrt(pow(deltaA_pT1bg[iii],2)+pow(deltaA_pT1yg[iii],2));//total asym err, pT bin 1, eta > 0
					errA_pT1l[iii]= .5*sqrt(pow(deltaA_pT1bl[iii],2)+pow(deltaA_pT1yl[iii],2));//total asym err, pT bin 1, eta < 0
					errA_pT2g[iii]= .5*sqrt(pow(deltaA_pT2bg[iii],2)+pow(deltaA_pT2yg[iii],2));
					errA_pT2l[iii]= .5*sqrt(pow(deltaA_pT2bl[iii],2)+pow(deltaA_pT2yl[iii],2));
					errA_pT3g[iii]= .5*sqrt(pow(deltaA_pT3bg[iii],2)+pow(deltaA_pT3yg[iii],2));
					errA_pT3l[iii]= .5*sqrt(pow(deltaA_pT3bl[iii],2)+pow(deltaA_pT3yl[iii],2));
					errA_pT4g[iii]= .5*sqrt(pow(deltaA_pT4bg[iii],2)+pow(deltaA_pT4yg[iii],2));
					errA_pT4l[iii]= .5*sqrt(pow(deltaA_pT4bl[iii],2)+pow(deltaA_pT4yl[iii],2));
					errA_pT5g[iii]= .5*sqrt(pow(deltaA_pT5bg[iii],2)+pow(deltaA_pT5yg[iii],2));
					errA_pT5l[iii]= .5*sqrt(pow(deltaA_pT5bl[iii],2)+pow(deltaA_pT5yl[iii],2));
				}	
		
			
  Output << "********************   FINAL ASYMMETRY , CONE > 0.7, pT Binning ****************" << endl;

 			Output << "<<<<<<<<<<<<<<<   pT Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" <<endl;
	Output << "double  pT1[9]= {"<<pT1[0]<<", "<<pT1[1]<<", "<<pT1[2]<<", "<<pT1[3]<<", "<<pT1[4]<< ", "<<pT1[5]<<", "<<pT1[6]<<", "<<pT1[7]<<", "<<pT1[8]<<"}"<<endl;	
	Output << "double  pT2[9]= {"<<pT2[0]<<", "<<pT2[1]<<", "<<pT2[2]<<", "<<pT2[3]<<", "<<pT2[4]<< ", "<<pT2[5]<<", "<<pT2[6]<<", "<<pT2[7]<<", "<<pT2[8]<<"}"<<endl;	
	Output << "double  pT3[9]= {"<<pT3[0]<<", "<<pT3[1]<<", "<<pT3[2]<<", "<<pT3[3]<<", "<<pT3[4]<< ", "<<pT3[5]<<", "<<pT3[6]<<", "<<pT3[7]<<", "<<pT3[8]<<"}"<<endl;	
	Output << "double  pT4[9]= {"<<pT4[0]<<", "<<pT4[1]<<", "<<pT4[2]<<", "<<pT4[3]<<", "<<pT4[4]<< ", "<<pT4[5]<<", "<<pT4[6]<<", "<<pT4[7]<<", "<<pT4[8]<<"}"<<endl;	
	Output << "double  pT5[9]= {"<<pT5[0]<<", "<<pT5[1]<<", "<<pT5[2]<<", "<<pT5[3]<<", "<<pT5[4]<< ", "<<pT5[5]<<", "<<pT5[6]<<", "<<pT5[7]<<", "<<pT5[8]<<"}"<<endl;	
			
			Output << "<<<<<<<< Average Polarization Values >>>>>>>>" <<endl;
	Output << "BLUE <P>: "<< avgPolB << ", Err_P: "<<rmsB<< endl;
	Output << "YELLOW <P>: "<< avgPolY << ", Err_P: "<<rmsY<< endl;



		Output << "<<<<<<<<<<<       Eta < 0, Standard Average A_{UT}   >>>>>>>>>>>>>>>>>"<<endl;
  Output<< "double A_pTgt1[9] ={"<<avgA_pT1l[0]<<", "<<avgA_pT1l[1]<<", "<<avgA_pT1l[2]<<", "<<avgA_pT1l[3]<<", "<<avgA_pT1l[4]<<", "<<avgA_pT1l[5]<<", "<<avgA_pT1l[6]<<", "<<avgA_pT1l[7]<<", "<<avgA_pT1l[8]<<"}"<<endl;
  Output<< "double Aerr_pTgt1[9]={"<<errA_pT1l[0]<<", "<<errA_pT1l[1]<<", "<<errA_pT1l[2]<<", "<<errA_pT1l[3]<<", "<<errA_pT1l[4]<<", "<<errA_pT1l[5]<<", "<<errA_pT1l[6]<<", "<<errA_pT1l[7]<<", "<<errA_pT1l[8]<<"}"<<endl;
  Output<< "double A_pTgt2[9]={"<<avgA_pT2l[0]<<", "<<avgA_pT2l[1]<<", "<<avgA_pT2l[2]<<", "<<avgA_pT2l[3]<<", "<<avgA_pT2l[4]<<", "<<avgA_pT2l[5]<<", "<<avgA_pT2l[6]<<", "<<avgA_pT2l[7]<<", "<<avgA_pT2l[8]<<"}"<<endl;
  Output<< "double Aerr_pTgt2[9]={"<<errA_pT2l[0]<<", "<<errA_pT2l[1]<<", "<<errA_pT2l[2]<<", "<<errA_pT2l[3]<<", "<<errA_pT2l[4]<<", "<<errA_pT2l[5]<<", "<<errA_pT2l[6]<<", "<<errA_pT2l[7]<<", "<<errA_pT2l[8]<<"}"<<endl;
  Output<< "double A_pTgt3[9] ={"<<avgA_pT3l[0]<<", "<<avgA_pT3l[1]<<", "<<avgA_pT3l[2]<<", "<<avgA_pT3l[3]<<", "<<avgA_pT3l[4]<<", "<<avgA_pT3l[5]<<", "<<avgA_pT3l[6]<<", "<<avgA_pT3l[7]<<", "<<avgA_pT3l[8]<<"}"<<endl;
  Output<< "double Aerr_pTgt3[9] ={"<<errA_pT3l[0]<<", "<<errA_pT3l[1]<<", "<<errA_pT3l[2]<<", "<<errA_pT3l[3]<<", "<<errA_pT3l[4]<<", "<<errA_pT3l[5]<<", "<<errA_pT3l[6]<<", "<<errA_pT3l[7]<<", "<<errA_pT3l[8]<<"}"<<endl;
  Output<< "double A_pTgt4[9] ={"<<avgA_pT4l[0]<<", "<<avgA_pT4l[1]<<", "<<avgA_pT4l[2]<<", "<<avgA_pT4l[3]<<", "<<avgA_pT4l[4]<<", "<<avgA_pT4l[5]<<", "<<avgA_pT4l[6]<<", "<<avgA_pT4l[7]<<", "<<avgA_pT4l[8]<<"}"<<endl;
  Output<< "double Aerr_pTgt4[9] ={"<<errA_pT4l[0]<<", "<<errA_pT4l[1]<<", "<<errA_pT4l[2]<<", "<<errA_pT4l[3]<<", "<<errA_pT4l[4]<<", "<<errA_pT4l[5]<<", "<<errA_pT4l[6]<<", "<<errA_pT4l[7]<<", "<<errA_pT4l[8]<<"}"<<endl;
  Output<< "double A_pTgt5[9]={"<<avgA_pT5l[0]<<", "<<avgA_pT5l[1]<<", "<<avgA_pT5l[2]<<", "<<avgA_pT5l[3]<<", "<<avgA_pT5l[4]<<", "<<avgA_pT5l[5]<<", "<<avgA_pT5l[6]<<", "<<avgA_pT5l[7]<<", "<<avgA_pT5l[8]<<"}"<<endl;
  Output<< "double Aerr_pTgt5[9] ={"<<errA_pT5l[0]<<", "<<errA_pT5l[1]<<", "<<errA_pT5l[2]<<", "<<errA_pT5l[3]<<", "<<errA_pT5l[4]<<", "<<errA_pT5l[5]<<", "<<errA_pT5l[6]<<", "<<errA_pT5l[7]<<", "<<errA_pT5l[8]<<"}"<<endl;
	 	
		Output << "<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>"<<endl;
  Output<< "double A_pTlt1[9] ={"<<avgA_pT1g[0]<<", "<<avgA_pT1g[1]<<", "<<avgA_pT1g[2]<<", "<<avgA_pT1g[3]<<", "<<avgA_pT1g[4]<<", "<<avgA_pT1g[5]<<", "<<avgA_pT1g[6]<<", "<<avgA_pT1g[7]<<", "<<avgA_pT1g[8]<<"}"<<endl;
  Output<< "double Aerr_pTlt1[9] ={"<<errA_pT1g[0]<<", "<<errA_pT1g[1]<<", "<<errA_pT1g[2]<<", "<<errA_pT1g[3]<<", "<<errA_pT1g[4]<<", "<<errA_pT1g[5]<<", "<<errA_pT1g[6]<<", "<<errA_pT1g[7]<<", "<<errA_pT1g[8]<<"}"<<endl;
  Output<< "double A_pTlt2[9] ={"<<avgA_pT2g[0]<<", "<<avgA_pT2g[1]<<", "<<avgA_pT2g[2]<<", "<<avgA_pT2g[3]<<", "<<avgA_pT2g[4]<<", "<<avgA_pT2g[5]<<", "<<avgA_pT2g[6]<<", "<<avgA_pT2g[7]<<", "<<avgA_pT2g[8]<<"}"<<endl;
  Output<< "double Aerr_pTlt2[9] ={"<<errA_pT2g[0]<<", "<<errA_pT2g[1]<<", "<<errA_pT2g[2]<<", "<<errA_pT2g[3]<<", "<<errA_pT2g[4]<<", "<<errA_pT2g[5]<<", "<<errA_pT2g[6]<<", "<<errA_pT2g[7]<<", "<<errA_pT2g[8]<<"}"<<endl;
  Output<< "double A_pTlt3[9] ={"<<avgA_pT3g[0]<<", "<<avgA_pT3g[1]<<", "<<avgA_pT3g[2]<<", "<<avgA_pT3g[3]<<", "<<avgA_pT3g[4]<<", "<<avgA_pT3g[5]<<", "<<avgA_pT3g[6]<<", "<<avgA_pT3g[7]<<", "<<avgA_pT3g[8]<<"}"<<endl;
  Output<< "double Aerr_pTlt3[9] ={"<<errA_pT3g[0]<<", "<<errA_pT3g[1]<<", "<<errA_pT3g[2]<<", "<<errA_pT3g[3]<<", "<<errA_pT3g[4]<<", "<<errA_pT3g[5]<<", "<<errA_pT3g[6]<<", "<<errA_pT3g[7]<<", "<<errA_pT3g[8]<<"}"<<endl;
  Output<< "double A_pTlt4[9] ={"<<avgA_pT4g[0]<<", "<<avgA_pT4g[1]<<", "<<avgA_pT4g[2]<<", "<<avgA_pT4g[3]<<", "<<avgA_pT4g[4]<<", "<<avgA_pT4g[5]<<", "<<avgA_pT4g[6]<<", "<<avgA_pT4g[7]<<", "<<avgA_pT4g[8]<<"}"<<endl;
  Output<< "double Aerr_pTlt4[9] ={"<<errA_pT4g[0]<<", "<<errA_pT4g[1]<<", "<<errA_pT4g[2]<<", "<<errA_pT4g[3]<<", "<<errA_pT4g[4]<<", "<<errA_pT4g[5]<<", "<<errA_pT4g[6]<<", "<<errA_pT4g[7]<<", "<<errA_pT4g[8]<<"}"<<endl;
  Output<< "double A_pTlt5[9] ={"<<avgA_pT5g[0]<<", "<<avgA_pT5g[1]<<", "<<avgA_pT5g[2]<<", "<<avgA_pT5g[3]<<", "<<avgA_pT5g[4]<<", "<<avgA_pT5g[5]<<", "<<avgA_pT5g[6]<<", "<<avgA_pT5g[7]<<", "<<avgA_pT5g[8]<<"}"<<endl;
  Output<< "double Aerr_pTlt5[9] ={"<<errA_pT5g[0]<<", "<<errA_pT5g[1]<<", "<<errA_pT5g[2]<<", "<<errA_pT5g[3]<<", "<<errA_pT5g[4]<<", "<<errA_pT5g[5]<<", "<<errA_pT5g[6]<<", "<<errA_pT5g[7]<<", "<<errA_pT5g[8]<<"}"<<endl;
	 	

	 	
	double  ratioMGt1[9]={0.9543,0.9543,0.9543,0.9543,0.9543,0.9543,0.9543,0.9543,0.9543};
 	double  ratioMGt2[9]={1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001};
 	double  ratioMGt3[9]={0.9893,0.9893,0.9893,0.9893,0.9893,0.9893,0.9893,0.9893,0.9893};
 	double  ratioMGt4[9]={1.064,1.064,1.064,1.064,1.064,1.064,1.064,1.064,1.064};
 	double  ratioMGt5[9]={1.063,1.063,1.063,1.063,1.063,1.063,1.063,1.063,1.063};
 
 	double  ratioMLt1[9]={1.056,1.056,1.056,1.056,1.056,1.056,1.056,1.056,1.056};
 	double  ratioMLt2[9]={1.044,1.044,1.044,1.044,1.044,1.044,1.044,1.044,1.044};
 	double  ratioMLt3[9]={0.9292,0.9292,0.9292,0.9292,0.9292,0.9292,0.9292,0.9292,0.9292};
 	double  ratioMLt4[9]={1.049,1.049,1.049,1.049,1.049,1.049,1.049,1.049,1.049};
 	double  ratioMLt5[9]={1.134,1.134,1.134,1.134,1.134,1.134,1.134,1.134,1.134};


	/*//preliminary PID fractions-------
        double pairpurityPtGt[5]={0.67307,0.817963,0.794931,0.765595,0.793264};//prelim PID fractions
        double pairpurityPtLt[5]={0.660449,0.8359,0.808052,0.739741,0.80209};
        double pairpurityMGt[5]={0.778802,0.792258,0.788996,0.779743,0.726311};
        double pairpurityMLt[5]={0.793783,0.804098,0.820309,0.790227,0.73937};
        double pairpurityEta[9]={0.777056,0.798529,0.800554,0.803539,0.802013,0.798089,0.79153,0.784865,0.760883};
	//--------------------------------
	*/
 //final pid fractions after TOF cut---------------------------------
 //pion
 double pionPairPtGt[5]={0.929242,0.911722,0.900156,0.882439,0.853673};
 double pionPairPtLt[5]={0.945092,0.925782,0.912807,0.895107,0.866101};
 double pionPairMGt[5]={0.878427,0.887055,0.895579,0.878562,0.849027};
 double pionPairMLt[5]={0.892918,0.901121,0.907527,0.891572,0.863809};
 //kaon+proton
 double kpPairPtGt[5]={0.00123106,0.00188561,0.00230285,0.00293531,0.00377407};
 double kpPairPtLt[5]={0.000727659,0.00131544,0.00173617,0.00229945,0.00303562};
 double kpPairMGt[5]={0.00298283,0.00287837,0.00237468,0.00333972,0.0051258};
 double kpPairMLt[5]={0.00223246,0.00218268,0.00184161,0.00263674,0.00410077};
 // pair purity in eta bins
 double pionPairEta[9]={0.862209,0.901796,0.902833,0.894983,0.886446,0.890468,0.896209,0.891722,0.844428};
 double kpPairEta[9]={0.00365175,0.00220525,0.00219847,0.00249915,0.00277542,0.00263735,0.00246155,0.00265875,0.00498614};
 double electronPairEta[9]={0.000116222,1.11912e-05,8.49341e-06,1.53427e-05,3.21867e-05,2.39963e-05,1.36267e-05,1.67752e-05,0.000106381};
 //only pion fractions for plotting
 double pairpurityPtGt[5] = {0.929242,0.911722,0.900156,0.882439,0.853673};
 double pairpurityPtLt[5] = {0.945092,0.925782,0.912807,0.895107,0.866101};
 double pairpurityMGt[5]  = {0.878427,0.887055,0.895579,0.878562,0.849027};
 double pairpurityMLt[5]  = {0.892918,0.901121,0.907527,0.891572,0.863809};
 double pairpurityEta[9]  = {0.862209,0.901796,0.902833,0.894983,0.886446,0.890468,0.896209,0.891722,0.844428};
 
 //-------------------------------------------------------------------------
 
        double biasMgt1[9]={0};double biasMgt2[9]={0};double biasMgt3[9]={0};double biasMgt4[9]={0};double biasMgt5[9]={0};
        double biasMlt1[9]={0};double biasMlt2[9]={0};double biasMlt3[9]={0};double biasMlt4[9]={0};double biasMlt5[9]={0};

        double pidMgt1[9]={0};double pidMgt2[9]={0};double pidMgt3[9]={0};double pidMgt4[9]={0};double pidMgt5[9]={0};
        double pidMlt1[9]={0};double pidMlt2[9]={0};double pidMlt3[9]={0};double pidMlt4[9]={0};double pidMlt5[9]={0};  

        double combSysMgt1[9]={0};double combSysMgt2[9]={0};double combSysMgt3[9]={0};double combSysMgt4[9]={0};double combSysMgt5[9]={0};
        double combSysMlt1[9]={0};double combSysMlt2[9]={0};double combSysMlt3[9]={0};double combSysMlt4[9]={0};double combSysMlt5[9]={0};

	for (int pbin=0; pbin<9; pbin++)
        {
        biasMgt1[pbin]=(1-ratioMGt1[pbin])*TMath::Max(WavgA_pT1g[pbin],WerrA_pT1g[pbin]);
        biasMgt2[pbin]=(1-ratioMGt2[pbin])*TMath::Max(WavgA_pT2g[pbin],WerrA_pT2g[pbin]);
        biasMgt3[pbin]=(1-ratioMGt3[pbin])*TMath::Max(WavgA_pT3g[pbin],WerrA_pT3g[pbin]);
        biasMgt4[pbin]=(1-ratioMGt4[pbin])*TMath::Max(WavgA_pT4g[pbin],WerrA_pT4g[pbin]);
        biasMgt5[pbin]=(1-ratioMGt5[pbin])*TMath::Max(WavgA_pT5g[pbin],WerrA_pT5g[pbin]);
                                      
        biasMlt1[pbin]=(1-ratioMLt1[pbin])*TMath::Max(WavgA_pT1l[pbin],WerrA_pT1l[pbin]);
        biasMlt2[pbin]=(1-ratioMLt2[pbin])*TMath::Max(WavgA_pT2l[pbin],WerrA_pT2l[pbin]);
        biasMlt3[pbin]=(1-ratioMLt3[pbin])*TMath::Max(WavgA_pT3l[pbin],WerrA_pT3l[pbin]);
        biasMlt4[pbin]=(1-ratioMLt4[pbin])*TMath::Max(WavgA_pT4l[pbin],WerrA_pT4l[pbin]);
        biasMlt5[pbin]=(1-ratioMLt5[pbin])*TMath::Max(WavgA_pT5l[pbin],WerrA_pT5l[pbin]);
        
        pidMgt1[pbin] =(1-pairpurityMGt[0])*TMath::Max(WavgA_pT1g[pbin],WerrA_pT1g[pbin]);
        pidMgt2[pbin] =(1-pairpurityMGt[1])*TMath::Max(WavgA_pT2g[pbin],WerrA_pT2g[pbin]);
        pidMgt3[pbin] =(1-pairpurityMGt[2])*TMath::Max(WavgA_pT3g[pbin],WerrA_pT3g[pbin]);
        pidMgt4[pbin] =(1-pairpurityMGt[3])*TMath::Max(WavgA_pT4g[pbin],WerrA_pT4g[pbin]);
        pidMgt5[pbin] =(1-pairpurityMGt[4])*TMath::Max(WavgA_pT5g[pbin],WerrA_pT5g[pbin]);
                                                                                                
        pidMlt1[pbin] =(1-pairpurityMLt[0])*TMath::Max(WavgA_pT1l[pbin],WerrA_pT1l[pbin]);
        pidMlt2[pbin] =(1-pairpurityMLt[1])*TMath::Max(WavgA_pT2l[pbin],WerrA_pT2l[pbin]);
        pidMlt3[pbin] =(1-pairpurityMLt[2])*TMath::Max(WavgA_pT3l[pbin],WerrA_pT3l[pbin]);
        pidMlt4[pbin] =(1-pairpurityMLt[3])*TMath::Max(WavgA_pT4l[pbin],WerrA_pT4l[pbin]);
        pidMlt5[pbin] =(1-pairpurityMLt[4])*TMath::Max(WavgA_pT5l[pbin],WerrA_pT5l[pbin]);

        combSysMgt1[pbin]=sqrt(pow(biasMgt1[pbin],2) + pow(pidMgt1[pbin],2));
        combSysMgt2[pbin]=sqrt(pow(biasMgt2[pbin],2) + pow(pidMgt2[pbin],2));
        combSysMgt3[pbin]=sqrt(pow(biasMgt3[pbin],2) + pow(pidMgt3[pbin],2));
        combSysMgt4[pbin]=sqrt(pow(biasMgt4[pbin],2) + pow(pidMgt4[pbin],2));
        combSysMgt5[pbin]=sqrt(pow(biasMgt5[pbin],2) + pow(pidMgt5[pbin],2));
        combSysMlt1[pbin]=sqrt(pow(biasMlt1[pbin],2) + pow(pidMlt1[pbin],2));
        combSysMlt2[pbin]=sqrt(pow(biasMlt2[pbin],2) + pow(pidMlt2[pbin],2));
        combSysMlt3[pbin]=sqrt(pow(biasMlt3[pbin],2) + pow(pidMlt3[pbin],2));
        combSysMlt4[pbin]=sqrt(pow(biasMlt4[pbin],2) + pow(pidMlt4[pbin],2));
        combSysMlt5[pbin]=sqrt(pow(biasMlt5[pbin],2) + pow(pidMlt5[pbin],2));
        }

		Output << "<<<<<<<<<<<       Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>"<<endl;
  Output<< "double A_pTlt1[9] ={"<<WavgA_pT1l[0]<<", "<<WavgA_pT1l[1]<<", "<<WavgA_pT1l[2]<<", "<<WavgA_pT1l[3]<<", "<<WavgA_pT1l[4]<<", "<<WavgA_pT1l[5]<<", "<<WavgA_pT1l[6]<<", "<<WavgA_pT1l[7]<<", "<<WavgA_pT1l[8]<<"}"<<endl;
  Output<< "double AWerr_pTlt1[9]={"<<WerrA_pT1l[0]<<", "<<WerrA_pT1l[1]<<", "<<WerrA_pT1l[2]<<", "<<WerrA_pT1l[3]<<", "<<WerrA_pT1l[4]<<", "<<WerrA_pT1l[5]<<", "<<WerrA_pT1l[6]<<", "<<WerrA_pT1l[7]<<", "<<WerrA_pT1l[8]<<"}"<<endl;
  Output<< "double Sys_pTlt1[9]={"<<biasMlt1[0]<<", "<<biasMlt1[1]<<", "<<biasMlt1[2]<<", "<<biasMlt1[3]<<", "<<biasMlt1[4]<<", "<<biasMlt1[5]<<", "<<biasMlt1[6]<<", "<<biasMlt1[7]<<", "<<biasMlt1[8]<<"}"<<endl;
  Output<< "double A_pTlt2[9]={"<<WavgA_pT2l[0]<<", "<<WavgA_pT2l[1]<<", "<<WavgA_pT2l[2]<<", "<<WavgA_pT2l[3]<<", "<<WavgA_pT2l[4]<<", "<<WavgA_pT2l[5]<<", "<<WavgA_pT2l[6]<<", "<<WavgA_pT2l[7]<<", "<<WavgA_pT2l[8]<<"}"<<endl;
  Output<< "double AWerr_pTlt2[9]={"<<WerrA_pT2l[0]<<", "<<WerrA_pT2l[1]<<", "<<WerrA_pT2l[2]<<", "<<WerrA_pT2l[3]<<", "<<WerrA_pT2l[4]<<", "<<WerrA_pT2l[5]<<", "<<WerrA_pT2l[6]<<", "<<WerrA_pT2l[7]<<", "<<WerrA_pT2l[8]<<"}"<<endl;
  Output<< "double Sys_pTlt2[9]={"<<biasMlt2[0]<<", "<<biasMlt2[1]<<", "<<biasMlt2[2]<<", "<<biasMlt2[3]<<", "<<biasMlt2[4]<<", "<<biasMlt2[5]<<", "<<biasMlt2[6]<<", "<<biasMlt2[7]<<", "<<biasMlt2[8]<<"}"<<endl;
  Output<< "double A_pTlt3[9] ={"<<WavgA_pT3l[0]<<", "<<WavgA_pT3l[1]<<", "<<WavgA_pT3l[2]<<", "<<WavgA_pT3l[3]<<", "<<WavgA_pT3l[4]<<", "<<WavgA_pT3l[5]<<", "<<WavgA_pT3l[6]<<", "<<WavgA_pT3l[7]<<", "<<WavgA_pT3l[8]<<"}"<<endl;
  Output<< "double AWerr_pTlt3[9] ={"<<WerrA_pT3l[0]<<", "<<WerrA_pT3l[1]<<", "<<WerrA_pT3l[2]<<", "<<WerrA_pT3l[3]<<", "<<WerrA_pT3l[4]<<", "<<WerrA_pT3l[5]<<", "<<WerrA_pT3l[6]<<", "<<WerrA_pT3l[7]<<", "<<WerrA_pT3l[8]<<"}"<<endl;
  Output<< "double Sys_pTlt3[9]={"<<biasMlt3[0]<<", "<<biasMlt3[1]<<", "<<biasMlt3[2]<<", "<<biasMlt3[3]<<", "<<biasMlt3[4]<<", "<<biasMlt3[5]<<", "<<biasMlt3[6]<<", "<<biasMlt3[7]<<", "<<biasMlt3[8]<<"}"<<endl;
  Output<< "double A_pTlt4[9] ={"<<WavgA_pT4l[0]<<", "<<WavgA_pT4l[1]<<", "<<WavgA_pT4l[2]<<", "<<WavgA_pT4l[3]<<", "<<WavgA_pT4l[4]<<", "<<WavgA_pT4l[5]<<", "<<WavgA_pT4l[6]<<", "<<WavgA_pT4l[7]<<", "<<WavgA_pT4l[8]<<"}"<<endl;
  Output<< "double AWerr_pTlt4[9] ={"<<WerrA_pT4l[0]<<", "<<WerrA_pT4l[1]<<", "<<WerrA_pT4l[2]<<", "<<WerrA_pT4l[3]<<", "<<WerrA_pT4l[4]<<", "<<WerrA_pT4l[5]<<", "<<WerrA_pT4l[6]<<", "<<WerrA_pT4l[7]<<", "<<WerrA_pT4l[8]<<"}"<<endl;
  Output<< "double Sys_pTlt4[9]={"<<biasMlt4[0]<<", "<<biasMlt4[1]<<", "<<biasMlt4[2]<<", "<<biasMlt4[3]<<", "<<biasMlt4[4]<<", "<<biasMlt4[5]<<", "<<biasMlt4[6]<<", "<<biasMlt4[7]<<", "<<biasMlt4[8]<<"}"<<endl;
  Output<< "double A_pTlt5[9]={"<<WavgA_pT5l[0]<<", "<<WavgA_pT5l[1]<<", "<<WavgA_pT5l[2]<<", "<<WavgA_pT5l[3]<<", "<<WavgA_pT5l[4]<<", "<<WavgA_pT5l[5]<<", "<<WavgA_pT5l[6]<<", "<<WavgA_pT5l[7]<<", "<<WavgA_pT5l[8]<<"}"<<endl;
  Output<< "double AWerr_pTlt5[9] ={"<<WerrA_pT5l[0]<<", "<<WerrA_pT5l[1]<<", "<<WerrA_pT5l[2]<<", "<<WerrA_pT5l[3]<<", "<<WerrA_pT5l[4]<<", "<<WerrA_pT5l[5]<<", "<<WerrA_pT5l[6]<<", "<<WerrA_pT5l[7]<<", "<<WerrA_pT5l[8]<<"}"<<endl;
  Output<< "double Sys_pTlt5[9]={"<<biasMlt5[0]<<", "<<biasMlt5[1]<<", "<<biasMlt5[2]<<", "<<biasMlt5[3]<<", "<<biasMlt5[4]<<", "<<biasMlt5[5]<<", "<<biasMlt5[6]<<", "<<biasMlt5[7]<<", "<<biasMlt5[8]<<"}"<<endl;
	 	
		Output << "<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>"<<endl;
  Output<< "double A_pTgt1[9] ={"<<WavgA_pT1g[0]<<", "<<WavgA_pT1g[1]<<", "<<WavgA_pT1g[2]<<", "<<WavgA_pT1g[3]<<", "<<WavgA_pT1g[4]<<", "<<WavgA_pT1g[5]<<", "<<WavgA_pT1g[6]<<", "<<WavgA_pT1g[7]<<", "<<WavgA_pT1g[8]<<"}"<<endl;
  Output<< "double AWerr_pTgt1[9] ={"<<WerrA_pT1g[0]<<", "<<WerrA_pT1g[1]<<", "<<WerrA_pT1g[2]<<", "<<WerrA_pT1g[3]<<", "<<WerrA_pT1g[4]<<", "<<WerrA_pT1g[5]<<", "<<WerrA_pT1g[6]<<", "<<WerrA_pT1g[7]<<", "<<WerrA_pT1g[8]<<"}"<<endl;
  Output<< "double Sys_pTgt1[9]={"<<biasMgt1[0]<<", "<<biasMgt1[1]<<", "<<biasMgt1[2]<<", "<<biasMgt1[3]<<", "<<biasMgt1[4]<<", "<<biasMgt1[5]<<", "<<biasMgt1[6]<<", "<<biasMgt1[7]<<", "<<biasMgt1[8]<<"}"<<endl;
  Output<< "double A_pTgt2[9] ={"<<WavgA_pT2g[0]<<", "<<WavgA_pT2g[1]<<", "<<WavgA_pT2g[2]<<", "<<WavgA_pT2g[3]<<", "<<WavgA_pT2g[4]<<", "<<WavgA_pT2g[5]<<", "<<WavgA_pT2g[6]<<", "<<WavgA_pT2g[7]<<", "<<WavgA_pT2g[8]<<"}"<<endl;
  Output<< "double AWerr_pTgt2[9] ={"<<WerrA_pT2g[0]<<", "<<WerrA_pT2g[1]<<", "<<WerrA_pT2g[2]<<", "<<WerrA_pT2g[3]<<", "<<WerrA_pT2g[4]<<", "<<WerrA_pT2g[5]<<", "<<WerrA_pT2g[6]<<", "<<WerrA_pT2g[7]<<", "<<WerrA_pT2g[8]<<"}"<<endl;
  Output<< "double Sys_pTgt2[9]={"<<biasMgt2[0]<<", "<<biasMgt2[1]<<", "<<biasMgt2[2]<<", "<<biasMgt2[3]<<", "<<biasMgt2[4]<<", "<<biasMgt2[5]<<", "<<biasMgt2[6]<<", "<<biasMgt2[7]<<", "<<biasMgt2[8]<<"}"<<endl;
  Output<< "double A_pTgt3[9] ={"<<WavgA_pT3g[0]<<", "<<WavgA_pT3g[1]<<", "<<WavgA_pT3g[2]<<", "<<WavgA_pT3g[3]<<", "<<WavgA_pT3g[4]<<", "<<WavgA_pT3g[5]<<", "<<WavgA_pT3g[6]<<", "<<WavgA_pT3g[7]<<", "<<WavgA_pT3g[8]<<"}"<<endl;
  Output<< "double AWerr_pTgt3[9] ={"<<WerrA_pT3g[0]<<", "<<WerrA_pT3g[1]<<", "<<WerrA_pT3g[2]<<", "<<WerrA_pT3g[3]<<", "<<WerrA_pT3g[4]<<", "<<WerrA_pT3g[5]<<", "<<WerrA_pT3g[6]<<", "<<WerrA_pT3g[7]<<", "<<WerrA_pT3g[8]<<"}"<<endl;
  Output<< "double Sys_pTgt3[9]={"<<biasMgt3[0]<<", "<<biasMgt3[1]<<", "<<biasMgt3[2]<<", "<<biasMgt3[3]<<", "<<biasMgt3[4]<<", "<<biasMgt3[5]<<", "<<biasMgt3[6]<<", "<<biasMgt3[7]<<", "<<biasMgt3[8]<<"}"<<endl;
  Output<< "double A_pTgt4[9] ={"<<WavgA_pT4g[0]<<", "<<WavgA_pT4g[1]<<", "<<WavgA_pT4g[2]<<", "<<WavgA_pT4g[3]<<", "<<WavgA_pT4g[4]<<", "<<WavgA_pT4g[5]<<", "<<WavgA_pT4g[6]<<", "<<WavgA_pT4g[7]<<", "<<WavgA_pT4g[8]<<"}"<<endl;
  Output<< "double AWerr_pTgt4[9] ={"<<WerrA_pT4g[0]<<", "<<WerrA_pT4g[1]<<", "<<WerrA_pT4g[2]<<", "<<WerrA_pT4g[3]<<", "<<WerrA_pT4g[4]<<", "<<WerrA_pT4g[5]<<", "<<WerrA_pT4g[6]<<", "<<WerrA_pT4g[7]<<", "<<WerrA_pT4g[8]<<"}"<<endl;
  Output<< "double Sys_pTgt4[9]={"<<biasMgt4[0]<<", "<<biasMgt4[1]<<", "<<biasMgt4[2]<<", "<<biasMgt4[3]<<", "<<biasMgt4[4]<<", "<<biasMgt4[5]<<", "<<biasMgt4[6]<<", "<<biasMgt4[7]<<", "<<biasMgt4[8]<<"}"<<endl;
  Output<< "double A_pTgt5[9] ={"<<WavgA_pT5g[0]<<", "<<WavgA_pT5g[1]<<", "<<WavgA_pT5g[2]<<", "<<WavgA_pT5g[3]<<", "<<WavgA_pT5g[4]<<", "<<WavgA_pT5g[5]<<", "<<WavgA_pT5g[6]<<", "<<WavgA_pT5g[7]<<", "<<WavgA_pT5g[8]<<"}"<<endl;
  Output<< "double AWerr_pTgt5[9] ={"<<WerrA_pT5g[0]<<", "<<WerrA_pT5g[1]<<", "<<WerrA_pT5g[2]<<", "<<WerrA_pT5g[3]<<", "<<WerrA_pT5g[4]<<", "<<WerrA_pT5g[5]<<", "<<WerrA_pT5g[6]<<", "<<WerrA_pT5g[7]<<", "<<WerrA_pT5g[8]<<"}"<<endl;
  Output<< "double Sys_pTgt5[9]={"<<biasMgt5[0]<<", "<<biasMgt5[1]<<", "<<biasMgt5[2]<<", "<<biasMgt5[3]<<", "<<biasMgt5[4]<<", "<<biasMgt5[5]<<", "<<biasMgt5[6]<<", "<<biasMgt5[7]<<", "<<biasMgt5[8]<<"}"<<endl;



	TGraphErrors *gr_AvsMgts[5];
        TGraphErrors *gr_AvsMlts[5];
        TGraphErrors *gr_AvsMgt[5];
        TGraphErrors *gr_AvsMlt[5];
        double pTAvg[9]={0};
        double A_Mgt[9]={0};
        double Aerr_Mgt[9]={0};
        double A_Mlt[9]={0};
        double Aerr_Mlt[9]={0};
	double pidSysGt[9]={0};
        double pidSysLt[9]={0};
        double trigSysGt[9]={0};
        double trigSysLt[9]={0};
        TLegend *leg;
        TLine *line[5];
        TLatex tex[5];

        double finalSysgt[9]={0};
        double finalSyslt[9]={0};
        double errx[9]={0};
        double errxps[9]={0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15};
        double errxts[9]={0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15};
	
	//print values in table format
	  ofstream ftable;
        ftable.open("ifftable_AutVspTTPConly.txt");
        ftable<<"//values are printed in the following order:\n";
        ftable<<"//1.Bin range 2.<M> 3.<pT> 4.A_UT 5.Sigma_Stat 6.Sigma_PID 7.Sigma_TriggerBias 8.Sigma_Total \n";
        ftable<<"// 9.A_UT 10.Sigma_Stat 11.Sigma_PID 12.Sigma_TriggerBias 13.Sigma_Total "<<endl;
        ftable<<"//(Columns 4-8 for eta>0 and columns 9-13 for eta<0)"<<endl;
        ftable<<"\n"<<endl;


	 TCanvas *myCanAvg = new TCanvas("myCanAvg","myCanAvg",150,10,1100,700);
  gStyle -> SetOptStat(0);
  myCanAvg->Divide(3,2,0,0);
  for(int pad=0;pad<6; pad++ ){
          myCanAvg->cd(pad+1);
          setPad(pad);
          if(pad==5){
                  TLatex text;
                  text.SetTextSize(0.08);
                  text.DrawLatex(0.1,0.7,"#font[22]{#color[2]{STAR Preliminary 2015}}");
                  text.SetTextSize(0.06);
                  text.DrawLatex(0.1,0.6,"#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 200 GeV}");
                  text.DrawLatex(0.1,0.6,"#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 200 GeV}");
                  text.DrawLatex(0.1,0.5,"#font[22]{#pm 3% scale uncertainty from beam}");
                  text.DrawLatex(0.1,0.4,"#font[22]{polarization (not shown)}");
          }else{
                  gPad->SetGrid(0,0);
                  if(pad==0)for(int mbin=0; mbin<9; mbin++){pTAvg[mbin]=pT1[mbin]; A_Mgt[mbin]=WavgA_pT1g[mbin]; Aerr_Mgt[mbin]=WerrA_pT1g[mbin]; A_Mlt[mbin]=WavgA_pT1l[mbin]; Aerr_Mlt[mbin]=WerrA_pT1l[mbin];finalSysgt[mbin]=combSysMgt1[mbin]; finalSyslt[mbin]=combSysMlt1[mbin];pidSysGt[mbin]=pidMgt1[mbin]; trigSysGt[mbin]=biasMgt1[mbin]; pidSysLt[mbin]=pidMlt1[mbin]; trigSysLt[mbin]=biasMlt1[mbin];}
                  if(pad==1)for(int mbin=0; mbin<9; mbin++){pTAvg[mbin]=pT2[mbin]; A_Mgt[mbin]=WavgA_pT2g[mbin]; Aerr_Mgt[mbin]=WerrA_pT2g[mbin]; A_Mlt[mbin]=WavgA_pT2l[mbin]; Aerr_Mlt[mbin]=WerrA_pT2l[mbin];finalSysgt[mbin]=combSysMgt2[mbin]; finalSyslt[mbin]=combSysMlt2[mbin];pidSysGt[mbin]=pidMgt2[mbin]; trigSysGt[mbin]=biasMgt2[mbin]; pidSysLt[mbin]=pidMlt2[mbin]; trigSysLt[mbin]=biasMlt2[mbin];}
                  if(pad==2)for(int mbin=0; mbin<9; mbin++){pTAvg[mbin]=pT3[mbin]; A_Mgt[mbin]=WavgA_pT3g[mbin]; Aerr_Mgt[mbin]=WerrA_pT3g[mbin]; A_Mlt[mbin]=WavgA_pT3l[mbin]; Aerr_Mlt[mbin]=WerrA_pT3l[mbin];finalSysgt[mbin]=combSysMgt3[mbin]; finalSyslt[mbin]=combSysMlt3[mbin];pidSysGt[mbin]=pidMgt3[mbin]; trigSysGt[mbin]=biasMgt3[mbin]; pidSysLt[mbin]=pidMlt3[mbin]; trigSysLt[mbin]=biasMlt3[mbin];}
                  if(pad==3)for(int mbin=0; mbin<9; mbin++){pTAvg[mbin]=pT4[mbin]; A_Mgt[mbin]=WavgA_pT4g[mbin]; Aerr_Mgt[mbin]=WerrA_pT4g[mbin]; A_Mlt[mbin]=WavgA_pT4l[mbin]; Aerr_Mlt[mbin]=WerrA_pT4l[mbin];finalSysgt[mbin]=combSysMgt4[mbin]; finalSyslt[mbin]=combSysMlt4[mbin];pidSysGt[mbin]=pidMgt4[mbin]; trigSysGt[mbin]=biasMgt4[mbin]; pidSysLt[mbin]=pidMlt4[mbin]; trigSysLt[mbin]=biasMlt4[mbin];}
                  if(pad==4)for(int mbin=0; mbin<9; mbin++){pTAvg[mbin]=pT5[mbin]; A_Mgt[mbin]=WavgA_pT5g[mbin]; Aerr_Mgt[mbin]=WerrA_pT5g[mbin]; A_Mlt[mbin]=WavgA_pT5l[mbin]; Aerr_Mlt[mbin]=WerrA_pT5l[mbin];finalSysgt[mbin]=combSysMgt5[mbin]; finalSyslt[mbin]=combSysMlt5[mbin];pidSysGt[mbin]=pidMgt5[mbin]; trigSysGt[mbin]=biasMgt5[mbin]; pidSysLt[mbin]=pidMlt5[mbin]; trigSysLt[mbin]=biasMlt5[mbin];}

		  ftable<<M[pad]<<" - "<<M[pad+1]<<" & ";
                  for(int i=0; i<9; i++){
                  if(i==0){ftable<<avgMinvGt[pad]<<" & "<<pTAvg[i]<<" & "<< A_Mgt[i]<<" & "<<Aerr_Mgt[i]<< " & "<<pidSysGt[i]<<" & "<<trigSysGt[i]<<" & " <<finalSysgt[i]<<" & "<<A_Mlt[i]<<" & "<<Aerr_Mlt[i]<< " & "<<pidSysLt[i]<<" & "<<trigSysLt[i]<<" & "<< finalSyslt[i]<<" \\\\ "<<endl;
                  }else{
                  ftable<<"\t & "<<avgMinvGt[pad]<<" & "<<pTAvg[i]<<" & "<< A_Mgt[i]<<" & "<<Aerr_Mgt[i]<< " & "<<pidSysGt[i]<<" & "<<trigSysGt[i]<<" & " <<finalSysgt[i]<<" & "<<A_Mlt[i]<<" & "<<Aerr_Mlt[i]<< " & "<<pidSysLt[i]<<" & "<<trigSysLt[i]<<" & "<< finalSyslt[i]<<" \\\\ "<<endl;
		  }
		 }

                  gr_AvsMgts[pad] = new TGraphErrors(9,pTAvg,A_Mgt,errxts,trigSysGt);
                  gr_AvsMgts[pad]->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}} ");
                  gr_AvsMgts[pad]->SetTitle("");
                  gr_AvsMgts[pad]->GetYaxis()->SetTitleOffset(1.2);
                  gr_AvsMgts[pad]->GetYaxis()->SetTitleSize(0.06);
                  gr_AvsMgts[pad]->GetYaxis()->SetLabelSize(0.06);
                  gr_AvsMgts[pad]->GetYaxis()->SetLabelFont(22);
                  gr_AvsMgts[pad]->GetXaxis()->SetTitle("#font[22]{p_{T}^{#pi^{+}#pi^{-}} (GeV/c)}   ");
                  gr_AvsMgts[pad]->GetXaxis()->SetTitleSize(0.06);
                  gr_AvsMgts[pad]->GetXaxis()->SetLabelSize(0.06);
                  gr_AvsMgts[pad]->GetXaxis()-> SetLabelFont(22);
		  gr_AvsMgts[pad]->GetXaxis()->SetLimits(2.5,11.5);
		  gr_AvsMgts[pad]->GetYaxis()->SetRangeUser(-0.015, 0.075);
                  gr_AvsMgts[pad]->SetFillStyle(0000);
                  gr_AvsMgts[pad]->SetLineColor(1);
                  gr_AvsMgts[pad]->SetLineWidth(1);
                  gr_AvsMgts[pad]->GetXaxis()->SetNdivisions(505);
                  gr_AvsMgts[pad]->GetYaxis()->SetNdivisions(505);
                  gr_AvsMgts[pad]->Draw("AP2");
                  if(pad==2){
                          gr_AvsMgts[pad]->GetYaxis()->SetTitle("");
                          gr_AvsMgts[pad]->GetYaxis()->SetLabelSize(0);
                  }
                  if(pad==1){
                          gr_AvsMgts[pad]->GetYaxis()->SetTitle("");
                          gr_AvsMgts[pad]->GetXaxis()->SetTitle("");
                          gr_AvsMgts[pad]->GetXaxis()->SetLabelSize(0);

                  }
                  if(pad==0){
                          gr_AvsMgts[pad]->GetXaxis()->SetTitle("");
                          gr_AvsMgts[pad]->GetXaxis()->SetLabelSize(0.0);

                  }

		  gr_AvsMgt[pad] = new TGraphErrors(9,pTAvg,A_Mgt,errx,Aerr_Mgt);//Backward average 
                  gr_AvsMgt[pad]-> SetMarkerStyle(20);
                  gr_AvsMgt[pad]-> SetMarkerColor(2);
                  gr_AvsMgt[pad]-> SetLineColor(2);
                  gr_AvsMgt[pad]-> SetLineWidth(1);
                  gr_AvsMgt[pad]-> Draw("same P");

                  gr_AvsMlts[pad] = new TGraphErrors(9,pTAvg,A_Mlt,errxts,trigSysLt);//Backward average 
                  gr_AvsMlts[pad]-> SetFillStyle(0000);
                  gr_AvsMlts[pad]-> SetLineColor(1);
                  gr_AvsMlts[pad]-> SetLineWidth(1);
                  gr_AvsMlts[pad]-> Draw("same 2");
                  
                  gr_AvsMlt[pad] = new TGraphErrors(9,pTAvg,A_Mlt,errx,Aerr_Mlt);//Backward average 
                  gr_AvsMlt[pad]-> SetMarkerStyle(20);
                  gr_AvsMlt[pad]-> SetMarkerColor(4);
                  gr_AvsMlt[pad]-> SetLineColor(4);
                  gr_AvsMlt[pad]-> SetLineWidth(1);
                  gr_AvsMlt[pad]-> Draw("same P");

                  gPad->Update();

                  line[pad]=  new TLine(myCanAvg->cd(pad+1)->GetUxmin(),0.,myCanAvg->cd(pad+1)->GetUxmax(),0.);
                  line[pad]->SetLineStyle(2);
                  line[pad]->SetLineWidth(1);
                  line[pad]->Draw();

                  tex[pad].SetTextSize(0.055);
                  tex[pad].SetTextAlign(13);
                  tex[pad].DrawLatex(3,0.067,Form("#font[22]{#color[1]{< M_{inv} > = %3.2lf GeV/c^{2}}}",avgMinvGt[pad]));//in %3.2lf, 3 sets the number of digits and 2 sets the numbers after decimal

                  if(pad==0){
                          leg = new TLegend(0.7,.61, 0.97, 0.81);
                          leg->AddEntry(gr_AvsMgt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} > 0}", "lp");
                          leg->AddEntry(gr_AvsMlt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} < 0}", "lp");
                          leg->AddEntry(gr_AvsMgts[pad], " #font[22]{ Syst. Error }", "f");
                          leg->SetTextSize(0.05);
                          leg->Draw();
                  }
                  gPad->Update();
          }

  }
myCanAvg->SaveAs("AsymVsPtTOFcutTPConly.pdf");


	//---------------------------------//

	//Forward Asymmetry BLUE and YELLOW
		TCanvas *myCan = new TCanvas("myCan","myCan",150,10,990,660);
		gStyle -> SetOptStat(0);
		//gStyle -> SetTitleX(.35);
		gStyle->SetLegendBorderSize(0);
		myCan->SetGrid(0,0);
		myCan -> Divide(3,2);
		myCan -> cd(1);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		//gPad->SetBottomMargin(0.005);
		//gPad->SetTopMargin(0.1);
		pt1 = new TGraphErrors(9,pT1,A_pT1bg,0,deltaA_pT1bg);
		pt1->GetYaxis()-> SetTitle("A_{UT}");
		pt1->SetTitle("");
		//pt1->GetYaxis()-> CenterTitle();
		//pt1->GetYaxis()-> SetTitleSize(.06);
	  	pt1->GetYaxis()->SetTitleOffset(1.66);
	  	//pt1->GetYaxis()->SetLabelSize(0.06);
	  	//pt1->GetXaxis()->SetLabelSize(0);
		pt1->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		//pt1-> GetXaxis()-> CenterTitle();
		//pt1-> SetTitle(" 2.85 < p_{T} < 3.65)");
		pt1-> SetMarkerStyle(20);
		pt1-> SetMarkerColor(4);
		pt1-> GetXaxis()->SetLimits(2,11);
		pt1-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt1->Draw("AP");
		myCan->Update();
		pt1s= new TGraphErrors(9,pT1,A_pT1yg,0,deltaA_pT1yg);
		pt1s->SetMarkerStyle(20);
		pt1s->SetMarkerColor(2);
		pt1s->Draw("same P");
		TLegend *leg1 = new TLegend(0.2,.7, 0.4, 0.8);
		leg1->AddEntry(pt1, "#eta > 0 , BLUE", "lp");
		leg1->AddEntry(pt1s, "#eta > 0, YELLOW", "lp");
		leg1->SetTextSize(0.04);
		leg1->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[0]));


		myCan->Update();
		//line1 = new TLine(0.3,0.,1.45,0);
 		TLine *line1=  new TLine(myCan->cd(1)->GetUxmin(),0.,myCan->cd(1)->GetUxmax(),0.);
		line1->SetLineStyle(2);
		line1->Draw();
		
		myCan -> cd(2);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		//gPad->SetBottomMargin(0.005);
		//gPad->SetTopMargin(0.1);
		pt2 = new TGraphErrors(9,pT2,A_pT2bg,0,deltaA_pT2bg);
		pt2->GetYaxis()-> SetTitle("A_{UT}");
		//pt2->GetYaxis()-> CenterTitle();
	  	pt2->GetYaxis()->SetTitleOffset(1.5);
		pt2->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		//pt2->GetXaxis()-> CenterTitle();
	  	//pt2->GetXaxis()->SetLabelSize(0);
	  	//pt2->GetYaxis()->SetLabelSize(0);
		//pt2-> SetTitle("3.65 < p_{T} < 4.20");
		pt2->SetTitle("");
		pt2-> SetMarkerStyle(20);
		pt2-> SetMarkerColor(4);
		pt2-> GetXaxis()->SetLimits(2,11);
		pt2-> GetYaxis()->SetRangeUser(-0.015, .065);
		pt2->Draw("AP"); 
		myCan->Update();
		pt2s= new TGraphErrors(9,pT2,A_pT2yg,0,deltaA_pT2yg);
		pt2s->SetMarkerStyle(20);
		pt2s->SetMarkerColor(2);
		pt2s->Draw("same P");
		TLegend *leg2 = new TLegend(0.2,0.7, 0.4, 0.8);
		leg2->AddEntry(pt2, "#eta > 0, BLUE", "lp");
		leg2->AddEntry(pt2s, "#eta > 0, YELLOW", "lp");
		leg2->SetTextSize(0.04);
		leg2->Draw(); 
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[1]));
		myCan -> Update();
		  TLine *line2=  new TLine(myCan->cd(2)->GetUxmin(),0.,myCan->cd(2)->GetUxmax(),0.);
 		line2->SetLineStyle(2);
		line2->Draw();

		myCan -> cd (3);
		gPad->SetGrid(0,0); 
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		//gPad->SetRightMargin(0.1);
		//gPad->SetBottomMargin(0.05);
		//gPad->SetTopMargin(0.1);
		pt3 = new TGraphErrors(9,pT3,A_pT3bg,0,deltaA_pT3bg);
		pt3->GetYaxis()-> SetTitle("A_{UT}");
		//pt3->GetYaxis()->CenterTitle();
	  	pt3->GetYaxis()->SetTitleOffset(1.5);
	  	//pt3->GetXaxis()->SetLabelSize(0);
	  	//pt3->GetYaxis()->SetLabelSize(0);
		pt3->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt3->SetTitle("");
		//pt3->GetXaxis()->CenterTitle();
		//pt3-> SetTitle(" 4.20 < p_{T} < 4.95");
		pt3-> SetMarkerStyle(20);
		pt3-> SetMarkerColor(4);
		pt3-> GetXaxis()->SetLimits(2,11);
		pt3-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt3->Draw("AP");
		myCan->Update();
		pt3s = new TGraphErrors(9, pT3, A_pT3yg, 0, deltaA_pT3yg);
		pt3s->SetMarkerStyle(20);
		pt3s->SetMarkerColor(2);
		pt3s->Draw("same P");
		TLegend *leg3 = new TLegend(0.2,0.7, 0.4, 0.8);
		leg3->AddEntry(pt3, "#eta > 0, BLUE", "lp");
		leg3->AddEntry(pt3s, "#eta > 0, YELLOW", "lp");
		leg3->SetTextSize(0.04);
		leg3->Draw(); 
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[2]));
		myCan -> Update();
		TLine *line3=  new TLine(myCan->cd(3)->GetUxmin(),0.,myCan->cd(3)->GetUxmax(),0.);
 		line3->SetLineStyle(2);
		line3->Draw();


		myCan -> cd (4); 
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		//gPad->SetBottomMargin(0.15);
		//gPad->SetTopMargin(0.005);
		pt4= new TGraphErrors(9,pT4,A_pT4bg,0,deltaA_pT4bg);
		pt4->GetYaxis()-> SetTitle("A_{UT}");
		//pt4->GetYaxis()-> CenterTitle();
	  	pt4->GetYaxis()->SetTitleOffset(1.5);
		pt4->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		//pt4->GetXaxis()-> CenterTitle();
		//pt4->GetXaxis()-> SetLabelSize(0);
		pt4->SetTitle("");
		//pt4->SetTitle("4.95 < p_{T} < 6.00");
		pt4-> SetMarkerStyle(20);
		pt4-> SetMarkerColor(4);
		pt4-> GetXaxis()->SetLimits(2,11);
		pt4-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		//pt4->GetYaxis()-> SetTitleSize(.05);
		//pt4->GetXaxis()-> SetTitleSize(.05);
                pt4->GetYaxis()->SetTitleOffset(1.6);
                //pt4->GetYaxis()->SetLabelSize(0.05);
                //pt4->GetXaxis()->SetLabelSize(0.05);

		pt4->Draw("AP");
		myCan->Update();
		pt4s = new TGraphErrors(9, pT4, A_pT4yg, 0, deltaA_pT4yg);
		pt4s->SetMarkerStyle(20);
		pt4s->SetMarkerColor(2);
		pt4s->Draw("same P");
		TLegend *leg4 = new TLegend(0.2,0.7, 0.4, 0.8);
		leg4->AddEntry(pt4, "#eta > 0, BLUE", "lp");
		leg4->AddEntry(pt4s, "#eta > 0, YELLOW", "lp");
		leg4->SetTextSize(0.04);
		leg4->Draw(); 
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[3]));
		myCan -> Update();
		TLine *line4=  new TLine(myCan->cd(4)->GetUxmin(),0.,myCan->cd(4)->GetUxmax(),0.);
 		line4->SetLineStyle(2);
		line4->Draw();
	
		myCan -> cd (5); 
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		//gPad->SetBottomMargin(0.15);
		//gPad->SetTopMargin(0.01);
		pt5 = new TGraphErrors(9,pT5,A_pT5bg,0,deltaA_pT5bg);
		pt5->GetYaxis()-> SetTitle("A_{UT}");
		//pt5->GetYaxis()-> CenterTitle();
		pt5->SetTitle("");
	  	pt5->GetYaxis()->SetTitleOffset(1.5);
		pt5->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		//pt5->GetXaxis()-> CenterTitle();
	  	//pt5->GetXaxis()->SetLabelSize(0.05);
	  	//pt5->GetXaxis()->SetTitleSize(0.05);
		//pt5-> SetTitle(" 6.00< p_{T} <20)");
		//pt5->GetYaxis()-> SetTitleSize(.06);
                //pt5->GetYaxis()->SetTitleOffset(1.5);
                //pt5->GetYaxis()->SetLabelSize(0.05);
		pt5-> SetMarkerStyle(20);
		pt5-> SetMarkerColor(4);
		pt5-> GetXaxis()->SetLimits(2,11);
		pt5-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt5->Draw("AP");
		myCan->Update();
		pt5s = new TGraphErrors(9, pT5, A_pT5yg, 0, deltaA_pT5yg);
		pt5s->SetMarkerStyle(20);
		pt5s->SetMarkerColor(2);
		pt5s->Draw("same P");
		TLegend *leg5 = new TLegend(0.2,0.7, 0.4, 0.8);
		leg5->AddEntry(pt5, "#eta > 0, BLUE", "lp");
		leg5->AddEntry(pt5s, "#eta > 0, YELLOW", "lp");
		leg5->SetTextSize(0.04);
		leg5->Draw(); 
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[4]));
		myCan -> Update();
		TLine *line5=  new TLine(myCan->cd(5)->GetUxmin(),0.,myCan->cd(5)->GetUxmax(),0.);
 		line5->SetLineStyle(2);
		line5->Draw();
	
	// Average Asymmetry 
	
		 TCanvas *myCanA = new TCanvas("myCanA","myCanA",150,10,990,660);
		gStyle -> SetOptStat(0);
		gStyle->SetLegendBorderSize(0);
		myCanA->SetGrid(0,0);
		myCanA -> Divide(3,2);
		myCanA -> cd(1);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt1AF = new TGraphErrors(9,pT1,avgA_pT1g,0,errA_pT1g);//Forward average
		pt1AF->GetYaxis()-> SetTitle("A_{UT}");
		pt1AF->SetTitle("");
	  	pt1AF->GetYaxis()->SetTitleOffset(1.66);
		pt1AF->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt1AF-> SetMarkerStyle(20);
		pt1AF-> SetMarkerColor(2);
		pt1AF-> GetXaxis()->SetLimits(2,11);
		pt1AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt1AF->Draw("AP");
		pt1AB = new TGraphErrors(9,pT1,avgA_pT1l,0,errA_pT1l);//Backward average 
		pt1AB-> SetMarkerStyle(20);
		pt1AB-> SetMarkerColor(4);
		pt1AB->Draw("same P");
		myCanA->Update();
		TLegend *leg1A = new TLegend(0.2,0.7, 0.4, 0.8);
		leg1A->AddEntry(pt1AF, "<A_{UT}>, #eta > 0", "lp");
		leg1A->AddEntry(pt1AB, "<A_{UT}>, #eta < 0", "lp");
		leg1A->SetTextSize(0.04);
		leg1A->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[0]));


		myCanA->Update();
 		TLine *line1A=  new TLine(myCanA->cd(1)->GetUxmin(),0.,myCanA->cd(1)->GetUxmax(),0.);
		line1A->SetLineStyle(2);
		line1A->Draw();
	
		myCanA -> cd(2);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt2AF = new TGraphErrors(9,pT2,avgA_pT2g,0,errA_pT2g);//Forward average
		pt2AF->GetYaxis()-> SetTitle("A_{UT}");
		pt2AF->SetTitle("");
	  	pt2AF->GetYaxis()->SetTitleOffset(1.66);
		pt2AF->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt2AF-> SetMarkerStyle(20);
		pt2AF-> SetMarkerColor(2);
		pt2AF-> GetXaxis()->SetLimits(2,11);
		pt2AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt2AF->Draw("AP");
		pt2AB = new TGraphErrors(9,pT2,avgA_pT2l,0,errA_pT2l);//Backward average
		pt2AB-> SetMarkerStyle(20);
		pt2AB-> SetMarkerColor(4);
		pt2AB->Draw("same P");
		gPad->Update();
		TLegend *leg2A = new TLegend(0.2,0.7, 0.4, 0.8);
		leg2A->AddEntry(pt2AF, "<A_{UT}>, #eta > 0", "lp");
		leg2A->AddEntry(pt2AB, "<A_{UT}>, #eta < 0", "lp");
		leg2A->SetTextSize(0.04);
		leg2A->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[1]));
		gPad->Update();
 		TLine *line2A=  new TLine(myCanA->cd(2)->GetUxmin(),0.,myCanA->cd(2)->GetUxmax(),0.);
		line2A->SetLineStyle(2);
		line2A->Draw();
		gPad->Update();
	
		myCanA -> cd(3);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt3AF = new TGraphErrors(9,pT3,avgA_pT3g,0,errA_pT3g);
		pt3AF->GetYaxis()-> SetTitle("A_{UT}");
		pt3AF->SetTitle("");
	  	pt3AF->GetYaxis()->SetTitleOffset(1.66);
		pt3AF->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c");
		pt3AF-> SetMarkerStyle(20);
		pt3AF-> SetMarkerColor(2);
		pt3AF-> GetXaxis()->SetLimits(2,11);
		pt3AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt3AF->Draw("AP");
		pt3AB = new TGraphErrors(9,pT3,avgA_pT3l,0,errA_pT3l);
		pt3AB-> SetMarkerStyle(20);
		pt3AB-> SetMarkerColor(4);
		pt3AB->Draw("same P");
		gPad->Update();
		TLegend *leg3A = new TLegend(0.2,0.7, 0.4, 0.8);
		leg3A->AddEntry(pt3AF, "<A_{UT}>, #eta > 0", "lp");
		leg3A->AddEntry(pt3AB, "<A_{UT}>, #eta < 0", "lp");
		leg3A->SetTextSize(0.04);
		leg3A->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[2]));
		gPad->Update();
 		TLine *line3A=  new TLine(myCanA->cd(3)->GetUxmin(),0.,myCanA->cd(3)->GetUxmax(),0.);
		line3A->SetLineStyle(2);
		line3A->Draw();
		gPad->Update();

		myCanA -> cd(4);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt4AF = new TGraphErrors(9,pT4,avgA_pT4g,0,errA_pT4g);
		pt4AF->GetYaxis()-> SetTitle("A_{UT}");
		pt4AF->SetTitle("");
	  	pt4AF->GetYaxis()->SetTitleOffset(1.66);
		pt4AF->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt4AF-> SetMarkerStyle(20);
		pt4AF-> SetMarkerColor(2);
		pt4AF-> GetXaxis()->SetLimits(2,11);
		pt4AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt4AF->Draw("AP");
		pt4AB = new TGraphErrors(9,pT4,avgA_pT4l,0,errA_pT4l);
		pt4AB-> SetMarkerStyle(20);
		pt4AB-> SetMarkerColor(4);
		pt4AB->Draw("same P");
		gPad->Update();
		TLegend *leg4A = new TLegend(0.2,0.7, 0.4, 0.8);
		leg4A->AddEntry(pt4AF, "<A_{UT}>, #eta > 0", "lp");
		leg4A->AddEntry(pt4AB, "<A_{UT}>, #eta < 0", "lp");
		leg4A->SetTextSize(0.04);
		leg4A->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[3]));
		gPad->Update();

 		TLine *line4A=  new TLine(myCanA->cd(4)->GetUxmin(),0.,myCanA->cd(4)->GetUxmax(),0.);
		line4A->SetLineStyle(2);
		line4A->Draw();
		gPad->Update();

		myCanA -> cd(5);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt5AF = new TGraphErrors(9,pT5,avgA_pT5g,0,errA_pT5g);
		pt5AF->GetYaxis()-> SetTitle("A_{UT}");
		pt5AF->SetTitle("");
	  	pt5AF->GetYaxis()->SetTitleOffset(1.66);
		pt5AF->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt5AF-> SetMarkerStyle(20);
		pt5AF-> SetMarkerColor(2);
		pt5AF-> GetXaxis()->SetLimits(2,11);
		pt5AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt5AF->Draw("AP");
		pt5AB= new TGraphErrors(9,pT5,avgA_pT5l,0,errA_pT5l);
		pt5AB-> SetMarkerStyle(20);
		pt5AB-> SetMarkerColor(4);
		pt5AB-> Draw("same P");
		gPad->Update();
		TLegend *leg5A = new TLegend(0.2,0.7, 0.4, 0.8);
		leg5A->AddEntry(pt5AF, "<A_{UT}>, #eta > 0", "lp");
		leg5A->AddEntry(pt5AB, "<A_{UT}>, #eta < 0", "lp");
		leg5A->SetTextSize(0.04);
		leg5A->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[4]));
		myCanA->Update();

 		TLine *line5A=  new TLine(myCanA->cd(5)->GetUxmin(),0.,myCanA->cd(5)->GetUxmax(),0.);
		line5A->SetLineStyle(2);
		line5A->Draw();
		gPad->Update();

	myCan->SaveAs("AsymVspT_Cone.7_ForwardTPConly.pdf");
	myCanA->SaveAs("AsymVspT_Cone.7_AverageTPConly.pdf");

	//Backward Asymmetry BLUE and YELLOW
	
		 TCanvas *cB = new TCanvas("cB","cB",150,10,990,660);
		gStyle -> SetOptStat(0);
		gStyle->SetLegendBorderSize(0);
		cB->SetGrid(0,0);
		cB -> Divide(3,2);
		cB -> cd(1);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt1Bl = new TGraphErrors(9,pT1,A_pT1bl,0,deltaA_pT1bl);
		pt1Bl->GetYaxis()-> SetTitle("A_{UT}");
		pt1Bl->SetTitle("");
	  	pt1Bl->GetYaxis()->SetTitleOffset(1.66);
		pt1Bl->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt1Bl-> SetMarkerStyle(20);
		pt1Bl-> SetMarkerColor(4);
		pt1Bl-> GetXaxis()->SetLimits(2,11);
		pt1Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt1Bl->Draw("AP");
		pt1Yl = new TGraphErrors(9,pT1,A_pT1yl,0,deltaA_pT1yl);
		pt1Yl-> SetMarkerStyle(20);
		pt1Yl-> SetMarkerColor(2);
		pt1Yl-> Draw("same P");
		gPad->Update();
		TLegend *lBl = new TLegend(0.2,0.7, 0.4, 0.8);
		lBl->AddEntry(pt1Bl, "BLUE, #eta < 0", "lp");
		lBl->AddEntry(pt1Yl, "YELLOW, #eta < 0", "lp");
		lBl->SetTextSize(0.04);
		lBl->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[0]));
 		TLine *line1A=  new TLine(cB->cd(1)->GetUxmin(),0.,cB->cd(1)->GetUxmax(),0.);
		line1A->SetLineStyle(2);
		line1A->Draw();
		gPad->Update();
	
		cB -> cd(2);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt2Bl = new TGraphErrors(9,pT2,A_pT2bl,0,deltaA_pT2bl);
		pt2Bl->GetYaxis()-> SetTitle("A_{UT}");
		pt2Bl->SetTitle("");
	  	pt2Bl->GetYaxis()->SetTitleOffset(1.66);
		pt2Bl->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt2Bl-> SetMarkerStyle(20);
		pt2Bl-> SetMarkerColor(4);
		pt2Bl-> GetXaxis()->SetLimits(2,11);
		pt2Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt2Bl->Draw("AP");
		pt2Yl = new TGraphErrors(9,pT2, A_pT2yl,0,deltaA_pT2yl);
		pt2Yl-> SetMarkerStyle(20);
		pt2Yl-> SetMarkerColor(2);
		pt2Yl-> Draw("same P");
		gPad->Update();
		TLegend *lBl = new TLegend(0.2,0.7, 0.4, 0.8);
		lBl->AddEntry(pt2Bl, "BLUE, #eta < 0", "lp");
		lBl->AddEntry(pt2Yl, "YELLOW, #eta < 0", "lp");
		lBl->SetTextSize(0.04);
		lBl->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[1]));
 		TLine *line1A=  new TLine(cB->cd(2)->GetUxmin(),0.,cB->cd(2)->GetUxmax(),0.);
		line1A->SetLineStyle(2);
		line1A->Draw();
		gPad->Update();
		
		cB -> cd(3);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt3Bl = new TGraphErrors(9,pT3,A_pT3bl,0,deltaA_pT3bl);
		pt3Bl->GetYaxis()-> SetTitle("A_{UT}");
		pt3Bl->SetTitle("");
	  	pt3Bl->GetYaxis()->SetTitleOffset(1.66);
		pt3Bl->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt3Bl-> SetMarkerStyle(20);
		pt3Bl-> SetMarkerColor(4);
		pt3Bl-> GetXaxis()->SetLimits(2,11);
		pt3Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt3Bl->Draw("AP");
		pt3Yl = new TGraphErrors(9,pT3,A_pT3yl,0,deltaA_pT3yl);
		pt3Yl-> SetMarkerStyle(20);
		pt3Yl-> SetMarkerColor(2);
		pt3Yl-> Draw("same P");
		gPad->Update();
		TLegend *lBl = new TLegend(0.2,0.7, 0.4, 0.8);
		lBl->AddEntry(pt3Bl, "BLUE, #eta < 0", "lp");
		lBl->AddEntry(pt3Yl, "YELLOW, #eta < 0", "lp");
		lBl->SetTextSize(0.04);
		lBl->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[2]));
 		TLine *line1A=  new TLine(cB->cd(3)->GetUxmin(),0.,cB->cd(3)->GetUxmax(),0.);
		line1A->SetLineStyle(2);
		line1A->Draw();
		gPad->Update();

		cB -> cd(4);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt4Bl = new TGraphErrors(9,pT4,A_pT4bl,0,deltaA_pT4bl);
		pt4Bl->GetYaxis()-> SetTitle("A_{UT}");
		pt4Bl->SetTitle("");
	  	pt4Bl->GetYaxis()->SetTitleOffset(1.66);
		pt4Bl->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt4Bl-> SetMarkerStyle(20);
		pt4Bl-> SetMarkerColor(4);
		pt4Bl-> GetXaxis()->SetLimits(2,11);
		pt4Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt4Bl->Draw("AP");
		pt4Yl = new TGraphErrors(9,pT4,A_pT4yl,0,deltaA_pT4yl);
		pt4Yl-> SetMarkerStyle(20);
		pt4Yl-> SetMarkerColor(2);
		pt4Yl-> Draw("same P");
		gPad->Update();
		TLegend *lBl = new TLegend(0.2,0.7, 0.4, 0.8);
		lBl->AddEntry(pt4Bl, "BLUE, #eta < 0", "lp");
		lBl->AddEntry(pt4Yl, "YELLOW, #eta < 0", "lp");
		lBl->SetTextSize(0.04);
		lBl->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c^{2}}}",avg_Minv[3]));
 		TLine *line1A=  new TLine(cB->cd(4)->GetUxmin(),0.,cB->cd(4)->GetUxmax(),0.);
		line1A->SetLineStyle(2);
		line1A->Draw();
		gPad->Update();

		cB -> cd(5);
		gPad->SetGrid(0,0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.01);
		pt5Bl = new TGraphErrors(9,pT5,A_pT5bl,0,deltaA_pT5bl);
		pt5Bl->GetYaxis()-> SetTitle("A_{UT}");
		pt5Bl->SetTitle("");
	  	pt5Bl->GetYaxis()->SetTitleOffset(1.66);
		pt5Bl->GetXaxis()->SetTitle("p_{T}^{#pi^{+}#pi^{-}}(GeV/c)");
		pt5Bl-> SetMarkerStyle(20);
		pt5Bl-> SetMarkerColor(4);
		pt5Bl-> GetXaxis()->SetLimits(2,11);
		pt5Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
		pt5Bl->Draw("AP");
		pt5Yl = new TGraphErrors(9,pT5,A_pT5yl,0,deltaA_pT5yl);
		pt5Yl-> SetMarkerStyle(20);
		pt5Yl-> SetMarkerColor(2);
		pt5Yl-> Draw("same P");
		gPad->Update();
		TLegend *lBl = new TLegend(0.2,0.7, 0.4, 0.8);
		lBl->AddEntry(pt5Bl, "BLUE, #eta < 0", "lp");
		lBl->AddEntry(pt5Yl, "YELLOW, #eta < 0", "lp");
		lBl->SetTextSize(0.04);
		lBl->Draw();
	
		TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(2.2,0.062,Form("#font[12]{#color[2]{< M_{inv} > = %g GeV/c}}",avg_Minv[4]));
 		TLine *line1A=  new TLine(cB->cd(5)->GetUxmin(),0.,cB->cd(5)->GetUxmax(),0.);
		line1A->SetLineStyle(2);
		line1A->Draw();
		gPad->Update();
	cB->Update();
	cB->SaveAs("AsymmetryVspT_Cone.7_BackwardTPConly.pdf");
	
chi2f->Write();
chi2f->Close();


}//main

//set pad for money plot

void setPad(int i){

   if(i==2){
   gPad->SetTopMargin(0.137);
   gPad->SetLeftMargin(0.02);
   gPad->SetBottomMargin(0.15);
   gPad->SetRightMargin(0.02);
   gPad->SetPad(.6634,0.4197,0.97,.9897);
   gPad->SetFillStyle(4000);
   gPad->Update();
   }else if(i==5){
   gPad->SetTopMargin(0.4);
   gPad->SetLeftMargin(0.01);
   gPad->SetBottomMargin(0.15);
   gPad->SetRightMargin(0.02);
   gPad->SetPad(.6634,0.01,0.97,.49);
   gPad->SetFillStyle(4000);
   gPad->SetFrameFillColor(0);
   gPad->SetFrameFillStyle(0);
   gPad->SetFrameBorderMode(0);

   gPad->Update();
   }else if(i==1){
  gPad->SetTopMargin(0.16);
  gPad->SetLeftMargin(0.0);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(.37,0.5,0.67,.99);
  gPad->Update();
  }else if(i==0){
  gPad->SetTopMargin(0.16);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(0.,0.5,0.37,.99);
  }else if(i==3){
  gPad->SetTopMargin(0.0);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(0.0,0.02052,0.37,.5052);
  }else if(i==4){
  gPad->SetTopMargin(0.0);
  gPad->SetLeftMargin(0.0);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.0035);
  gPad->SetPad(.37,0.02052,0.67,.5052);
  }
}
                                  
