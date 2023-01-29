//This code is modefied to calculate asymmetry for individual beams and averaged for total asymmetry
#include<iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

void integratedAsym_pT_9bins(const char *ifile="/gpfs01/star/pwg/pokhrel/IFF_TREES/Trees_NewBatch/ntuples_final.root")//interactive 
//oid test_pT_Cone3(const char *ifile)//schedular 
{
	TH1D *hpolB = new TH1D("hpolB","", 100, 0, 1);	
	TH1D *hpolY = new TH1D("hpolY","", 100, 0, 1);	
	//chi2/ndf distribution for sine fit 
	TFile *f = new TFile(ifile);
	//TFile *f = new TFile("/gpfs01/star/pwg/pokhrel/IFF_TREES/AllTrig/ntuple/ntuple_allData.root");
	//TFile *f = new TFile("/star/data05/scratch/pokhrel/IFF_TREES/ntuple_test10.root");
	TTree *ntuple1JP = (TTree*)f->Get("ntuple1JP"); //get trees 
	TTree *ntuple2JP = (TTree*)f->Get("ntuple2JP");
	TTree *ntuple4JP = (TTree*)f->Get("ntuple4JP");
	TTree *ntuple5JP = (TTree*)f->Get("ntuple5JP");
	TTree *ntuple6 = (TTree*)f->Get("ntuple6");
	//define variables to hold trees variables 

	float eta_pair;
	float PhiRS;
	float fillNum;
	float runNum;
	float fspinconfig;
	float cone;
	float pT_pair;
	float Minv;
	float PhiRSB;
	float PhiRSY;
	float fitPts_min_pair;
	float ftrigger;
	float polB_corr, polY_corr;
	double pi = 3.14159265359;


	double Npairs[9];
	double pTpairs[9];
	double Mpairs[9];
	double NgtB[9];
	double MgtB[9];
	double pTgtB[9];
	double NltB[9];
	double MltB[9];
	double pTltB[9];
	double NgtY[9];
	double MgtY[9];
	double pTgtY[9];
	double NltY[9];
	double MltY[9];
	double pTltY[9];
	for(int k=0; k<9; k++) 
	{	
		Npairs[k]=0.;
		pTpairs[k]=0.;
		Mpairs[k]=0.;
		NgtB[k]=0;
		MgtB[k]=0;
		pTgtB[k]=0;
		NltB[k]=0;
		MltB[k]=0;
		pTltB[k]=0;
		NgtY[k]=0;
		MgtY[k]=0;
		pTgtY[k]=0;
		NltY[k]=0;
		MltY[k]=0;
		pTltY[k]=0;
	}

	//get the values of variables from tree 


	ntuple1JP->SetBranchAddress("fspinconfig",&fspinconfig);
	ntuple1JP->SetBranchAddress("cone",&cone);
	ntuple2JP->SetBranchAddress("Minv",&Minv);
	ntuple2JP->SetBranchAddress("pT_pair",&pT_pair);
	ntuple2JP->SetBranchAddress("eta_pair",&eta_pair);
	ntuple4JP->SetBranchAddress("PhiRSB",&PhiRSB);
	ntuple4JP->SetBranchAddress("PhiRSY",&PhiRSY);
	ntuple5JP->SetBranchAddress("fitPts_min_pair",&fitPts_min_pair);
	ntuple6->SetBranchAddress("polB_corr", &polB_corr);
	ntuple6->SetBranchAddress("polY_corr", &polY_corr);

	int nentries = (int)ntuple1JP->GetEntries(); //entries of trees 
	int ne = (int)ntuple6->GetEntries(); //entries of trees 
		for(int p=0; p<ne; p++)
		//for(int p=0; p<10000; p++)
		{
			ntuple6->GetEntry(p);
			hpolB->Fill(polB_corr);			
			hpolY->Fill(polY_corr);			
		}

		
	 cout<<"Entries  = "<<nentries<<endl;
	//add friend to the tree. no root file to add friend means the tree is on the same root file 
	ntuple1JP->AddFriend("ntuple2JP"); 
	ntuple1JP->AddFriend("ntuple4JP");
	ntuple1JP->AddFriend("ntuple5JP");
	//ntuple1->AddFriend("ntuple6");
	Double_t pT[10]={2.8, 3.459, 3.7717, 4.10, 4.461, 4.8923, 5.4431, 6.22, 7.566, 15};//binned dodging K0 mass range (0.4876 – 0.5076) 
	//Double_t pT[6]={2.70,3.97,4.71,5.69,7.2,20};//binned for 0.05< cone < 0.3 dodging K0 mass range (0.4876 – 0.5076) 
	//Double_t M[6]={.28, .54, .75, 1.0, 1.4, 4 };//binning with dodiging K0 mass range (0.4876 – 0.5076)      
	
	double eta_range[9];

	//To store average polarization values from histograms
	double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;
		
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


		//cout << "Entries: " << endl;
		for (int i=0; i < nentries; i++ )
		//for (int i=0; i < 10000; i++ )
		{

			ntuple1JP -> GetEntry(i);

			if(cone>0.7)continue;
			if(Minv>4.)continue;
			if(fitPts_min_pair<15)continue;
			
			 if(Minv>=0.4876 && Minv<0.5076 ) continue;//dodge K0 mass range, doesn't cause asymmetry.
			//if(Minv>0.4876 || Minv<0.5076)  cout <<"Something is wrong!! "<< "Minv: "<< Minv << endl;
			//cout << Minv << endl;

			//cout << "selection cuts are good ...working on phi loop.... " << endl; 
			//-------------BLUE-------------------------------------
			//Phi
			for(int phi=0;phi<16;phi++)
			{
				if(PhiRSB>=(phi-8.)/8.*pi && PhiRSB<=(phi-7.)/8.*pi){
					for(int m=0;m<9;m++)
					{
						if(pT_pair>=pT[m] && pT_pair<pT[m+1])
						{
							Npairs[m]=Npairs[m]+1;
							pTpairs[m]=pTpairs[m]+pT_pair;
							Mpairs[m]=Mpairs[m]+Minv;
							if(eta_pair>0){
							NgtB[m]=NgtB[m]+1;
							pTgtB[m]=pTgtB[m]+pT_pair;
							MgtB[m]=MgtB[m]+Minv;
							}
							if(eta_pair<0){
							NltB[m]=NltB[m]+1;
							pTltB[m]=pTltB[m]+pT_pair;
							MltB[m]=MltB[m]+Minv;
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
			///// BLUE BEAM ENDS ////////
			//// YELLOW BEAM ////////////
			//Phi loop 
			for(int phi=0;phi<16;phi++)
			{
				if(PhiRSY>=(phi-8.)/8.*pi && PhiRSY<=(phi-7.)/8.*pi)
				{	
					for(int m=0;m<9;m++)
					{
						if(pT_pair>=pT[m] && pT_pair<pT[m+1])
						{
							Npairs[m]=Npairs[m]+1;
							pTpairs[m]=pTpairs[m]+pT_pair;
							Mpairs[m]=Mpairs[m]+Minv;
							if(eta_pair<0){
							NgtY[m]=NgtY[m]+1;
							pTgtY[m]=pTgtY[m]+pT_pair;
							MgtY[m]=MgtY[m]+Minv;
							}
							if(eta_pair>0){
							NltY[m]=NltY[m]+1;
							pTltY[m]=pTltY[m]+pT_pair;
							MltY[m]=MltY[m]+Minv;
							}
							if(fspinconfig==51 || fspinconfig==83)
							{
								if(eta_pair<0)//eta<0 is forward direction for YELLOW beam
								{
									NpTGtUpY[m*16+phi]++;
								}
								if(eta_pair>0)//eta>0 is backward for YELLOW beam
								{
									NpTLtUpY[m*16+phi]++;
								}
							}
							if(fspinconfig==53 || fspinconfig==85)
							{
								//eta_pair<0 for yellow is similar to eta_pair>0 for blue
								if(eta_pair<0)//eta_pair direction should be flipped for yellow beam 
								{
									NpTGtDnY[m*16+phi]++;
								}
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
	 //YELLOW BEAM ENDS

		//Get average polarizations and rms from histograms 
		avgPolB = hpolB->GetMean();
		avgPolY = hpolY->GetMean();
		rmsB    = hpolB->GetRMS();
		rmsY    = hpolY->GetRMS();
		//cout<<"Blue: "<<avgPolB<<" Yellow: "<<avgPolY<<endl;
		double avg_M[9]={0};
		double avg_pT[9]={0};
		double avg_MgtB[9]={0};
		double avg_MltB[9]={0};
		double avg_pTgtB[9]={0};
		double avg_pTltB[9]={0};
		double avg_MgtY[9]={0};
		double avg_MltY[9]={0};
		double avg_pTgtY[9]={0};
		double avg_pTltY[9]={0};
		double avg_Mgt[9]={0};
		double avg_Mlt[9]={0};
		double avg_pTgt[9]={0};
		double avg_pTlt[9]={0};
			for(int j=0; j<9; j++)
			{
				avg_M[j]=(double)Mpairs[j]/(double)Npairs[j];
				avg_pT[j]=(double)pTpairs[j]/(double)Npairs[j];
				avg_MgtB[j]=(double)MgtB[j]/(double)NgtB[j];
				avg_pTgtB[j]=(double)pTgtB[j]/(double)NgtB[j];
				avg_MltB[j]=(double)MltB[j]/(double)NltB[j];
				avg_pTltB[j]=(double)pTltB[j]/(double)NltB[j];
				avg_MgtY[j]=(double)MgtY[j]/(double)NgtY[j];
				avg_pTgtY[j]=(double)pTgtY[j]/(double)NgtY[j];
				avg_MltY[j]=(double)MltY[j]/(double)NltY[j];
				avg_pTltY[j]=(double)pTltY[j]/(double)NltY[j];
				
				avg_Mgt[j]=.5*(avg_MgtB[j]+avg_MgtY[j]);
				avg_pTgt[j]=.5*(avg_pTgtB[j]+avg_pTgtY[j]);
				avg_Mlt[j]=.5*(avg_MltB[j]+avg_MltY[j]);
				avg_pTlt[j]=.5*(avg_pTltB[j]+avg_pTltY[j]);


			}


		double dAdB,dAdC,dAdD,dAdE,dAdP;//variables for error calculation
		//Use rms for polarization error from distribution!
		double dP_B = rmsB; // Polarization errors blue
		double dP_Y = rmsY; // Polarization errors yellow 
		double dA[165];// Asymmetry amplitude from sine fit 
		double B,C,D,E; 
		double a,b;
		double Asym[165];
		double BE, DC;

		// Asymmetry calculation for BLUE eta >0
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NpTGtUpB[m*16+ang];
					C = NpTGtUpB[m*16+ang+8];
					D = NpTGtDnB[m*16+ang];
					E = NpTGtDnB[m*16+ang+8];
					BE= (double)(B*E);
                                        DC= (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);
				}
				if(ang>7)
				{
					B = NpTGtUpB[m*16+ang];
					C = NpTGtUpB[m*16+ang-8];
					D = NpTGtDnB[m*16+ang];
					E = NpTGtDnB[m*16+ang-8];		
					BE= (double)(B*E);
                                        DC= (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);
				}
			Asym[ang]=(1./avgPolB)*((a-b)/(a+b));
                        dAdB = (1./avgPolB)*E*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdE = (1./avgPolB)*B*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdD = (-1./avgPolB)*C*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdC = (-1./avgPolB)*D*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdP = (-1./(avgPolB*avgPolB))*(sqrt(BE)-sqrt(D*C))/(sqrt(BE)+sqrt(DC));
                        dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_B)**2);

	}//angle loop
			double Abg[9];// Here Abg refers to Asymmetry for Blue in eta >0
			double deltaAbg[9];
			double chi2Ndf[9];
			double chi2[9];
			double angle[16]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi,pi/16, 3*pi/16, 5*pi/16, 7*pi/16,9*pi/16,11*pi/16,13*pi/16,15*pi/16};
			double ex[16]={0};
			char name[600];
			char title[600];
			double chiSquareCal[16]={0};
			double diff[16]={0};
			double fitVal[16]={0};
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *grt = new TGraphErrors(16,angle,Asym,ex,dA);
			grt->SetMarkerStyle(20);
			sprintf(title, "pT bin %i, BLUE #eta>0",m);
			grt->SetTitle(title);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,3.14159265359);
			fit->SetParameter(0,0.0001);
			grt->Fit(fit,"R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2[m] = fit->GetChisquare();
			cout<<"ChiSquare = "<<chi2Ndf[m]<<endl;
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->CenterTitle();
			grt->GetYaxis()->SetRangeUser(-0.09,0.09);
			
			Abg[m]=fit->GetParameter(0);
			deltaAbg[m]=fit->GetParError(0);
			
			for(int k=0; k<16; k++){
                        TF1 *fun = grt->GetFunction("fit");
                        fitVal[k]=fun->Eval(angle[k],0,0,0);//should give functional value at x=angle, y=0, z=0, t=0 
                        diff[k]=(Asym[k]-fitVal[k]);
                        chiSquareCal[m]+=(diff[k]*diff[k])/(dA[k]*dA[k]);
                        }

                        TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(0.5,0.082,Form("#font[12]{#color[2]{#chi^{2}_{cal} = %g}}",chiSquareCal[m]));


                        latex.SetTextAlign(13);
                        latex.SetTextSize(0.04);
                        latex.DrawLatex(0.5,0.068,"#font[12]{#color[2]{f(x)=[0]*sin[#phi_{RS}]}}");
                        gPad->Update();

			sprintf(name,"pT_bin%i_etaGt_Blue.png",m);
			//c1->Print(name);
		}//invariant mass loop
	
		//Asymmetry for BLUE , eta_pair<0 (Backward asymmetry)
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NpTLtUpB[m*16+ang];
					C = NpTLtUpB[m*16+ang+8];
					D = NpTLtDnB[m*16+ang];
					E = NpTLtDnB[m*16+ang+8];
					BE= (double)(B*E);
                                        DC= (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);
				}
				if(ang>7)
				{
					B = NpTLtUpB[m*16+ang];
					C = NpTLtUpB[m*16+ang-8];
					D = NpTLtDnB[m*16+ang];
					E = NpTLtDnB[m*16+ang-8];
					BE= (double)(B*E);
                                        DC= (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);		
				}
			Asym[ang]=(1./avgPolB)*((a-b)/(a+b));
                        dAdB = (1./avgPolB)*E*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdE = (1./avgPolB)*B*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdD = (-1./avgPolB)*C*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdC = (-1./avgPolB)*D*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdP = (-1./(avgPolB*avgPolB))*(sqrt(BE)-sqrt(D*C))/(sqrt(BE)+sqrt(DC));
                        dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_B)**2);
			}//angle loop
			double Abl[9];// Here Abl refers to Asymmetry for Blue in eta <0
			double deltaAbl[9];
			double chi2Ndf[9];
			double chi2[9];
			double angle[16]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi,pi/16, 3*pi/16, 5*pi/16, 7*pi/16,9*pi/16,11*pi/16,13*pi/16,15*pi/16};
			double ex[16]={0};
			char name[600];
			char title[600];
			double chiSquareCal[16]={0};
			double diff[16]={0};
			double fitVal[16]={0};
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *grbl = new TGraphErrors(16,angle,Asym,ex,dA);
			grbl->SetMarkerStyle(20);
			sprintf(title, "pT bin %i,BLUE #eta<0",m);
			grbl->SetTitle(title);
			grbl->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,3.14159265359);
			//fit->SetParameter(0,0.0001);
			grbl->Fit(fit,"R");
			grbl->SetMarkerColor(4);
			grbl->GetXaxis()->SetTitle("#Phi_{RS}");
			grbl->GetXaxis()->SetTitleOffset(1);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2[m] = fit->GetChisquare();
			cout<<"ChiSquare = "<<chi2Ndf[m]<<endl;
			grbl->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grbl->GetYaxis()->CenterTitle();
			grbl->GetYaxis()->SetRangeUser(-0.09,0.09);
		
			Abl[m]=fit->GetParameter(0);
			deltaAbl[m]=fit->GetParError(0);
		
			for(int k=0; k<16; k++){
                        TF1 *fun = grbl->GetFunction("fit");
                        fitVal[k]=fun->Eval(angle[k],0,0,0);//should give functional value at x=angle, y=0, z=0, t=0 
                        diff[k]=(Asym[k]-fitVal[k]);
                        chiSquareCal[m]+=(diff[k]*diff[k])/(dA[k]*dA[k]);
                        }

                        TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(0.5,0.082,Form("#font[12]{#color[2]{#chi^{2}_{cal} = %g}}",chiSquareCal[m]));


                        latex.SetTextAlign(13);
                        latex.SetTextSize(0.04);
                        latex.DrawLatex(0.5,0.068,"#font[12]{#color[2]{f(x)=[0]*sin[#phi_{RS}]}}");
                        gPad->Update();
		
			sprintf(name,"pT_bin%i_etaLt_Blue.png",m);
			//c1->Print(name);
		}//invariant mass loop
	
	
		//Asymmetry for YELLOW beam , eta > 0  
		for(int m=0;m<9;m++)
		{	
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NpTGtUpY[m*16+ang];
					C = NpTGtUpY[m*16+ang+8];
					D = NpTGtDnY[m*16+ang];
					E = NpTGtDnY[m*16+ang+8];
                                        BE= (double)(B*E);
                                        DC= (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);

				}
				if(ang>7)
				{
					B = NpTGtUpY[m*16+ang];
					C = NpTGtUpY[m*16+ang-8];
					D = NpTGtDnY[m*16+ang];
					E = NpTGtDnY[m*16+ang-8];
                                        BE= (double)(B*E);
                                        DC= (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);
				}
			Asym[ang]=(1./avgPolY)*((a-b)/(a+b));
                        dAdB = (1./avgPolY)*E*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdE = (1./avgPolY)*B*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdD = (-1./avgPolY)*C*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdC = (-1./avgPolY)*D*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdP = (-1./(avgPolY*avgPolY))*(sqrt(BE)-sqrt(D*C))/(sqrt(BE)+sqrt(DC));
                        dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_Y)**2);
			}
			double Ayg[9];
			double deltaAyg[9];
			double chi2Ndf[9];
			double chi2[9];
			double angle[16]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi,pi/16, 3*pi/16, 5*pi/16, 7*pi/16,9*pi/16,11*pi/16,13*pi/16,15*pi/16};
			double ex[16]={0};
			char name[600];
			char title[600];
			double chiSquareCal[16]={0};
			double diff[16]={0};
			double fitVal[16]={0};
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *gr = new TGraphErrors(16,angle,Asym,ex,dA);
			gr->SetMarkerStyle(20);
			sprintf(title, "pT bin %i, YELLOW #eta>0", m);
			gr->SetTitle(title);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,3.14159265359);
			fit->SetParameter(0,0.0001);
			gr->Fit(fit,"R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2[m] = fit->GetChisquare();
			cout<<"ChiSquare= "<<chi2Ndf[m]<<endl;
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle();
			gr->GetYaxis()->SetRangeUser(-0.09,0.09);
		
			Ayg[m]=fit->GetParameter(0);
			deltaAyg[m]=fit->GetParError(0);
		
			for(int k=0; k<16; k++){
                        TF1 *fun = gr->GetFunction("fit");
                        fitVal[k]=fun->Eval(angle[k],0,0,0);//should give functional value at x=angle, y=0, z=0, t=0 
                        diff[k]=(Asym[k]-fitVal[k]);
                        chiSquareCal[m]+=(diff[k]*diff[k])/(dA[k]*dA[k]);
                        }

                        TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(0.5,0.082,Form("#font[12]{#color[2]{#chi^{2}_{cal} = %g}}",chiSquareCal[m]));


                        latex.SetTextAlign(13);
                        latex.SetTextSize(0.04);
                        latex.DrawLatex(0.5,0.068,"#font[12]{#color[2]{f(x)=[0]*sin[#phi_{RS}]}}");
                        gPad->Update();
			sprintf(name,"pT_bin%i_etaGt_Yellow.png",m);
			//c1->Print(name);

		}//pT loop
	
		//Asymmetry for Yellow beam , eta<0 (backward asymmetry)
		for(int m=0;m<9;m++)
		{	
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NpTLtUpY[m*16+ang];
					C = NpTLtUpY[m*16+ang+8];
					D = NpTLtDnY[m*16+ang];
					E = NpTLtDnY[m*16+ang+8];
					BE = (double)(B*E);
                                        DC = (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);
				}
				if(ang>7)
				{
					B = NpTLtUpY[m*16+ang];
					C = NpTLtUpY[m*16+ang-8];
					D = NpTLtDnY[m*16+ang];
					E = NpTLtDnY[m*16+ang-8];
					BE = (double)(B*E);
                                        DC = (double)(D*C);
                                        a = sqrt(BE);
                                        b = sqrt(DC);
				}
			Asym[ang]=(1./avgPolY)*((a-b)/(a+b));
                        dAdB = (1./avgPolY)*E*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdE = (1./avgPolY)*B*sqrt(DC)/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
                        dAdD = (-1./avgPolY)*C*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdC = (-1./avgPolY)*D*sqrt(BE)/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
                        dAdP = (-1./(avgPolY*avgPolY))*(sqrt(BE)-sqrt(D*C))/(sqrt(BE)+sqrt(DC));
                        dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_Y)**2);	
			}
			double Ayl[9];
			double deltaAyl[9];
			double chi2Ndf[9];
			double chi2[9];
			double angle[16]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi,pi/16, 3*pi/16, 5*pi/16, 7*pi/16,9*pi/16,11*pi/16,13*pi/16,15*pi/16};
			double ex[16]={0};
			char name[600];
			char title[600];
			double chiSquareCal[16]={0};
			double diff[16]={0};
			double fitVal[16]={0};
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *gryl = new TGraphErrors(16,angle,Asym,ex,dA);
			gryl->SetMarkerStyle(20);
			sprintf(title,"pT bin %i, YELLOW #eta<0",m);
			gryl->SetTitle(title);
			gryl->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,3.14159265359);
			fit->SetParameter(0,0.0001);
			gryl->Fit(fit,"R");
			gryl->SetMarkerColor(4);
			gryl->GetXaxis()->SetTitle("#Phi_{RS}");
			gryl->GetXaxis()->SetTitleOffset(1);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2[m] = fit->GetChisquare();
			cout<<"ChiSquare= "<<chi2Ndf[m]<<endl;
			gryl->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gryl->GetYaxis()->CenterTitle();
			gryl->GetYaxis()->SetRangeUser(-0.09,0.09);
		
			Ayl[m]=fit->GetParameter(0);
			deltaAyl[m]=fit->GetParError(0);
		
			for(int k=0; k<16; k++){
                        TF1 *fun = gryl->GetFunction("fit");
                        fitVal[k]=fun->Eval(angle[k],0,0,0);//should give functional value at x=angle, y=0, z=0, t=0 
                        diff[k]=(Asym[k]-fitVal[k]);
                        chiSquareCal[m]+=(diff[k]*diff[k])/(dA[k]*dA[k]);
                        }

                        TLatex latex;
                        latex.SetTextSize(0.04);
                        latex.SetTextAlign(13);
                        latex.DrawLatex(0.5,0.082,Form("#font[12]{#color[2]{#chi^{2}_{cal} = %g}}",chiSquareCal[m]));


                        latex.SetTextAlign(13);
                        latex.SetTextSize(0.04);
                        latex.DrawLatex(0.5,0.068,"#font[12]{#color[2]{f(x)=[0]*sin[#phi_{RS}]}}");
                        gPad->Update();
			sprintf(name,"pT_bin%i_etaLt_Yellow.png",m);
			//c1->Print(name);

		}//pT loop
		
			//for Blue beam eta > 0
			double	A_pTbg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_pTbg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
			//for Blue beam eta < 0
			double	A_pTbl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4],  Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_pTbl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
			// for Yellow beam eta >0
			double A_pTyg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4],  Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pTyg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};
			
			double A_pTyl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pTyl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		
	

			//calculate average asymmetry  
			double avgA_pTg[9] ={0};
			double avgA_pTl[9] ={0};
				for (Int_t ii=0; ii<9; ii++)
				{
					avgA_pTg[ii]= (A_pTbg[ii]+A_pTyg[ii])/2.;
					avgA_pTl[ii]= (A_pTbl[ii]+A_pTyl[ii])/2.;
				}	
			//total asymmetry error calculation
			//since two asymmetry(BLUE+YELLOW) are averaged, error propagation formaula is used for total asymmetry error 
			//if function, f = (a+b)/2, Then Error, df =1/2 √{(∂f/∂a)^2*(∆a)^2+(∂f/∂b)^2*(∆b)^2}
			double errA_pTg[9] ={0};
			double errA_pTl[9] ={0};
				for (Int_t iii=0; iii<9; iii++)
				{
					errA_pTg[iii]= .5*sqrt(pow(deltaA_pTbg[iii],2)+pow(deltaA_pTyg[iii],2));//total asym err, pT bin 1, eta > 0
					errA_pTl[iii]= .5*sqrt(pow(deltaA_pTbl[iii],2)+pow(deltaA_pTyl[iii],2));//total asym err, pT bin 1, eta < 0
				}	
		cout << "/////////////////////Average Asymmetry (pT) Values/////////////////////"<< endl;
               	 cout << "Avg_MgtB[9] = {" <<avg_MgtB[0]<<", "<<avg_MgtB[1]<<", "<<avg_MgtB[2]<<", "<<avg_MgtB[3]<<", "<<avg_MgtB[4]<<", "<<avg_MgtB[5]<<", "<<avg_MgtB[6]<<", "<<avg_MgtB[7]<<", "<<avg_MgtB[8]<<"}"<< endl;
                cout << "Avg_MltB[9] = {" <<avg_MltB[0]<<", "<<avg_MltB[1]<<", "<<avg_MltB[2]<<", "<<avg_MltB[3]<<", "<<avg_MltB[4]<<", "<<avg_MltB[5]<<", "<<avg_MltB[6]<<", "<<avg_MltB[7]<<", "<<avg_MltB[8]<<"}"<< endl;
                cout << "Avg_pTgtB[9] = {"<<avg_pTgtB[0]<<", "<<avg_pTgtB[1]<<", "<<avg_pTgtB[2]<<", "<<avg_pTgtB[3]<<", "<<avg_pTgtB[4]<<", "<<avg_pTgtB[5]<<", "<<avg_pTgtB[6]<<", "<<avg_pTgtB[7]<<", "<<avg_pTgtB[8]<<"}"<< endl;
                cout << "Avg_pTltB[9] = {"<<avg_pTltB[0]<<", "<<avg_pTltB[1]<<", "<<avg_pTltB[2]<<", "<<avg_pTltB[3]<<", "<<avg_pTltB[4]<<", "<<avg_pTltB[5]<<", "<<avg_pTltB[6]<<", "<<avg_pTltB[7]<<", "<<avg_pTltB[8]<<"}"<< endl;
                cout << "AsymBlueF[9]={"<<Abg[0]<<", "<< Abg[1]<<", "<< Abg[2]<<", "<< Abg[3]<<", "<< Abg[4]<<", "<<Abg[5]<<", "<<Abg[6]<<", "<<Abg[7]<<", "<<Abg[8]<<"};"<<endl; 
                cout << "errBlueF[9]={"<<deltaAbg[0]<<", "<< deltaAbg[1]<<", "<< deltaAbg[2]<<", "<< deltaAbg[3]<<", "<< deltaAbg[4]<<", "<<deltaAbg[5]<<", "<<deltaAbg[6]<<", "<<deltaAbg[7]<<", "<<deltaAbg[8]<<"}"<<endl; 
                cout << "AsymBlueB[9]={"<<Abl[0]<<", "<< Abl[1]<<", "<< Abl[2]<<", "<< Abl[3]<<", "<< Abl[4]<<", "<<Abl[5]<<", "<<Abl[6]<<", "<<Abl[7]<<", "<<Abl[8]<<"}"<<endl; 
                cout << "errBlueB[9]={"<<deltaAbl[0]<<", "<< deltaAbl[1]<<", "<< deltaAbl[2]<<", "<< deltaAbl[3]<<", "<< deltaAbl[4]<<", "<<deltaAbl[5]<<", "<<deltaAbl[6]<<", "<<deltaAbl[7]<<", "<<deltaAbl[8]<<"}"<<endl; 
                cout << "AsymYellowF[9]={"<<Ayg[0]<<", "<< Ayg[1]<<", "<< Ayg[2]<<", "<< Ayg[3]<<", "<< Ayg[4]<<", "<<Ayg[5]<<", "<<Ayg[6]<<", "<<Ayg[7]<<", "<<Ayg[8]<<"}"<<endl; 
                cout << "errYellowF[9]={"<<deltaAyg[0]<<", "<< deltaAyg[1]<<", "<< deltaAyg[2]<<", "<< deltaAyg[3]<<", "<< deltaAyg[4]<<", "<<deltaAyg[5]<<", "<<deltaAyg[6]<<", "<<deltaAyg[7]<<", "<<deltaAyg[8]<<"}"<<endl; 
                cout << "AsymYellowB[9]={"<<Ayl[0]<<", "<< Ayl[1]<<", "<< Ayl[2]<<", "<< Ayl[3]<<", "<< Ayl[4]<<", "<<Ayl[5]<<", "<<Ayl[6]<<", "<<Ayl[7]<<", "<<Ayl[8]<<"}"<<endl; 
                cout << "errYellowB[9]={"<<deltaAyl[0]<<", "<< deltaAyl[1]<<", "<< deltaAyl[2]<<", "<< deltaAyl[3]<<", "<< deltaAyl[4]<<", "<<deltaAyl[5]<<", "<<deltaAyl[6]<<", "<<deltaAyl[7]<<", "<<deltaAyl[8]<<"}"<<endl; 
		cout << "avgAsympTF[9]={"<<avgA_pTg[0]<<", "<<avgA_pTg[1]<< ", "<<avgA_pTg[2]<< ", "<<avgA_pTg[3]<<", "<<avgA_pTg[4]<<", "<<avgA_pTg[5]<<", "<<avgA_pTg[6]<<", "<<avgA_pTg[7]<<", "<<avgA_pTg[8]<<"};"<<endl;
                cout << "avgErrpTF[9]={"<<errA_pTg[0]<<", "<<errA_pTg[1]<< ", "<<errA_pTg[2]<< ", "<<errA_pTg[3]<<", "<<errA_pTg[4]<<", "<<errA_pTg[5]<<", "<<errA_pTg[6]<<", "<<errA_pTg[7]<<", "<<errA_pTg[8]<<"};"<<endl;
                cout << "avgAsympTB[9]={"<<avgA_pTl[0]<<", "<<avgA_pTl[1]<< ", "<<avgA_pTl[2]<< ", "<<avgA_pTl[3]<<", "<<avgA_pTl[4]<<", "<<avgA_pTl[5]<<", "<<avgA_pTl[6]<<", "<<avgA_pTl[7]<<", "<<avgA_pTl[8]<<"};"<<endl;
                cout << "avgErrpTB[9]={"<<errA_pTl[0]<<", "<<errA_pTl[1]<< ", "<<errA_pTl[2]<< ", "<<errA_pTl[3]<<", "<<errA_pTl[4]<<", "<<errA_pTl[5]<<", "<<errA_pTl[6]<<", "<<errA_pTl[7]<<", "<<errA_pTl[8]<<"};"<<endl;
                cout << "/////////////////////////////////////////////////////////////////"<< endl;





		//Data from Run2006 analysis paper

		double asym06f[9]={.0092, .011, .029, .0071, .053};
                double err06f[9]={.010, .0086, .0082, .011, .011};
                double asym06b[9]={-.022, -.007, .018, .021, .0088};
                double err06b[9]={.010, .0087, .0082, .011, .011};
                double minv06f[9]={.45, .51, .59, .67, .81};
                double pt06f[9]={3.62, 4.49, 5.70, 7.18, 10.45};
                double minv06b[9]={.38, .41, .47, .55, .68};
                double pt06b[9]={3.61, 4.49, 5.68, 7.17, 10.32};
		

		// eta > 0 
		TCanvas *myCan = new TCanvas("myCan","myCan",500,550);
		gStyle -> SetOptStat(0);
		gStyle -> SetTitle("");
		gStyle->SetLegendBorderSize(0);
		myCan -> Divide(2,1);
		myCan -> cd(1);
		gPad->SetPad(0.0,.35,1.,1.);
		gPad->SetGrid(0,0);
		gPad->SetBottomMargin(0.01);
		double err[9]={0};
		//forward asymmetry for BLUE
		TGraphErrors *gr1 = new TGraphErrors(9,avg_pTgtB,A_pTbg,err,deltaA_pTbg);
		gr1->SetTitle("");
		gr1->GetYaxis()-> SetTitle("A_{UT}");
		gr1->GetYaxis()-> CenterTitle();
	  	gr1->GetYaxis()->SetTitleOffset(1.15);
	  	gr1->GetYaxis()->SetLabelSize(0.045);
	  	gr1->GetXaxis()->SetLabelSize(0);
	  	gr1->GetYaxis()->SetTitleSize(0.045);
		gr1-> SetMarkerStyle(21);
		gr1-> SetMarkerColor(kBlue);
		gr1-> GetXaxis()->SetRangeUser(0.,11);
		gr1-> GetYaxis()->SetRangeUser(-0.03, 0.12);
		gr1->Draw("AP");
		myCan->Update();
 		//TLine *line1=  new TLine(myCan->cd(1)->GetUxmin(),0.,myCan->cd(1)->GetUxmax(),0.);
 		TLine *line1=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
		line1->SetLineStyle(2);
		line1->Draw();
		//Forward asymmetry for yellow
		TGraphErrors *gr3 = new TGraphErrors(9,avg_pTgtY,A_pTyg,err,deltaA_pTyg);
		gr3-> SetMarkerStyle(25);
		gr3-> SetMarkerColor(2);
		gr3->Draw("same P");
		TLegend *leg1 = new TLegend(0.2,.65, 0.5, 0.85);
		leg1->AddEntry("", "#eta > 0, Cone < 0.7 ", "");
		leg1->AddEntry(gr1, "Blue beam ", "lp");
		leg1->AddEntry(gr3, "Yellow beam ", "lp");
		//leg1->AddEntry(gr4, "< A_{UT} > ", "lp");
		leg1->Draw();
                gPad->Update();


		myCan->cd(2);
		gPad->SetPad(0.0,0.0,1.,.35);
		gPad->SetGrid(0,0);
		gPad->SetTopMargin(.02);
		gPad->SetBottomMargin(0.2);
		TGraphErrors *gr2= new TGraphErrors(9,avg_pTgtB,avg_MgtB,err,errA_pTg);
		gr2->SetTitle("");
		gr2->GetYaxis()-> SetTitle("< M^{#pi^{+}#pi^{-}_{inv}} > (GeV/c^{2})");
		gr2->GetYaxis()-> CenterTitle();
	  	gr2->GetYaxis()->SetTitleOffset(0.6);
	  	gr2->GetYaxis()->SetTitleSize(.08);
	  	gr2->GetXaxis()->SetRangeUser(0,11);
	  	gr2->GetYaxis()->SetRangeUser(0,1);
	  	gr2->GetXaxis()->SetTitleSize(.08);
	  	gr2->GetXaxis()->SetTitleOffset(1.2);
	  	gr2->GetXaxis()->SetLabelOffset(.01);
	  	gr2->GetYaxis()->SetLabelOffset(.01);
	  	gr2->GetYaxis()->SetLabelSize(0.08);
		gr2->GetXaxis()-> SetTitle("pT^{#pi^{+}#pi^{-}} (GeV/c)");
	  	gr2->GetXaxis()->SetLabelSize(0.08);
		gr2->GetXaxis()-> CenterTitle();
		gr2->SetMarkerStyle(21);
		gr2->SetMarkerSize(1);
		gr2->SetMarkerColor(kBlue);
		gr2->Draw("AP");
		TGraphErrors *grYgt= new TGraphErrors(9,avg_pTgtY,avg_MgtY,err,errA_pTl);
		grYgt->SetMarkerStyle(25);
		grYgt->SetMarkerColor(2);
		grYgt->Draw("same P");
		myCan->Update();

//myCan--> Forward asymmetry for BLUE and YELLOW -----DONE



//myCan2 ----> Bavkward asymmetry for BLUE ans YELLOW		
		TCanvas *myCan2 = new TCanvas("myCan2","myCan2",500,550);
		gStyle -> SetOptStat(0);
		gStyle -> SetTitle("");
		gStyle->SetLegendBorderSize(0);
		myCan2 -> Divide(2,1);
		myCan2 -> cd(1);
		gPad->SetPad(0.0,.35,1.,1.);
		gPad->SetGrid(0,0);
		gPad->SetBottomMargin(0.01);
		double err[9]={0};
		TGraphErrors *gr3 = new TGraphErrors(9,avg_pTltB,A_pTbl,err,deltaA_pTbl);
		gr3->SetTitle("");
		gr3->GetYaxis()-> SetTitle("A_{UT}");
		gr3->GetYaxis()-> CenterTitle();
	  	gr3->GetYaxis()->SetTitleOffset(1.15);
	  	gr3->GetYaxis()->SetLabelSize(0.045);
	  	gr3->GetYaxis()->SetTitleSize(0.045);
		gr3-> SetMarkerStyle(21);
		gr3-> SetMarkerColor(kBlue);
		gr3-> GetXaxis()->SetRangeUser(0.,11);
                gr3-> GetXaxis()->SetLabelSize(0);
		gr3-> GetYaxis()->SetRangeUser(-0.03, 0.12);
		gr3->Draw("AP");
		myCan2->Update();
 		//TLine *line1=  new TLine(myCan->cd(1)->GetUxmin(),0.,myCan->cd(1)->GetUxmax(),0.);
 		TLine *line2=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
		line2->SetLineStyle(2);
		line2->Draw();
		
		TGraphErrors *gr4 = new TGraphErrors(9,avg_pTltY,A_pTyl,err,deltaA_pTyl);
		gr4-> SetMarkerStyle(25);
		gr4-> SetMarkerColor(2);
		gr4->Draw("same P");
		TLegend *leg2 = new TLegend(0.2,.65, 0.5, 0.85);
		leg2->AddEntry("", "#eta < 0, Cone < 0.7 ", "");
		leg2->AddEntry(gr3, "Blue beam ", "lp");
		leg2->AddEntry(gr4, "Yellow beam ", "lp");
		//leg2->AddEntry(gr5, "< A_{UT} > ", "lp");
		leg2->Draw();
                gPad->Update();
		

		myCan2->cd(2);
		gPad->SetPad(0.0,0.0,1.,.35);
		gPad->SetGrid(0,0);
		gPad->SetTopMargin(.02);
		gPad->SetBottomMargin(0.2);
		TGraphErrors *gr6= new TGraphErrors(9,avg_pTltB,avg_MltB,err,err);
		gr6->SetTitle("");
		gr6->GetYaxis()-> SetTitle("<M^{#pi^{+}#pi^{-}}_{inv} > (GeV/c^{2})");
		gr6->GetYaxis()-> CenterTitle();
	  	gr6->GetYaxis()->SetTitleOffset(0.6);
	  	gr6->GetYaxis()->SetTitleSize(.08);
	  	gr6->GetXaxis()->SetRangeUser(0,11);
	  	gr6->GetYaxis()->SetRangeUser(0,1);
	  	gr6->GetXaxis()->SetTitleSize(.08);
	  	gr6->GetXaxis()->SetTitleOffset(1.2);
	  	gr6->GetXaxis()->SetLabelOffset(.01);
	  	gr6->GetYaxis()->SetLabelOffset(.01);
	  	gr6->GetYaxis()->SetLabelSize(0.08);
		gr6->GetXaxis()-> SetTitle("pT^{#pi^{+}#pi^{-}} (GeV/c)");
	  	gr6->GetXaxis()->SetLabelSize(0.08);
		gr6->GetXaxis()-> CenterTitle();
		gr6->SetMarkerStyle(21);
		gr6->SetMarkerSize(1);
		gr6->SetMarkerColor(kBlue);
		gr6->Draw("AP");
		TGraphErrors *gr6ltY= new TGraphErrors(9,avg_pTltY,avg_MltY,err,err);
		gr6ltY->SetMarkerStyle(25);
                gr6ltY->SetMarkerColor(2);
                gr6ltY->Draw("same P");
		
		myCan2->Update();
//myCan2---> Done drawing Backward asymmetry --------------->>>	


	
		//Plot average asymmetry
		TCanvas *myCan3 = new TCanvas("myCan3","myCan3",500,550);
                gStyle -> SetOptStat(0);
                gStyle -> SetTitle("");
                gStyle->SetLegendBorderSize(0);
                myCan3 -> Divide(2,1);
                myCan3 -> cd(1);
                gPad->SetPad(0.0,.35,1.,1.);
                gPad->SetGrid(0,0);
                gPad->SetBottomMargin(0.01);
                TGraphErrors *grAg = new TGraphErrors(9,avg_pTgt,avgA_pTg,err,errA_pTg);
                grAg->SetTitle("");
                grAg->GetYaxis()-> SetTitle("A_{UT}");
                grAg->GetYaxis()-> CenterTitle();
                grAg->GetYaxis()->SetTitleOffset(1.15);
                grAg->GetYaxis()->SetLabelSize(0.045);
                grAg->GetYaxis()->SetTitleSize(0.045);
                grAg-> SetMarkerStyle(21);
                grAg-> SetMarkerColor(kRed);
                grAg-> GetXaxis()->SetRangeUser(0.,11);
                grAg-> GetXaxis()->SetLabelSize(0);
                grAg-> GetYaxis()->SetRangeUser(-0.03, 0.12);
                grAg->Draw("AP");
                TGraphErrors *grAl = new TGraphErrors(9,avg_pTlt,avgA_pTl,err,errA_pTl);
		grAl-> SetMarkerStyle(20);
                grAl-> SetMarkerColor(4);
		grAl->Draw("same P");
		myCan3->Update();
		TLine *lineA=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
                lineA->SetLineStyle(2);
                lineA->Draw();
                TLegend *legA = new TLegend(0.19,.6, 0.5, 0.85);
                legA->AddEntry("", "p#uparrow + p, #sqrt{s} = 200 GeV", "");
                legA->AddEntry("", "Cone < 0.7", "");
                legA->AddEntry(grAg, "< A_{UT} > , #eta > 0 ", "lp");
                legA->AddEntry(grAl, "< A_{UT} > , #eta < 0", "lp");
                legA->SetTextSize(.04);
                legA->Draw();
		myCan3->Update();
		

		myCan3->cd(2);
                gPad->SetPad(0.0,0.0,1.,.35);
                gPad->SetGrid(0,0);
                gPad->SetTopMargin(.02);
                gPad->SetBottomMargin(0.2);
                TGraphErrors *grB= new TGraphErrors(9,avg_pTgt,avg_Mgt,err,err);
                grB->SetTitle("");
                grB->GetYaxis()-> SetTitle("< M^{#pi^{+}#pi^{-}}_{inv} > (GeV/c^{2})");
                grB->GetYaxis()-> CenterTitle();
                grB->GetYaxis()->SetTitleOffset(0.6);
                grB->GetYaxis()->SetTitleSize(.08);
                grB->GetXaxis()->SetRangeUser(0,11);
                grB->GetYaxis()->SetRangeUser(0,1);
                grB->GetXaxis()->SetTitleSize(.08);
                grB->GetXaxis()->SetTitleOffset(1.1);
                grB->GetXaxis()->SetLabelOffset(.01);
                grB->GetYaxis()->SetLabelOffset(.01);
                grB->GetYaxis()->SetLabelSize(0.08);
                grB->GetXaxis()-> SetTitle("p_{T}^{#pi^{+}#pi^{-}} (GeV/c)");
                grB->GetXaxis()->SetLabelSize(0.08);
                grB->GetXaxis()-> CenterTitle();
                grB->SetMarkerStyle(21);
                grB->SetMarkerColor(kRed);
                grB->SetMarkerSize(1.);
                grB->Draw("AP");
		gPad->Update();
                TGraphErrors *grLt= new TGraphErrors(9,avg_pTlt,avg_Mlt,err,err);
		grLt->SetMarkerStyle(25);
                grLt->SetMarkerColor(2);
                grLt->SetMarkerSize(1.);
                grLt->Draw("same P");
		myCan3->Update();

//Done--------->Drawing average asymmetry in forward and backward----------->>>		

    /* 
//Draw comparison plot in forward direction---------------->>>>
		TCanvas *myCanGt = new TCanvas("myCanGt","myCanGt",550,500);
                gStyle -> SetOptStat(0);
                gStyle -> SetTitle("");
                gStyle->SetLegendBorderSize(0);
                gPad->SetGrid(0,0);
                gPad->SetLeftMargin(0.1);
                TGraphErrors *grAgt = new TGraphErrors(5,avg_pTgt,avgA_pTg,err,errA_pTg);
                grAgt->SetTitle("");
                grAgt->GetYaxis()-> SetTitle("A_{UT}");
                grAgt->GetXaxis()-> SetTitle("p_{T}^{#pi^{+}#pi^{-}}");
                grAgt->GetYaxis()-> CenterTitle();
                grAgt->GetXaxis()-> CenterTitle();
                grAgt->GetYaxis()->SetTitleOffset(1.12);
                grAgt->GetXaxis()->SetTitleOffset(.9);
                grAgt->GetYaxis()->SetLabelSize(0.035);
                grAgt-> SetMarkerStyle(21);
                grAgt-> SetMarkerColor(kRed);
                grAgt-> GetXaxis()->SetRangeUser(0.,11);
                grAgt-> GetXaxis()->SetLabelSize(0.035);
                grAgt-> GetYaxis()->SetRangeUser(-0.03, 0.12);
                grAgt->Draw("AP");
	myCanGt->Update();
		TGraphErrors *gr06f = new TGraphErrors(5,pt06f, asym06f,err,err06f);
                gr06f-> SetMarkerStyle(25);
                gr06f-> SetMarkerColor(kRed);
                gr06f->Draw("same P");

                TLegend *leg06 = new TLegend(0.14,.64, 0.4, 0.85);
                leg06->AddEntry("", "p#uparrow + p, #sqrt{s} = 200 GeV ", " ");
                leg06->AddEntry("", "#eta > 0", " ");
                leg06->AddEntry(grAgt, "< A_{UT} >, 2015", "lp");
                leg06->AddEntry(gr06f, "< A_{UT} >, 2006 ", "lp");
    		leg06->SetTextSize(0.03);
	        leg06->Draw();
		gPad->Update();
		 TLine *line06f=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
                line06f->SetLineStyle(2);
                line06f->Draw();
	myCanGt->Update();

//Comparison for FORWARD asymmetry done-----------------///
//Draw comparison plot for BACKWARD asymmetry-------------->>>>
		TCanvas *myCanLt = new TCanvas("myCanLt","myCanLt",550,500);
                gStyle -> SetOptStat(0);
                gStyle -> SetTitle("");
                gStyle->SetLegendBorderSize(0);
                gPad->SetGrid(0,0);
                gPad->SetLeftMargin(0.1);
                TGraphErrors *grAlt = new TGraphErrors(5,avg_pTlt,avgA_pTl,err,errA_pTl);
                grAlt->SetTitle("");
                grAlt->GetYaxis()-> SetTitle("A_{UT}");
                grAlt->GetXaxis()-> SetTitle("p_{T}^{#pi^{+}#pi^{-}} (GeV/c)");
                grAlt->GetYaxis()-> CenterTitle();
                grAlt->GetXaxis()-> CenterTitle();
                grAlt->GetYaxis()->SetTitleOffset(1.12);
                grAlt->GetXaxis()->SetTitleOffset(.9);
                grAlt->GetYaxis()->SetLabelSize(0.035);
                grAlt-> SetMarkerStyle(21);
                grAlt-> SetMarkerColor(kRed);
                grAlt-> GetXaxis()->SetRangeUser(0.,11);
                grAlt-> GetXaxis()->SetLabelSize(0.035);
                grAlt-> GetYaxis()->SetRangeUser(-0.035, 0.12);
                grAlt->Draw("AP");

                TGraphErrors *gr06b = new TGraphErrors(5,pt06b, asym06b,err,err06b);
                gr06b-> SetMarkerStyle(25);
                gr06b-> SetMarkerColor(2);
                gr06b->Draw("same P");

                TLegend *leg07 = new TLegend(0.14,.64, 0.4, 0.85);
                leg07->AddEntry("", "p#uparrow + p, #sqrt{s} = 200 GeV ", " ");
                leg07->AddEntry("", "#eta < 0", " ");
                leg07->AddEntry(grAlt, "< A_{UT} >, 2015", "lp");
                leg07->AddEntry(gr06b, "< A_{UT} >, 2006 ", "lp");
    		leg07->SetTextSize(0.03);
	        leg07->Draw();
		gPad->Update();
		 TLine *line07b=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
                line07b->SetLineStyle(2);
                line07b->Draw();
           myCanLt->Update();

*/
	myCan->SaveAs("AsymVspT_Forward.pdf");
	myCan2->SaveAs("AsymVspT_Backward.pdf");
	myCan3->SaveAs("AsymVspT_Average.pdf");
	//myCanGt->SaveAs("AsymVspT_ConeLt3_Run15VSRun06GT.png");
	//myCanLt->SaveAs("AsymVspT_ConeLt3_Run15VSRun06LT.pn	i
}//main
