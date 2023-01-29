//This code is modefied to calculate asymmetry for individual beams and averaged for total asymmetry
#include <iostream>
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

ofstream Output;

void integratedAsym_Minv_9bins(const char *ifile="/star/u/pokhrel/GPFS/IFF_TREES/StartLessTOF/iffNtuplesFinalA.root")
{
	TH1D *hpolB = new TH1D("hpolB","", 100, 0, 1);	
	TH1D *hpolY = new TH1D("hpolY","", 100, 0, 1);	
	TFile *f = new TFile(ifile);

	TTree *ntuple1tpc = (TTree*)f->Get("ntuple1tpc"); //get trees 
	TTree *ntuple2tpc = (TTree*)f->Get("ntuple2tpc");
	TTree *ntuple4tpc = (TTree*)f->Get("ntuple4tpc");
	TTree *ntuple5tpc = (TTree*)f->Get("ntuple5tpc");
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

	Output.open("DataYieldsTPConly.txt");
	ofstream out;
	out.open("intAsymMinvTPConly.txt");


	//get the values of variables from tree 
	

	ntuple1tpc->SetBranchAddress("fspinconfig",&fspinconfig);
	ntuple1tpc->SetBranchAddress("cone",&cone);
	ntuple2tpc->SetBranchAddress("Minv",&Minv);
	ntuple2tpc->SetBranchAddress("pT_pair",&pT_pair);
	ntuple2tpc->SetBranchAddress("eta_pair",&eta_pair);
	ntuple4tpc->SetBranchAddress("PhiRSB",&PhiRSB);
	ntuple4tpc->SetBranchAddress("PhiRSY",&PhiRSY);
	ntuple5tpc->SetBranchAddress("fitPts_min_pair",&fitPts_min_pair);
	ntuple6->SetBranchAddress("polB_corr", &polB_corr);
	ntuple6->SetBranchAddress("polY_corr", &polY_corr);

	int nentries = (int)ntuple1tpc->GetEntries(); //entries of trees 
	int ne = (int)ntuple6->GetEntries(); //entries of trees 
	//add friend to the tree. no root file to add friend means the tree is on the same root file 
	ntuple1tpc->AddFriend("ntuple2tpc"); 
	ntuple1tpc->AddFriend("ntuple4tpc");
	ntuple1tpc->AddFriend("ntuple5tpc");
	
	//Binning for Cone < 0.7
	Double_t M[10]={0.20, 0.4035, 0.517, 0.6124, 0.7108, 0.8031, 0.9553, 1.1224, 1.379, 4};//binning for 0.05< cone < .3 with dodiging K0 mass range (0.4876 – 0.5076)      

	//Binning for Cone<0.3
	//Double_t M[10]={0.20, 0.3529, 0.4067, 0.4572, 0.5227, 0.5727, 0.6328, 0.7113, 0.8172, 4.000};//binning for 0.05< cone < .3 with dodiging K0 mass range (0.4876 – 0.5076)      
	

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

		for(int p=0; p<ne; p++)
		//for(int p=0; p<100000; p++)
		{
			ntuple6->GetEntry(p);	
			hpolB->Fill(polB_corr);			
			hpolY->Fill(polY_corr);			
		}
		
		for (int i=0; i < nentries; i++ )
		//for (int i=0; i < 10000; i++ )
		{

			ntuple1tpc -> GetEntry(i);

			if(cone>0.7) continue;
			if(Minv>4.) continue;
			if(fitPts_min_pair<15)continue;
			
			if( Minv>=0.4876 && Minv<0.5076 ) continue;//dodge K0 mass range, doesn't cause asymmetry.
			//-------------BLUE-------------------------------------
			//Phi
			for(int phi=0;phi<16;phi++)
			{
				if(PhiRSB>=(phi-8.)/8.*pi && PhiRSB<=(phi-7.)/8.*pi){
					for(int m=0;m<9;m++)
					{
						if(Minv>=M[m] && Minv<M[m+1])
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
									NMinvGtUpB[m*16+phi]++;
								}
								if(eta_pair<0)
								{
									NMinvLtUpB[m*16+phi]++;
								}
							}
							if(fspinconfig==83 || fspinconfig==85)
							{
								if(eta_pair>0)
								{
									NMinvGtDnB[m*16+phi]++;
								}
								if(eta_pair<0)
								{
									NMinvLtDnB[m*16+phi]++;
								}
							}
						}
					}//Minv loop

				}//phi-range loop 
			}//phi loop
			/**************** BLUE BEAM ENDS *****************************////////
			//// YELLOW BEAM ////////////
			//Phi loop 
			for(int phi=0;phi<16;phi++)
			{
				if(PhiRSY>=(phi-8.)/8.*pi && PhiRSY<=(phi-7.)/8.*pi)
				{	
					//Minv loop
					//for(int m=0;m<5;m++)
					for(int m=0;m<9;m++)
					{
						if(Minv>=M[m] && Minv<M[m+1])
						{
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
								//if(eta_pair>0)
								if(eta_pair<0)
								{
									NMinvGtUpY[m*16+phi]++;
								}
								//if(eta_pair<0)
								if(eta_pair>0)
								{
									NMinvLtUpY[m*16+phi]++;
								}
							}
							if(fspinconfig==53 || fspinconfig==85)
							{
								//if(eta_pair>0)
								//eta_pair<0 for yellow is similar to eta_pair>0 for blue
								if(eta_pair<0)//eta_pair direction should be flipped for yellow beam 
								{
									NMinvGtDnY[m*16+phi]++;
								}
								if(eta_pair>0)
								{
									NMinvLtDnY[m*16+phi]++;
								}
							}
						}
					}//Minv loop 

				}//PhiRS range loop
			}//phi loop
		}//entries  for loop 
	 //YELLOW BEAM ENDS

 	//Get average polarizations and rms from histograms 
 		//avgPolB=0.57534; rmsB=0.0370551;
 		//avgPolY=0.58560; rmsY=0.0386434;


		//Get average polarizations and rms from histograms 
		avgPolB = hpolB->GetMean();
		avgPolY = hpolY->GetMean();
		rmsB    = hpolB->GetRMS();
		rmsY    = hpolY->GetRMS();
		
	
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
			
			Output<<"//Hadron Yields  for Blue beam"<<endl;
			Output<<"M BinRange\t"<<"\teta_pair>0\t"<<"\teta_pair<0"<<endl;
		for(int j=0; j<9; j++){
			 cout<<"Bin: "<<j<<" M_gt = "<<avg_Mgt[j]<<endl;
			Output<<M[j]<<"-"<<M[j+1]<<"\t"<<"\t"<<NgtB[j]<<"\t"<<NltB[j]<<endl;
			}
			Output<<"//Hadron Yields  for Yellow beam"<<endl;
			Output<<"M BinRange\t"<<"\teta_pair>0\t"<<"\teta_pair<0"<<endl;
		for(int j=0; j<9; j++){
			Output<<M[j]<<"-"<<M[j+1]<<"\t\t"<<NgtY[j]<<"\t\t"<<NltY[j]<<endl;
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
					B = NMinvGtUpB[m*16+ang];
					C = NMinvGtUpB[m*16+ang+8];
					D = NMinvGtDnB[m*16+ang];
					E = NMinvGtDnB[m*16+ang+8];
					BE = (double)(B*E);
					DC = (double)(D*C);
					a = sqrt(BE);
					b = sqrt(DC);
				}
				if(ang>7)
				{
					B = NMinvGtUpB[m*16+ang];
					C = NMinvGtUpB[m*16+ang-8];
					D = NMinvGtDnB[m*16+ang];
					E = NMinvGtDnB[m*16+ang-8];		
					BE = (double)(B*E);
					DC = (double)(D*C);
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
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			double chiSquareCal[8]={0};
			double diff[8]={0};
			double fitVal[8]={0};
			char name[600];
			char title[600];
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *grt = new TGraphErrors(8,angle,Asym,ex,dA);
			sprintf(title,"Minv bin %i #eta>0, BLUE beam",m);
			grt->SetTitle(title);
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0);
			//fit->SetParameter(0,0.0001);
			grt->Fit(fit,"R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->CenterTitle();
			grt->GetYaxis()->SetRangeUser(-0.09,0.09);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2[m] = fit->GetChisquare();
			cout<<"ChiSquare = "<<chi2Ndf[m]<<endl;

			Abg[m]=fit->GetParameter(0);
			deltaAbg[m]=fit->GetParError(0);

			for(int k=0; k<8; k++){
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
			sprintf(name,"Minv_bin%i_etaGt_Blue.png",m);
			c1->SaveAs( name);
		}//invariant mass loop
		
		//Asymmetry for BLUE , eta_pair<0 (Backward asymmetry)
		for(int m=0;m<9;m++)
		{
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NMinvLtUpB[m*16+ang];
					C = NMinvLtUpB[m*16+ang+8];
					D = NMinvLtDnB[m*16+ang];
					E = NMinvLtDnB[m*16+ang+8];
					BE = (double)(B*E);
					DC = (double)(D*C);
					a=sqrt(BE);
					b=sqrt(DC);
				}
				if(ang>7)
				{
					B = NMinvLtUpB[m*16+ang];
					C = NMinvLtUpB[m*16+ang-8];
					D = NMinvLtDnB[m*16+ang];
					E = NMinvLtDnB[m*16+ang-8];		
					BE = (double)(B*E);
					DC = (double)(D*C);
					a=sqrt(BE);
					b=sqrt(DC);
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
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double chi2[9];
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[16]={0};
			double chiSquareCal[8]={0};
			double diff[8]={0};
			double fitVal[8]={0};
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *grbl = new TGraphErrors(8,angle,Asym,ex,dA);
			grbl->SetMarkerStyle(20);
			sprintf(title,"Minv bin %i #eta<0, BLUE beam",m);
                        grbl->SetTitle(title);
			grbl->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0);
			//fit->SetParameter(0,0.0001);
			grbl->Fit(fit,"R");
			grbl->SetMarkerColor(4);
			grbl->GetXaxis()->SetTitle("#Phi_{RS}");
			grbl->GetXaxis()->SetTitleOffset(1);
			chi2Ndf[m] = fit->GetChisquare()/fit->GetNDF();
			chi2[m] = fit->GetChisquare();
			cout<<"ChiSquareNdf = "<<chi2Ndf[m]<<endl;
			grbl->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grbl->GetYaxis()->CenterTitle();
			grbl->GetYaxis()->SetRangeUser(-0.09,0.09);
		
			Abl[m]=fit->GetParameter(0);
			deltaAbl[m]=fit->GetParError(0);
		
			for(int k=0; k<8; k++){
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
			
			sprintf(name,"Minv_bin%i_etaLt_Blue.png",m);
			c1->SaveAs(name);
		}//invariant mass loop
	
		
		//Asymmetry for YELLOW beam , eta > 0  
		for(int m=0;m<9;m++)
		{	
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NMinvGtUpY[m*16+ang];
					C = NMinvGtUpY[m*16+ang+8];
					D = NMinvGtDnY[m*16+ang];
					E = NMinvGtDnY[m*16+ang+8];
					BE= (double)(B*E); 
					DC= (double)(D*C);
					a = sqrt(BE);
					b = sqrt(DC);
				}
				if(ang>7)
				{
					B = NMinvGtUpY[m*16+ang];
					C = NMinvGtUpY[m*16+ang-8];
					D = NMinvGtDnY[m*16+ang];
					E = NMinvGtDnY[m*16+ang-8];
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
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			double chiSquareCal[8]={0};
			double diff[8]={0};
			double fitVal[8]={0};
			char name[600];
			char title[600];
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *gr = new TGraphErrors(8,angle,Asym,ex,dA);
			sprintf(title,"Minv bin %i #eta>0, YELLOW beam",m);
			gr->SetTitle(title);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
			//fit->SetParameter(0,0.0001);
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
	
			for(int k=0; k<8; k++){
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
	
			sprintf(name,"Minv_bin%i_etaGt_Yellow.png",m);
			c1->Print(name);

		}//Minv loop
		
		//Asymmetry for Yellow beam , eta<0 (backward asymmetry)
		for(int m=0;m<9;m++)
		{	
			for(int ang=0;ang<16;ang++)
			{
				if(ang<8)
				{
					B = NMinvLtUpY[m*16+ang];
					C = NMinvLtUpY[m*16+ang+8];
					D = NMinvLtDnY[m*16+ang];
					E = NMinvLtDnY[m*16+ang+8];
					BE= (double)(B*E); 
					DC= (double)(D*C);
					a = sqrt(BE);
					b = sqrt(DC);
				}
				if(ang>7)
				{
					B = NMinvLtUpY[m*16+ang];
					C = NMinvLtUpY[m*16+ang-8];
					D = NMinvLtDnY[m*16+ang];
					E = NMinvLtDnY[m*16+ang-8];
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
			double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
			double ex[8]={0};
			double chiSquareCal[8]={0};
			double diff[8]={0};
			double fitVal[8]={0};
			char name[600];
			char title[600];
			gStyle->SetOptDate(0);
			gStyle->SetStatX(0.9);
			TGraphErrors *gryl = new TGraphErrors(8,angle,Asym,ex,dA);
			gryl->SetMarkerStyle(20);
			sprintf(title,"Minv bin %i #eta<0, YELLOW beam",m);
			gryl->SetTitle(title);
			gryl->Draw("AP");
			TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0);
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
		
			for(int k=0; k<8; k++){
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


			sprintf(name,"Minv_bin%i_etaLt_yellow.png",m);
			c1->Print(name);

		}//Minv loop
		
			//for Blue beam eta > 0
			double	A_Mbg[9]={Abg[0],Abg[1],Abg[2],Abg[3],Abg[4],Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_Mbg[9]={deltaAbg[0],deltaAbg[1],deltaAbg[2],deltaAbg[3],deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
			//for Blue beam eta < 0
			double	A_Mbl[9]={Abl[0],Abl[1],Abl[2],Abl[3],Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_Mbl[9]={deltaAbl[0],deltaAbl[1],deltaAbl[2],deltaAbl[3],deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
			// for Yellow beam eta >0
			double A_Myg[9]={Ayg[0],Ayg[1],Ayg[2],Ayg[3],Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_Myg[9]={deltaAyg[0],deltaAyg[1],deltaAyg[2],deltaAyg[3],deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};
			
			double A_Myl[9]={Ayl[0],Ayl[1],Ayl[2],Ayl[3],Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_Myl[9]={deltaAyl[0],deltaAyl[1],deltaAyl[2],deltaAyl[3],deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		
	

			//calculate average asymmetry  
			double avgA_Mg[9] ={0};
			double avgA_Ml[9] ={0};
			double WavgA_Mg[9] ={0};
			double WavgA_Ml[9] ={0};
				for (Int_t ii=0; ii<9; ii++)
				{
					avgA_Mg[ii]= (A_Mbg[ii]+A_Myg[ii])/2.;
					avgA_Ml[ii]= (A_Mbl[ii]+A_Myl[ii])/2.;
					WavgA_Mg[ii]= (A_Mbg[ii]*(1/pow(deltaA_Mbg[ii],2))+A_Myg[ii]*(1/pow(deltaA_Myg[ii],2)))/((1/pow(deltaA_Mbg[ii],2))+(1/pow(deltaA_Myg[ii],2)));
					WavgA_Ml[ii]= (A_Mbl[ii]*(1/pow(deltaA_Mbl[ii],2))+A_Myl[ii]*(1/pow(deltaA_Myl[ii],2)))/((1/pow(deltaA_Mbl[ii],2))+(1/pow(deltaA_Myl[ii],2)));
				}	
			//total asymmetry error calculation
			//since two asymmetry(BLUE+YELLOW) are averaged, error propagation formaula is used for total asymmetry error 
			//if function, f = (a+b)/2, Then Error, df =1/2 √{(∂f/∂a)^2*(∆a)^2+(∂f/∂b)^2*(∆b)^2}
			double errA_Mg[9] ={0};
			double errA_Ml[9] ={0};
			double WerrA_Mg[9] ={0};
			double WerrA_Ml[9] ={0};
				for (Int_t iii=0; iii<9; iii++)
				{
					errA_Mg[iii]= .5*sqrt(pow(deltaA_Mbg[iii],2)+pow(deltaA_Myg[iii],2));//total asym err, pT bin 1, eta > 0
					errA_Ml[iii]= .5*sqrt(pow(deltaA_Mbl[iii],2)+pow(deltaA_Myl[iii],2));//total asym err, pT bin 1, eta < 0
					WerrA_Mg[iii]= 1/sqrt((1/pow(deltaA_Mbg[iii],2))+(1/pow(deltaA_Myg[iii],2)));
					WerrA_Ml[iii]= 1/sqrt((1/pow(deltaA_Mbl[iii],2))+(1/pow(deltaA_Myl[iii],2)));
				
				}	
		out << "/////////////////////Average Asymmetry Values : Minv Binning/////////////////////"<< endl;	
                out << "double Avg_MgtB[9] = {" <<avg_MgtB[0]<<", "<<avg_MgtB[1]<<", "<<avg_MgtB[2]<<", "<<avg_MgtB[3]<<", "<<avg_MgtB[4]<<", "<<avg_MgtB[5]<<", "<<avg_MgtB[6]<<", "<<avg_MgtB[7]<<", "<<avg_MgtB[8]<<"}"<< endl;
                out << "double Avg_MltB[9] = {" <<avg_MltB[0]<<", "<<avg_MltB[1]<<", "<<avg_MltB[2]<<", "<<avg_MltB[3]<<", "<<avg_MltB[4]<<", "<<avg_MltB[5]<<", "<<avg_MltB[6]<<", "<<avg_MltB[7]<<", "<<avg_MltB[8]<<"}"<< endl;
                out << "double Avg_pTgtB[9] = {"<<avg_pTgtB[0]<<", "<<avg_pTgtB[1]<<", "<<avg_pTgtB[2]<<", "<<avg_pTgtB[3]<<", "<<avg_pTgtB[4]<<", "<<avg_pTgtB[5]<<", "<<avg_pTgtB[6]<<", "<<avg_pTgtB[7]<<", "<<avg_pTgtB[8]<<"}"<< endl;
                out << "double Avg_pTltB[9] = {"<<avg_pTltB[0]<<", "<<avg_pTltB[1]<<", "<<avg_pTltB[2]<<", "<<avg_pTltB[3]<<", "<<avg_pTltB[4]<<", "<<avg_pTltB[5]<<", "<<avg_pTltB[6]<<", "<<avg_pTltB[7]<<", "<<avg_pTltB[8]<<"}"<< endl;
		out << "double AsymBlueF[9]={"<<Abg[0]<<", "<< Abg[1]<<", "<< Abg[2]<<", "<< Abg[3]<<", "<< Abg[4]<<", "<<Abg[5]<<", "<<Abg[6]<<", "<<Abg[7]<<", "<<Abg[8]<<"};"<<endl; 
		out << "double errBlueF[9]={"<<deltaAbg[0]<<", "<< deltaAbg[1]<<", "<< deltaAbg[2]<<", "<< deltaAbg[3]<<", "<< deltaAbg[4]<<", "<<deltaAbg[5]<<", "<<deltaAbg[6]<<", "<<deltaAbg[7]<<", "<<deltaAbg[8]<<"}"<<endl; 
		out << "double AsymBlueB[9]={"<<Abl[0]<<", "<< Abl[1]<<", "<< Abl[2]<<", "<< Abl[3]<<", "<< Abl[4]<<"}"<<", "<<Abl[5]<<", "<<Abl[6]<<", "<<Abl[7]<<", "<<Abl[8]<<"}"<<endl; 
		out << "double errBlueB[9]={"<<deltaAbl[0]<<", "<< deltaAbl[1]<<", "<< deltaAbl[2]<<", "<< deltaAbl[3]<<", "<< deltaAbl[4]<<", "<<deltaAbl[5]<<", "<<deltaAbl[6]<<", "<<deltaAbl[7]<<", "<<deltaAbl[8]<<"}"<<endl; 
		out << "double AsymYellowF[9]={"<<Ayg[0]<<", "<< Ayg[1]<<", "<< Ayg[2]<<", "<< Ayg[3]<<", "<< Ayg[4]<<", "<<Ayg[5]<<", "<<Ayg[6]<<", "<<Ayg[7]<<", "<<Ayg[8]<<"}"<<endl; 
		out << "double errYellowF[9]={"<<deltaAyg[0]<<", "<< deltaAyg[1]<<", "<< deltaAyg[2]<<", "<< deltaAyg[3]<<", "<< deltaAyg[4]<<", "<<deltaAyg[5]<<", "<<deltaAyg[6]<<", "<<deltaAyg[7]<<", "<<deltaAyg[8]<<"}"<<endl; 
		out << "double AsymYellowB[9]={"<<Ayl[0]<<", "<< Ayl[1]<<", "<< Ayl[2]<<", "<< Ayl[3]<<", "<< Ayl[4]<<", "<<Ayl[5]<<", "<<Ayl[6]<<", "<<Ayl[7]<<", "<<Ayl[8]<<"}"<<endl; 
		out << "double errYellowB[9]={"<<deltaAyl[0]<<", "<< deltaAyl[1]<<", "<< deltaAyl[2]<<", "<< deltaAyl[3]<<", "<< deltaAyl[4]<<", "<<deltaAyl[5]<<", "<<deltaAyl[6]<<", "<<deltaAyl[7]<<", "<<deltaAyl[8]<<"}"<<endl; 
		out << "double avgAsymF[9]={"<<avgA_Mg[0]<<", "<<avgA_Mg[1]<< ","<<avgA_Mg[2]<< ","<<avgA_Mg[3]<<", "<<avgA_Mg[4]<<", "<<avgA_Mg[5]<<", "<<avgA_Mg[6]<<", "<<avgA_Mg[7]<<", "<<avgA_Mg[8]<<"}"<<endl; 
		out << "double avgErrF[9]={"<<errA_Mg[0]<<", "<<errA_Mg[1]<< ","<<errA_Mg[2]<< ","<<errA_Mg[3]<<", "<<errA_Mg[4]<<", "<<errA_Mg[5]<<", "<<errA_Mg[6]<<", "<<errA_Mg[7]<<", "<<errA_Mg[8]<<"}"<<endl; 
		out << "double avgAsymB[9]={"<<avgA_Ml[0]<<", "<<avgA_Ml[1]<< ","<<avgA_Ml[2]<< ","<<avgA_Ml[3]<<", "<<avgA_Ml[4]<<", "<<avgA_Ml[5]<<", "<<avgA_Ml[6]<<", "<<avgA_Ml[7]<<", "<<avgA_Ml[8]<<"}"<<endl; 
		out << "double avgErrB[9]={"<<errA_Ml[0]<<", "<<errA_Ml[1]<< ","<<errA_Ml[2]<< ","<<errA_Ml[3]<<", "<<errA_Ml[4]<<", "<<errA_Ml[5]<<", "<<errA_Ml[6]<<", "<<errA_Ml[7]<<", "<<errA_Ml[8]<<"}"<<endl; 
		out << "/////////////////////////////////////////////////////////////////"<< endl;

double pionPairPtGt[5]={0.995685,0.99413,0.992819,0.981137,0.951615};
double pionPairPtLt[5]={0.997217,0.995865,0.994613,0.982959,0.954391};
double pionPairMGt[5]={0.978901,0.98193,0.981563,0.980261,0.971465};
double pionPairMLt[5]={0.981487,0.984854,0.98391,0.982889,0.973322};
//kaon+proton
double kpPairPtGt[5]={0,0,0,2.08895e-05,0.000131916};
double kpPairPtLt[5]={0,0,0,1.90933e-05,0.000122819};
double kpPairMGt[5]={1.88864e-05,3.60047e-05,3.40103e-05,4.00047e-05,7.8545e-05};
double kpPairMLt[5]={1.49303e-05,2.6513e-05,2.64161e-05,3.2977e-05,7.32954e-05};
//// pair purity in eta bins
double pionPairEta[9]={0.97593,0.98328,0.982627,0.978905,0.975671,0.978282,0.983113,0.983608,0.976228};
double kpPairEta[9]={3.55615e-05,2.98475e-05,3.96205e-05,5.86822e-05,7.06504e-05,5.12718e-05,2.85773e-05,2.04411e-05,2.54495e-05};
double electronPairEta[9]={3.69162e-05,8.40101e-06,5.82411e-06,8.56687e-06,1.40317
//preliminary pion PID fractions
//double pairPurityMGt[9]={0.768775,0.815674,0.79902,0.796715,0.825616,0.805815,0.785091,0.758806,0.718494};
//double pairpurityMLt[9]={0.781753,0.826413,0.81047,0.809187,0.834432,0.81532,0.795783,0.772011,0.732798};
//-----------------------------



//bias ratio
double trigBiasGt[9]={1.07026, 1.01152, 1.03012, 1.00929, 0.987516, 1.14202, 1.05634, 1.24891, 1.08074};
double trigBiasGtErr[9]={0.0522475, 0.0436093, 0.0512523, 0.0500996, 0.0469354, 0.0518734, 0.0542341, 0.0689986, 0.0443317};
double trigBiasLt[9]={1.15514, 1.10721, 0.886907, 0.908908, 1.0177, 1.00483, 1.15098, 1.0318, 1.16321};
double trigBiasLtErr[9]={0.0561346, 0.051315, 0.0472823, 0.0488303, 0.0450848, 0.052088, 0.0680179, 0.0571762, 0.0591917};
 //calculate systematic error
 double pidSysGt[9]={0};
 double biasSysGt[9]={0};
 double combSysGt[9]={0};
 double pidSysLt[9]={0};
 double biasSysLt[9]={0};
 double combSysLt[9]={0};
//print values in table format
  ofstream ftable;
        ftable.open("ifftable_AutVsMinvIntPtTPConly.txt");
        ftable<<"//values are printed in the following order:\n";
        ftable<<"//1.Bin range 2.<p_T> 3.<M> 4.A_UT 5.Sigma_Stat 6.Sigma_PID 7.Sigma_TriggerBias 8.Sigma_Total \n";
        ftable<<"// 9.A_UT 10.Sigma_Stat 11.Sigma_PID 12.Sigma_TriggerBias 13.Sigma_Total "<<endl;
        ftable<<"//(Columns 4-8 for eta>0 and columns 9-13 for eta<0)"<<endl;
        ftable<<"\n"<<endl;




 for(int i=0; i<9; i++){
         pidSysGt[i]=(1-pairPurityMGt[i])*TMath::Max(WavgA_Mg[i],WerrA_Mg[i]);
         biasSysGt[i]=(1-trigBiasGt[i])*TMath::Max(WavgA_Mg[i],WerrA_Mg[i]);
         combSysGt[i]=sqrt(pow(pidSysGt[i],2)+pow(biasSysGt[i],2));
         pidSysLt[i]=(1-pairPurityMLt[i])*TMath::Max(WavgA_Ml[i],WerrA_Ml[i]);
         biasSysLt[i]=(1-trigBiasGt[i])*TMath::Max(WavgA_Ml[i],WerrA_Ml[i]);
         combSysLt[i]=sqrt(pow(pidSysLt[i],2)+pow(biasSysLt[i],2));

 	ftable<<M[i]<<" - "<<M[i+1]<< " & "<<avg_pTgt[i]<<" & "<<avg_Mgt[i]<<" & "<< WavgA_Mg[i]<<" & "<<WerrA_Mg[i]<< " & "<<pidSysGt[i]<<" & "<<biasSysGt[i]<<" & " <<combSysGt[i]<< " & "<<WavgA_Ml[i]<<" & "<<WerrA_Ml[i]<< " & "<<pidSysLt[i]<<" & "<<biasSysLt[i]<<" & " <<combSysLt[i]<<" \\\\ "<<endl;

 }
	double avg_pT7=5.25;
        double run06M[5]={0.36, 0.50, 0.69, 0.88, 1.19};
        double run06asym[5]={0.0054, 0.018, 0.023, 0.070, 0.039};
        double run06err[5]={0.010, 0.0068, 0.0081, 0.013, 0.020};

     //theory from Radici et. al. @√s = 200 GeV

      //Mh     lower    upper
      // 0.36   0.0027   0.0076
      //  0.5    0.0077   0.0132
      //  0.69   0.0173   0.0286
      //  0.88   0.0262   0.0438
      //  1.19   0.0200   0.0360
       
      double thM200[5]={0.36, 0.5,  0.69, 0.88, 1.19};
      double thEl200[5]={0.0027, 0.0077, 0.0173, 0.0262, 0.02};
      double thEh200[5]={0.0076, 0.0132, 0.0286, 0.0438, 0.0360};

     double err[9]={0};
     double errxps[9]={0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015};
     double errxps11[7]={0.015,0.015,0.015,0.015,0.015,0.015,0.015};

      double thA200[5]={0};
      double thE200[5]={0};
      double thEtest[5]={0};
      for(int i=0; i<5; i++){
              thA200[i]=0.5*(thEl200[i]+thEh200[i]);
              thE200[i]=thEh200[i]-thA200[i];
              thEtest[i]=(-1)*thEl200[i]+thA200[i];
              //cout<<thE200[i]<<" == "<<thEtest[i]<<endl;
      }

        TCanvas *myCanF = new TCanvas("myCanF","myCanF",700,500);
        gStyle -> SetOptStat(0);
        gStyle -> SetOptDate(0);
        gStyle -> SetTitle("");
        gStyle->SetLegendBorderSize(0);
        gPad->SetGrid(0,0);
        gPad->SetTopMargin(0.05);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);

        TGraphErrors *th200 = new TGraphErrors(5,thM200,thA200,0,thE200);
        th200->SetTitle("");
        th200->SetFillStyle(1001);
        th200->SetFillColorAlpha(39,0.8);

        th200->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})  }}");
        th200->SetTitle("");
        th200->GetYaxis()->SetLabelSize(0.04);
        th200->GetYaxis()->SetLabelFont(22);
        th200->GetYaxis()->SetTitleSize(0.04);
        th200->GetYaxis()->SetTitleOffset(1.2);
        th200->GetYaxis()->SetRangeUser(-0.01,0.099);
        th200->GetXaxis()->SetLimits(0.21, 1.9);
        th200->GetXaxis()->SetTitle("#font[22]{M^{#pi^{+}#pi^{-}}_{inv}(GeV/c^{2})}");
        th200->GetXaxis()->SetTitleSize(0.04);
        th200-> GetXaxis()->SetLabelSize(0.04);
        th200-> GetXaxis()->SetLabelFont(22);
        th200->GetXaxis()->SetTitleOffset(1.1);
        th200->Draw("AE3");

     TGraphErrors *Run06= new TGraphErrors(5,run06M,run06asym,err,run06err);
     Run06-> SetMarkerStyle(20);
     Run06-> SetMarkerColor(4);
     Run06-> SetLineColor(4);
     Run06-> SetMarkerSize(1.2);
     Run06-> GetXaxis()->SetLimits(0.2,2.5);
     Run06->Draw("same P");

     TGraphErrors *Run15s = new TGraphErrors(9,avg_Mgt,WavgA_Mg,errxps,combSys);
     Run15s->SetFillStyle(0);
     Run15s->SetLineColor(2);
     Run15s->Draw("same 2");

     TGraphErrors *Run15 = new TGraphErrors(9,avg_Mgt,WavgA_Mg,err,WerrA_Mg);
     Run15-> SetMarkerStyle(20);
     Run15-> SetMarkerSize(1.2);
     Run15-> SetMarkerColor(kRed);
     Run15-> SetLineColor(kRed);
     Run15-> GetXaxis()->SetLimits(0.2,2.5);
     Run15->Draw("P same");
     myCanF->Update();

     TLine *line1=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
     line1->SetLineStyle(2);
     line1->Draw();

     TLegend *legend=new TLegend(0.68,0.56,0.93, 0.94);
     legend->AddEntry(th200, "Radici et. al.","f");
     legend->AddEntry(Run15, "Run 15, Cone < 0.7","lp");
     legend->AddEntry("", "#LT p_{T} #GT = 5.25 GeV/c","");
     legend->AddEntry(Run06, "Run 06, Cone < 0.3","lp");
     legend->AddEntry("", "#LT p_{T} #GT = 6 GeV/c","");
     legend->AddEntry(Run15s, " Syst. Error","f");
     legend->AddEntry("", "#eta^{#pi^{+}#pi^{-}} > 0 "," ");
     legend->SetTextSize(0.036);
     legend->SetTextFont(22);
     legend->Draw();

     TLatex latex;
     latex.SetTextFont(22);
     latex.SetTextSize(0.028);
     latex.SetTextAlign(13);//align at top
     latex.SetTextSize(0.055);
     latex.DrawLatex(0.27,0.096,"#color[2]{#font[22]{STAR Preliminary 2015}}");
     latex.SetTextSize(0.04);
     latex.DrawLatex(0.27,0.0895,"#color[1]{#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X} at #sqrt{s} = 200 GeV}");
     latex.SetTextSize(0.03);
     latex.DrawLatex(0.6,-0.003,"#color[1]{#font[22]{#pm 3% scale uncertainty from beam polarization (not shown)}}");
     myCanF->Update();
     myCanF->SaveAs("asymVsMinvIntpTTPConly.pdf");
//------------------------
}//main
