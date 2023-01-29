//////////////////////////////////////////////////////////////
//Edited : Jan 6 , 2020 Babu 				    //	
// This code produce histograms for TSSA as a function of   //
// Eta.							    //
//////////////////////////////////////////////////////////////
#include <iostream> 
//#include "std_lib_facilities.h" 
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace std;

ofstream Output;

//void asymmetryVsEta_9bin(const char *ifile)
void asymmetryVsEta_9bin(const char *ifile="/star/u/pokhrel/GPFS/IFF_TREES/StartLessTOF/iffNtuplesFinalA.root")
{
	gStyle->SetLegendBorderSize(0);

	TH1D *hpolB = new TH1D("hpolB", "", 100, 0, 1);
	TH1D *hpolY = new TH1D("hpolY", "", 100, 0, 1);
	TFile *hist4EtaBins=new TFile("hist4EtaBins.root","RECREATE");
	TH1F *hchi2NdfB=new TH1F("hchi2NdfB","", 100, 0, 50);
	TH1F *hchi2NdfY=new TH1F("hchi2NdfY","", 100, 0, 50);
	//pT bin boundaries (boundaries are set to ensure roughlyequal stat in all bins)
	Double_t pT[10]= {2.800, 3.470,  3.790, 4.120, 4.490, 4.938,  5.505, 6.300, 7.660, 15.00 };
	Double_t M[10]= {0.250, 0.403, 0.516,   0.612,  0.711,  0.803,  0.921,  1.070,  1.286,  4.000}; 
	//eta bin boundaries 
	double eta_range[10]={-1., -0.668, -0.469, -0.281, -0.096, 0.089, 0.275, 0.470, 0.675, 1.};

	//histograms for pT-pair and Minv in different eta-bins(To check the average pT and Minv in each eta bins. Since pT and Minv are integrated over eta bins, these quantities supposed to similar in each eta bins. Carl's suggestion spin PWG-2020/11/04)
	const int nBeam=2;
	const char *beam[2]={"B","Y"};
	const int nBins=9;
	TH1D *pTpair[2][9];
	TH1D *minv[2][9];
	TH1D *heta[2][9];
	
	for(int i=0; i<2; i++){
	  for(int ii=0; ii<9; ii++){
	    pTpair[i][ii]= new TH1D(Form("%s_pTpair_EtaBin%i",beam[i], ii), "", 100, 0,15);
	    minv[i][ii]= new TH1D(Form("%s_Minv_EtaBin%i",beam[i], ii), "", 100, 0,4);
	    heta[i][ii]=new TH1D(Form("%s_EtaBin%i",beam[i], ii),"", 100, eta_range[ii], eta_range[ii+1]);	
	  }
	}
	cout<<"initialization is good..."<<endl;
	TFile *f = new TFile(ifile);
	Output.open("AsymVsEtaTOFonly.txt");

	//get the trees  
	TTree *ntuple1tof = (TTree*)f->Get("ntuple1tof"); 
	TTree *ntuple2tof = (TTree*)f->Get("ntuple2tof");
	TTree *ntuple4tof = (TTree*)f->Get("ntuple4tof");
	TTree *ntuple5tof = (TTree*)f->Get("ntuple5tof");
	TTree *ntuple6 = (TTree*)f->Get("ntuple6");

	//Define variables to hold leaves content 
	float eta_pair;   
	float PhiRS;
	float fspinconfig;
	float cone;
	float pT_pair;
	float Minv;
	float trigger;
	float PhiRSB;
	float PhiRSY;
	float fitPts_min_pair;
	float polB_corr, polY_corr;
	//To store average polarization values from histograms
	double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;
	double pi = 3.14159265359;

	//variables for eta bin average
	 double Npairs[9]={0};
	double pTpairs[9]={0};
	 double etapairs[9]={0};
	 double Mpairs[9]={0};
	 double avg_eta[9]={0};
	 double avg_Minv[9]={0};
	 double avg_pT[9]={0};

	//Get the leaves content and store on variables 
	//structure:: tree -> SetBanchAddress("variable to store leaf content", &Leaf)
	ntuple1tof->SetBranchAddress("fspinconfig",&fspinconfig);
	ntuple1tof->SetBranchAddress("cone",&cone);
	ntuple2tof->SetBranchAddress("Minv",&Minv);
	ntuple2tof->SetBranchAddress("pT_pair",&pT_pair);
	ntuple2tof->SetBranchAddress("eta_pair",&eta_pair);
	ntuple4tof->SetBranchAddress("PhiRSB",&PhiRSB);
	ntuple4tof->SetBranchAddress("PhiRSY",&PhiRSY);
	ntuple5tof->SetBranchAddress("fitPts_min_pair",&fitPts_min_pair);
	ntuple6->SetBranchAddress("polB_corr", &polB_corr);
	ntuple6->SetBranchAddress("polY_corr", &polY_corr);


	//Get number of entries in the tree 
	int  nentries = (int)ntuple1tof->GetEntries();
	cout<<"nentries = "<<nentries<<endl;
	//add friend to the tree 
	//structure :: tree_where_to_add -> AddFriend("tree_name_to_add", "root_file_where_the_tree_t0_add_is"
	//if no file name is given it is assumed that the friend tree is on the same root file 
	ntuple1tof->AddFriend("ntuple2tof");
	ntuple1tof->AddFriend("ntuple4tof");
	ntuple1tof->AddFriend("ntuple5tof");
	int NpTGtUpB[175]={0};
	int NpTLtUpB[175]={0};
	int NpTGtDnB[175]={0};
	int NpTLtDnB[175]={0};
	int NMinvGtUpB[175]={0};
	int NMinvLtUpB[175]={0};
	int NMinvGtDnB[175]={0};
	int NMinvLtDnB[175]={0};
	int NetaUpB[175]={0};
	int NetaDnB[175]={0};
	int NpTGtUpY[175]={0};
	int NpTLtUpY[175]={0};
	int NpTGtDnY[175]={0};
	int NpTLtDnY[175]={0};
	int NMinvGtUpY[175]={0};
	int NMinvLtUpY[175]={0};
	int NMinvGtDnY[175]={0};
	int NMinvLtDnY[175]={0};
	int NetaUpY[175]={0};
	int NetaDnY[175]={0};
	int NpTGtUpT[175]={0};
	int NpTLtUpT[175]={0};
	int NpTGtDnT[175]={0};
	int NpTLtDnT[175]={0};
	int NMinvGtUpT[175]={0};
	int NMinvLtUpT[175]={0};
	int NMinvGtDnT[175]={0};
	int NMinvLtDnT[175]={0};
	int NetaUpT[175]={0};
	int NetaDnT[175]={0};


	int ne = (int )ntuple6->GetEntries(); //entries of trees 
/*	for(int p=0; p<ne; p++)
	{
		ntuple6->GetEntry(p);
		hpolB->Fill(polB_corr);
		hpolY->Fill(polY_corr);
	}
*/	
	 //this doubles the code running time. Make use of direct values. 

	for(int j=0;j<nentries;j++)
	//for(int j=0;j<10000;j++)
	{
		ntuple1tof->GetEntry(j);
		if(cone>0.7)continue;
		if(Minv>4.)continue;
		if(fitPts_min_pair<15)continue;
		if(Minv>0.4876 && Minv<0.5076) continue;//dodge K0 mass range, doesn't cause asymmetry.
		double total_pT = 0;
		double total_N = 0;
		total_pT = total_pT + pT_pair;
		total_N = total_N +1;

		//-------------BLUE-------------------------------------
		//Phi
		for(int phi=0;phi<16;phi++)
		{
			if(PhiRSB>=(phi-8.)/8.*pi && PhiRSB<=(phi-7.)/8.*pi)
			{
				//Eta
				for(int eta=0;eta<9;eta++)
				{
					if(eta_pair>=eta_range[eta] && eta_pair<eta_range[eta+1])
					{
						Npairs[eta]=Npairs[eta]+1;
						pTpairs[eta]=pTpairs[eta]+pT_pair;
						etapairs[eta]=etapairs[eta]+eta_pair;
						Mpairs[eta]=Mpairs[eta]+Minv;
						//fill histograms for pT-pair and Minv 
						pTpair[0][eta]->Fill(pT_pair);
						minv[0][eta]->Fill(Minv);
						heta[0][eta]->Fill(eta_pair);
						if(fspinconfig==51 || fspinconfig==53)  
						{
							NetaUpB[eta*16+phi]++;   
						}
						if(fspinconfig==83 || fspinconfig==85)  
						{
							NetaDnB[eta*16+phi]++;
						}

					}
				}//Eta loop  
			}//Phi loop
		}


		//--------------end BLUE---------------------------------------

		//-----------------YELLOW---------------------------------------
		//Phi
		for(int phi=0;phi<16;phi++)
		{
			if(PhiRSY>=(phi-8.)/8.*pi && PhiRSY<=(phi-7.)/8.*pi) 
			{
				//Eta
				for(int eta=0;eta<9;eta++)
				{
						double etaPair=(-1)*eta_pair;
						//cout<<"etaPair: "<<eta_pair<<", "<<etaPair<<endl;
						if(etaPair>=eta_range[eta] && etaPair<eta_range[eta+1]){
						  //cout<<etaPair<<" "<<pT_pair<<" "<<Minv<<endl;
						  pTpair[1][eta]->Fill(pT_pair);
						  minv[1][eta]->Fill(Minv);
						heta[1][eta]->Fill(etaPair);
						}
					if(eta_pair>=eta_range[eta] && eta_pair<eta_range[eta+1])
					{
						if(fspinconfig==51 || fspinconfig==83)
						{
							if(eta==0)NetaUpY[(eta+8)*16+phi]++;
							if(eta==1)NetaUpY[(eta+6)*16+phi]++;
							if(eta==2)NetaUpY[(eta+4)*16+phi]++;
							if(eta==3)NetaUpY[(eta+2)*16+phi]++;
							if(eta==4)NetaUpY[(eta)*16+phi]++;
							if(eta==5)NetaUpY[(eta-2)*16+phi]++;
							if(eta==6)NetaUpY[(eta-4)*16+phi]++;
							if(eta==7)NetaUpY[(eta-6)*16+phi]++;
							if(eta==8)NetaUpY[(eta-8)*16+phi]++;
						}
						if(fspinconfig==53 || fspinconfig==85)
						{
							if(eta==0)NetaDnY[(eta+8)*16+phi]++;
							if(eta==1)NetaDnY[(eta+6)*16+phi]++;
							if(eta==2)NetaDnY[(eta+4)*16+phi]++;
							if(eta==3)NetaDnY[(eta+2)*16+phi]++;
							if(eta==4)NetaDnY[(eta)*16+phi]++;
							if(eta==5)NetaDnY[(eta-2)*16+phi]++;
							if(eta==6)NetaDnY[(eta-4)*16+phi]++;
							if(eta==7)NetaDnY[(eta-6)*16+phi]++;
							if(eta==8)NetaDnY[(eta-8)*16+phi]++;
						}

					} pTpair[0][eta]->Fill(pT_pair);
                                                minv[0][eta]->Fill(Minv);
				}//Eta loop
			}
		}//phi loop
		//---------------end YELLOW--------------------------------------


	}

Output<<"*******************  A_{UT vs #eta , cone < 0.7 , trigger: 1, 2 *************}"<<endl;
	for(int k=0;k<9;k++)
	{	Output<<"Eta Bin: "<<k<<endl;
		avg_eta[k] = (double)etapairs[k]/(double)Npairs[k];		
		avg_Minv[k]=(double)Mpairs[k]/(double)Npairs[k];
		avg_pT[k]  = (double)pTpairs[k]/(double)Npairs[k];	
		Output<<"no. of pairs["<<k<<"]="<<Npairs[k]<<endl;
		Output<<"pT of pairs["<<k<<"]="<<pTpairs[k]<<endl;
		Output<<"eta of pairs["<<k<<"]="<<etapairs[k]<<endl;
		Output<<"Minv of pairs["<<k<<"]="<<Mpairs[k]<<endl;
		Output<<"<pT>["<<k<<"]="<<pTpairs[k]/Npairs[k]<<endl;
		Output<<"<eta>["<<k<<"]="<<etapairs[k]/Npairs[k]<<endl;
		Output<<"<Minv>["<<k<<"]="<<Mpairs[k]/Npairs[k]<<endl;
	//outfile<<"avg_eta["<<k<<"] = "<<avg_eta[k]<< ",  avg_Minv["<<k<<"] = "<<avg_Minv[k]<< ", avg_pT["<<k<<"] = "<<avg_pT[k]<<endl;


	}



	//-----------------end total counts----------------------------- 

	
/*	//Get polarization mean and rms from histograms 
	avgPolB = hpolB->GetMean();
	avgPolY = hpolY->GetMean();
	rmsB    = hpolB->GetRMS();
	rmsY    = hpolY->GetRMS();
*/	
	//polarization values (using direct values runs the code fast.)
	// Blue Beam:    Avg Pol: 0.5753,  RMS: 0.03706
	// Yellow Beam:  Avg Pol: 0.5856,  RMS: 0.03864
	avgPolB = 0.5753;
	avgPolY = 0.5856;
	rmsB    = 0.03706;
	rmsY    = 0.03864;
	
	
	double a,b;
	double dAdB,dAdC,dAdD,dAdE,dAdP;
	double dP_B = rmsB;//use rms from total distribution for now. Should do independent calculatioon for blue and yellow then averaged.
	double dP_Y = rmsY;
	double dA[175],Asym[175];
	double B,C,D,E;
	double pi = 3.14159265359;
	double BE, DC;
	TCanvas *c1=new TCanvas("c1","c1",600,600);
	c1->Print("FitPlots_etaTOFonly.pdf(");
	//asymmetry for BLUE eta > 0
	for(int eta=0;eta<9;eta++)
	{
		for(int ang=0;ang<16;ang++)
		{
			if(ang<8)
			{  
				B = NetaUpB[eta*16+ang];
				C = NetaUpB[eta*16+ang+8];
				D = NetaDnB[eta*16+ang];
				E = NetaDnB[eta*16+ang+8];
				BE = (double)(B*E);
				DC = (double)(D*C);
				a = sqrt(BE);
				b = sqrt(DC);
			}
			if(ang>7)
			{
				B = NetaUpB[eta*16+ang];
				C = NetaUpB[eta*16+ang-8];
				D = NetaDnB[eta*16+ang];
				E = NetaDnB[eta*16+ang-8];
				BE = (double)(B*E);
				DC = (double)(D*C);
				a = sqrt(BE);
				b = sqrt(DC);
			}
			Asym[ang]=(1./avgPolB)*((a-b)/(a+b));
			dAdB = (1./avgPolB)*(E*sqrt(DC))/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));     
			dAdE = (1./avgPolB)*(B*sqrt(DC))/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
			dAdD = (-1./avgPolB)*(C*sqrt(BE))/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
			dAdC = (-1./avgPolB)*(D*sqrt(BE))/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
			dAdP = (-1./(avgPolB*avgPolB))*(sqrt(BE)-sqrt(DC))/(sqrt(BE)+sqrt(DC));
			dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_B)**2); 

		//cout << "***************** eta bin: " << eta << ", phi bin: "<< ang << "  **********************"<< endl; 
		//cout <<"B: "<< B << ", C: "<< C << ", D: "<< D << ", E: "<< E << ", a: "<< a << ", b: "<< b <<", BE: "<<BE << ", DC:  "<<DC<<endl;
		//cout << "dAdB: "<< dAdB <<  ", dAdE: "<< dAdE << ", dAdD: "<< dAdD << ", dAdC: "<< dAdC << ", dAdP:  "<< dAdP << ", Asym["<<ang << "]: "<<Asym[ang]<< ", dA["<<ang<<"]: "<< dA[ang]<<", avgPolB: "<< avgPolB << endl;
	
		}//angle loop

		double AB[9];
		double deltaAB[9];
		char name[600];
		char title[600];
		double chi2Ndf[9]; 
		double chi2[9]; 
		double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};

		double ex[8]={0};
		gStyle->SetOptDate(0);
		grt = new TGraphErrors(8,angle,Asym,ex,dA);
		grt->SetMarkerStyle(20);
		sprintf(title,"eta bin %i, BLUE", eta);
		grt->SetTitle(title);
		grt->Draw("AP");
		TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
		fit->SetParameter(0,0.0001);
		grt->Fit(fit,"R");
		grt->SetMarkerColor(4);
		grt->GetXaxis()->SetTitle("#Phi_{RS}");
		grt->GetYaxis()->SetRangeUser(-0.06, 0.06);
		grt->GetXaxis()->SetTitleOffset(1);
		chi2Ndf[eta] = fit->GetChisquare()/fit->GetNDF();
		chi2[eta] = fit->GetChisquare();
		 hchi2NdfB->Fill(chi2Ndf[eta]);

		cout<<"chisquare: "<<chi2Ndf[eta]<<endl;
		grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
		grt->GetYaxis()->CenterTitle(kTRUE);
		AB[eta]=fit->GetParameter(0);
		deltaAB[eta]=fit->GetParError(0);

	c1->Print("FitPlots_etaTOFonly.pdf");
		
	}//eta loop

	// asymmetry for yellow beam 
	for(int eta=0;eta<9;eta++)
	{
		for(int ang=0;ang<16;ang++)
		{
			if(ang<8)
			{
				B = NetaUpY[eta*16+ang];
				C = NetaUpY[eta*16+ang+8];
				D = NetaDnY[eta*16+ang];
				E = NetaDnY[eta*16+ang+8];
				BE =(double)(B*E);
				DC =(double)(D*C);
				a = sqrt(BE);
				b = sqrt(DC);

			}
			if(ang>7)
			{
				B = NetaUpY[eta*16+ang];
				C = NetaUpY[eta*16+ang-8];
				D = NetaDnY[eta*16+ang];
				E = NetaDnY[eta*16+ang-8];
				BE =(double)(B*E);
				DC = (double)(D*C);
				a = sqrt(BE);
				b = sqrt(DC);
			}
			Asym[ang]=(1./avgPolB)*((a-b)/(a+b));
			dAdB = (1./avgPolB)*(E*sqrt(DC))/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));     
			dAdE = (1./avgPolB)*(B*sqrt(DC))/(sqrt(BE)*((sqrt(BE)+sqrt(DC))**2));
			dAdD = (-1./avgPolB)*(C*sqrt(BE))/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
			dAdC = (-1./avgPolB)*(D*sqrt(BE))/(sqrt(DC)*((sqrt(BE)+sqrt(DC))**2));
			dAdP = (-1./(avgPolB*avgPolB))*(sqrt(BE)-sqrt(DC))/(sqrt(BE)+sqrt(DC));
			dA[ang] = sqrt((fabs(dAdB)*sqrt(B))**2 + (fabs(dAdC)*sqrt(C))**2 + (fabs(dAdD)*sqrt(D))**2 + (fabs(dAdE)*sqrt(E))**2 + (fabs(dAdP)*dP_Y)**2); 

		//cout << "***************** eta bin: " << eta << ", phi bin: "<< ang << " Yellow Beam **********************"<< endl; 
		//cout <<"B: "<< B << ", C: "<< C << ", D: "<< D << ", E: "<< E << ", a: "<< a << ", b: "<< b <<", BE: "<<BE << ", DC:  "<<DC<<endl;
		//cout << "dAdB: "<< dAdB <<  ", dAdE: "<< dAdE << ", dAdD: "<< dAdD << ", dAdC: "<< dAdC << ", dAdP:  "<< dAdP << ", Asym["<<ang << "]: "<<Asym[ang]<< ", dA["<<ang<<"]: "<< dA[ang]<<", avgPolB: "<< avgPolB << endl;
	
		}//angle loop
		double AY[9];
		double deltaAY[9];
		char name[600];
		char title[600];
		double chi2Ndf[9]; 
		double chi2[9]; 
		double angle[8]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
		double ex[8]={0};
		gStyle->SetOptDate(0);
		grt = new TGraphErrors(8,angle,Asym,ex,dA);
		sprintf(title,"eta bin %i, YELLOW",eta);
		grt->SetTitle(title);
		grt->SetMarkerStyle(20);
		grt->Draw("AP");
		TF1 *fit = new TF1("fit","[0]*sin(x)",-3.14159265359,0.);
		fit->SetParameter(0,0.0001);
		grt->Fit(fit,"R");
		grt->SetMarkerColor(4);
		grt->GetXaxis()->SetTitle("#Phi_{RS}");
		grt->GetXaxis()->SetTitleOffset(1);
		grt->GetYaxis()->SetRangeUser(-0.06, 0.06);
		chi2Ndf[eta] = fit->GetChisquare()/fit->GetNDF();
		chi2[eta] = fit->GetChisquare();
		 hchi2NdfY->Fill(chi2Ndf[eta]);
		cout<<"chisquare: "<<chi2Ndf[eta]<<endl;
		grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
		grt->GetYaxis()->CenterTitle(kTRUE);
		AY[eta]=fit->GetParameter(0);
		deltaAY[eta]=fit->GetParError(0);
	c1->Print("FitPlots_etaTOFonly.pdf");

	}//eta loop
	c1->Print("FitPlots_etaTOFonly.pdf)");
	//average x and z values in 9 eta bins from MC simulation.


	double A_B[9]={AB[0],AB[1],AB[2],AB[3],AB[4], AB[5],AB[6], AB[7], AB[8]};
	double deltaA_B[9]={deltaAB[0],deltaAB[1],deltaAB[2],deltaAB[3],deltaAB[4], deltaAB[5], deltaAB[6], deltaAB[7], deltaAB[8]};
	double A_Y[9]={AY[0],AY[1],AY[2],AY[3],AY[4], AY[5],AY[6], AY[7], AY[8]};
	double deltaA_Y[9]={deltaAY[0],deltaAY[1],deltaAY[2],deltaAY[3],deltaAY[4], deltaAY[5], deltaAY[6], deltaAY[7], deltaAY[8]};
	double avgA[9]={0};  double WavgA[9]={0};
	double errA[9]={0};  double WerrA[9]={0};
	//calculate avaerage asymmetry and total error
	for (Int_t l = 0; l<9; l++)
	{
		avgA[l] = (A_B[l]+A_Y[l])/2.; 
		errA[l] = .5*sqrt(pow(deltaA_B[l],2)+pow(deltaA_Y[l],2));

		 WavgA[l]= (A_B[l]*(1/pow(deltaA_B[l],2))+A_Y[l]*(1/pow(deltaA_Y[l],2)))/((1/pow(deltaA_B[l],2))+(1/pow(deltaA_Y[l],2)));
		 WerrA[l]=1/sqrt((1/pow(deltaA_B[l],2))+(1/pow(deltaA_Y[l],2)));

		//cout << "avgA["<<l<<"]="<< avgA[l]<< endl;
		//cout << "errA["<<l<<"]="<< errA[l]<< endl;
	}
ofstream aout;
        aout.open("AutVsEtaIntTofOnly.txt");
        for(int i=0; i<9; i++){
                aout<<WavgA[i]<<" "<<WerrA[i]<<endl;
        }



	Output<<"********** BLUE Asymmetry **************"<<endl;
	Output<<"double A_B[9]={"<<AB[0]<<","<<AB[1]<<","<<AB[2]<<","<<AB[3]<<","<<AB[4]<<","<< AB[5]<<","<<AB[6]<<","<< AB[7]<<","<< AB[8]<<"};"<<endl;
	Output<<"double deltaA_B[9]={"<<deltaAB[0]<<","<<deltaAB[1]<<","<<deltaAB[2]<<","<<deltaAB[3]<<","<<deltaAB[4]<<","<< deltaAB[5]<<","<<deltaAB[6]<<","<< deltaAB[7]<<","<< deltaAB[8]<<"};"<<endl;
	Output<<"****************************************"<<endl;
	
	Output<<"********** YELLOW Asymmetry **************"<<endl;
	Output<<"double A_Y[9]={"<<AY[0]<<","<<AY[1]<<","<<AY[2]<<","<<AY[3]<<","<<AY[4]<<","<< AY[5]<<","<<AY[6]<<","<< AY[7]<<","<< AY[8]<<"};"<<endl;
	Output<<"double delta A_Y[9]={"<<deltaAY[0]<<","<<deltaAY[1]<<","<<deltaAY[2]<<","<<deltaAY[3]<<","<<deltaAY[4]<<","<< deltaAY[5]<<","<<deltaAY[6]<<","<< deltaAY[7]<<","<< deltaAY[8]<<"};"<<endl;
	Output<<"****************************************"<<endl;
	
	Output << "*******Average Asymmetry  Values **************"<< endl;
                Output << "avgEta[9]={"<<avg_eta[0]<<", "<<avg_eta[1]<< ","<<avg_eta[2]<< ","<<avg_eta[3]<<","<<avg_eta[4]<<avg_eta[5]<<", "<<avg_eta[6]<<","<<avg_eta[7]<<", "<<avg_eta[8]<<"};"<<endl;
                Output << "avgAsymEta[9]={"<<avgA[0]<<", "<<avgA[1]<< ","<<avgA[2]<< ","<<avgA[3]<<","<<avgA[4]<<", "<< avgA[5]<<", "<<avgA[6]<<", "<<avgA[7]<<", "<<avgA[8]<<"};"<<endl;
                Output << "avgErrEta[9]={"<<errA[0]<<", "<<errA[1]<< ","<<errA[2]<< ","<<errA[3]<<","<<errA[4]<<", "<<errA[5]<<", "<<errA[6]<<", "<<errA[7]<< ", "<<errA[8]<<"};"<<endl;
                Output << "*************************************************"<< endl;
	Output << "*******Weighted Average Asymmetry  Values **************"<< endl;
                Output << "WavgAsymEta[9]={"<<WavgA[0]<<", "<<WavgA[1]<< ","<<WavgA[2]<< ","<<WavgA[3]<<","<<WavgA[4]<<", "<< WavgA[5]<<", "<<WavgA[6]<<", "<<WavgA[7]<<", "<<WavgA[8]<<"};"<<endl;
                Output << "avgErrEta[9]={"<<WerrA[0]<<", "<<WerrA[1]<< ","<<WerrA[2]<< ","<<WerrA[3]<<","<<WerrA[4]<<", "<<WerrA[5]<<", "<<WerrA[6]<<", "<<WerrA[7]<< ", "<<WerrA[8]<<"};"<<endl;



	// Run2006 asymmetry values for comparison
	double asym06[4]={-0.00020, .0075, .013, .031};
        double err06[4]={.0060, .0061, .0061, .0060};	
        double eta06[4]={-0.75, -0.26, 0.26, 0.75};	
	//money plot -----------------//
	//trigger bias 
	double ratioEta[9]={1.184,1.184,1.184,1.184,1.184,1.184,1.184,1.184,1.184};
 	double ratioEtaErr[9]={0.0626708,0.0531546,0.0418376,0.0523667,0.0605381,0.0520935,0.0694973,0.0488803,0.041424};
	//average x and x values in eta bins
	double avg_x[9]={0.134072,0.117516,0.116531,0.126561,0.142683,0.157495,0.170704,0.179955,0.219244};
	double avg_xerr[9]={0.00361079,0.00585883,0.00509625,0.00289428,0.00393939,0.00262938,0.00264746,0.00321909,0.00204874};
	double avg_z[9]={0.491561,0.498538,0.494753,0.469221,0.450428,0.466991,0.481748,0.499678,0.480493};
	double avg_zerr[9]={0.00540925,0.00908226,0.00869527,0.00471911,0.00423186,0.00473377,0.00446456,0.00558819,0.00330105};
	/*//preliminary pion fractions
	double pairpurityPtGt[5]={0.67307,0.817963,0.794931,0.765595,0.793264};
	double pairpurityPtLt[5]={0.660449,0.8359,0.808052,0.739741,0.80209};
	double pairpurityMGt[5]={0.524506,0.767845,0.81264,0.827867,0.842939};
	double pairpurityMLt[5]={0.528519,0.779716,0.83881,0.83958,0.848677};
	double pairpurityEta[9]={0.777056,0.798529,0.800554,0.803539,0.802013,0.798089,0.79153,0.784865,0.760883};	
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
     
      
	//systematics
	double pidEta[9]={0};
	double biasEta[9]={0};
	double combSys[9]={0};
	for (int pbin=0; pbin<9; pbin++){
   	biasEta[pbin]=(1-ratioEta[pbin])*TMath::Max(WavgA[pbin],WerrA[pbin]);
   	pidEta[pbin]=(1-pairpurityEta[pbin])*TMath::Max(WavgA[pbin],WerrA[pbin]);
   	combSys[pbin]=sqrt(pow(biasEta[pbin],2) + pow(pidEta[pbin],2));
 	}
	        //print values in tabular format
	       ofstream ftable;
	       ftable.open("ifftable_AutVsEtaTOFonly.txt");
	       ftable<<"//values are printed in the following order:\n";
	       ftable<<"//1.eta Bin range 2.<eta> 3.<pt> 4.<M> 5.<x> 6.<z> 7.A_UT 8.Sigma_AUT 9.Sigma_PID 10.Sigma_TriggerBias 11.Sigma_Total "<<endl;
	       ftable<<"\n"<<endl;
	       for(int i=0; i<9; i++){ 
	       ftable<<eta_range[i]<<" - "<<eta_range[i+1]<<" & "<<avg_eta[i]<<" & "<<avg_pT[i]<<" & "<<avg_Minv[i]<<" & "<<avg_x[i]<<" & "<<avg_z[i]<< " & "<<WavgA[i]<<" & "<< WerrA[i]<<" & "<<pidEta[i]<<" & "<<biasEta[i]<<" & "<< combSys[i]<<" \\\\ "<<endl; 
	       	}


 
	TGraphErrors *syst;
	TGraphErrors *data;
	TGraphErrors *gX;
	TGraphErrors *gZ;
	const int n=9;	
	double errxs[9]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
	TCanvas *myCanA = new TCanvas("myCanA","myCanA",600,500);
 	gStyle -> SetOptStat(0);
 	myCanA -> Divide(2,1);
	for(int pad=0; pad<2; pad++){
		myCanA->cd(pad+1);
		if(pad==0){
			setPad(pad);
			syst = new TGraphErrors(n,avg_eta,WavgA,errxs,combSys);//Forward average
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


			data = new TGraphErrors(n, avg_eta, WavgA, 0, WerrA);//Forward average
			data-> SetMarkerStyle(20);
			data-> SetMarkerColor(2);
			data-> SetLineColor(2);
			data->Draw("P same ");


			myCanA->Update();
			TLine *lineA=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
			lineA->SetLineStyle(2);
			lineA->Draw();

			TLegend *leg = new TLegend(0.18,0.64, 0.38, 0.71);
			leg->AddEntry(syst," #font[22]{Syst. Error}","f");
			leg->SetTextSize(.06);
			leg->Draw();

			TLatex text;
			text.SetTextSize(0.07);
			text.SetTextFont(22);
			text.DrawLatex(-0.905, 0.032,"#font[22]{#color[2]{STAR Preliminary 2015}}");
			text.SetTextSize(0.06);
			text.DrawLatex(-0.905, 0.0285, "#font[22]{ p^{#uparrow}+ p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 200 GeV}");
			text.DrawLatex(-0.77, -0.0022, "#scale[0.8]{#pm 3% scale uncertainty from beam polarization (not shown)}");

		}else if(pad==1){
		setPad(pad);
		gZ = new TGraphErrors(n,avg_eta,avg_z,0,avg_zerr);
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
		gZ->GetXaxis()->SetTitleFont(22);
		gZ->GetXaxis()->SetLabelFont(22);
		gZ->GetXaxis()->SetLabelSize(0.12);
		gZ->GetXaxis()->SetTitleSize(0.15);
		gZ->GetYaxis()->SetLabelSize(0.13);
		gZ->GetYaxis()->SetLabelFont(22);
		gZ->Draw("AP");

		gPad->Update();
		gX = new TGraphErrors(n,avg_eta,avg_x,0,avg_xerr);
		gX->SetMarkerStyle(20);
		gX->SetMarkerColor(4);
		gX->Draw("P same");

		TLegend *leg1=new TLegend(0.3,0.53,0.5,0.65);
		leg1->SetNColumns(2);
		leg1->AddEntry(gZ,"#font[22]{#color[2]{< z >}}", "lp");
		leg1->AddEntry(gX,"#font[22]{#color[4]{< x >}}", "lp");
		leg1->SetTextSize(0.12);
		leg1->Draw();
		gPad->Update();
		}
	

	}
	myCanA->SaveAs("asymVsEtaTOFcutTOFonly.pdf");	
	//----------------------------//
	/*		
	gStyle->SetOptDate(0);
	gROOT->ForceStyle(0);
//Blue and yelllow asymmetry 
	TCanvas *myCan = new TCanvas("myCan", "myCan",500,450);
	gStyle -> SetTitleX(.35);
	gStyle->SetLegendBorderSize(0);
	myCan -> cd();
	myCan->SetGrid(0,0);
	gPad->SetLeftMargin(.12);

	gr1 = new TGraphErrors(9,avg_eta,A_B,0,deltaA_B);
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerSize(1);
	gr1->SetTitle("A_{UT} vs #eta");
	gr1->GetXaxis()->SetRangeUser(-1.,1.);
	gr1->GetXaxis()->SetLabelSize(0.05);
	gr1->GetYaxis()->SetRangeUser(-.01,.04);
	gr1->SetMarkerColor(kBlue);
	gr1->GetYaxis()-> SetTitle("A_{UT}");
	gr1->GetXaxis()-> SetTitle("#eta^{#pi^{+}#pi^{-}}");
	gr1->GetXaxis()-> CenterTitle();
	gr1->GetYaxis()-> CenterTitle();
	gr1->GetYaxis()->SetTitleOffset(1.43);
	gr1->GetXaxis()->SetTitleOffset(1.);
	gr1->GetYaxis()->SetLabelSize(0.035);
	gr1->GetXaxis()->SetLabelSize(0.035);
	gr1->GetYaxis()->SetTitleSize(0.04);
	gr1->GetXaxis()->SetTitleSize(0.04);
	gr1-> SetMarkerStyle(20);
	gr1->Draw("AP");
	myCan->Update();
	gr2 = new TGraphErrors(9,avg_eta,A_Y,0,deltaA_Y);
	gr2-> SetMarkerStyle(20);
	gr2-> SetMarkerColor(kRed);
	gr2->Draw("same P");	
	myCan->Update();	
	TLine *line1=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
	line1->SetLineStyle(2);
	line1->Draw();
	TLegend *leg1 = new TLegend(0.2,.68, 0.5, 0.88);
	leg1->AddEntry("", " Cone < 0.7 ", "");
	leg1->AddEntry(gr1, " Blue Beam  ", "lp");
	leg1->AddEntry(gr2, " Yellow Beam  ", "lp");
	leg1->SetTextSize(0.04);
	leg1->Draw();
	myCan->Update();
//average of yellow and blue asymmetry
	TCanvas *avg = new TCanvas("avg", "avg",500,450);
        gStyle -> SetTitle("");
	gStyle->SetTitleX(.4);
        gStyle->SetLegendBorderSize(0);
        avg -> cd();
        avg->SetGrid(0,0);
	gPad->SetLeftMargin(.12);
	grA = new TGraphErrors(9,avg_eta,avgA,0,errA);
	grA->SetMarkerStyle(20);
	grA->SetMarkerColor(2);
	grA->SetTitle(" < A_{UT} > vs #eta");
	grA->GetXaxis()->SetRangeUser(-1.,1.);
	grA->GetXaxis()->SetLabelSize(0.05);
	grA->GetYaxis()->SetRangeUser(-.01,.04);
	grA->SetMarkerColor(kRed);
	grA->GetYaxis()-> SetTitle("< A_{UT} >");
	grA->GetXaxis()-> SetTitle("#eta^{#pi^{+}#pi^{-}}");
	grA->GetXaxis()-> CenterTitle();
	grA->GetYaxis()-> CenterTitle();
	grA->GetYaxis()->SetTitleOffset(1.43);
	grA->GetXaxis()->SetTitleOffset(1.);
	grA->GetYaxis()->SetLabelSize(0.035);
	grA->GetXaxis()->SetLabelSize(0.035);
	grA->GetYaxis()->SetTitleSize(0.04);
	grA->GetXaxis()->SetTitleSize(0.04);
	grA->Draw("AP");
avg->Update();
 	avg->Update();
        TLine *avgL=  new TLine(gPad->GetUxmin(),0.,gPad->GetUxmax(),0.);
        avgL->SetLineStyle(2);
        avgL->Draw();
	TLegend *leg2 = new TLegend(0.2,.68, 0.5, 0.88);
	leg2->AddEntry("", "Run15,p#uparrow + p, #sqrt{s}= 200 GeV ", "");
	leg2->AddEntry("", " Cone < 0.7 ", "");
	leg2->AddEntry(grA, " < A_{UT} >  ", "lp");
	leg2->SetTextSize(0.04);
	leg2->Draw();
avg->Update();
avg->SaveAs("AsymVsEta_ConeLt7_9bin_AverageTOFonly.pdf");

	myCan->SaveAs("AsymVsEta_ConeLt7_9binTOFonly.pdf");
*/
hist4EtaBins->Write();
hist4EtaBins->Close();
}

void setPad(int i){
	if(i==0){
		gPad->SetPad(0.0,0.35,1.,1.);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.00);
		gPad->SetGrid(0,0);

		gPad->Update();
	}else if(i==1){
		gPad->SetPad(0.0,0.0,1,0.35);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.3);
		gPad->SetTopMargin(0.0);
		gPad->SetGrid(0,0);

		gPad->Update();
	}
}

