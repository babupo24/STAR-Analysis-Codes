// This code is modefied to calculate asymmetry for individual beams and averaged for total asymmetry
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace std;
// function to set pad for plotting
void setPad(int i);
ofstream Output;
double calErr(double aTpc, double eTpc, double aTof, double eTof);
double relDiff(double aTpc, double aTof);
double *getTofAsym(TFile *file);
double *getTofErr(TFile *file);

void asymmetryVsinvMass_averaged(const char *ifile = "/star/u/pokhrel/GPFS/IFF_TREES/StartLessTOF/iffNtuplesFinalA.root")
{
	// Histograms to get mean polarization for yellow and blue beam
	TH1D *hpolB = new TH1D("hpolB", "", 100, 0, 1);
	TH1D *hpolY = new TH1D("hpolY", "", 100, 0, 1);

	TFile *chi2f = new TFile("chi2.root", "RECREATE");
	TH1F *chi2NdfBGt = new TH1F("chi2NdfBGt", "Blue eta > 0", 100, 0, 50);
	TH1F *chi2NdfBLt = new TH1F("chi2NdfBLt", "Blue eta < 0", 100, 0, 50);
	TH1F *chi2NdfYGt = new TH1F("chi2NdfYGt", "Yellow eta > 0", 100, 0, 50);
	TH1F *chi2NdfYLt = new TH1F("chi2NdfYLt", "Yellow eta < 0", 100, 0, 50);

	// pT distribution for each pT bin for average pT
	TH1D *hpT_pair[5];
	TH1D *hpT_pairGtB[5];
	TH1D *hpT_pairLtB[5];
	TH1D *hpT_pairGtY[5];
	TH1D *hpT_pairLtY[5];
	for (Int_t n = 0; n < 5; n++)
	{
		hpT_pair[n] = new TH1D(Form("pT_pair_pbin_%i", n), "", 100, 0, 20);
	}
	for (Int_t n = 0; n < 5; n++)
	{
		hpT_pairGtB[n] = new TH1D(Form("pT_pairGtB_pbin_%i", n), "", 100, 0, 20);
	}
	for (Int_t n = 0; n < 5; n++)
	{
		hpT_pairLtB[n] = new TH1D(Form("pT_pairLtB_pbin_%i", n), "", 100, 0, 20);
	}
	for (Int_t n = 0; n < 5; n++)
	{
		hpT_pairGtY[n] = new TH1D(Form("pT_pairGtY_pbin_%i", n), "", 100, 0, 20);
	}
	for (Int_t n = 0; n < 5; n++)
	{
		hpT_pairLtY[n] = new TH1D(Form("pT_pairLtY_pbin_%i", n), "", 100, 0, 20);
	}

	TFile *f = new TFile(ifile);
	TTree *ntuple1f = (TTree *)f->Get("ntuple1f"); // get trees
	TTree *ntuple2f = (TTree *)f->Get("ntuple2f");
	TTree *ntuple4f = (TTree *)f->Get("ntuple4f");
	TTree *ntuple5f = (TTree *)f->Get("ntuple5f");
	TTree *ntuple6 = (TTree *)f->Get("ntuple6");
	// define variables to hold trees variables

	Output.open("AsymVsMinvTPCOrTOF.txt");

	float eta_pair;
	float PhiRS;
	float fspinconfig;
	float cone;
	float pT_pair;
	float Minv;
	float PhiRSB;
	float PhiRSY;
	float fitPts_min_pair;
	float polB_corr, polY_corr;
	double pi = 3.14159265359;

	double NpairsBGt[5][9] = {0};
	double NpairsBLt[5][9] = {0};
	double pTpairsBGt[5][9] = {0};
	double pTpairsBLt[5][9] = {0};
	double MpairsBGt[5][9] = {0};
	double MpairsBLt[5][9] = {0};
	double NpairsYGt[5][9] = {0};
	double NpairsYLt[5][9] = {0};
	double pTpairsYGt[5][9] = {0};
	double pTpairsYLt[5][9] = {0};
	double MpairsYGt[5][9] = {0};
	double MpairsYLt[5][9] = {0};
	double Npairs[9];
	double pTpairs[9];
	double Mpairs[9];
	double Minv1[9]; // for invariant mass bin average
	double Minv2[9];
	double Minv3[9];
	double Minv4[9];
	double Minv5[9];
	for (int k = 0; k < 9; k++)
	{
		Npairs[k] = 0.;
		pTpairs[k] = 0.;
		Mpairs[k] = 0.;
		Minv1[k] = 0;
		Minv2[k] = 0;
		Minv3[k] = 0;
		Minv4[k] = 0;
		Minv5[k] = 0;
	}

	// get the values of variables from tree

	ntuple1f->SetBranchAddress("fspinconfig", &fspinconfig);
	ntuple1f->SetBranchAddress("cone", &cone);
	ntuple2f->SetBranchAddress("Minv", &Minv);
	ntuple2f->SetBranchAddress("pT_pair", &pT_pair);
	ntuple2f->SetBranchAddress("eta_pair", &eta_pair);
	ntuple4f->SetBranchAddress("PhiRSB", &PhiRSB);
	ntuple4f->SetBranchAddress("PhiRSY", &PhiRSY);
	ntuple5f->SetBranchAddress("fitPts_min_pair", &fitPts_min_pair);
	ntuple6->SetBranchAddress("polB_corr", &polB_corr);
	ntuple6->SetBranchAddress("polY_corr", &polY_corr);

	// To store average polarization values from histograms
	double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;

	Int_t nentries = (Int_t)ntuple1f->GetEntries(); // entries of trees
	// add friend to the tree. no root file to add friend means the tree is on the same root file
	ntuple1f->AddFriend("ntuple2f");
	ntuple1f->AddFriend("ntuple4f");
	ntuple1f->AddFriend("ntuple5f");
	// double pT[10]={2.80, 3.47,  3.79, 4.12, 4.49, 4.938,  5.505, 6.300, 7.66, 15 };
	// double eta_range[10]={-1.200,   -0.668,  -0.469, -0.281, -0.096, 0.089, 0.275, 0.47,   0.675,  1.2 };
	Double_t avg_pT_pair[5] = {0};
	Double_t avg_pT_pairGtB[5] = {0};
	Double_t avg_pT_pairLtB[5] = {0};
	Double_t avg_pT_pairGtY[5] = {0};
	Double_t avg_pT_pairLtY[5] = {0};
	// double pT[6]={2.80, 3.727, 4.343, 5.157, 6.529, 15.00};
	// double M[10]={0.250, 0.403, 0.516, 0.612, 0.711, 0.803, 0.921, 1.070, 1.286, 4.000};
	double pT[6] = {2.80, 3.71, 4.303, 5.084, 6.404, 15.00};
	// Invariant mass binning for each pT-Bin
	double M1[10] = {0.20, 0.3891, 0.4792, 0.5852, 0.6742, 0.7623, 0.8566, 0.9689, 1.098, 4.000};
	double M2[10] = {0.20, 0.3900, 0.4786, 0.5858, 0.6792, 0.7719, 0.8799, 1.019, 1.207, 4.000};
	double M3[10] = {0.20, 0.3973, 0.4900, 0.6004, 0.6985, 0.7921, 0.912, 1.075, 1.304, 4.000};
	double M4[10] = {0.200, 0.4096, 0.524, 0.6246, 0.7246, 0.8218, 0.9521, 1.1418, 1.418, 4.000};
	double M5[10] = {0.20, 0.4397, 0.5648, 0.6778, 0.7744, 0.8823, 1.031, 1.2582, 1.6275, 4.00};

	int ne = (int)ntuple6->GetEntries();
	for (int pol = 0; pol < ne; pol++)
	{
		ntuple6->GetEntry(pol);
		hpolB->Fill(polB_corr);
		hpolY->Fill(polY_corr);
	}

	TCanvas *fitCanvas[5];
	for (Int_t pbin = 0; pbin < 5; pbin++)
	{
		int NpTGtUpB[165] = {0};
		int NpTLtUpB[165] = {0};
		int NpTGtDnB[165] = {0};
		int NpTLtDnB[165] = {0};
		int NMinvGtUpB[165] = {0};
		int NMinvLtUpB[165] = {0};
		int NMinvGtDnB[165] = {0};
		int NMinvLtDnB[165] = {0};
		int NetaUpB[165] = {0};
		int NetaDnB[165] = {0};
		int NpTGtUpY[165] = {0};
		int NpTLtUpY[165] = {0};
		int NpTGtDnY[165] = {0};
		int NpTLtDnY[165] = {0};
		int NMinvGtUpY[165] = {0};
		int NMinvLtUpY[165] = {0};
		int NMinvGtDnY[165] = {0};
		int NMinvLtDnY[165] = {0};
		int NetaUpY[165] = {0};
		int NetaDnY[165] = {0};
		int NpTGtUpT[165] = {0};
		int NpTLtUpT[165] = {0};
		int NpTGtDnT[165] = {0};
		int NpTLtDnT[165] = {0};
		int NMinvGtUpT[165] = {0};
		int NMinvLtUpT[165] = {0};
		int NMinvGtDnT[165] = {0};
		int NMinvLtDnT[165] = {0};
		int NetaUpT[165] = {0};
		int NetaDnT[165] = {0};

		double M[10] = {0};
		// choose correspondig invariant mass bin for different  pT-Bin
		if (pbin == 0)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M1[k];
			}
		if (pbin == 1)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M2[k];
			}
		if (pbin == 2)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M3[k];
			}
		if (pbin == 3)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M4[k];
			}
		if (pbin == 4)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M5[k];
			}

		// cout << "Pbin: "<< pbin<<", M[10]={"<<M[0]<<", "<< M[1]<<", "<< M[2]<<", "<< M[3]<<", "<< M[4]<<", "<< M[5]<<", "<< M[6]<<", "<< M[7]<<", "<< M[8]<<", "<< M[9]<<"}"<<endl;;
		cout << "Entries: " << nentries << endl;
		for (Int_t i = 0; i < nentries; i++)
		// for (Int_t i=0; i < 10000; i++ )
		{

			ntuple1f->GetEntry(i);
			if (pT_pair < pT[pbin] || pT_pair >= pT[pbin + 1])
				continue;
			if (cone > 0.7)
				continue;
			if (Minv > 4.)
				continue;
			if (Minv > 0.4876 && Minv < 0.5076)
				continue; // dodge K0 mass range, doesn't cause asymmetry. Implemention was not working before.This is fine now.
			if (fitPts_min_pair < 15)
				continue;

			hpT_pair[pbin]->Fill(pT_pair);
			// cout << "polarizztion values   "<< polB_corr << "  "<<polY_corr <<endl;
			// cout << "selection cuts are good ...working on phi loop.... " << endl;

			// cout <<cone<<endl;
			//-------------BLUE-------------------------------------
			// Phi
			for (int phi = 0; phi < 16; phi++)
			{
				if (PhiRSB >= (phi - 8.) / 8. * pi && PhiRSB <= (phi - 7.) / 8. * pi)
				{
					// Invariant mass loop
					for (int m = 0; m < 9; m++)
					{
						if (Minv >= M[m] && Minv < M[m + 1])
						{
							Npairs[m] = Npairs[m] + 1;
							pTpairs[m] = pTpairs[m] + pT_pair;
							Mpairs[m] = Mpairs[m] + Minv;
							if (eta_pair > 0)
							{
								hpT_pairGtB[pbin]->Fill(pT_pair);
								NpairsBGt[pbin][m] = NpairsBGt[pbin][m] + 1;
								pTpairsBGt[pbin][m] = pTpairsBGt[pbin][m] + pT_pair;
								MpairsBGt[pbin][m] = MpairsBGt[pbin][m] + Minv;
							}
							if (eta_pair < 0)
							{
								hpT_pairLtB[pbin]->Fill(pT_pair);
								NpairsBLt[pbin][m] = NpairsBLt[pbin][m] + 1;
								pTpairsBLt[pbin][m] = pTpairsBLt[pbin][m] + pT_pair;
								MpairsBLt[pbin][m] = MpairsBLt[pbin][m] + Minv;
							}

							if (fspinconfig == 51 || fspinconfig == 53)
							{
								if (eta_pair > 0)
								{
									NMinvGtUpB[m * 16 + phi]++;
								}
								if (eta_pair < 0)
								{
									NMinvLtUpB[m * 16 + phi]++;
								}
							}
							if (fspinconfig == 83 || fspinconfig == 85)
							{
								if (eta_pair > 0)
								{
									NMinvGtDnB[m * 16 + phi]++;
								}
								if (eta_pair < 0)
								{
									NMinvLtDnB[m * 16 + phi]++;
								}
							}
						}
					} // Minv loop

				}															// phi-range loop
			}																// phi loop
			/**************** BLUE BEAM ENDS *****************************/ ///////

			////***************** YELLOW BEAM *****************************////////////
			// Phi loop
			for (int phi = 0; phi < 16; phi++)
			{
				if (PhiRSY >= (phi - 8.) / 8. * pi && PhiRSY <= (phi - 7.) / 8. * pi)
				{
					// Minv loop
					for (int m = 0; m < 9; m++)
					{
						if (Minv >= M[m] && Minv < M[m + 1])
						{
							if (eta_pair < 0)
							{
								hpT_pairGtY[pbin]->Fill(pT_pair);
								NpairsYGt[pbin][m] = NpairsYGt[pbin][m] + 1;
								pTpairsYGt[pbin][m] = pTpairsYGt[pbin][m] + pT_pair;
								MpairsYGt[pbin][m] = MpairsYGt[pbin][m] + Minv;
							}
							if (eta_pair > 0)
							{
								hpT_pairLtY[pbin]->Fill(pT_pair);
								NpairsYLt[pbin][m] = NpairsYLt[pbin][m] + 1;
								pTpairsYLt[pbin][m] = pTpairsYLt[pbin][m] + pT_pair;
								MpairsYLt[pbin][m] = MpairsYLt[pbin][m] + Minv;
							}
							if (fspinconfig == 51 || fspinconfig == 83)
							{
								// if(eta_pair>0)
								if (eta_pair < 0)
								{
									NMinvGtUpY[m * 16 + phi]++;
								}
								// if(eta_pair<0)
								if (eta_pair > 0)
								{
									NMinvLtUpY[m * 16 + phi]++;
								}
							}
							if (fspinconfig == 53 || fspinconfig == 85)
							{
								// if(eta_pair>0)
								if (eta_pair < 0)
								{
									NMinvGtDnY[m * 16 + phi]++;
								}
								// if(eta_pair<0)
								if (eta_pair > 0)
								{
									NMinvLtDnY[m * 16 + phi]++;
								}
							}
						}
					}													 // Minv loop
				}														 // PhiRS range loop
			}															 // phi loop
		}																 // entries  for loop
		/******************** YELLOW BEAM ENDS ************************/ /////

		// Get average polarizations and rms from histograms
		// avgPolB=0.57534; rmsB=0.0370551;
		// avgPolY=0.58560; rmsY=0.0386434;

		avgPolB = hpolB->GetMean();
		avgPolY = hpolY->GetMean();
		rmsB = hpolB->GetRMS();
		rmsY = hpolY->GetRMS();

		// BLUE <P>: 0.57534, Err_P: 0.0370551
		// YELLOW <P>: 0.5856, Err_P: 0.0386434

		// Store Minv bin averages for each pT bin in an array for plotting purpose
		/*if(pbin==0){for(int k0=0;k0<9;k0++){Minv1[k0]=Mpairs[k0]/(double)Npairs[k0];}}
		if(pbin==1){for(int k1=0;k1<9;k1++){Minv2[k1]=Mpairs[k1]/(double)Npairs[k1];}}
		if(pbin==2){for(int k2=0;k2<9;k2++){Minv3[k2]=Mpairs[k2]/(double)Npairs[k2];}}
		if(pbin==3){for(int k3=0;k3<9;k3++){Minv4[k3]=Mpairs[k3]/(double)Npairs[k3];}}
		if(pbin==4){for(int k4=0;k4<9;k4++){Minv5[k4]=Mpairs[k4]/(double)Npairs[k4];}}
		*/
		if (pbin == 0)
		{
			for (int k0 = 0; k0 < 9; k0++)
			{
				Minv1[k0] = 0.5 * ((MpairsBGt[0][k0] / (double)NpairsBGt[0][k0]) + (MpairsYGt[0][k0] / (double)NpairsYGt[0][k0]));
			}
		}
		if (pbin == 1)
		{
			for (int k0 = 0; k0 < 9; k0++)
			{
				Minv2[k0] = 0.5 * ((MpairsBGt[1][k0] / (double)NpairsBGt[1][k0]) + (MpairsYGt[1][k0] / (double)NpairsYGt[1][k0]));
			}
		}
		if (pbin == 2)
		{
			for (int k0 = 0; k0 < 9; k0++)
			{
				Minv3[k0] = 0.5 * ((MpairsBGt[2][k0] / (double)NpairsBGt[2][k0]) + (MpairsYGt[2][k0] / (double)NpairsYGt[2][k0]));
			}
		}
		if (pbin == 3)
		{
			for (int k0 = 0; k0 < 9; k0++)
			{
				Minv4[k0] = 0.5 * ((MpairsBGt[3][k0] / (double)NpairsBGt[3][k0]) + (MpairsYGt[3][k0] / (double)NpairsYGt[3][k0]));
			}
		}
		if (pbin == 4)
		{
			for (int k0 = 0; k0 < 9; k0++)
			{
				Minv5[k0] = 0.5 * ((MpairsBGt[4][k0] / (double)NpairsBGt[4][k0]) + (MpairsYGt[4][k0] / (double)NpairsYGt[4][k0]));
			}
		}
		// Asymmetry calculation for BLUE eta >0
		double dAdB, dAdC, dAdD, dAdE, dAdP; // variables for error calculation
		// Use rms for polarization error from distribution!
		double dP_B = rmsB; // Polarization errors blue
		double dP_Y = rmsY; // Polarization errors yellow
		double dA[165];		// Asymmetry amplitude from sine fit
		double B, C, D, E;
		double a, b;
		double Asym[165];
		double pi = 3.14159265359;

		fitCanvas[pbin] = new TCanvas(Form("fitCanvas_pTBin%i", pbin), "", 650, 550);
		fitCanvas[pbin]->Print(Form("FitPlots_pTbin_%i_AsymVsMinvTPCOrTOF.pdf(", pbin));

		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(NMinvGtUpB[m * 16 + ang] * NMinvGtDnB[m * 16 + ang + 8]);
					b = sqrt(NMinvGtDnB[m * 16 + ang] * NMinvGtUpB[m * 16 + ang + 8]);
					B = NMinvGtUpB[m * 16 + ang];
					C = NMinvGtUpB[m * 16 + ang + 8];
					D = NMinvGtDnB[m * 16 + ang];
					E = NMinvGtDnB[m * 16 + ang + 8];
				}
				if (ang > 7)
				{
					a = sqrt(NMinvGtUpB[m * 16 + ang] * NMinvGtDnB[m * 16 + ang - 8]);
					b = sqrt(NMinvGtDnB[m * 16 + ang] * NMinvGtUpB[m * 16 + ang - 8]);
					B = NMinvGtUpB[m * 16 + ang];
					C = NMinvGtUpB[m * 16 + ang - 8];
					D = NMinvGtDnB[m * 16 + ang];
					E = NMinvGtDnB[m * 16 + ang - 8];
				}
				Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b));
				// cout << "asymmetry for blue beam for bin "<< ang <<" --> "<< Asym[ang] << endl;
				dAdB = (1. / avgPolB) * E * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdE = (1. / avgPolB) * B * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdD = (-1. / avgPolB) * C * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdC = (-1. / avgPolB) * D * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
				dA[ang] = sqrt((fabs(dAdB) * sqrt(B)) * *2 + (fabs(dAdC) * sqrt(C)) * *2 + (fabs(dAdD) * sqrt(D)) * *2 + (fabs(dAdE) * sqrt(E)) * *2 + (fabs(dAdP) * dP_B) * *2);
				// cout <<"polB:"<<avgPolB<<", polY: "<<avgPolY<< ", a="<<a<<", b="<<b<<", B="<<B<<", C="<<C<<", D="<<D<<", E="<<E<<", Asym["<<ang<<"]="<<Asym[ang]<<", dA["<<ang<<"]="<<dA[ang]<<endl;
			} // angle loop

			double Abg[9]; // Here Abg refers to Asymmetry for Blue in eta >0
			double deltaAbg[9];
			char name[600];
			char title[600];
			double chi2Ndf[9];
			// double angle[16]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi,1./16.*pi,3./16.*pi,5./16.*pi,7./16.*pi,9./16.*pi,11./16.*pi,13./16.*pi,15./16.*pi};

			// fitting must be done for the one hemisphere (comments from collaboration meeting 2021/03/03)
			double angle[8] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi};
			// double ex[16]={0};
			double ex[8] = {0};
			gStyle->SetOptDate(0);
			// grt = new TGraphErrors(16,angle,Asym,ex,dA);
			grt = new TGraphErrors(8, angle, Asym, ex, dA);
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0.);
			fit->SetParameter(0, 0.0001);
			grt->Fit(fit, "R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			sprintf(title, "pT bin %i, Minv bin %i, BLUE, #eta > 0", pbin, m);
			grt->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			chi2NdfBGt->Fill(chi2Ndf[m]);
			// cout<<"ChiSquare = "<<chi2Ndf[m]<<endl;
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->CenterTitle();
			grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Abg[m] = fit->GetParameter(0);
			deltaAbg[m] = fit->GetParError(0);
			sprintf(name, "pT%i_Minvbin%i_etaGt_BLUETPCOrTOF.pdf", pbin, m);
			// c1->SaveAs(name);
			fitCanvas[pbin]->Print(Form("FitPlots_pTbin_%i_AsymVsMinvTPCOrTOF.pdf", pbin));

		} // invariant mass loop

		// YELLOW eta > 0
		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(NMinvGtUpY[m * 16 + ang] * NMinvGtDnY[m * 16 + ang + 8]);
					b = sqrt(NMinvGtDnY[m * 16 + ang] * NMinvGtUpY[m * 16 + ang + 8]);
					B = NMinvGtUpY[m * 16 + ang];
					C = NMinvGtUpY[m * 16 + ang + 8];
					D = NMinvGtDnY[m * 16 + ang];
					E = NMinvGtDnY[m * 16 + ang + 8];
				}
				if (ang > 7)
				{
					a = sqrt(NMinvGtUpY[m * 16 + ang] * NMinvGtDnY[m * 16 + ang - 8]);
					b = sqrt(NMinvGtDnY[m * 16 + ang] * NMinvGtUpY[m * 16 + ang - 8]);
					B = NMinvGtUpY[m * 16 + ang];
					C = NMinvGtUpY[m * 16 + ang - 8];
					D = NMinvGtDnY[m * 16 + ang];
					E = NMinvGtDnY[m * 16 + ang - 8];
				}
				Asym[ang] = (1. / avgPolY) * ((a - b) / (a + b));

				dAdB = (1. / avgPolY) * E * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdE = (1. / avgPolY) * B * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdD = (-1. / avgPolY) * C * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdC = (-1. / avgPolY) * D * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdP = (-1. / (avgPolY * avgPolY)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
				dA[ang] = sqrt((fabs(dAdB) * sqrt(B)) * *2 + (fabs(dAdC) * sqrt(C)) * *2 + (fabs(dAdD) * sqrt(D)) * *2 + (fabs(dAdE) * sqrt(E)) * *2 + (fabs(dAdP) * dP_Y) * *2);
			}
			double Ayg[9];
			double deltaAyg[9];
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi};
			double ex[8] = {0};
			gStyle->SetOptDate(0);
			gr = new TGraphErrors(8, angle, Asym, ex, dA);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0.);
			fit->SetParameter(0, 0.0001);
			gr->Fit(fit, "R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			sprintf(title, "pT bin %i, Minv bin %i, YELLOW, #eta > 0", pbin, m);
			gr->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			chi2NdfYGt->Fill(chi2Ndf[m]);
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle(kTRUE);
			gr->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Ayg[m] = fit->GetParameter(0);
			deltaAyg[m] = fit->GetParError(0);
			sprintf(name, "pT%i_Minvbin%i_etaGt_YELLOWTPCOrTOF.pdf", pbin, m);
			// c1->SaveAs(name);
			fitCanvas[pbin]->Print(Form("FitPlots_pTbin_%i_AsymVsMinvTPCOrTOF.pdf", pbin));

		} // Minv loop
		  // BLUE beam in eta < 0
		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(NMinvLtUpB[m * 16 + ang] * NMinvLtDnB[m * 16 + ang + 8]);
					b = sqrt(NMinvLtDnB[m * 16 + ang] * NMinvLtUpB[m * 16 + ang + 8]);
					B = NMinvLtUpB[m * 16 + ang];
					C = NMinvLtUpB[m * 16 + ang + 8];
					D = NMinvLtDnB[m * 16 + ang];
					E = NMinvLtDnB[m * 16 + ang + 8];
				}
				if (ang > 7)
				{
					a = sqrt(NMinvLtUpB[m * 16 + ang] * NMinvLtDnB[m * 16 + ang - 8]);
					b = sqrt(NMinvLtDnB[m * 16 + ang] * NMinvLtUpB[m * 16 + ang - 8]);
					B = NMinvLtUpB[m * 16 + ang];
					C = NMinvLtUpB[m * 16 + ang - 8];
					D = NMinvLtDnB[m * 16 + ang];
					E = NMinvLtDnB[m * 16 + ang - 8];
				}
				Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b));
				dAdB = (1. / avgPolB) * E * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdE = (1. / avgPolB) * B * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdD = (-1. / avgPolB) * C * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdC = (-1. / avgPolB) * D * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
				dA[ang] = sqrt((fabs(dAdB) * sqrt(B)) * *2 + (fabs(dAdC) * sqrt(C)) * *2 + (fabs(dAdD) * sqrt(D)) * *2 + (fabs(dAdE) * sqrt(E)) * *2 + (fabs(dAdP) * dP_B) * *2);
			} // angle loop

			double Abl[9];
			double deltaAbl[9]; /// bl refers to blue beam eta less than 0
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi};
			double ex[8] = {0};
			gStyle->SetOptDate(0);
			grt = new TGraphErrors(8, angle, Asym, ex, dA);
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0.);
			fit->SetParameter(0, 0.0001);
			grt->Fit(fit, "R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			sprintf(title, "pT bin %i, Minv bin %i, BLUE, #eta < 0", pbin, m);
			grt->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			chi2NdfBLt->Fill(chi2Ndf[m]);
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->SetTitleOffset(1);
			grt->GetYaxis()->CenterTitle(kTRUE);
			grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Abl[m] = fit->GetParameter(0);
			deltaAbl[m] = fit->GetParError(0);
			sprintf(name, "pT%i_Minvbin%i_etaLt_BLUETPCOrTOF.pdf", pbin, m);
			// c1->Print(name);
			fitCanvas[pbin]->Print(Form("FitPlots_pTbin_%i_AsymVsMinvTPCOrTOF.pdf", pbin));

		} // invariant mass loop

		// YELLOW Beam Eta<0
		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(NMinvLtUpY[m * 16 + ang] * NMinvLtDnY[m * 16 + ang + 8]);
					b = sqrt(NMinvLtDnY[m * 16 + ang] * NMinvLtUpY[m * 16 + ang + 8]);
					B = NMinvLtUpY[m * 16 + ang];
					C = NMinvLtUpY[m * 16 + ang + 8];
					D = NMinvLtDnY[m * 16 + ang];
					E = NMinvLtDnY[m * 16 + ang + 8];
				}
				if (ang > 7)
				{
					a = sqrt(NMinvLtUpY[m * 16 + ang] * NMinvLtDnY[m * 16 + ang - 8]);
					b = sqrt(NMinvLtDnY[m * 16 + ang] * NMinvLtUpY[m * 16 + ang - 8]);
					B = NMinvLtUpY[m * 16 + ang];
					C = NMinvLtUpY[m * 16 + ang - 8];
					D = NMinvLtDnY[m * 16 + ang];
					E = NMinvLtDnY[m * 16 + ang - 8];
				}
				Asym[ang] = (1. / avgPolY) * ((a - b) / (a + b));

				dAdB = (1. / avgPolY) * E * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdE = (1. / avgPolY) * B * sqrt(D * C) / (sqrt(B * E) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdD = (-1. / avgPolY) * C * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdC = (-1. / avgPolY) * D * sqrt(B * E) / (sqrt(D * C) * ((sqrt(B * E) + sqrt(D * C)) * *2));
				dAdP = (-1. / (avgPolY * avgPolY)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
				dA[ang] = sqrt((fabs(dAdB) * sqrt(B)) * *2 + (fabs(dAdC) * sqrt(C)) * *2 + (fabs(dAdD) * sqrt(D)) * *2 + (fabs(dAdE) * sqrt(E)) * *2 + (fabs(dAdP) * dP_Y) * *2);
			}
			double Ayl[9];
			double deltaAyl[9];
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[8] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi};
			double ex[8] = {0};
			gStyle->SetOptDate(0);
			gr = new TGraphErrors(8, angle, Asym, ex, dA);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0.);
			fit->SetParameter(0, 0.0001);
			gr->Fit(fit, "R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			sprintf(title, "pT bin %i, Minv bin %i, YELLOW, #eta < 0", pbin, m);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			chi2NdfYLt->Fill(chi2Ndf[m]);
			gr->SetTitle(title);
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle(kTRUE);
			gr->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Ayl[m] = fit->GetParameter(0);
			deltaAyl[m] = fit->GetParError(0);
			sprintf(name, "pTbin%i_Minvbin%i_etaLt_YELLOWTPCOrTOF.pdf", pbin, m);
			fitCanvas[pbin]->Print(Form("FitPlots_pTbin_%i_AsymVsMinvTPCOrTOF.pdf", pbin));

		} // Minv loop
		fitCanvas[pbin]->Print(Form("FitPlots_pTbin_%i_AsymVsMinvTPCOrTOF.pdf)", pbin));

		// Yellow beam asymmetry ends

		if (pbin == 0)
		{ // for Blue beam eta > 0
			double A_pT1bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_pT1bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
			// for Blue beam eta < 0
			double A_pT1bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_pT1bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
			// for Yellow beam eta >0
			double A_pT1yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT1yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT1yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT1yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}

		if (pbin == 1)
		{ // for Blue beam eta > 0
			double A_pT2bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_pT2bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT2bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_pT2bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT2yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT2yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT2yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT2yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		if (pbin == 2)
		{
			double A_pT3bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_pT3bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT3bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_pT3bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT3yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT3yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT3yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT3yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		if (pbin == 3)
		{
			double A_pT4bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_pT4bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT4bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_pT4bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT4yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT4yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT4yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT4yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}
		if (pbin == 4)
		{
			double A_pT5bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
			double deltaA_pT5bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};

			double A_pT5bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
			double deltaA_pT5bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};

			double A_pT5yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
			double deltaA_pT5yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};

			double A_pT5yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			double deltaA_pT5yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyg[7], deltaAyg[8]};
		}

	} // pbin bin loop

	for (Int_t n = 0; n < 5; n++)
	{
		avg_pT_pair[n] = hpT_pair[n]->GetMean();
		avg_pT_pairGtB[n] = hpT_pairGtB[n]->GetMean();
		avg_pT_pairLtB[n] = hpT_pairLtB[n]->GetMean();
		avg_pT_pairGtY[n] = hpT_pairGtY[n]->GetMean();
		avg_pT_pairLtY[n] = hpT_pairLtY[n]->GetMean();
	} // cout << "avg_pT_pair["<<n<<"]="<<avg_pT_pair[n]<< endl;}

	double MinvBGt[5][9] = {0};
	double MinvBLt[5][9] = {0};
	double MinvYGt[5][9] = {0};
	double MinvYLt[5][9] = {0};
	double MinvAvgGt[5][9] = {0};
	double MinvAvgLt[5][9] = {0};
	double pTBGt[5][9] = {0};
	double pTBLt[5][9] = {0};
	double pTYGt[5][9] = {0};
	double pTYLt[5][9] = {0};
	double pTAvgGt[5][9] = {0};
	double pTAvgLt[5][9] = {0};
	double avgpTpairGt[5] = {0};
	double avgpTpairLt[5] = {0};
	for (int pbin = 0; pbin < 5; pbin++)
	{
		avgpTpairGt[pbin] = 0.5 * (avg_pT_pairGtB[pbin] + avg_pT_pairGtY[pbin]);
		avgpTpairLt[pbin] = 0.5 * (avg_pT_pairLtB[pbin] + avg_pT_pairLtY[pbin]);

		for (int mbin = 0; mbin < 9; mbin++)
		{
			MinvBGt[pbin][mbin] = MpairsBGt[pbin][mbin] / (double)NpairsBGt[pbin][mbin];
			pTBGt[pbin][mbin] = pTpairsBGt[pbin][mbin] / (double)NpairsBGt[pbin][mbin];
			MinvYGt[pbin][mbin] = MpairsYGt[pbin][mbin] / (double)NpairsYGt[pbin][mbin];
			pTYGt[pbin][mbin] = pTpairsYGt[pbin][mbin] / (double)NpairsYGt[pbin][mbin];
			MinvAvgGt[pbin][mbin] = 0.5 * (MinvBGt[pbin][mbin] + MinvYGt[pbin][mbin]);
			pTAvgGt[pbin][mbin] = 0.5 * (pTBGt[pbin][mbin] + pTYGt[pbin][mbin]);

			MinvBLt[pbin][mbin] = MpairsBLt[pbin][mbin] / (double)NpairsBLt[pbin][mbin];
			pTBLt[pbin][mbin] = pTpairsBLt[pbin][mbin] / (double)NpairsBLt[pbin][mbin];
			MinvYLt[pbin][mbin] = MpairsYLt[pbin][mbin] / (double)NpairsYLt[pbin][mbin];
			pTYLt[pbin][mbin] = pTpairsYLt[pbin][mbin] / (double)NpairsYLt[pbin][mbin];
			MinvAvgLt[pbin][mbin] = 0.5 * (MinvBLt[pbin][mbin] + MinvYLt[pbin][mbin]);
			pTAvgLt[pbin][mbin] = 0.5 * (pTBLt[pbin][mbin] + pTYLt[pbin][mbin]);
			// cout<<"pT BGt: "<< pTBGt[pbin][mbin]<<"  pT YGt: "<< pTYGt[pbin][mbin]<<endl;
		}
	}

	Output << "avgpTpairGtB[5]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << avg_pT_pairGtB[pbin] << ",";
		if (pbin == 4)
			Output << "};" << endl;
	}
	Output << "avgpTpairGtY[5]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << avg_pT_pairGtY[pbin] << ",";
		if (pbin == 4)
			Output << "};" << endl;
	}
	Output << "avgpTpairLtB[5]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << avg_pT_pairLtB[pbin] << ",";
		if (pbin == 4)
			Output << "};" << endl;
	}
	Output << "avgpTpairLtY[5]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << avg_pT_pairLtY[pbin] << ",";
		if (pbin == 4)
			Output << "};" << endl;
	}

	Output << "avgpTpairGt[5]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << avgpTpairGt[pbin] << ",";
		if (pbin == 4)
			Output << "};" << endl;
	}
	Output << "avgpTpairLt[5]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << avgpTpairLt[pbin] << ",";
		if (pbin == 4)
			Output << "};" << endl;
	}

	Output << "AvgMinvGt[5][9]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << "{";
		for (int mbin = 0; mbin < 9; mbin++)
		{
			Output << MinvAvgGt[pbin][mbin] << ",";
			if (mbin == 8)
				Output << "};" << endl;
		}
	}
	Output << "};" << endl;
	Output << "AvgpTGt[5][9]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << "{";
		for (int mbin = 0; mbin < 9; mbin++)
		{
			Output << pTAvgGt[pbin][mbin] << ",";
			if (mbin == 8)
				Output << "};" << endl;
		}
	}
	Output << "};" << endl;
	Output << "AvgMinvLt[5][9]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << "{";
		for (int mbin = 0; mbin < 9; mbin++)
		{
			Output << MinvAvgLt[pbin][mbin] << ",";
			if (mbin == 8)
				Output << "};" << endl;
		}
	}
	Output << "};" << endl;
	Output << "AvgpTLt[5][9]={";
	for (int pbin = 0; pbin < 5; pbin++)
	{
		Output << "{";
		for (int mbin = 0; mbin < 9; mbin++)
		{
			Output << pTAvgLt[pbin][mbin] << ",";
			if (mbin == 8)
				Output << "};" << endl;
		}
	}
	Output << "};" << endl;
	// calculate average asymmetry

	// simple average asymmetry between BLUe and YELLOW
	double avgA_pT1g[9] = {0};
	double errA_pT1g[9] = {0};
	double avgA_pT1l[9] = {0};
	double errA_pT1l[9] = {0};
	double avgA_pT2g[9] = {0};
	double errA_pT2g[9] = {0};
	double avgA_pT2l[9] = {0};
	double errA_pT2l[9] = {0};
	double avgA_pT3g[9] = {0};
	double errA_pT3g[9] = {0};
	double avgA_pT3l[9] = {0};
	double errA_pT3l[9] = {0};
	double avgA_pT4g[9] = {0};
	double errA_pT4g[9] = {0};
	double avgA_pT4l[9] = {0};
	double errA_pT4l[9] = {0};
	double avgA_pT5g[9] = {0};
	double errA_pT5g[9] = {0};
	double avgA_pT5l[9] = {0};
	double errA_pT5l[9] = {0};
	// weighted average asymmetry between BLUE and YELLOW with number of π^+π^- events in each bin.
	double WavgA_pT1g[9] = {0};
	double WerrA_pT1g[9] = {0};
	double WavgA_pT1l[9] = {0};
	double WerrA_pT1l[9] = {0};
	double WavgA_pT2g[9] = {0};
	double WerrA_pT2g[9] = {0};
	double WavgA_pT2l[9] = {0};
	double WerrA_pT2l[9] = {0};
	double WavgA_pT3g[9] = {0};
	double WerrA_pT3g[9] = {0};
	double WavgA_pT3l[9] = {0};
	double WerrA_pT3l[9] = {0};
	double WavgA_pT4g[9] = {0};
	double WerrA_pT4g[9] = {0};
	double WavgA_pT4l[9] = {0};
	double WerrA_pT4l[9] = {0};
	double WavgA_pT5g[9] = {0};
	double WerrA_pT5g[9] = {0};
	double WavgA_pT5l[9] = {0};
	double WerrA_pT5l[9] = {0};

	for (Int_t ii = 0; ii < 9; ii++)
	{
		// simple averages
		avgA_pT1g[ii] = (A_pT1bg[ii] + A_pT1yg[ii]) / 2.;
		avgA_pT1l[ii] = (A_pT1bl[ii] + A_pT1yl[ii]) / 2.;
		avgA_pT2g[ii] = (A_pT2bg[ii] + A_pT2yg[ii]) / 2.;
		avgA_pT2l[ii] = (A_pT2bl[ii] + A_pT2yl[ii]) / 2.;
		avgA_pT3g[ii] = (A_pT3bg[ii] + A_pT3yg[ii]) / 2.;
		avgA_pT3l[ii] = (A_pT3bl[ii] + A_pT3yl[ii]) / 2.;
		avgA_pT4g[ii] = (A_pT4bg[ii] + A_pT4yg[ii]) / 2.;
		avgA_pT4l[ii] = (A_pT4bl[ii] + A_pT4yl[ii]) / 2.;
		avgA_pT5g[ii] = (A_pT5bg[ii] + A_pT5yg[ii]) / 2.;
		avgA_pT5l[ii] = (A_pT5bl[ii] + A_pT5yl[ii]) / 2.;

		// total asymmetry error calculation
		// since two asymmetry(BLUE+YELLOW) are averaged, error propagation formaula is used for total asymmetry error
		// if function, f = (a+b)/2, Then Error, df =1/2 √{(∂f/∂a)^2*(∆a)^2+(∂f/∂b)^2*(∆b)^2}

		errA_pT1g[ii] = .5 * sqrt(pow(deltaA_pT1bg[ii], 2) + pow(deltaA_pT1yg[ii], 2)); // total asym err, pT bin 1, eta > 0
		errA_pT1l[ii] = .5 * sqrt(pow(deltaA_pT1bl[ii], 2) + pow(deltaA_pT1yl[ii], 2)); // total asym err, pT bin 1, eta < 0
		errA_pT2g[ii] = .5 * sqrt(pow(deltaA_pT2bg[ii], 2) + pow(deltaA_pT2yg[ii], 2));
		errA_pT2l[ii] = .5 * sqrt(pow(deltaA_pT2bl[ii], 2) + pow(deltaA_pT2yl[ii], 2));
		errA_pT3g[ii] = .5 * sqrt(pow(deltaA_pT3bg[ii], 2) + pow(deltaA_pT3yg[ii], 2));
		errA_pT3l[ii] = .5 * sqrt(pow(deltaA_pT3bl[ii], 2) + pow(deltaA_pT3yl[ii], 2));
		errA_pT4g[ii] = .5 * sqrt(pow(deltaA_pT4bg[ii], 2) + pow(deltaA_pT4yg[ii], 2));
		errA_pT4l[ii] = .5 * sqrt(pow(deltaA_pT4bl[ii], 2) + pow(deltaA_pT4yl[ii], 2));
		errA_pT5g[ii] = .5 * sqrt(pow(deltaA_pT5bg[ii], 2) + pow(deltaA_pT5yg[ii], 2));
		errA_pT5l[ii] = .5 * sqrt(pow(deltaA_pT5bl[ii], 2) + pow(deltaA_pT5yl[ii], 2));

		// weighted averages
		WavgA_pT1g[ii] = (A_pT1bg[ii] * (1 / pow(deltaA_pT1bg[ii], 2)) + A_pT1yg[ii] * (1 / pow(deltaA_pT1yg[ii], 2))) / ((1 / pow(deltaA_pT1bg[ii], 2)) + (1 / pow(deltaA_pT1yg[ii], 2)));
		WavgA_pT2g[ii] = (A_pT2bg[ii] * (1 / pow(deltaA_pT2bg[ii], 2)) + A_pT2yg[ii] * (1 / pow(deltaA_pT2yg[ii], 2))) / ((1 / pow(deltaA_pT2bg[ii], 2)) + (1 / pow(deltaA_pT2yg[ii], 2)));
		WavgA_pT3g[ii] = (A_pT3bg[ii] * (1 / pow(deltaA_pT3bg[ii], 2)) + A_pT3yg[ii] * (1 / pow(deltaA_pT3yg[ii], 2))) / ((1 / pow(deltaA_pT3bg[ii], 2)) + (1 / pow(deltaA_pT3yg[ii], 2)));
		WavgA_pT4g[ii] = (A_pT4bg[ii] * (1 / pow(deltaA_pT4bg[ii], 2)) + A_pT4yg[ii] * (1 / pow(deltaA_pT4yg[ii], 2))) / ((1 / pow(deltaA_pT4bg[ii], 2)) + (1 / pow(deltaA_pT4yg[ii], 2)));
		WavgA_pT5g[ii] = (A_pT5bg[ii] * (1 / pow(deltaA_pT5bg[ii], 2)) + A_pT5yg[ii] * (1 / pow(deltaA_pT5yg[ii], 2))) / ((1 / pow(deltaA_pT5bg[ii], 2)) + (1 / pow(deltaA_pT5yg[ii], 2)));

		WerrA_pT1g[ii] = 1 / sqrt((1 / pow(deltaA_pT1bg[ii], 2)) + (1 / pow(deltaA_pT1yg[ii], 2))); // total asym err, pT bin 1, eta > 0
		WerrA_pT2g[ii] = 1 / sqrt((1 / pow(deltaA_pT2bg[ii], 2)) + (1 / pow(deltaA_pT2yg[ii], 2)));
		WerrA_pT3g[ii] = 1 / sqrt((1 / pow(deltaA_pT3bg[ii], 2)) + (1 / pow(deltaA_pT3yg[ii], 2)));
		WerrA_pT4g[ii] = 1 / sqrt((1 / pow(deltaA_pT4bg[ii], 2)) + (1 / pow(deltaA_pT4yg[ii], 2)));
		WerrA_pT5g[ii] = 1 / sqrt((1 / pow(deltaA_pT5bg[ii], 2)) + (1 / pow(deltaA_pT5yg[ii], 2)));

		WavgA_pT1l[ii] = (A_pT1bl[ii] * (1 / pow(deltaA_pT1bl[ii], 2)) + A_pT1yl[ii] * (1 / pow(deltaA_pT1yl[ii], 2))) / ((1 / pow(deltaA_pT1bl[ii], 2)) + (1 / pow(deltaA_pT1yl[ii], 2)));
		WavgA_pT2l[ii] = (A_pT2bl[ii] * (1 / pow(deltaA_pT2bl[ii], 2)) + A_pT2yl[ii] * (1 / pow(deltaA_pT2yl[ii], 2))) / ((1 / pow(deltaA_pT2bl[ii], 2)) + (1 / pow(deltaA_pT2yl[ii], 2)));
		WavgA_pT3l[ii] = (A_pT3bl[ii] * (1 / pow(deltaA_pT3bl[ii], 2)) + A_pT3yl[ii] * (1 / pow(deltaA_pT3yl[ii], 2))) / ((1 / pow(deltaA_pT3bl[ii], 2)) + (1 / pow(deltaA_pT3yl[ii], 2)));
		WavgA_pT4l[ii] = (A_pT4bl[ii] * (1 / pow(deltaA_pT4bl[ii], 2)) + A_pT4yl[ii] * (1 / pow(deltaA_pT4yl[ii], 2))) / ((1 / pow(deltaA_pT4bl[ii], 2)) + (1 / pow(deltaA_pT4yl[ii], 2)));
		WavgA_pT5l[ii] = (A_pT5bl[ii] * (1 / pow(deltaA_pT5bl[ii], 2)) + A_pT5yl[ii] * (1 / pow(deltaA_pT5yl[ii], 2))) / ((1 / pow(deltaA_pT5bl[ii], 2)) + (1 / pow(deltaA_pT5yl[ii], 2)));

		WerrA_pT1l[ii] = 1 / sqrt((1 / pow(deltaA_pT1bl[ii], 2)) + (1 / pow(deltaA_pT1yl[ii], 2))); // total asym err, pT bin 1, eta < 0
		WerrA_pT2l[ii] = 1 / sqrt((1 / pow(deltaA_pT2bl[ii], 2)) + (1 / pow(deltaA_pT2yl[ii], 2)));
		WerrA_pT3l[ii] = 1 / sqrt((1 / pow(deltaA_pT3bl[ii], 2)) + (1 / pow(deltaA_pT3yl[ii], 2)));
		WerrA_pT4l[ii] = 1 / sqrt((1 / pow(deltaA_pT4bl[ii], 2)) + (1 / pow(deltaA_pT4yl[ii], 2)));
		WerrA_pT5l[ii] = 1 / sqrt((1 / pow(deltaA_pT5bl[ii], 2)) + (1 / pow(deltaA_pT5yl[ii], 2)));
	}

	Output << "********************   FINAL ASYMMETRY , CONE < 0.7 Minv Binning ****************" << endl;

	Output << "<<<<<<<<<<<<<<<   Minv Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
	Output << "double  Minv1[9]= {" << Minv1[0] << ", " << Minv1[1] << ", " << Minv1[2] << ", " << Minv1[3] << ", " << Minv1[4] << ", " << Minv1[5] << ", " << Minv1[6] << ", " << Minv1[7] << ", " << Minv1[8] << "}" << endl;
	Output << "double  Minv2[9]= {" << Minv2[0] << ", " << Minv2[1] << ", " << Minv2[2] << ", " << Minv2[3] << ", " << Minv2[4] << ", " << Minv2[5] << ", " << Minv2[6] << ", " << Minv2[7] << ", " << Minv2[8] << "}" << endl;
	Output << "double  Minv3[9]= {" << Minv3[0] << ", " << Minv3[1] << ", " << Minv3[2] << ", " << Minv3[3] << ", " << Minv3[4] << ", " << Minv3[5] << ", " << Minv3[6] << ", " << Minv3[7] << ", " << Minv3[8] << "}" << endl;
	Output << "double  Minv4[9]= {" << Minv4[0] << ", " << Minv4[1] << ", " << Minv4[2] << ", " << Minv4[3] << ", " << Minv4[4] << ", " << Minv4[5] << ", " << Minv4[6] << ", " << Minv4[7] << ", " << Minv4[8] << "}" << endl;
	Output << "double  Minv5[9]= {" << Minv5[0] << ", " << Minv5[1] << ", " << Minv5[2] << ", " << Minv5[3] << ", " << Minv5[4] << ", " << Minv5[5] << ", " << Minv5[6] << ", " << Minv5[7] << ", " << Minv5[8] << "}" << endl;

	Output << "<<<<<<<< Average Polarization Values >>>>>>>>" << endl;
	Output << "BLUE <P>: " << avgPolB << ", Err_P: " << rmsB << endl;
	Output << "YELLOW <P>: " << avgPolY << ", Err_P: " << rmsY << endl;

	Output << "<<<<<<<<<<<       Eta < 0, Standard Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "double  A_pT1Lt[9] ={" << avgA_pT1l[0] << ", " << avgA_pT1l[1] << ", " << avgA_pT1l[2] << ", " << avgA_pT1l[3] << ", " << avgA_pT1l[4] << ", " << avgA_pT1l[5] << ", " << avgA_pT1l[6] << ", " << avgA_pT1l[7] << ", " << avgA_pT1l[8] << "}" << endl;
	Output << "double  Aerr_pT1Lt[9]={" << errA_pT1l[0] << ", " << errA_pT1l[1] << ", " << errA_pT1l[2] << ", " << errA_pT1l[3] << ", " << errA_pT1l[4] << ", " << errA_pT1l[5] << ", " << errA_pT1l[6] << ", " << errA_pT1l[7] << ", " << errA_pT1l[8] << "}" << endl;
	Output << "double  A_pT2Lt[9]={" << avgA_pT2l[0] << ", " << avgA_pT2l[1] << ", " << avgA_pT2l[2] << ", " << avgA_pT2l[3] << ", " << avgA_pT2l[4] << ", " << avgA_pT2l[5] << ", " << avgA_pT2l[6] << ", " << avgA_pT2l[7] << ", " << avgA_pT2l[8] << "}" << endl;
	Output << "double  Aerr_pT2Lt[9] ={" << errA_pT2l[0] << ", " << errA_pT2l[1] << ", " << errA_pT2l[2] << ", " << errA_pT2l[3] << ", " << errA_pT2l[4] << ", " << errA_pT2l[5] << ", " << errA_pT2l[6] << ", " << errA_pT2l[7] << ", " << errA_pT2l[8] << "}" << endl;
	Output << "double A_pT3Lt[9] ={" << avgA_pT3l[0] << ", " << avgA_pT3l[1] << ", " << avgA_pT3l[2] << ", " << avgA_pT3l[3] << ", " << avgA_pT3l[4] << ", " << avgA_pT3l[5] << ", " << avgA_pT3l[6] << ", " << avgA_pT3l[7] << ", " << avgA_pT3l[8] << "}" << endl;
	Output << "double Aerr_pT3Lt[9] ={" << errA_pT3l[0] << ", " << errA_pT3l[1] << ", " << errA_pT3l[2] << ", " << errA_pT3l[3] << ", " << errA_pT3l[4] << ", " << errA_pT3l[5] << ", " << errA_pT3l[6] << ", " << errA_pT3l[7] << ", " << errA_pT3l[8] << "}" << endl;
	Output << "double  A_pT4Lt[9] ={" << avgA_pT4l[0] << ", " << avgA_pT4l[1] << ", " << avgA_pT4l[2] << ", " << avgA_pT4l[3] << ", " << avgA_pT4l[4] << ", " << avgA_pT4l[5] << ", " << avgA_pT4l[6] << ", " << avgA_pT4l[7] << ", " << avgA_pT4l[8] << "}" << endl;
	Output << "double Aerr_pT4Lt[9] ={" << errA_pT4l[0] << ", " << errA_pT4l[1] << ", " << errA_pT4l[2] << ", " << errA_pT4l[3] << ", " << errA_pT4l[4] << ", " << errA_pT4l[5] << ", " << errA_pT4l[6] << ", " << errA_pT4l[7] << ", " << errA_pT4l[8] << "}" << endl;
	Output << "double  A_pT5Lt[9] ={" << avgA_pT5l[0] << ", " << avgA_pT5l[1] << ", " << avgA_pT5l[2] << ", " << avgA_pT5l[3] << ", " << avgA_pT5l[4] << ", " << avgA_pT5l[5] << ", " << avgA_pT5l[6] << ", " << avgA_pT5l[7] << ", " << avgA_pT5l[8] << "}" << endl;
	Output << "double Aerr_pT5Lt[9] ={" << errA_pT5l[0] << ", " << errA_pT5l[1] << ", " << errA_pT5l[2] << ", " << errA_pT5l[3] << ", " << errA_pT5l[4] << ", " << errA_pT5l[5] << ", " << errA_pT5l[6] << ", " << errA_pT5l[7] << ", " << errA_pT5l[8] << "}" << endl;

	Output << "<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "double  A_pT1Gt[9] ={" << avgA_pT1g[0] << ", " << avgA_pT1g[1] << ", " << avgA_pT1g[2] << ", " << avgA_pT1g[3] << ", " << avgA_pT1g[4] << ", " << avgA_pT1g[5] << ", " << avgA_pT1g[6] << ", " << avgA_pT1g[7] << ", " << avgA_pT1g[8] << "}" << endl;
	Output << "double  Aerr_pT1Gt[9] ={" << errA_pT1g[0] << ", " << errA_pT1g[1] << ", " << errA_pT1g[2] << ", " << errA_pT1g[3] << ", " << errA_pT1g[4] << ", " << errA_pT1g[5] << ", " << errA_pT1g[6] << ", " << errA_pT1g[7] << ", " << errA_pT1g[8] << "}" << endl;
	Output << "double  A_pT2Gt[9] ={" << avgA_pT2g[0] << ", " << avgA_pT2g[1] << ", " << avgA_pT2g[2] << ", " << avgA_pT2g[3] << ", " << avgA_pT2g[4] << ", " << avgA_pT2g[5] << ", " << avgA_pT2g[6] << ", " << avgA_pT2g[7] << ", " << avgA_pT2g[8] << "}" << endl;
	Output << "double  Aerr_pT2Gt[9] ={" << errA_pT2g[0] << ", " << errA_pT2g[1] << ", " << errA_pT2g[2] << ", " << errA_pT2g[3] << ", " << errA_pT2g[4] << ", " << errA_pT2g[5] << ", " << errA_pT2g[6] << ", " << errA_pT2g[7] << ", " << errA_pT2g[8] << "}" << endl;
	Output << "double  A_pT3Gt[9] ={" << avgA_pT3g[0] << ", " << avgA_pT3g[1] << ", " << avgA_pT3g[2] << ", " << avgA_pT3g[3] << ", " << avgA_pT3g[4] << ", " << avgA_pT3g[5] << ", " << avgA_pT3g[6] << ", " << avgA_pT3g[7] << ", " << avgA_pT3g[8] << "}" << endl;
	Output << "double  Aerr_pT3Gt[9] ={" << errA_pT3g[0] << ", " << errA_pT3g[1] << ", " << errA_pT3g[2] << ", " << errA_pT3g[3] << ", " << errA_pT3g[4] << ", " << errA_pT3g[5] << ", " << errA_pT3g[6] << ", " << errA_pT3g[7] << ", " << errA_pT3g[8] << "}" << endl;
	Output << "double  A_pT4Gt[9] ={" << avgA_pT4g[0] << ", " << avgA_pT4g[1] << ", " << avgA_pT4g[2] << ", " << avgA_pT4g[3] << ", " << avgA_pT4g[4] << ", " << avgA_pT4g[5] << ", " << avgA_pT4g[6] << ", " << avgA_pT4g[7] << ", " << avgA_pT4g[8] << "}" << endl;
	Output << "double  Aerr_pT4Gt[9] ={" << errA_pT4g[0] << ", " << errA_pT4g[1] << ", " << errA_pT4g[2] << ", " << errA_pT4g[3] << ", " << errA_pT4g[4] << ", " << errA_pT4g[5] << ", " << errA_pT4g[6] << ", " << errA_pT4g[7] << ", " << errA_pT4g[8] << "}" << endl;
	Output << "double  A_pT5Gt[9] ={" << avgA_pT5g[0] << ", " << avgA_pT5g[1] << ", " << avgA_pT5g[2] << ", " << avgA_pT5g[3] << ", " << avgA_pT5g[4] << ", " << avgA_pT5g[5] << ", " << avgA_pT5g[6] << ", " << avgA_pT5g[7] << ", " << avgA_pT5g[8] << "}" << endl;
	Output << "double  Aerr_pT5Gt[9] ={" << errA_pT5g[0] << ", " << errA_pT5g[1] << ", " << errA_pT5g[2] << ", " << errA_pT5g[3] << ", " << errA_pT5g[4] << ", " << errA_pT5g[5] << ", " << errA_pT5g[6] << ", " << errA_pT5g[7] << ", " << errA_pT5g[8] << "}" << endl;

	Output << "<<<<<<<<<<<       Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "double  A_pT1Lt[9] ={" << WavgA_pT1l[0] << ", " << WavgA_pT1l[1] << ", " << WavgA_pT1l[2] << ", " << WavgA_pT1l[3] << ", " << WavgA_pT1l[4] << ", " << WavgA_pT1l[5] << ", " << WavgA_pT1l[6] << ", " << WavgA_pT1l[7] << ", " << WavgA_pT1l[8] << "}" << endl;
	Output << "double  AWerr_pT1Lt[9]={" << WerrA_pT1l[0] << ", " << WerrA_pT1l[1] << ", " << WerrA_pT1l[2] << ", " << WerrA_pT1l[3] << ", " << WerrA_pT1l[4] << ", " << WerrA_pT1l[5] << ", " << WerrA_pT1l[6] << ", " << WerrA_pT1l[7] << ", " << WerrA_pT1l[8] << "}" << endl;
	Output << "double  A_pT2Lt[9]={" << WavgA_pT2l[0] << ", " << WavgA_pT2l[1] << ", " << WavgA_pT2l[2] << ", " << WavgA_pT2l[3] << ", " << WavgA_pT2l[4] << ", " << WavgA_pT2l[5] << ", " << WavgA_pT2l[6] << ", " << WavgA_pT2l[7] << ", " << WavgA_pT2l[8] << "}" << endl;
	Output << "double  AWerr_pT2Lt[9] ={" << WerrA_pT2l[0] << ", " << WerrA_pT2l[1] << ", " << WerrA_pT2l[2] << ", " << WerrA_pT2l[3] << ", " << WerrA_pT2l[4] << ", " << WerrA_pT2l[5] << ", " << WerrA_pT2l[6] << ", " << WerrA_pT2l[7] << ", " << WerrA_pT2l[8] << "}" << endl;
	Output << "double A_pT3Lt[9] ={" << WavgA_pT3l[0] << ", " << WavgA_pT3l[1] << ", " << WavgA_pT3l[2] << ", " << WavgA_pT3l[3] << ", " << WavgA_pT3l[4] << ", " << WavgA_pT3l[5] << ", " << WavgA_pT3l[6] << ", " << WavgA_pT3l[7] << ", " << WavgA_pT3l[8] << "}" << endl;
	Output << "double AWerr_pT3Lt[9] ={" << WerrA_pT3l[0] << ", " << WerrA_pT3l[1] << ", " << WerrA_pT3l[2] << ", " << WerrA_pT3l[3] << ", " << WerrA_pT3l[4] << ", " << WerrA_pT3l[5] << ", " << WerrA_pT3l[6] << ", " << WerrA_pT3l[7] << ", " << WerrA_pT3l[8] << "}" << endl;
	Output << "double  A_pT4Lt[9] ={" << WavgA_pT4l[0] << ", " << WavgA_pT4l[1] << ", " << WavgA_pT4l[2] << ", " << WavgA_pT4l[3] << ", " << WavgA_pT4l[4] << ", " << WavgA_pT4l[5] << ", " << WavgA_pT4l[6] << ", " << WavgA_pT4l[7] << ", " << WavgA_pT4l[8] << "}" << endl;
	Output << "double AWerr_pT4Lt[9] ={" << WerrA_pT4l[0] << ", " << WerrA_pT4l[1] << ", " << WerrA_pT4l[2] << ", " << WerrA_pT4l[3] << ", " << WerrA_pT4l[4] << ", " << WerrA_pT4l[5] << ", " << WerrA_pT4l[6] << ", " << WerrA_pT4l[7] << ", " << WerrA_pT4l[8] << "}" << endl;
	Output << "double  A_pT5Lt[9] ={" << WavgA_pT5l[0] << ", " << WavgA_pT5l[1] << ", " << WavgA_pT5l[2] << ", " << WavgA_pT5l[3] << ", " << WavgA_pT5l[4] << ", " << WavgA_pT5l[5] << ", " << WavgA_pT5l[6] << ", " << WavgA_pT5l[7] << ", " << WavgA_pT5l[8] << "}" << endl;
	Output << "double AWerr_pT5Lt[9] ={" << WerrA_pT5l[0] << ", " << WerrA_pT5l[1] << ", " << WerrA_pT5l[2] << ", " << WerrA_pT5l[3] << ", " << WerrA_pT5l[4] << ", " << WerrA_pT5l[5] << ", " << WerrA_pT5l[6] << ", " << WerrA_pT5l[7] << ", " << WerrA_pT5l[8] << "}" << endl;

	Output << "<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "double  A_pT1Gt[9] ={" << WavgA_pT1g[0] << ", " << WavgA_pT1g[1] << ", " << WavgA_pT1g[2] << ", " << WavgA_pT1g[3] << ", " << WavgA_pT1g[4] << ", " << WavgA_pT1g[5] << ", " << WavgA_pT1g[6] << ", " << WavgA_pT1g[7] << ", " << WavgA_pT1g[8] << "}" << endl;
	Output << "double  AWerr_pT1Gt[9] ={" << WerrA_pT1g[0] << ", " << WerrA_pT1g[1] << ", " << WerrA_pT1g[2] << ", " << WerrA_pT1g[3] << ", " << WerrA_pT1g[4] << ", " << WerrA_pT1g[5] << ", " << WerrA_pT1g[6] << ", " << WerrA_pT1g[7] << ", " << WerrA_pT1g[8] << "}" << endl;
	Output << "double  A_pT2Gt[9] ={" << WavgA_pT2g[0] << ", " << WavgA_pT2g[1] << ", " << WavgA_pT2g[2] << ", " << WavgA_pT2g[3] << ", " << WavgA_pT2g[4] << ", " << WavgA_pT2g[5] << ", " << WavgA_pT2g[6] << ", " << WavgA_pT2g[7] << ", " << WavgA_pT2g[8] << "}" << endl;
	Output << "double  AWerr_pT2Gt[9] ={" << WerrA_pT2g[0] << ", " << WerrA_pT2g[1] << ", " << WerrA_pT2g[2] << ", " << WerrA_pT2g[3] << ", " << WerrA_pT2g[4] << ", " << WerrA_pT2g[5] << ", " << WerrA_pT2g[6] << ", " << WerrA_pT2g[7] << ", " << WerrA_pT2g[8] << "}" << endl;
	Output << "double  A_pT3Gt[9] ={" << WavgA_pT3g[0] << ", " << WavgA_pT3g[1] << ", " << WavgA_pT3g[2] << ", " << WavgA_pT3g[3] << ", " << WavgA_pT3g[4] << ", " << WavgA_pT3g[5] << ", " << WavgA_pT3g[6] << ", " << WavgA_pT3g[7] << ", " << WavgA_pT3g[8] << "}" << endl;
	Output << "double  AWerr_pT3Gt[9] ={" << WerrA_pT3g[0] << ", " << WerrA_pT3g[1] << ", " << WerrA_pT3g[2] << ", " << WerrA_pT3g[3] << ", " << WerrA_pT3g[4] << ", " << WerrA_pT3g[5] << ", " << WerrA_pT3g[6] << ", " << WerrA_pT3g[7] << ", " << WerrA_pT3g[8] << "}" << endl;
	Output << "double  A_pT4Gt[9] ={" << WavgA_pT4g[0] << ", " << WavgA_pT4g[1] << ", " << WavgA_pT4g[2] << ", " << WavgA_pT4g[3] << ", " << WavgA_pT4g[4] << ", " << WavgA_pT4g[5] << ", " << WavgA_pT4g[6] << ", " << WavgA_pT4g[7] << ", " << WavgA_pT4g[8] << "}" << endl;
	Output << "double  AWerr_pT4Gt[9] ={" << WerrA_pT4g[0] << ", " << WerrA_pT4g[1] << ", " << WerrA_pT4g[2] << ", " << WerrA_pT4g[3] << ", " << WerrA_pT4g[4] << ", " << WerrA_pT4g[5] << ", " << WerrA_pT4g[6] << ", " << WerrA_pT4g[7] << ", " << WerrA_pT4g[8] << "}" << endl;
	Output << "double  A_pT5Gt[9] ={" << WavgA_pT5g[0] << ", " << WavgA_pT5g[1] << ", " << WavgA_pT5g[2] << ", " << WavgA_pT5g[3] << ", " << WavgA_pT5g[4] << ", " << WavgA_pT5g[5] << ", " << WavgA_pT5g[6] << ", " << WavgA_pT5g[7] << ", " << WavgA_pT5g[8] << "}" << endl;
	Output << "double  AWerr_pT5Gt[9] ={" << WerrA_pT5g[0] << ", " << WerrA_pT5g[1] << ", " << WerrA_pT5g[2] << ", " << WerrA_pT5g[3] << ", " << WerrA_pT5g[4] << ", " << WerrA_pT5g[5] << ", " << WerrA_pT5g[6] << ", " << WerrA_pT5g[7] << ", " << WerrA_pT5g[8] << "}" << endl;

	// trigger bias ratio
	// 5 invariant mass bins eta > 0
	double ratioMGt1[9] = {0.9006, 0.9006, 0.9006, 0.9006, 0.9006, 0.9006, 0.9006, 0.9006, 0.9006}; // trigger bias
	double ratioMGt2[9] = {0.9156, 0.9156, 0.9156, 0.9156, 0.9156, 0.9156, 0.9156, 0.9156, 0.9156};
	double ratioMGt3[9] = {1.108, 1.108, 1.108, 1.108, 1.108, 1.108, 1.108, 1.108, 1.108};
	double ratioMGt4[9] = {1.083, 1.083, 1.083, 1.083, 1.083, 1.083, 1.083, 1.083, 1.083};
	double ratioMGt5[9] = {1.032, 1.032, 1.032, 1.032, 1.032, 1.032, 1.032, 1.032, 1.032};

	double ratioMLt1[9] = {0.9051, 0.9051, 0.9051, 0.9051, 0.9051, 0.9051, 0.9051, 0.9051, 0.9051}; // trigger bias
	double ratioMLt2[9] = {1.103, 1.103, 1.103, 1.103, 1.103, 1.103, 1.103, 1.103, 1.103};
	double ratioMLt3[9] = {1.085, 1.085, 1.085, 1.085, 1.085, 1.085, 1.085, 1.085, 1.085};
	double ratioMLt4[9] = {1.042, 1.042, 1.042, 1.042, 1.042, 1.042, 1.042, 1.042, 1.042};
	double ratioMLt5[9] = {1.058, 1.058, 1.058, 1.058, 1.058, 1.058, 1.058, 1.058, 1.058};

	/*//particle fractions for preliminary------------------------
	double pairpurityPtGt[5]={0.67307,0.817963,0.794931,0.765595,0.793264};//prelim PID fractions
	double pairpurityPtLt[5]={0.660449,0.8359,0.808052,0.739741,0.80209};
	double pairpurityMGt[5]={0.778802,0.792258,0.788996,0.779743,0.726311};
	double pairpurityMLt[5]={0.793783,0.804098,0.820309,0.790227,0.73937};
	double pairpurityEta[9]={0.777056,0.798529,0.800554,0.803539,0.802013,0.798089,0.79153,0.784865,0.760883};
	//-----------------------------------------------------------
	*/
	// final pid fractions after TOF cut---------------------------------
	// pion
	double pionPairPtGt[5] = {0.929242, 0.911722, 0.900156, 0.882439, 0.853673};
	double pionPairPtLt[5] = {0.945092, 0.925782, 0.912807, 0.895107, 0.866101};
	double pionPairMGt[5] = {0.878427, 0.887055, 0.895579, 0.878562, 0.849027};
	double pionPairMLt[5] = {0.892918, 0.901121, 0.907527, 0.891572, 0.863809};
	// kaon+proton
	double kpPairPtGt[5] = {0.00123106, 0.00188561, 0.00230285, 0.00293531, 0.00377407};
	double kpPairPtLt[5] = {0.000727659, 0.00131544, 0.00173617, 0.00229945, 0.00303562};
	double kpPairMGt[5] = {0.00298283, 0.00287837, 0.00237468, 0.00333972, 0.0051258};
	double kpPairMLt[5] = {0.00223246, 0.00218268, 0.00184161, 0.00263674, 0.00410077};
	// pair purity in eta bins
	double pionPairEta[9] = {0.862209, 0.901796, 0.902833, 0.894983, 0.886446, 0.890468, 0.896209, 0.891722, 0.844428};
	double kpPairEta[9] = {0.00365175, 0.00220525, 0.00219847, 0.00249915, 0.00277542, 0.00263735, 0.00246155, 0.00265875, 0.00498614};
	double electronPairEta[9] = {0.000116222, 1.11912e-05, 8.49341e-06, 1.53427e-05, 3.21867e-05, 2.39963e-05, 1.36267e-05, 1.67752e-05, 0.000106381};
	// only pion fractions for plotting
	double pairpurityPtGt[5] = {0.929242, 0.911722, 0.900156, 0.882439, 0.853673};
	double pairpurityPtLt[5] = {0.945092, 0.925782, 0.912807, 0.895107, 0.866101};
	double pairpurityMGt[5] = {0.878427, 0.887055, 0.895579, 0.878562, 0.849027};
	double pairpurityMLt[5] = {0.892918, 0.901121, 0.907527, 0.891572, 0.863809};
	double pairpurityEta[9] = {0.862209, 0.901796, 0.902833, 0.894983, 0.886446, 0.890468, 0.896209, 0.891722, 0.844428};

	//-------------------------------------------------------------------------
	/*//PID fractions  TOF only-------------------------------------------------
	///////// Signal fractions /////////////////
	//pion
	double pionPairPtGt[5]={0.995685,0.99413,0.992819,0.981137,0.951615};
	double pionPairPtLt[5]={0.997217,0.995865,0.994613,0.982959,0.954391};
	double pionPairMGt[5]={0.978901,0.98193,0.981563,0.980261,0.971465};
	double pionPairMLt[5]={0.981487,0.984854,0.98391,0.982889,0.973322};
	//kaon+proton
	double kpPairPtGt[5]={0,0,0,2.08895e-05,0.000131916};
	double kpPairPtLt[5]={0,0,0,1.90933e-05,0.000122819};
	double kpPairMGt[5]={1.88864e-05,3.60047e-05,3.40103e-05,4.00047e-05,7.8545e-05};
	double kpPairMLt[5]={1.49303e-05,2.6513e-05,2.64161e-05,3.2977e-05,7.32954e-05};
	// pair purity in eta bins
	double pionPairEta[9]={0.97593,0.98328,0.982627,0.978905,0.975671,0.978282,0.983113,0.983608,0.976228};
	double kpPairEta[9]={3.55615e-05,2.98475e-05,3.96205e-05,5.86822e-05,7.06504e-05,5.12718e-05,2.85773e-05,2.04411e-05,2.54495e-05};
	double electronPairEta[9]={3.69162e-05,8.40101e-06,5.82411e-06,8.56687e-06,1.40317e-05,1.3517e-05,9.81749e-06,1.35601e-05,4.67372e-05};
	*/
	//-----------------------------

	double biasMgt1[9] = {0};
	double biasMgt2[9] = {0};
	double biasMgt3[9] = {0};
	double biasMgt4[9] = {0};
	double biasMgt5[9] = {0};
	double biasMlt1[9] = {0};
	double biasMlt2[9] = {0};
	double biasMlt3[9] = {0};
	double biasMlt4[9] = {0};
	double biasMlt5[9] = {0};

	double pidMgt1[9] = {0};
	double pidMgt2[9] = {0};
	double pidMgt3[9] = {0};
	double pidMgt4[9] = {0};
	double pidMgt5[9] = {0};
	double pidMlt1[9] = {0};
	double pidMlt2[9] = {0};
	double pidMlt3[9] = {0};
	double pidMlt4[9] = {0};
	double pidMlt5[9] = {0};

	double combSysMgt1[9] = {0};
	double combSysMgt2[9] = {0};
	double combSysMgt3[9] = {0};
	double combSysMgt4[9] = {0};
	double combSysMgt5[9] = {0};
	double combSysMlt1[9] = {0};
	double combSysMlt2[9] = {0};
	double combSysMlt3[9] = {0};
	double combSysMlt4[9] = {0};
	double combSysMlt5[9] = {0};

	for (int pbin = 0; pbin < 9; pbin++)
	{
		biasMgt1[pbin] = (1 - ratioMGt1[pbin]) * TMath::Max(WavgA_pT1g[pbin], WerrA_pT1g[pbin]);
		biasMgt2[pbin] = (1 - ratioMGt2[pbin]) * TMath::Max(WavgA_pT2g[pbin], WerrA_pT2g[pbin]);
		biasMgt3[pbin] = (1 - ratioMGt3[pbin]) * TMath::Max(WavgA_pT3g[pbin], WerrA_pT3g[pbin]);
		biasMgt4[pbin] = (1 - ratioMGt4[pbin]) * TMath::Max(WavgA_pT4g[pbin], WerrA_pT4g[pbin]);
		biasMgt5[pbin] = (1 - ratioMGt5[pbin]) * TMath::Max(WavgA_pT5g[pbin], WerrA_pT5g[pbin]);

		biasMlt1[pbin] = (1 - ratioMLt1[pbin]) * TMath::Max(WavgA_pT1l[pbin], WerrA_pT1l[pbin]);
		biasMlt2[pbin] = (1 - ratioMLt2[pbin]) * TMath::Max(WavgA_pT2l[pbin], WerrA_pT2l[pbin]);
		biasMlt3[pbin] = (1 - ratioMLt3[pbin]) * TMath::Max(WavgA_pT3l[pbin], WerrA_pT3l[pbin]);
		biasMlt4[pbin] = (1 - ratioMLt4[pbin]) * TMath::Max(WavgA_pT4l[pbin], WerrA_pT4l[pbin]);
		biasMlt5[pbin] = (1 - ratioMLt5[pbin]) * TMath::Max(WavgA_pT5l[pbin], WerrA_pT5l[pbin]);

		pidMgt1[pbin] = (1 - pairpurityMGt[0]) * TMath::Max(WavgA_pT1g[pbin], WerrA_pT1g[pbin]);
		pidMgt2[pbin] = (1 - pairpurityMGt[1]) * TMath::Max(WavgA_pT2g[pbin], WerrA_pT2g[pbin]);
		pidMgt3[pbin] = (1 - pairpurityMGt[2]) * TMath::Max(WavgA_pT3g[pbin], WerrA_pT3g[pbin]);
		pidMgt4[pbin] = (1 - pairpurityMGt[3]) * TMath::Max(WavgA_pT4g[pbin], WerrA_pT4g[pbin]);
		pidMgt5[pbin] = (1 - pairpurityMGt[4]) * TMath::Max(WavgA_pT5g[pbin], WerrA_pT5g[pbin]);

		pidMlt1[pbin] = (1 - pairpurityMLt[0]) * TMath::Max(WavgA_pT1l[pbin], WerrA_pT1l[pbin]);
		pidMlt2[pbin] = (1 - pairpurityMLt[1]) * TMath::Max(WavgA_pT2l[pbin], WerrA_pT2l[pbin]);
		pidMlt3[pbin] = (1 - pairpurityMLt[2]) * TMath::Max(WavgA_pT3l[pbin], WerrA_pT3l[pbin]);
		pidMlt4[pbin] = (1 - pairpurityMLt[3]) * TMath::Max(WavgA_pT4l[pbin], WerrA_pT4l[pbin]);
		pidMlt5[pbin] = (1 - pairpurityMLt[4]) * TMath::Max(WavgA_pT5l[pbin], WerrA_pT5l[pbin]);

		combSysMgt1[pbin] = sqrt(pow(biasMgt1[pbin], 2) + pow(pidMgt1[pbin], 2));
		combSysMgt2[pbin] = sqrt(pow(biasMgt2[pbin], 2) + pow(pidMgt2[pbin], 2));
		combSysMgt3[pbin] = sqrt(pow(biasMgt3[pbin], 2) + pow(pidMgt3[pbin], 2));
		combSysMgt4[pbin] = sqrt(pow(biasMgt4[pbin], 2) + pow(pidMgt4[pbin], 2));
		combSysMgt5[pbin] = sqrt(pow(biasMgt5[pbin], 2) + pow(pidMgt5[pbin], 2));
		combSysMlt1[pbin] = sqrt(pow(biasMlt1[pbin], 2) + pow(pidMlt1[pbin], 2));
		combSysMlt2[pbin] = sqrt(pow(biasMlt2[pbin], 2) + pow(pidMlt2[pbin], 2));
		combSysMlt3[pbin] = sqrt(pow(biasMlt3[pbin], 2) + pow(pidMlt3[pbin], 2));
		combSysMlt4[pbin] = sqrt(pow(biasMlt4[pbin], 2) + pow(pidMlt4[pbin], 2));
		combSysMlt5[pbin] = sqrt(pow(biasMlt5[pbin], 2) + pow(pidMlt5[pbin], 2));
	}

	TGraphErrors *gr_AvsMgts[5];
	TGraphErrors *gr_AvsMlts[5];
	TGraphErrors *gr_AvsMgt[5];
	TGraphErrors *gr_AvsMlt[5];
	double minvAvg[9] = {0};
	double A_Mgt[9] = {0};
	double Aerr_Mgt[9] = {0};
	double A_Mlt[9] = {0};
	double Aerr_Mlt[9] = {0};
	double pidSysGt[9] = {0};
	double pidSysLt[9] = {0};
	double trigSysGt[9] = {0};
	double trigSysLt[9] = {0};

	TLegend *leg;
	TLine *line[5];
	TLatex tex[5];

	double finalSysgt[9] = {0};
	double finalSyslt[9] = {0};
	double errx[9] = {0};
	double errxps[9] = {0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035};
	double errxts[9] = {0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035};
	// print values in tabular format
	ofstream ftable;
	ftable.open("ifftable_AutVsMinvTPCOrTOF.txt");
	ftable << "//values are printed in the following order:\n";
	ftable << "//1.Bin range 2.<p_T> 3.<M> 4.A_UT 5.Sigma_Stat 6.Sigma_PID 7.Sigma_TriggerBias 8.Sigma_Total \n";
	ftable << "// 9.A_UT 10.Sigma_Stat 11.Sigma_PID 12.Sigma_TriggerBias 13.Sigma_Total " << endl;
	ftable << "//(Columns 4-8 for eta>0 and columns 9-13 for eta<0)" << endl;
	ftable << "\n"
		   << endl;

	TCanvas *myCanAvg = new TCanvas("myCanAvg", "myCanAvg", 150, 10, 1100, 700);
	gStyle->SetOptStat(0);
	myCanAvg->Divide(3, 2, 0, 0);
	for (int pad = 0; pad < 6; pad++)
	{
		myCanAvg->cd(pad + 1);
		setPad(pad);
		if (pad == 5)
		{
			TLatex text;
			text.SetTextSize(0.08);
			text.DrawLatex(0.1, 0.7, "#font[22]{#color[2]{STAR Preliminary 2015}}");
			text.SetTextSize(0.06);
			text.DrawLatex(0.1, 0.6, "#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 200 GeV}");
			text.DrawLatex(0.1, 0.6, "#font[22]{p^{#uparrow} + p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 200 GeV}");
			text.DrawLatex(0.1, 0.5, "#font[22]{#pm 3% scale uncertainty from beam}");
			text.DrawLatex(0.1, 0.4, "#font[22]{polarization (not shown)}");
		}
		else
		{
			gPad->SetGrid(0, 0);
			if (pad == 0)
				for (int mbin = 0; mbin < 9; mbin++)
				{
					minvAvg[mbin] = Minv1[mbin];
					A_Mgt[mbin] = WavgA_pT1g[mbin];
					Aerr_Mgt[mbin] = WerrA_pT1g[mbin];
					A_Mlt[mbin] = WavgA_pT1l[mbin];
					Aerr_Mlt[mbin] = WerrA_pT1l[mbin];
					finalSysgt[mbin] = combSysMgt1[mbin];
					finalSyslt[mbin] = combSysMlt1[mbin];
					pidSysGt[mbin] = pidMgt1[mbin];
					trigSysGt[mbin] = biasMgt1[mbin];
					pidSysLt[mbin] = pidMlt1[mbin];
					trigSysLt[mbin] = biasMlt1[mbin];
				}
			if (pad == 1)
				for (int mbin = 0; mbin < 9; mbin++)
				{
					minvAvg[mbin] = Minv2[mbin];
					A_Mgt[mbin] = WavgA_pT2g[mbin];
					Aerr_Mgt[mbin] = WerrA_pT2g[mbin];
					A_Mlt[mbin] = WavgA_pT2l[mbin];
					Aerr_Mlt[mbin] = WerrA_pT2l[mbin];
					finalSysgt[mbin] = combSysMgt2[mbin];
					finalSyslt[mbin] = combSysMlt2[mbin];
					pidSysGt[mbin] = pidMgt2[mbin];
					trigSysGt[mbin] = biasMgt2[mbin];
					pidSysLt[mbin] = pidMlt2[mbin];
					trigSysLt[mbin] = biasMlt2[mbin];
				}
			if (pad == 2)
				for (int mbin = 0; mbin < 9; mbin++)
				{
					minvAvg[mbin] = Minv3[mbin];
					A_Mgt[mbin] = WavgA_pT3g[mbin];
					Aerr_Mgt[mbin] = WerrA_pT3g[mbin];
					A_Mlt[mbin] = WavgA_pT3l[mbin];
					Aerr_Mlt[mbin] = WerrA_pT3l[mbin];
					finalSysgt[mbin] = combSysMgt3[mbin];
					finalSyslt[mbin] = combSysMlt3[mbin];
					pidSysGt[mbin] = pidMgt3[mbin];
					trigSysGt[mbin] = biasMgt3[mbin];
					pidSysLt[mbin] = pidMlt3[mbin];
					trigSysLt[mbin] = biasMlt3[mbin];
				}
			if (pad == 3)
				for (int mbin = 0; mbin < 9; mbin++)
				{
					minvAvg[mbin] = Minv4[mbin];
					A_Mgt[mbin] = WavgA_pT4g[mbin];
					Aerr_Mgt[mbin] = WerrA_pT4g[mbin];
					A_Mlt[mbin] = WavgA_pT4l[mbin];
					Aerr_Mlt[mbin] = WerrA_pT4l[mbin];
					finalSysgt[mbin] = combSysMgt4[mbin];
					finalSyslt[mbin] = combSysMlt4[mbin];
					pidSysGt[mbin] = pidMgt4[mbin];
					trigSysGt[mbin] = biasMgt4[mbin];
					pidSysLt[mbin] = pidMlt4[mbin];
					trigSysLt[mbin] = biasMlt4[mbin];
				}
			if (pad == 4)
				for (int mbin = 0; mbin < 9; mbin++)
				{
					minvAvg[mbin] = Minv5[mbin];
					A_Mgt[mbin] = WavgA_pT5g[mbin];
					Aerr_Mgt[mbin] = WerrA_pT5g[mbin];
					A_Mlt[mbin] = WavgA_pT5l[mbin];
					Aerr_Mlt[mbin] = WerrA_pT5l[mbin];
					finalSysgt[mbin] = combSysMgt5[mbin];
					finalSyslt[mbin] = combSysMlt5[mbin];
					pidSysGt[mbin] = pidMgt5[mbin];
					trigSysGt[mbin] = biasMgt5[mbin];
					pidSysLt[mbin] = pidMlt5[mbin];
					trigSysLt[mbin] = biasMlt5[mbin];
				}

			ftable << pT[pad] << " - " << pT[pad + 1] << " & ";
			for (int i = 0; i < 9; i++)
			{
				if (i == 0)
				{
					ftable << avgpTpairGt[pad] << " & " << minvAvg[i] << " & " << A_Mgt[i] << " & " << Aerr_Mgt[i] << " & " << pidSysGt[i] << " & " << trigSysGt[i] << " & " << finalSysgt[i] << " & " << A_Mlt[i] << " & " << Aerr_Mlt[i] << " & " << pidSysLt[i] << " & " << trigSysLt[i] << " & " << finalSyslt[i] << " \\\\ " << endl;
				}
				else
				{
					ftable << "\t & " << minvAvg[i] << " & " << A_Mgt[i] << " & " << Aerr_Mgt[i] << " & " << pidSysGt[i] << " & " << trigSysGt[i] << " & " << finalSysgt[i] << " & " << A_Mlt[i] << " & " << Aerr_Mlt[i] << " & " << pidSysLt[i] << " & " << trigSysLt[i] << " & " << finalSyslt[i] << " \\\\ " << endl;
				}
			}
			gr_AvsMgts[pad] = new TGraphErrors(9, minvAvg, A_Mgt, errxts, finalSysgt);
			gr_AvsMgts[pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}} ");
			gr_AvsMgts[pad]->SetTitle("");
			gr_AvsMgts[pad]->GetYaxis()->SetTitleOffset(1.2);
			gr_AvsMgts[pad]->GetYaxis()->SetTitleSize(0.06);
			gr_AvsMgts[pad]->GetYaxis()->SetLabelSize(0.06);
			gr_AvsMgts[pad]->GetYaxis()->SetLabelFont(22);
			gr_AvsMgts[pad]->GetXaxis()->SetTitle("#font[22]{M_{#pi^{+}#pi^{-}} (GeV/c^{2})}   ");
			gr_AvsMgts[pad]->GetXaxis()->SetTitleSize(0.06);
			gr_AvsMgts[pad]->GetXaxis()->SetLabelSize(0.06);
			gr_AvsMgts[pad]->GetXaxis()->SetLabelFont(22);
			gr_AvsMgts[pad]->GetXaxis()->SetLimits(0.2, 2.3);
			gr_AvsMgts[pad]->GetYaxis()->SetRangeUser(-0.015, 0.07);
			gr_AvsMgts[pad]->SetFillStyle(0000);
			gr_AvsMgts[pad]->SetLineColor(1);
			gr_AvsMgts[pad]->SetLineWidth(1);
			gr_AvsMgts[pad]->GetXaxis()->SetNdivisions(505);
			gr_AvsMgts[pad]->GetYaxis()->SetNdivisions(505);
			gr_AvsMgts[pad]->Draw("AP2");
			if (pad == 2)
			{
				gr_AvsMgts[pad]->GetYaxis()->SetTitle("");
				gr_AvsMgts[pad]->GetYaxis()->SetLabelSize(0);
			}
			if (pad == 1)
			{
				gr_AvsMgts[pad]->GetYaxis()->SetTitle("");
				gr_AvsMgts[pad]->GetXaxis()->SetTitle("");
				gr_AvsMgts[pad]->GetXaxis()->SetLabelSize(0);
			}
			if (pad == 0)
			{
				gr_AvsMgts[pad]->GetXaxis()->SetTitle("");
				gr_AvsMgts[pad]->GetXaxis()->SetLabelSize(0.0);
			}

			gr_AvsMgt[pad] = new TGraphErrors(9, minvAvg, A_Mgt, errx, Aerr_Mgt); // Backward average
			gr_AvsMgt[pad]->SetMarkerStyle(20);
			gr_AvsMgt[pad]->SetMarkerColor(2);
			gr_AvsMgt[pad]->SetLineColor(2);
			gr_AvsMgt[pad]->SetLineWidth(1);
			gr_AvsMgt[pad]->Draw("same P");

			gr_AvsMlts[pad] = new TGraphErrors(9, minvAvg, A_Mlt, errxts, finalSyslt); // Backward average
			gr_AvsMlts[pad]->SetFillStyle(0000);
			gr_AvsMlts[pad]->SetLineColor(1);
			gr_AvsMlts[pad]->SetLineWidth(1);
			gr_AvsMlts[pad]->Draw("same 2");

			gr_AvsMlt[pad] = new TGraphErrors(9, minvAvg, A_Mlt, errx, Aerr_Mlt); // Backward average
			gr_AvsMlt[pad]->SetMarkerStyle(20);
			gr_AvsMlt[pad]->SetMarkerColor(4);
			gr_AvsMlt[pad]->SetLineColor(4);
			gr_AvsMlt[pad]->SetLineWidth(1);
			gr_AvsMlt[pad]->Draw("same P");

			gPad->Update();

			line[pad] = new TLine(myCanAvg->cd(pad + 1)->GetUxmin(), 0., myCanAvg->cd(pad + 1)->GetUxmax(), 0.);
			line[pad]->SetLineStyle(2);
			line[pad]->SetLineWidth(1);
			line[pad]->Draw();

			tex[pad].SetTextSize(0.055);
			tex[pad].SetTextAlign(13);
			tex[pad].DrawLatex(1.3, 0.065, Form("#font[22]{#color[1]{< p_{T} > = %3.2lf GeV/c}}", avgpTpairGt[pad])); // in %3.2lf, 3 sets the number of digits and 2 sets the numbers after decimal

			if (pad == 0)
			{
				leg = new TLegend(0.24, .6, 0.6, 0.8);
				leg->AddEntry(gr_AvsMgt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} > 0}", "lp");
				leg->AddEntry(gr_AvsMlt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} < 0}", "lp");
				leg->AddEntry(gr_AvsMgts[pad], " #font[22]{ Syst. Error }", "f");
				leg->SetTextSize(0.05);
				leg->Draw();
			}
			gPad->Update();
		}
	}
	myCanAvg->SaveAs("AsymVsMinvTOFcutTPCOrTOF.pdf");

	/*	//Forward Asymmetry BLUE and YELLOW
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
			pt1 = new TGraphErrors(9,Minv1,A_pT1bg,0,deltaA_pT1bg);
			pt1->GetYaxis()-> SetTitle("A_{UT}");
			pt1->SetTitle("");
			//pt1->GetYaxis()-> CenterTitle();
			//pt1->GetYaxis()-> SetTitleSize(.06);
			pt1->GetYaxis()->SetTitleOffset(1.66);
			//pt1->GetYaxis()->SetLabelSize(0.06);
			//pt1->GetXaxis()->SetLabelSize(0);
			pt1->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			//pt1-> GetXaxis()-> CenterTitle();
			//pt1-> SetTitle(" 2.85 < p_{T} < 3.65)");
			pt1-> SetMarkerStyle(20);
			pt1-> SetMarkerColor(4);
			pt1-> GetXaxis()->SetLimits(0.2,2.3);
			pt1-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt1->Draw("AP");
			myCan->Update();
			pt1s= new TGraphErrors(9,Minv1,A_pT1yg,0,deltaA_pT1yg);
			pt1s->SetMarkerStyle(20);
			pt1s->SetMarkerColor(2);
			pt1s->Draw("same P");
			TLegend *leg1 = new TLegend(0.7,.75, 0.8, 0.85);
			leg1->AddEntry(pt1, "#eta > 0 , BLUE", "lp");
			leg1->AddEntry(pt1s, "#eta > 0, YELLOW", "lp");
			leg1->SetTextSize(0.04);
			leg1->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[0]));


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
			pt2 = new TGraphErrors(9,Minv2,A_pT2bg,0,deltaA_pT2bg);
			pt2->GetYaxis()-> SetTitle("A_{UT}");
			//pt2->GetYaxis()-> CenterTitle();
			pt2->GetYaxis()->SetTitleOffset(1.5);
			pt2->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			//pt2->GetXaxis()-> CenterTitle();
			//pt2->GetXaxis()->SetLabelSize(0);
			//pt2->GetYaxis()->SetLabelSize(0);
			//pt2-> SetTitle("3.65 < p_{T} < 4.20");
			pt2->SetTitle("");
			pt2-> SetMarkerStyle(20);
			pt2-> SetMarkerColor(4);
			pt2-> GetXaxis()->SetLimits(0.2,2.3);
			pt2-> GetYaxis()->SetRangeUser(-0.015, .065);
			pt2->Draw("AP");
			myCan->Update();
			pt2s= new TGraphErrors(9,Minv2,A_pT2yg,0,deltaA_pT2yg);
			pt2s->SetMarkerStyle(20);
			pt2s->SetMarkerColor(2);
			pt2s->Draw("same P");
			TLegend *leg2 = new TLegend(0.65,0.75, 0.8, 0.85);
			leg2->AddEntry(pt2, "#eta > 0, BLUE", "lp");
			leg2->AddEntry(pt2s, "#eta > 0, YELLOW", "lp");
			leg2->SetTextSize(0.04);
			leg2->Draw();
			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[1]));
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
			pt3 = new TGraphErrors(9,Minv3,A_pT3bg,0,deltaA_pT3bg);
			pt3->GetYaxis()-> SetTitle("A_{UT}");
			//pt3->GetYaxis()->CenterTitle();
			pt3->GetYaxis()->SetTitleOffset(1.5);
			//pt3->GetXaxis()->SetLabelSize(0);
			//pt3->GetYaxis()->SetLabelSize(0);
			pt3->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt3->SetTitle("");
			//pt3->GetXaxis()->CenterTitle();
			//pt3-> SetTitle(" 4.20 < p_{T} < 4.95");
			pt3-> SetMarkerStyle(20);
			pt3-> SetMarkerColor(4);
			pt3-> GetXaxis()->SetLimits(0.2,2.3);
			pt3-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt3->Draw("AP");
			myCan->Update();
			pt3s = new TGraphErrors(9, Minv3, A_pT3yg, 0, deltaA_pT3yg);
			pt3s->SetMarkerStyle(20);
			pt3s->SetMarkerColor(2);
			pt3s->Draw("same P");
			TLegend *leg3 = new TLegend(0.65,0.75, 0.8, 0.85);
			leg3->AddEntry(pt3, "#eta > 0, BLUE", "lp");
			leg3->AddEntry(pt3s, "#eta > 0, YELLOW", "lp");
			leg3->SetTextSize(0.04);
			leg3->Draw();
			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[2]));
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
			pt4= new TGraphErrors(9,Minv4,A_pT4bg,0,deltaA_pT4bg);
			pt4->GetYaxis()-> SetTitle("A_{UT}");
			//pt4->GetYaxis()-> CenterTitle();
			pt4->GetYaxis()->SetTitleOffset(1.5);
			pt4->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			//pt4->GetXaxis()-> CenterTitle();
			//pt4->GetXaxis()-> SetLabelSize(0);
			pt4->SetTitle("");
			//pt4->SetTitle("4.95 < p_{T} < 6.00");
			pt4-> SetMarkerStyle(20);
			pt4-> SetMarkerColor(4);
			pt4-> GetXaxis()->SetLimits(0.2,2.);
			pt4-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			//pt4->GetYaxis()-> SetTitleSize(.05);
			//pt4->GetXaxis()-> SetTitleSize(.05);
					pt4->GetYaxis()->SetTitleOffset(1.6);
					//pt4->GetYaxis()->SetLabelSize(0.05);
					//pt4->GetXaxis()->SetLabelSize(0.05);

			pt4->Draw("AP");
			myCan->Update();
			pt4s = new TGraphErrors(9, Minv4, A_pT4yg, 0, deltaA_pT4yg);
			pt4s->SetMarkerStyle(20);
			pt4s->SetMarkerColor(2);
			pt4s->Draw("same P");
			TLegend *leg4 = new TLegend(0.65,0.75, 0.8, 0.85);
			leg4->AddEntry(pt4, "#eta > 0, BLUE", "lp");
			leg4->AddEntry(pt4s, "#eta > 0, YELLOW", "lp");
			leg4->SetTextSize(0.04);
			leg4->Draw();
			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[3]));
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
			pt5 = new TGraphErrors(9,Minv5,A_pT5bg,0,deltaA_pT5bg);
			pt5->GetYaxis()-> SetTitle("A_{UT}");
			pt5->GetYaxis()-> CenterTitle();
			pt5->SetTitle("");
			pt5->GetYaxis()->SetTitleOffset(1.5);
			pt5->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			//pt5->GetXaxis()-> CenterTitle();
			//pt5->GetXaxis()->SetLabelSize(0.05);
			//pt5->GetXaxis()->SetTitleSize(0.05);
			//pt5-> SetTitle(" 6.00< p_{T} <20)");
			//pt5->GetYaxis()-> SetTitleSize(.06);
					//pt5->GetYaxis()->SetTitleOffset(1.5);
					//pt5->GetYaxis()->SetLabelSize(0.05);
			pt5-> SetMarkerStyle(20);
			pt5-> SetMarkerColor(4);
			pt5-> GetXaxis()->SetLimits(0.2,2.2);
			pt5-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt5->Draw("AP");
			myCan->Update();
			pt5s = new TGraphErrors(9, Minv5, A_pT5yg, 0, deltaA_pT5yg);
			pt5s->SetMarkerStyle(20);
			pt5s->SetMarkerColor(2);
			pt5s->Draw("same P");
			TLegend *leg5 = new TLegend(0.65,0.75, 0.8, 0.85);
			leg5->AddEntry(pt5, "#eta > 0, BLUE", "lp");
			leg5->AddEntry(pt5s, "#eta > 0, YELLOW", "lp");
			leg5->SetTextSize(0.04);
			leg5->Draw();
			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[4]));
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
			pt1AF = new TGraphErrors(9,Minv1,avgA_pT1g,0,errA_pT1g);//Forward average
			pt1AF->GetYaxis()-> SetTitle("A_{UT}");
			pt1AF->SetTitle("");
			pt1AF->GetYaxis()->SetTitleOffset(1.66);
			pt1AF->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt1AF-> SetMarkerStyle(20);
			pt1AF-> SetMarkerColor(2);
			pt1AF-> GetXaxis()->SetLimits(0.2,2.3);
			pt1AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt1AF->Draw("AP");
			pt1AB = new TGraphErrors(9,Minv1,avgA_pT1l,0,errA_pT1l);//Backward average
			pt1AB-> SetMarkerStyle(20);
			pt1AB-> SetMarkerColor(4);
			pt1AB->Draw("same P");
			myCanA->Update();
			TLegend *leg1A = new TLegend(0.65,.77, 0.8, 0.87);
			leg1A->AddEntry(pt1AF, "#eta > 0", "lp");
			leg1A->AddEntry(pt1AB, "#eta < 0", "lp");
			leg1A->SetTextSize(0.04);
			leg1A->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[0]));


			myCanA->Update();
			TLine *line1A=  new TLine(myCanA->cd(1)->GetUxmin(),0.,myCanA->cd(1)->GetUxmax(),0.);
			line1A->SetLineStyle(2);
			line1A->Draw();

			myCanA -> cd(2);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt2AF = new TGraphErrors(9,Minv2,avgA_pT2g,0,errA_pT2g);//Forward average
			pt2AF->GetYaxis()-> SetTitle("A_{UT}");
			pt2AF->SetTitle("");
			pt2AF->GetYaxis()->SetTitleOffset(1.66);
			pt2AF->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt2AF-> SetMarkerStyle(20);
			pt2AF-> SetMarkerColor(2);
			pt2AF-> GetXaxis()->SetLimits(0.2,2.3);
			pt2AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt2AF->Draw("AP");
			pt2AB = new TGraphErrors(9,Minv2,avgA_pT2l,0,errA_pT2l);//Backward average
			pt2AB-> SetMarkerStyle(20);
			pt2AB-> SetMarkerColor(4);
			pt2AB->Draw("same P");
			gPad->Update();
			TLegend *leg2A = new TLegend(0.65,.77, 0.8, 0.87);
			leg2A->AddEntry(pt2AF, "#eta > 0", "lp");
			leg2A->AddEntry(pt2AB, "#eta < 0", "lp");
			leg2A->SetTextSize(0.04);
			leg2A->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[1]));
			gPad->Update();
			TLine *line2A=  new TLine(myCanA->cd(2)->GetUxmin(),0.,myCanA->cd(2)->GetUxmax(),0.);
			line2A->SetLineStyle(2);
			line2A->Draw();
			gPad->Update();

			myCanA -> cd(3);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt3AF = new TGraphErrors(9,Minv3,avgA_pT3g,0,errA_pT3g);
			pt3AF->GetYaxis()-> SetTitle("A_{UT}");
			pt3AF->SetTitle("");
			pt3AF->GetYaxis()->SetTitleOffset(1.66);
			pt3AF->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt3AF-> SetMarkerStyle(20);
			pt3AF-> SetMarkerColor(2);
			pt3AF-> GetXaxis()->SetLimits(0.2,2.3);
			pt3AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt3AF->Draw("AP");
			pt3AB = new TGraphErrors(9,Minv3,avgA_pT3l,0,errA_pT3l);
			pt3AB-> SetMarkerStyle(20);
			pt3AB-> SetMarkerColor(4);
			pt3AB->Draw("same P");
			gPad->Update();
			TLegend *leg3A = new TLegend(0.65,.77, 0.8, 0.87);
			leg3A->AddEntry(pt3AF, " #eta > 0", "lp");
			leg3A->AddEntry(pt3AB, " #eta < 0", "lp");
			leg3A->SetTextSize(0.04);
			leg3A->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[2]));
			gPad->Update();
			TLine *line3A=  new TLine(myCanA->cd(3)->GetUxmin(),0.,myCanA->cd(3)->GetUxmax(),0.);
			line3A->SetLineStyle(2);
			line3A->Draw();
			gPad->Update();

			myCanA -> cd(4);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt4AF = new TGraphErrors(9,Minv4,avgA_pT4g,0,errA_pT4g);
			pt4AF->GetYaxis()-> SetTitle("A_{UT}");
			pt4AF->SetTitle("");
			pt4AF->GetYaxis()->SetTitleOffset(1.66);
			pt4AF->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt4AF-> SetMarkerStyle(20);
			pt4AF-> SetMarkerColor(2);
			pt4AF-> GetXaxis()->SetLimits(0.2,2.3);
			pt4AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt4AF->Draw("AP");
			pt4AB = new TGraphErrors(9,Minv4,avgA_pT4l,0,errA_pT4l);
			pt4AB-> SetMarkerStyle(20);
			pt4AB-> SetMarkerColor(4);
			pt4AB->Draw("same P");
			gPad->Update();
			TLegend *leg4A = new TLegend(0.65,.77, 0.8, 0.87);
			leg4A->AddEntry(pt4AF, "#eta > 0", "lp");
			leg4A->AddEntry(pt4AB, "#eta < 0", "lp");
			leg4A->SetTextSize(0.04);
			leg4A->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[3]));
			gPad->Update();

			TLine *line4A=  new TLine(myCanA->cd(4)->GetUxmin(),0.,myCanA->cd(4)->GetUxmax(),0.);
			line4A->SetLineStyle(2);
			line4A->Draw();
			gPad->Update();

			myCanA -> cd(5);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt5AF = new TGraphErrors(9,Minv5,avgA_pT5g,0,errA_pT5g);
			pt5AF->GetYaxis()-> SetTitle("A_{UT}");
			pt5AF->SetTitle("");
			pt5AF->GetYaxis()->SetTitleOffset(1.66);
			pt5AF->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt5AF-> SetMarkerStyle(20);
			pt5AF-> SetMarkerColor(2);
			pt5AF-> GetXaxis()->SetLimits(0.2,2.3);
			pt5AF-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt5AF->Draw("AP");
			pt5AB= new TGraphErrors(9,Minv5,avgA_pT5l,0,errA_pT5l);
			pt5AB-> SetMarkerStyle(20);
			pt5AB-> SetMarkerColor(4);
			pt5AB-> Draw("same P");
			gPad->Update();
			TLegend *leg5A = new TLegend(0.65,.77, 0.8, 0.87);
			leg5A->AddEntry(pt5AF, "#eta > 0", "lp");
			leg5A->AddEntry(pt5AB, "#eta < 0", "lp");
			leg5A->SetTextSize(0.04);
			leg5A->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[4]));
			myCanA->Update();

			TLine *line5A=  new TLine(myCanA->cd(5)->GetUxmin(),0.,myCanA->cd(5)->GetUxmax(),0.);
			line5A->SetLineStyle(2);
			line5A->Draw();
			gPad->Update();

		myCan->SaveAs("AsymVsMinv_EtaGt0_ForwardTPCOrTOF.pdf");
		myCanA->SaveAs("AsymVsMinv_EtaGt0_AverageTPCOrTOF.pdf");

		myCan->SaveAs("AsymVsMinv_EtaGt0_ForwardTPCOrTOF.pdf");
		myCanA->SaveAs("AsymVsMinv_EtaGt0_AverageTPCOrTOF.pdf");
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
			pt1Bl = new TGraphErrors(9,Minv1,A_pT1bl,0,deltaA_pT1bl);
			pt1Bl->GetYaxis()-> SetTitle("A_{UT}");
			pt1Bl->SetTitle("");
			pt1Bl->GetYaxis()->SetTitleOffset(1.66);
			pt1Bl->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt1Bl-> SetMarkerStyle(20);
			pt1Bl-> SetMarkerColor(4);
			pt1Bl-> GetXaxis()->SetLimits(0.2,2.3);
			pt1Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt1Bl->Draw("AP");
			pt1Yl = new TGraphErrors(9,Minv1,A_pT1yl,0,deltaA_pT1yl);
			pt1Yl-> SetMarkerStyle(20);
			pt1Yl-> SetMarkerColor(2);
			pt1Yl-> Draw("same P");
			gPad->Update();
			TLegend *lBl = new TLegend(0.65,.77, 0.8, 0.87);
			lBl->AddEntry(pt1Bl, "BLUE, #eta < 0", "lp");
			lBl->AddEntry(pt1Yl, "YELLOW, #eta < 0", "lp");
			lBl->SetTextSize(0.04);
			lBl->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[0]));
			gPad->Update();
			TLine *line1A=  new TLine(cB->cd(1)->GetUxmin(),0.,cB->cd(1)->GetUxmax(),0.);
			line1A->SetLineStyle(2);
			line1A->Draw();
			gPad->Update();

			cB -> cd(2);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt2Bl = new TGraphErrors(9,Minv2,A_pT2bl,0,deltaA_pT2bl);
			pt2Bl->GetYaxis()-> SetTitle("A_{UT}");
			pt2Bl->SetTitle("");
			pt2Bl->GetYaxis()->SetTitleOffset(1.66);
			pt2Bl->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt2Bl-> SetMarkerStyle(20);
			pt2Bl-> SetMarkerColor(4);
			pt2Bl-> GetXaxis()->SetLimits(0.2,2.3);
			pt2Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt2Bl->Draw("AP");
			pt2Yl = new TGraphErrors(9,Minv2, A_pT2yl,0,deltaA_pT2yl);
			pt2Yl-> SetMarkerStyle(20);
			pt2Yl-> SetMarkerColor(2);
			pt2Yl-> Draw("same P");
			gPad->Update();
			TLegend *lBl = new TLegend(0.65,.77, 0.8, 0.87);
			lBl->AddEntry(pt2Bl, "BLUE, #eta < 0", "lp");
			lBl->AddEntry(pt2Yl, "YELLOW, #eta < 0", "lp");
			lBl->SetTextSize(0.04);
			lBl->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[1]));
			gPad->Update();
			TLine *line1A=  new TLine(cB->cd(2)->GetUxmin(),0.,cB->cd(2)->GetUxmax(),0.);
			line1A->SetLineStyle(2);
			line1A->Draw();
			gPad->Update();

			cB -> cd(3);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt3Bl = new TGraphErrors(9,Minv3,A_pT3bl,0,deltaA_pT3bl);
			pt3Bl->GetYaxis()-> SetTitle("A_{UT}");
			pt3Bl->SetTitle("");
			pt3Bl->GetYaxis()->SetTitleOffset(1.66);
			pt3Bl->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt3Bl-> SetMarkerStyle(20);
			pt3Bl-> SetMarkerColor(4);
			pt3Bl-> GetXaxis()->SetLimits(0.2,2.3);
			pt3Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt3Bl->Draw("AP");
			pt3Yl = new TGraphErrors(9,Minv3,A_pT3yl,0,deltaA_pT3yl);
			pt3Yl-> SetMarkerStyle(20);
			pt3Yl-> SetMarkerColor(2);
			pt3Yl-> Draw("same P");
			gPad->Update();
			TLegend *lBl = new TLegend(0.65,.77, 0.8, 0.87);
			lBl->AddEntry(pt3Bl, "BLUE, #eta < 0", "lp");
			lBl->AddEntry(pt3Yl, "YELLOW, #eta < 0", "lp");
			lBl->SetTextSize(0.04);
			lBl->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[2]));
			gPad->Update();
			TLine *line1A=  new TLine(cB->cd(3)->GetUxmin(),0.,cB->cd(3)->GetUxmax(),0.);
			line1A->SetLineStyle(2);
			line1A->Draw();
			gPad->Update();

			cB -> cd(4);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt4Bl = new TGraphErrors(9,Minv4,A_pT4bl,0,deltaA_pT4bl);
			pt4Bl->GetYaxis()-> SetTitle("A_{UT}");
			pt4Bl->SetTitle("");
			pt4Bl->GetYaxis()->SetTitleOffset(1.66);
			pt4Bl->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt4Bl-> SetMarkerStyle(20);
			pt4Bl-> SetMarkerColor(4);
			pt4Bl-> GetXaxis()->SetLimits(0.2,2.3);
			pt4Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt4Bl->Draw("AP");
			pt4Yl = new TGraphErrors(9,Minv4,A_pT4yl,0,deltaA_pT4yl);
			pt4Yl-> SetMarkerStyle(20);
			pt4Yl-> SetMarkerColor(2);
			pt4Yl-> Draw("same P");
			gPad->Update();
			TLegend *lBl = new TLegend(0.65,.77, 0.8, 0.87);
			lBl->AddEntry(pt4Bl, "BLUE, #eta < 0", "lp");
			lBl->AddEntry(pt4Yl, "YELLOW, #eta < 0", "lp");
			lBl->SetTextSize(0.04);
			lBl->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[3]));
			gPad->Update();
			TLine *line1A=  new TLine(cB->cd(4)->GetUxmin(),0.,cB->cd(4)->GetUxmax(),0.);
			line1A->SetLineStyle(2);
			line1A->Draw();
			gPad->Update();

			cB -> cd(5);
			gPad->SetGrid(0,0);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			pt5Bl = new TGraphErrors(9,Minv5,A_pT5bl,0,deltaA_pT5bl);
			pt5Bl->GetYaxis()-> SetTitle("A_{UT}");
			pt5Bl->SetTitle("");
			pt5Bl->GetYaxis()->SetTitleOffset(1.66);
			pt5Bl->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/c^{2})");
			pt5Bl-> SetMarkerStyle(20);
			pt5Bl-> SetMarkerColor(4);
			pt5Bl-> GetXaxis()->SetLimits(0.2,2.3);
			pt5Bl-> GetYaxis()->SetRangeUser(-0.015, 0.065);
			pt5Bl->Draw("AP");
			pt5Yl = new TGraphErrors(9,Minv5,A_pT5yl,0,deltaA_pT5yl);
			pt5Yl-> SetMarkerStyle(20);
			pt5Yl-> SetMarkerColor(2);
			pt5Yl-> Draw("same P");
			gPad->Update();
			TLegend *lBl = new TLegend(0.65,.77, 0.8, 0.87);
			lBl->AddEntry(pt5Bl, "BLUE, #eta < 0", "lp");
			lBl->AddEntry(pt5Yl, "YELLOW, #eta < 0", "lp");
			lBl->SetTextSize(0.04);
			lBl->Draw();

			TLatex latex;
							latex.SetTextSize(0.04);
							latex.SetTextAlign(13);
							latex.DrawLatex(0.25,0.062,Form("#font[12]{#color[2]{< p_{T} > = %g GeV/c}}",avgpTpairGt[4]));
			gPad->Update();
			TLine *line1A=  new TLine(cB->cd(5)->GetUxmin(),0.,cB->cd(5)->GetUxmax(),0.);
			line1A->SetLineStyle(2);
			line1A->Draw();
			gPad->Update();
		cB->Update();
		cB->SaveAs("AsymmetryVsMinv_BackwardTPCOrTOF.pdf");

	chi2f->Write();
	chi2f->Close();
	*/
} // main
// reltative difference between TpcOrTof and TofOnly data for PID systematic
double relDiff(double aTpc, double aTof)
{
	double reldiff;
	reldiff = (aTpc - aTof) / aTpc;
	return reldiff;
}

// relative error calculation
double calErr(double aTpc, double eTpc, double aTof, double eTof)
{
	double relErr;
	relErr = (1 / (aTpc * aTpc)) * sqrt(pow(aTof * eTpc, 2) + pow(aTpc * eTof, 2));
	return relErr;
}
// get Tof only asymmetry
double *getTofAsym(TFile *ifile)
{
	double asym, err;
	int i = 0;
	static double asymTof[9];
	static double errTof[9];
	if (!ifile.isOpen())
	{
		cout << "Tof asymmetry file not found!"
	}
	else
	{
		while (ifile >> asym >> err)
		{
			asymTof[i++] = asym;
			errTof[i++] = err;
		}
	}
	return asymTof;
}
// get Tof only error
double *getTofErr(TFile *ifile)
{
	double asym, err;
	int i = 0;
	static double asymTof[9];
	static double errTof[9];
	if (!ifile.isOpen())
	{
		cout << "Tof asymmetry file not found!"
	}
	else
	{
		while (ifile >> asym >> err)
		{
			asymTof[i++] = asym;
			errTof[i++] = err;
		}
	}
	return errTof;
}
// set pad for  canvas
void setPad(int i)
{

	if (i == 2)
	{
		gPad->SetTopMargin(0.137);
		gPad->SetLeftMargin(0.02);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.02);
		gPad->SetPad(.6634, 0.4197, 0.97, .9897);
		gPad->SetFillStyle(4000);
		gPad->Update();
	}
	else if (i == 5)
	{
		gPad->SetTopMargin(0.4);
		gPad->SetLeftMargin(0.01);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.02);
		gPad->SetPad(.6634, 0.01, 0.97, .49);
		gPad->SetFillStyle(4000);
		gPad->SetFrameFillColor(0);
		gPad->SetFrameFillStyle(0);
		gPad->SetFrameBorderMode(0);

		gPad->Update();
	}
	else if (i == 1)
	{
		gPad->SetTopMargin(0.16);
		gPad->SetLeftMargin(0.0);
		gPad->SetBottomMargin(0.01);
		gPad->SetRightMargin(0.0);
		gPad->SetPad(.37, 0.5, 0.67, .99);
		gPad->Update();
	}
	else if (i == 0)
	{
		gPad->SetTopMargin(0.16);
		gPad->SetLeftMargin(0.2);
		gPad->SetBottomMargin(0.01);
		gPad->SetRightMargin(0.0);
		gPad->SetPad(0., 0.5, 0.37, .99);
	}
	else if (i == 3)
	{
		gPad->SetTopMargin(0.0);
		gPad->SetLeftMargin(0.2);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.0);
		gPad->SetPad(0.0, 0.02052, 0.37, .5052);
	}
	else if (i == 4)
	{
		gPad->SetTopMargin(0.0);
		gPad->SetLeftMargin(0.0);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.0035);
		gPad->SetPad(.37, 0.02052, 0.67, .5052);
	}
}
