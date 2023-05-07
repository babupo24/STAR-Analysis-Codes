#include "Run12pp200Ana.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include <vector>
#include <iterator>
using namespace std;
ofstream Output;

void Run12pp200Ana::Loop()
{

	gROOT->Reset();

	Double_t p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1;
	Double_t p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
	Double_t cosPhiRB, sinPhiRB, cosPhiRY, sinPhiRY;
	Double_t PhiR_cosB, PhiR_sinB, PhiRB, PhiR_cosY, PhiR_sinY, PhiRY, phiRB, phiDiffB, phi_cos, phi_sin, phi_pair;
	Double_t msqr, msqr1, msqr2; //
	Double_t pi = 3.14159265359;
	Double_t Cone_Min = 0.02;
	Double_t Cone_Max7 = 0.7;
	Double_t Cone_Max6 = 0.6;
	Double_t Cone_Max5 = 0.5;
	Double_t Cone_Max4 = 0.4;
	Double_t Cone_Max3 = 0.3;
	Double_t Minv_Max = 4.0;	// GeV/c^2
	Double_t Minv_Min = 0.27;	// GeV/c^2
	Double_t ptPair_Max = 15.0; // GeV/c
	Double_t ptPair_Min = 1.0;	// GeV/c
	Double_t Eta_Cut = 1.0;		// GeV/c

	Double_t m_pion = 0.1396; // GeV

	const int nBins = 13;
	Double_t nBinsEdges[nBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

	// control plots
	hxsecMJP0Ctr = new TH1D("hxsecMJP0Ctr", "", nBins, nBinsEdges);
	hxsecMJP1Ctr = new TH1D("hxsecMJP1Ctr", "", nBins, nBinsEdges);
	hxsecMJP2Ctr = new TH1D("hxsecMJP2Ctr", "", nBins, nBinsEdges);

	// histograms for unfolding
	// Cone < 0.7, without Cone_CUT_MIN
	hMDatJP0 = new TH1D("hMDatJP0", "", 160, 0.27, 4.0);
	hMDatJP1 = new TH1D("hMDatJP1", "", 160, 0.27, 4.0);
	hMDatJP2 = new TH1D("hMDatJP2", "", 160, 0.27, 4.0);

	// Cone < 0.7, with Cone_CUT_MIN
	hMDatJP0_C7 = new TH1D("hMDatJP0_C7", "", 160, 0.27, 4.0);
	hMDatJP1_C7 = new TH1D("hMDatJP1_C7", "", 160, 0.27, 4.0);
	hMDatJP2_C7 = new TH1D("hMDatJP2_C7", "", 160, 0.27, 4.0);

	// Cone < 0.3, with Cone_CUT_MIN
	hMDatJP0_C3 = new TH1D("hMDatJP0_C3", "", 160, 0.27, 1.6);
	hMDatJP1_C3 = new TH1D("hMDatJP1_C3", "", 160, 0.27, 1.6);
	hMDatJP2_C3 = new TH1D("hMDatJP2_C3", "", 160, 0.27, 1.6);

	hMDatJP0Gt = new TH1D("hMDatJP0Gt", "", 160, 0.27, 4.0);
	hMDatJP1Gt = new TH1D("hMDatJP1Gt", "", 160, 0.27, 4.0);
	hMDatJP2Gt = new TH1D("hMDatJP2Gt", "", 160, 0.27, 4.0);

	hMDatJP0Lt = new TH1D("hMDatJP0Lt", "", 160, 0.27, 4.0);
	hMDatJP1Lt = new TH1D("hMDatJP1Lt", "", 160, 0.27, 4.0);
	hMDatJP2Lt = new TH1D("hMDatJP2Lt", "", 160, 0.27, 4.0);

	for (int i = 0; i < nBins; i++)
	{
		hMinv[i] = new TH1D(Form("hMinv_%i", i), "", 100, nBinsEdges[i], nBinsEdges[i + 1]);
		hPt[i] = new TH1D(Form("hPt_%i", i), "", 100, 0.2, 15.5);
		hEta[i] = new TH1D(Form("hEta_%i", i), "", 100, -1.2, 1.2);
		if (i < 8)
		{
			hMinv_C3[i] = new TH1D(Form("hMinv_C3%i", i), "", 100, nBinsEdges[i], nBinsEdges[i + 1]);
			hPt_C3[i] = new TH1D(Form("hPt_C3%i", i), "", 100, 0.2, 15.5);
			hEta_C3[i] = new TH1D(Form("hEta_C3%i", i), "", 100, -1.2, 1.2);
		}
	}

	if (fChain == 0)
		return;

	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	// Event Loop
	Double_t pairCount;

	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	// for (Long64_t jentry = 0; jentry < 10000; jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		if (fverRank < 1e6)
			continue; // Already applied at tree level
		if (fabs(fVZ) > 60.)
			continue;

		Bool_t isJP = kFALSE;
		Bool_t isJP0 = kFALSE;
		Bool_t isJP1 = kFALSE;
		Bool_t isJP2 = kFALSE;
		unsigned int triggerId = 0;
		for (int trig = 0; trig < ftrigger->size(); trig++)
		{
			triggerId = ftrigger->at(trig);
			if (triggerId == 370601 || triggerId == 370611 || triggerId == 370621)
				isJP = kTRUE;
			if (triggerId == 370601)
				isJP0 = kTRUE;
			if (triggerId == 370611)
				isJP1 = kTRUE;
			if (triggerId == 370621)
				isJP2 = kTRUE;
		}
		if (!isJP)
			continue; // only JP0, JP1, and JP2 triggers
		pairCount = 0;
		// Track Loop
		for (int tr = 0; tr < fmaxpar; tr++)
		{
			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];
			if (fpT[tr] > 0.5 && fpT[tr] < 15.0 && ffitPts[tr] > 15 && abs(feta[tr]) < 1. && fitPtsRatio_tr > .51 && fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2.0) // && fhitsdedx[tr] > 20)
			{

				// Track 2 loop
				for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
				{
					Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];
					if (fpT[tr2] > 0.5 && fpT[tr2] < 15. && ffitPts[tr2] > 15 && fabs(feta[tr2]) < 1. && fnSigmaPion[tr2] > -1.0 && fnSigmaPion[tr2] < 2.0 && fitPtsRatio_tr2 > .51) //&& fhitsdedx[tr2] > 20)
					{

						Double_t phiDiff = fphi[tr] - fphi[tr2];
						if (phiDiff > pi)
							phiDiff -= (2 * pi);
						if (phiDiff < ((-1) * pi))
							phiDiff += (2 * pi);
						// if (phiDiff > pi || phiDiff < -pi)
						//	continue;

						if (fcharge[tr] == fcharge[tr2])
							continue;

						cone = sqrt(pow(feta[tr] - feta[tr2], 2) + pow(phiDiff, 2));

						if (cone < Cone_Min || cone >= Cone_Max7)
							continue;
						// if (cone >= Cone_Max7)
						//	continue;

						if (fcharge[tr] > 0)
						{
							p1x = fpT[tr] * cos(fphi[tr]);
							p2x = fpT[tr2] * cos(fphi[tr2]);
							p1y = fpT[tr] * sin(fphi[tr]);
							p2y = fpT[tr2] * sin(fphi[tr2]);
							p1z = fpT[tr] * sinh(feta[tr]);
							p2z = fpT[tr2] * sinh(feta[tr2]);
							/*
							p1x = fpX[tr];
							p2x = fpX[tr2];
							p1y = fpY[tr];
							p2y = fpY[tr2];
							p1z = fpZ[tr];
							p2z = fpZ[tr2];
							*/
							p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
							p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
						}
						if (fcharge[tr] < 0)
						{
							p1x = fpT[tr2] * cos(fphi[tr2]);
							p2x = fpT[tr] * cos(fphi[tr]);
							p1y = fpT[tr2] * sin(fphi[tr2]);
							p2y = fpT[tr] * sin(fphi[tr]);
							p1z = fpT[tr2] * sinh(feta[tr2]);
							p2z = fpT[tr] * sinh(feta[tr]);
							/*
							p1x = fpX[tr2];
							p2x = fpX[tr];
							p1y = fpY[tr2];
							p2y = fpY[tr];
							p1z = fpZ[tr2];
							p2z = fpZ[tr];
							*/
							p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
							p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
						}

						// Components of sum vector
						psx = p1x + p2x;
						psy = p1y + p2y;
						psz = p1z + p2z;
						// Sum vector
						ps = sqrt(psx * psx + psy * psy + psz * psz);
						// difference of momentum
						Rx = (p1x - p2x);
						Ry = (p1y - p2y);
						Rz = (p1z - p2z);
						// relative momentum vector
						R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

						// calculate M,pt,eta
						E1 = sqrt(m_pion * m_pion + p1 * p1);
						E2 = sqrt(m_pion * m_pion + p2 * p2);
						Minv = sqrt(2 * m_pion * m_pion + 2 * (E1 * E2 - p1x * p2x - p1y * p2y - p1z * p2z));

						pT_pair = sqrt((p1x + p2x) * (p1x + p2x) + (p1y + p2y) * (p1y + p2y));

						eta_pair = TMath::ASinH((p1z + p2z) / pT_pair);

						if (Minv < Minv_Min || Minv > Minv_Max || pT_pair < ptPair_Min || pT_pair > ptPair_Max || fabs(eta_pair) > Eta_Cut)
							continue;

						phi_cos = acos((p1x + p2x) / pT_pair);
						phi_sin = asin((p1y + p2y) / pT_pair);
						if (phi_sin > 0)
							phi_pair = phi_cos;
						if (phi_sin < 0)
							phi_pair = -1. * phi_cos;

						// Vector definations for phiR calculation
						TVector3 v_blue(0, 0, 1);						// Blue beam
						TVector3 v_yell(0, 0, -1);						// yellow beam
						TVector3 pSum(psx, psy, psz);					// pion pair sum vector
						TVector3 pSumHat(psx / ps, psy / ps, psz / ps); // pion pair's momentum unit vector
						TVector3 Rc(Rx, Ry, Rz);						// pion pair's momemtum difference vector

						// blue beam
						cosPhiRB = (pSumHat.Cross(v_blue) * ((Double_t)1 / (pSumHat.Cross(v_blue).Mag()))).Dot((pSumHat.Cross(Rc)) * ((Double_t)1 / (pSumHat.Cross(Rc).Mag())));
						sinPhiRB = (v_blue.Cross(Rc)).Dot(pSumHat) * ((Double_t)1 / (pSumHat.Cross(v_blue).Mag())) * ((Double_t)1 / (pSumHat.Cross(Rc).Mag()));
						PhiR_cosB = acos(cosPhiRB);
						PhiR_sinB = asin(sinPhiRB);

						if (PhiR_sinB > 0)
							PhiRB = PhiR_cosB;
						if (PhiR_sinB < 0)
							PhiRB = (-1) * PhiR_cosB;

						// yellow beam
						cosPhiRY = (pSumHat.Cross(v_yell) * ((Double_t)1 / (pSumHat.Cross(v_yell).Mag()))).Dot((pSumHat.Cross(Rc)) * ((Double_t)1 / (pSumHat.Cross(Rc).Mag())));
						sinPhiRY = (v_yell.Cross(Rc)).Dot(pSumHat) * ((Double_t)1 / (pSumHat.Cross(v_yell).Mag())) * ((Double_t)1 / (pSumHat.Cross(Rc).Mag()));
						PhiR_cosY = acos(cosPhiRY);
						PhiR_sinY = asin(sinPhiRY);
						if (PhiR_sinY > 0)
							PhiRY = PhiR_cosY;
						if (PhiR_sinY < 0)
							PhiRY = -1. * PhiR_cosY;
						// phiR done-------

						if (isJP0)
						{
							hxsecMJP0Ctr->Fill(Minv);
							hMDatJP0->Fill(Minv);
							if (cone > Cone_Min && cone < Cone_Max7)
								hMDatJP0_C7->Fill(Minv);
							if (cone > Cone_Min && cone < Cone_Max3 && Minv < 1.6)
								hMDatJP0_C3->Fill(Minv);
						}
						if (isJP1)
						{
							hxsecMJP1Ctr->Fill(Minv);
							hMDatJP1->Fill(Minv);
							if (cone > Cone_Min && cone < Cone_Max7)
								hMDatJP1_C7->Fill(Minv);
							if (cone > Cone_Min && cone < Cone_Max3 && Minv < 1.6)
								hMDatJP1_C3->Fill(Minv);
						}
						if (isJP2)
						{
							hxsecMJP2Ctr->Fill(Minv);
							hMDatJP2->Fill(Minv);
							if (cone > Cone_Min && cone < Cone_Max7)
								hMDatJP2_C7->Fill(Minv);
							if (cone > Cone_Min && cone < Cone_Max3 && Minv < 1.6)
								hMDatJP2_C3->Fill(Minv);
						}

						for (int i = 0; i < 13; i++)
						{
							if (Minv >= nBinsEdges[i] && Minv < nBinsEdges[i + 1])
							{
								if (cone > Cone_Min)
								{
									hMinv[i]->Fill(Minv);
									hPt[i]->Fill(pT_pair);
									hEta[i]->Fill(eta_pair);
								}
								if (i < 8)
								{
									if (cone > Cone_Min && cone < Cone_Max3 && Minv < 1.6)
									{
										hMinv_C3[i]->Fill(Minv);
										hPt_C3[i]->Fill(pT_pair);
										hEta_C3[i]->Fill(eta_pair);
									}
								}
							}
						}

						pairCount++;
					} // track 2 cuts
				}	  // track 2 loop
			}		  // track 1 PID and selection cuts
		}			  // track loop

	} // event loop
} // Iff::Loop()

Run12pp200Ana::Run12pp200Ana(char *ifile)
{
	fChain = new TChain("ftree");
	fChain->Add(ifile);

	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	Init();
}

Run12pp200Ana::~Run12pp200Ana()
{
	if (!fChain)
		return;
	delete fChain->GetCurrentFile();
}

Int_t Run12pp200Ana::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain)
		return 0;
	return fChain->GetEntry(entry);
}
Long64_t Run12pp200Ana::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain)
		return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0)
		return centry;
	if (!fChain->InheritsFrom(TChain::Class()))
		return centry;
	TChain *chain = (TChain *)fChain;
	if (chain->GetTreeNumber() != fCurrent)
	{
		fCurrent = chain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void Run12pp200Ana::Init()
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers

	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("frefmult", &frefmult, &b_frefmult);
	fChain->SetBranchAddress("fmaxpar", &fmaxpar, &b_fmaxpar);
	fChain->SetBranchAddress("ffillNum", &ffillNum, &b_ffillNum);
	fChain->SetBranchAddress("frunNum", &frunNum, &b_frunNum);
	fChain->SetBranchAddress("ftrigger", &ftrigger);
	fChain->SetBranchAddress("fVZ", &fVZ, &b_fVZ);
	fChain->SetBranchAddress("fevTime", &fevTime, &b_fevTime);
	fChain->SetBranchAddress("fverRank", &fverRank, &b_fverRank);
	fChain->SetBranchAddress("fpT", fpT, &b_fpT);
	fChain->SetBranchAddress("fp", fp, &b_fp);
	fChain->SetBranchAddress("fpX", fpX, &b_fpX);
	fChain->SetBranchAddress("fpY", fpY, &b_fpY);
	fChain->SetBranchAddress("fpZ", fpZ, &b_fpZ);
	fChain->SetBranchAddress("feta", feta, &b_feta);
	fChain->SetBranchAddress("fphi", fphi, &b_fphi);
	fChain->SetBranchAddress("fcharge", fcharge, &b_fcharge);
	fChain->SetBranchAddress("fnSigmaPion", fnSigmaPion, &b_fnSigmaPion);
	fChain->SetBranchAddress("fnSigmaKaon", fnSigmaKaon, &b_fnSigmaKaon);
	fChain->SetBranchAddress("fnSigmaProton", fnSigmaProton, &b_fnSigmaProton);
	fChain->SetBranchAddress("fnSigmaElectron", fnSigmaElectron, &b_fnSigmaElectron);
	fChain->SetBranchAddress("fdEdx", fdEdx, &b_fdEdx);
	fChain->SetBranchAddress("fdca", fdca, &b_fdca);
	fChain->SetBranchAddress("ffitPts", ffitPts, &b_ffitPts);
	fChain->SetBranchAddress("ffitPtsPoss", ffitPtsPoss, &b_ffitPtsPoss);
	fChain->SetBranchAddress("fhitsdedx", fhitsdedx, &b_fhitsdedx);
	fChain->SetBranchAddress("fBetaToF", fBetaToF, &b_fBetaToF);
	fChain->SetBranchAddress("fvpdVz", &fvpdVz, &b_fvpdVz);
	Notify();
}

Bool_t Run12pp200Ana::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}
void Run12pp200Ana::Finish(char *ofile)
{
	TFile *fout = new TFile(ofile, "recreate");
	fout->cd();
	hMDatJP0->Write();
	hMDatJP1->Write();
	hMDatJP2->Write();

	hMDatJP0_C7->Write();
	hMDatJP1_C7->Write();
	hMDatJP2_C7->Write();

	hMDatJP0_C3->Write();
	hMDatJP1_C3->Write();
	hMDatJP2_C3->Write();

	hxsecMJP0Ctr->Write();
	hxsecMJP1Ctr->Write();
	hxsecMJP2Ctr->Write();
	for (int i = 0; i < 13; i++)
	{
		hMinv[i]->Write();
		hPt[i]->Write();
		hEta[i]->Write();
		if (i < 8)
		{
			hMinv_C3[i]->Write();
			hPt_C3[i]->Write();
			hEta_C3[i]->Write();
		}
	}

	fout->Close();
}

void Run12pp200Ana::Show(Long64_t entry)
{
	// eventTime-> Write();
	// pions->Write();
	//  Print contents of entry.
	//  If entry is not specified, print current entry
	if (!fChain)
		return;
	fChain->Show(entry);
}
Int_t Run12pp200Ana::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

ClassImp(Run12pp200Ana);
