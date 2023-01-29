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
using namespace std;
ofstream Output;

void Run12pp200Ana::Loop()
{
	gROOT->Reset();
	TH1::SetDefaultSumw2();
	Int_t nTrg = 3;
	Int_t nCh = 2;
	const char *tTrig[3] = {"JP0", "JP1", "JP2"};
	const char *tCh[2] = {"Pos", "Neg"};

	const int nBins = 16;
	Double_t nBinsEdges[nBins + 1] = {0.28, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};
	// hist in xsec-bins
	hxsecMJP0 = new TH1D("hxsecMJP0", "", nBins, nBinsEdges);
	hxsecMJP1 = new TH1D("hxsecMJP1", "", nBins, nBinsEdges);
	hxsecMJP2 = new TH1D("hxsecMJP2", "", nBins, nBinsEdges);

	// hist for combined trigger
	hVZJP = new TH1D("hVZJP", "", 100, -65.0, 65.0);
	hConeJP = new TH1D("hConeJP", "", 100, 0.0, 0.8);
	hMinvJP = new TH1D("hMinvJP", "", 100, 0.0, 4.0);
	hpTPairJP = new TH1D("hpTPairJP", "", 100, 0.0, 16.0);
	hEtaPairJP = new TH1D("hEtaPairJP", "", 100, -1.2, 1.2);

	htrackEtaJP = new TH1D("htrackEtaJP", "", 100, -1.2, 1.2);
	htrackpTJP = new TH1D("htrackpTJP", "", 100, 0.0, 16);
	htrackPhiJP = new TH1D("htrackPhiJP", "", 100, -3.3, 3.3);
	htracknSigmaPionJP = new TH1D("htracknSigmaPionJP", "", 100, -10, 10);

	hpionEtaJP = new TH1D("hpionEtaJP", "", 100, -1.2, 1.2);
	hpionpTJP = new TH1D("hpionpTJP", "", 100, 0.0, 16);
	hpionPhiJP = new TH1D("hpionPhiJP", "", 100, -3.3, 3.3);

	// Hists for JP0, JP1 and JP2
	for (int ntrg = 0; ntrg < nTrg; ntrg++)
	{
		hVZ[ntrg] = new TH1D(Form("hVZ_%s", tTrig[ntrg]), "", 100, -65.0, 65.0);
		hCone[ntrg] = new TH1D(Form("hCone_%s", tTrig[ntrg]), "", 100, 0.0, 0.8);
		hMinv[ntrg] = new TH1D(Form("hMinv_%s", tTrig[ntrg]), "", 100, 0.0, 4.0);
		hpTPair[ntrg] = new TH1D(Form("hpTPair_%s", tTrig[ntrg]), "", 100, 0.0, 16.0);
		hEtaPair[ntrg] = new TH1D(Form("hEtaPair_%s", tTrig[ntrg]), "", 100, -1.2, 1.2);

		htrackEta[ntrg] = new TH1D(Form("htrackEta_%s", tTrig[ntrg]), "", 100, -1.2, 1.2);
		htrackpT[ntrg] = new TH1D(Form("htrackpT_%s", tTrig[ntrg]), "", 100, 0.0, 16);
		htrackPhi[ntrg] = new TH1D(Form("htrackPhi_%s", tTrig[ntrg]), "", 100, -3.3, 3.3);
		htracknSigmaPion[ntrg] = new TH1D(Form("htracknSigmaPion_%s", tTrig[ntrg]), "", 100, -10, 10);
		hpionEta[ntrg] = new TH1D(Form("hpionEta_%s", tTrig[ntrg]), "", 100, -1.2, 1.2);
		hpionpT[ntrg] = new TH1D(Form("hpionpT_%s", tTrig[ntrg]), "", 100, 0.0, 16);
		hpionPhi[ntrg] = new TH1D(Form("hpionPhi_%s", tTrig[ntrg]), "", 100, -3.3, 3.3);
	}

	Double_t p1x = -999.0, p2x = -999.0, p1y = -999.0, p2y = -999.0, p1z = -999., p2z = -999., psx = -999., psy = -999., psz = -999., ps = -999., cone = -999., R = -999., Rx = -999., Ry = -999., Rz = -999., R1 = -999.;
	Double_t p1 = -999., p2 = -999., E1 = -999., E2 = -999., Minv = -999., pT_pair = -999., pT_min_pair = -999., eta_pair = -999., fitPts_min_pair = -999.;
	Double_t cosPhiRB = -999., sinPhiRB = -999., cosPhiRY = -999., sinPhiRY = -999.;
	Double_t PhiR_cosB = -999., PhiR_sinB = -999., PhiRB = -999., PhiR_cosY = -999., PhiR_sinY = -999., PhiRY = -999., phiRB = -999., phiDiffB = -999.;
	Double_t pi = 3.14159265359;
	Double_t CONE_CUT_MIN = 0.0;
	Double_t m_pion = 0.1396; // GeV

	if (fChain == 0)
		return;

	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	// Event Loop
	cout << "Event loop started with the number of events  = " << nentries << endl;
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		if (fverRank < 1e6 || fabs(fVZ) > 60)
			continue;

		Bool_t isJP = kFALSE;
		Bool_t isJP0 = kFALSE;
		Bool_t isJP1 = kFALSE;
		Bool_t isJP2 = kFALSE;
		unsigned int triggerId = 0;
		for (unsigned int trig = 0; trig < ftrigger->size(); trig++)
		{
			triggerId = ftrigger->at(trig);
			// cout << "Event No: " << jentry << "Trigger ID: " << triggerId << endl;
			if (triggerId == 370601 || triggerId == 370611 || triggerId == 370621)
				isJP = kTRUE;
			if (triggerId == 370601)
				isJP0 = kTRUE;
			if (triggerId == 370611)
				isJP1 = kTRUE;
			if (triggerId == 370621)
				isJP2 = kTRUE;
		}
		//if (!isJP)
		//	continue;

		if (isJP0)
			hVZ[0]->Fill(fVZ);
		if (isJP1)
			hVZ[1]->Fill(fVZ);
		if (isJP2)
			hVZ[2]->Fill(fVZ);

		if (isJP0 || isJP1 || isJP2)
			hVZJP->Fill(fVZ);

		// Track Loop
		for (int tr = 0; tr < fmaxpar; tr++)
		{
			Double_t fitPtsRatio_tr = 0.0;
			fitPtsRatio_tr = (Double_t)ffitPts[tr] / ffitPtsPoss[tr];

			// if (fpT[tr] > 0.5 && fpT[tr] < 15.0 && ffitPts[tr] > 15 && fabs(feta[tr]) < 1. && fitPtsRatio_tr > .51){
			if (fpT[tr] < 0.5 || fpT[tr] > 15.0 || ffitPts[tr] < 15 || fabs(feta[tr]) > 1. || fitPtsRatio_tr < .51 || !dcaCut(fpT[tr], fdca[tr]))
				continue;

			if (isJP2)
			{
				htrackEta[2]->Fill(feta[tr]);
				htrackPhi[2]->Fill(fphi[tr]);
				htrackpT[2]->Fill(fpT[tr]);
				htracknSigmaPion[2]->Fill(fnSigmaPion[tr]);
			}
			if (isJP1)
			{
				htrackEta[1]->Fill(feta[tr]);
				htrackPhi[1]->Fill(fphi[tr]);
				htrackpT[1]->Fill(fpT[tr]);
				htracknSigmaPion[1]->Fill(fnSigmaPion[tr]);
			}
			if (isJP0)
			{
				htrackEta[0]->Fill(feta[tr]);
				htrackPhi[0]->Fill(fphi[tr]);
				htrackpT[0]->Fill(fpT[tr]);
				htracknSigmaPion[0]->Fill(fnSigmaPion[tr]);
			}

			if (isJP0 || isJP1 || isJP2)
			{
				htrackEtaJP->Fill(feta[tr]);
				htrackPhiJP->Fill(fphi[tr]);
				htrackpTJP->Fill(fpT[tr]);
				htracknSigmaPionJP->Fill(fnSigmaPion[tr]);
			}

			if (!(fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2.0))
				continue;

			if (isJP2)
			{
				hpionEta[2]->Fill(feta[tr]);
				hpionPhi[2]->Fill(fphi[tr]);
				hpionpT[2]->Fill(fpT[tr]);
			}
			if (isJP1)
			{
				hpionEta[1]->Fill(feta[tr]);
				hpionPhi[1]->Fill(fphi[tr]);
				hpionpT[1]->Fill(fpT[tr]);
			}
			if (isJP0)
			{
				hpionEta[0]->Fill(feta[tr]);
				hpionPhi[0]->Fill(fphi[tr]);
				hpionpT[0]->Fill(fpT[tr]);
			}
			if (isJP0 || isJP1 || isJP2)
			{
				hpionEtaJP->Fill(feta[tr]);
				hpionPhiJP->Fill(fphi[tr]);
				hpionpTJP->Fill(fpT[tr]);
			}
			// Track 2 loop
			for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
			{
				Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];
				// if (fpT[tr2] > 0.5 && fpT[tr2] < 15.0 && ffitPts[tr2] > 15 && fabs(feta[tr2]) < 1. && fnSigmaPion[tr2] > -1.0 && fnSigmaPion[tr2] < 2.0 && fitPtsRatio_tr2 > .51){
				if (fpT[tr2] < 0.5 || fpT[tr2] > 15.0 || ffitPts[tr2] < 15 || fabs(feta[tr2]) > 1. || fnSigmaPion[tr2] < -1.0 || fnSigmaPion[tr2] > 2.0 || fitPtsRatio_tr2 < .51 || !dcaCut(fpT[tr2], fdca[tr2]))
					continue;

				Double_t phiDiff = fphi[tr] - fphi[tr2];
				if (phiDiff > pi)
					phiDiff -= (2 * pi);
				if (phiDiff < ((-1) * pi))
					phiDiff += (2 * pi);

				if (fcharge[tr] == fcharge[tr2])
					continue;

				cone = sqrt(pow(feta[tr] - feta[tr2], 2) + pow(phiDiff, 2)); // cone cut < 0.3
				if (cone >= 0.7)
					continue;

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

				psx = p1x + p2x;
				psy = p1y + p2y;
				psz = p1z + p2z;
				// Sum vector
				ps = sqrt(psx * psx + psy * psy + psz * psz);
				// Relative momentum
				Rx = (p1x - p2x);
				Ry = (p1y - p2y);
				Rz = (p1z - p2z);

				R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

				// calculate M,pt,eta
				E1 = sqrt(m_pion * m_pion + p1 * p1);
				E2 = sqrt(m_pion * m_pion + p2 * p2);
				Minv = sqrt(2 * m_pion * m_pion + 2 * (E1 * E2 - p1x * p2x - p1y * p2y - p1z * p2z));
				pT_pair = sqrt((p1x + p2x) * (p1x + p2x) + (p1y + p2y) * (p1y + p2y));
				eta_pair = TMath::ASinH((p1z + p2z) / pT_pair);

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

				if (Minv > 4.0 || fabs(eta_pair) > 1.0 || pT_pair > 15.0)
					continue;

				if (isJP0)
				{
					hxsecMJP0->Fill(Minv);
					hCone[0]->Fill(cone);
					hMinv[0]->Fill(Minv);
					hEtaPair[0]->Fill(eta_pair);
					hpTPair[0]->Fill(pT_pair);
				}
				if (isJP1)
				{
					hxsecMJP1->Fill(Minv);
					hCone[1]->Fill(cone);
					hMinv[1]->Fill(Minv);
					hEtaPair[1]->Fill(eta_pair);
					hpTPair[1]->Fill(pT_pair);
				}
				if (isJP2)
				{
					hxsecMJP2->Fill(Minv);
					hCone[2]->Fill(cone);
					hMinv[2]->Fill(Minv);
					hEtaPair[2]->Fill(eta_pair);
					hpTPair[2]->Fill(pT_pair);
				}

				if (isJP0 || isJP1 || isJP2)
				{
					hConeJP->Fill(cone);
					hMinvJP->Fill(Minv);
					hEtaPairJP->Fill(eta_pair);
					hpTPairJP->Fill(pT_pair);
				}
			} // track 2 loop
		}	  // track loop

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
	hVZJP->SetDirectory(fout);
	hVZJP->Write();

	hxsecMJP0->SetDirectory(fout);
	hxsecMJP0->Write();
	hxsecMJP1->SetDirectory(fout);
	hxsecMJP1->Write();
	hxsecMJP2->SetDirectory(fout);
	hxsecMJP2->Write();

	hConeJP->SetDirectory(fout);
	hConeJP->Write();
	hMinvJP->SetDirectory(fout);
	hMinvJP->Write();
	hpTPairJP->SetDirectory(fout);
	hpTPairJP->Write();
	hEtaPairJP->SetDirectory(fout);
	hEtaPairJP->Write();

	htrackEtaJP->SetDirectory(fout);
	htrackEtaJP->Write();
	htrackpTJP->SetDirectory(fout);
	htrackpTJP->Write();
	htrackPhiJP->SetDirectory(fout);
	htrackPhiJP->Write();
	htracknSigmaPionJP->SetDirectory(fout);
	htracknSigmaPionJP->Write();

	hpionEtaJP->SetDirectory(fout);
	hpionEtaJP->Write();
	hpionpTJP->SetDirectory(fout);
	hpionpTJP->Write();
	hpionPhiJP->SetDirectory(fout);
	hpionPhiJP->Write();

	for (int i = 0; i < 3; i++)
	{
		hVZ[i]->SetDirectory(fout);
		hVZ[i]->Write();
		hCone[i]->SetDirectory(fout);
		hCone[i]->Write();
		hMinv[i]->SetDirectory(fout);
		hMinv[i]->Write();
		hpTPair[i]->SetDirectory(fout);
		hpTPair[i]->Write();
		hEtaPair[i]->SetDirectory(fout);
		hEtaPair[i]->Write();

		htrackEta[i]->SetDirectory(fout);
		htrackEta[i]->Write();
		htrackpT[i]->SetDirectory(fout);
		htrackpT[i]->Write();
		htrackPhi[i]->SetDirectory(fout);
		htrackPhi[i]->Write();
		htracknSigmaPion[i]->SetDirectory(fout);
		htracknSigmaPion[i]->Write();

		hpionEta[i]->SetDirectory(fout);
		hpionEta[i]->Write();
		hpionpT[i]->SetDirectory(fout);
		hpionpT[i]->Write();
		hpionPhi[i]->SetDirectory(fout);
		hpionPhi[i]->Write();
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
Bool_t Run12pp200Ana::dcaCut(Double_t pT, Double_t dca)
{
	if ((pT < 0.5 && dca < 2.0) || (pT > 0.5 && pT < 1.5 && dca < ((-1) * pT + 2.5)) || (pT > 1.5 && dca < 1.0))
	{
		return true;
	}
	return false;
}

Int_t Run12pp200Ana::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

ClassImp(Run12pp200Ana);
