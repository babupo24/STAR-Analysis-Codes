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
	// initialize histograms
	Int_t nTrg = 3;
	Int_t nCh = 2;
	const char *tTrig[3] = {"JP0", "JP1", "JP2"};
	const char *tCh[2] = {"Pos", "Neg"};
	// final cross-section mass bins after resolution studies
	const int nBins = 16;
	Double_t nBinsEdges[nBins + 1] = {0.28, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};
	hpartonicPt = new TH1D("hpartonicPt", "", 100, 0.0, 100.0);
	// hist in xsec-bins
	hxsecMJP0 = new TH1D("hxsecMJP0", "", nBins, nBinsEdges);
	hxsecMJP1 = new TH1D("hxsecMJP1", "", nBins, nBinsEdges);
	hxsecMJP2 = new TH1D("hxsecMJP2", "", nBins, nBinsEdges);

	for (int ntrg = 0; ntrg < nTrg; ntrg++)
	{
		hVZ[ntrg] = new TH1D(Form("hVZ_%s", tTrig[ntrg]), "", 100, -65.0, 65.0);
		hCone[ntrg] = new TH1D(Form("hCone_%s", tTrig[ntrg]), "", 100, 0.0, 0.8);
		hMinv[ntrg] = new TH1D(Form("hMinv_%s", tTrig[ntrg]), "", 100, 0.0, 4.0);
		hpTPair[ntrg] = new TH1D(Form("hpTPair_%s", tTrig[ntrg]), "", 100, 0.0, 16.0);
		hEtaPair[ntrg] = new TH1D(Form("hEtaPair_%s", tTrig[ntrg]), "", 100, -1.2, 1.2);
		hPhiRY[ntrg] = new TH1D(Form("hPhiRY_%s", tTrig[ntrg]), "", 100, -3.2, 3.2);
		hPhiRB[ntrg] = new TH1D(Form("hPhiRB_%s", tTrig[ntrg]), "", 100, -3.2, 3.2);

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
	Double_t PhiR_cosB = -999., PhiR_sinB = -999., PhiRB = -999., PhiR_cosY = -999., PhiR_sinY = -999., PhiRY = -999., phiRB = -999., phiDiffB = -999., phi_pair = -999.0, phi_sin = -999.0, phi_cos = -999.0;
	Double_t pi = 3.14159265359;
	Double_t CONE_CUT_MIN = 0.0;
	Double_t m_pion = 0.1396; // GeV

	Double_t R1_p1 = -999, R2_p1 = -999, R1_p2 = -999, R2_p2 = -999, Rs_tr1 = -999, Rs_tr2 = -999, sourcePartEta = -999, sourcePartPhi = -999, Rs_pair = -999; // opening angles between parton and tracks and associated parton eta, phi

	if (fChain == 0)
		return;

	// calculate weight for each partonic pT bins
	const Int_t nptBin = 13;
	Double_t partPtRange[nptBin + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., -1};
	Double_t fudge[nptBin] = {1.228, 1.051, 1.014, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	Double_t crossSection[nptBin] = {9.0012e+00, 1.46253e+00, 3.54566e-01, 1.51622e-01, 2.49062e-02, 5.84527e-03, 2.30158e-03, 3.42755e-04, 4.57002e-05, 9.72535e-06, 4.69889e-07, 2.69202e-08, 1.43453e-09};
	// Double_t numEvents[nptBin] = {3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 2459112, 2459112, 1229556, 1229556}; // number of events produced
	Double_t numEvents[nptBin] = {3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 3688668, 2458337, 2459112, 1229556, 1229556}; // number of events produced

	// Double_t numEvents[nptBin] = {3318626, 3301413, 3291662, 3280010, 3282543, 3275693, 3276437, 3276795, 3272804, 2179660, 2183230, 1091927, 1090857}; // number of events in file after processing MuDsts
	//  Double_t numEvents[nptBin] = {1.84159e+06, 2.1678e+06, 2.31769e+06, 2.41826e+06, 2.5085e+06, 2.55509e+06, 2.5922e+06, 2.62385e+06, 2.63813e+06, 1.76531e+06, 1.7739e+06, 887798, 888318}; // number of good events processed for the analysis
	Double_t binLumi[nptBin] = {0}, binWt[nptBin] = {0};
	for (int mBin = 0; mBin < nptBin; mBin++)
	{
		binLumi[mBin] = (Double_t)(numEvents[mBin] * crossSection[12]) / (crossSection[mBin] * numEvents[12]); // bin luminosity
		binWt[mBin] = (Double_t)1. / (binLumi[mBin] * fudge[mBin]);											   // Bin weight for partonic pT weighting.
		cout << binWt[mBin] << endl;
	}

	Long64_t nentries = fChain->GetEntries();
	cout << "Event loop started with no of entries: " << nentries << endl;
	Long64_t nbytes = 0, nb = 0;
	// Event Loop
	Double_t nEventsProcessed = 0;
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	// for (Long64_t jentry = 0; jentry < 10000; jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		Double_t ptw = 0;
		Double_t extra_weight = 0;
		// weight from vertex shape correction
		// Double_t vz_weight = 0.9890 + (fVZ * 0.00110841) + (fVZ * fVZ * 1.79596e-05) + (fVZ * fVZ * fVZ * 5.71566e-07);
		Double_t vz_weight = 1.0;

		// soft reweighting from Dmitry
		double p[4] = {1.21733561, -0.32633577, 0.16913723, 0.82134143};
		extra_weight = (Double_t)1. / (1. + (p[0] + p[1] * (partonicPtBin - 2.) + p[2] * (partonicPtBin - 2.) * (partonicPtBin - 2.)) * exp(-p[3] * (partonicPtBin - 2.)));

		for (int j = 0; j < 13; j++)
		{
			if (j < 12 && partonicPtBin >= partPtRange[j] && partonicPtBin < partPtRange[j + 1])
			{
				ptw = binWt[j] * extra_weight * vz_weight;
				// ptw = binWt[j];
			}
			else if (j == 12 && partonicPtBin >= partPtRange[j])
			{
				ptw = binWt[j] * extra_weight * vz_weight;
				// ptw = binWt[j];
			}
		}

		if (fverRank < 1e6 || fabs(fVZ) > 60.)
			continue;

		hpartonicPt->Fill(partonicPtBin, ptw);

		Bool_t isJP = kFALSE;
		Bool_t isJP0 = kFALSE;
		Bool_t isJP1 = kFALSE;
		Bool_t isJP2 = kFALSE;
		int triggerId = 0;
		for (unsigned int trig = 0; trig < ftrigger->size(); trig++)
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
		// if (!isJP)
		//	continue;

		if (isJP2)
			hVZ[2]->Fill(fVZ, ptw);
		if (isJP1)
			hVZ[1]->Fill(fVZ, ptw);
		if (isJP0)
			hVZ[0]->Fill(fVZ, ptw);

		// Track Loop
		for (int tr = 0; tr < fmaxpar; tr++)
		{

			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];
			if (fpT[tr] < 0.5 || fpT[tr] > 15.0 || ffitPts[tr] < 15 || fabs(feta[tr]) > 1. || fitPtsRatio_tr < .51 || !dcaCut(fpT[tr], fdca[tr]))
				continue;

			if (isJP2)
			{
				htrackEta[2]->Fill(feta[tr], ptw);
				htrackPhi[2]->Fill(fphi[tr], ptw);
				htrackpT[2]->Fill(fpT[tr], ptw);
				htracknSigmaPion[2]->Fill(fnSigmaPion[tr], ptw);
			}
			if (isJP1)
			{
				htrackEta[1]->Fill(feta[tr], ptw);
				htrackPhi[1]->Fill(fphi[tr], ptw);
				htrackpT[1]->Fill(fpT[tr], ptw);
				htracknSigmaPion[1]->Fill(fnSigmaPion[tr], ptw);
			}
			if (isJP0)
			{
				htrackEta[0]->Fill(feta[tr], ptw);
				htrackPhi[0]->Fill(fphi[tr], ptw);
				htrackpT[0]->Fill(fpT[tr], ptw);
				htracknSigmaPion[0]->Fill(fnSigmaPion[tr], ptw);
			}
			if (!(fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2.0))
				continue;

			if (isJP2)
			{
				hpionEta[2]->Fill(feta[tr], ptw);
				hpionPhi[2]->Fill(fphi[tr], ptw);
				hpionpT[2]->Fill(fpT[tr], ptw);
			}
			if (isJP1)
			{
				hpionEta[1]->Fill(feta[tr], ptw);
				hpionPhi[1]->Fill(fphi[tr], ptw);
				hpionpT[1]->Fill(fpT[tr], ptw);
			}
			if (isJP0)
			{
				hpionEta[0]->Fill(feta[tr], ptw);
				hpionPhi[0]->Fill(fphi[tr], ptw);
				hpionpT[0]->Fill(fpT[tr], ptw);
			}

			//  Track 2 loop
			for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
			{
				Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];
				// if (fpT[tr2] < 0.5 || fpT[tr2] > 15.0 || ffitPts[tr2] < 15 || fabs(feta[tr2]) > 1. || fnSigmaPion[tr2] < -1.0 || fnSigmaPion[tr2] > 2.0 || fitPtsRatio_tr2 < .51 || !dcaCut(fpT[tr2], fdca[tr2]))
				if (fpT[tr2] < 0.5 || fpT[tr2] > 15.0 || ffitPts[tr2] < 15 || fabs(feta[tr2]) > 1. || fnSigmaPion[tr2] < -1.0 || fnSigmaPion[tr2] > 2.0 || fitPtsRatio_tr2 < .51 || !dcaCut(fpT[tr], fdca[tr]))
					continue;

				// if (fpT[tr2] > 0.5 && fpT[tr2] < 15. && ffitPts[tr2] > 15 && fabs(feta[tr2]) < 1. && fnSigmaPion[tr2] > -1.0 && fnSigmaPion[tr2] < 2.0 && fitPtsRatio_tr2 > .51 && dcaCut(fpT[tr2], fdca[tr2])){

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
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
					// cout << "Track 1 +ve, " << p1x << ", " << p1y << " , " << p1z << endl;
				}
				if (fcharge[tr] < 0)
				{
					p1x = fpT[tr2] * cos(fphi[tr2]);
					p2x = fpT[tr] * cos(fphi[tr]);
					p1y = fpT[tr2] * sin(fphi[tr2]);
					p2y = fpT[tr] * sin(fphi[tr]);
					p1z = fpT[tr2] * sinh(feta[tr2]);
					p2z = fpT[tr] * sinh(feta[tr]);
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
					// cout << "Track 1 -ve, " << p1x << ", " << p1y << " , " << p1z << endl;
				}
				// Components of sum vector
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
				// phi_cos = acos((fpT[tr] * cos(fphi[tr]) + fpT[tr2] * cos(fphi[tr2])) / pT_pair);
				// phi_sin = asin((fpT[tr] * sin(fphi[tr]) + fpT[tr2] * sin(fphi[tr2])) / pT_pair);
				phi_cos = acos((p1x + p2x) / pT_pair);
				phi_sin = asin((p1y + p2y) / pT_pair);
				// cout << "Cosine: " << (Double_t)acos((p1x + p2x) / pT_pair) << "  " << phi_cos << "\t Sin  ";
				// cout << (Double_t)asin((p1y + p2y) / pT_pair) << "  " << phi_sin << endl;
				if (phi_sin > 0)
					phi_pair = phi_cos;
				if (phi_sin < 0)
					phi_pair = -1. * phi_cos;

				Rs_pair = sqrt(pow(eta_pair - sourcePartEta, 2) + pow(phi_pair - sourcePartPhi, 2));
				// cout << "Rs_tr1 = " << Rs_tr1 << ", Rs_tr2 =  " << Rs_tr2 << ",  Rs_pair = " << Rs_pair << endl;
				// if (Rs_pair > 0.5 || Minv > 4.0 || fabs(eta_pair) > 1.0 || pT_pair > 15.0)
				// Rs cut is causing huge mismatch

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

				if (Minv > 4.0 || fabs(eta_pair) > 1.0 || pT_pair > 15.0)
					continue;

				if (isJP0)
				{
					hxsecMJP0->Fill(Minv, ptw);

					hCone[0]->Fill(cone, ptw);
					hMinv[0]->Fill(Minv, ptw);
					hEtaPair[0]->Fill(eta_pair, ptw);
					hpTPair[0]->Fill(pT_pair, ptw);
					hPhiRY[0]->Fill(PhiRY, ptw);
					hPhiRB[0]->Fill(PhiRB, ptw);
				}
				if (isJP1)
				{
					hxsecMJP1->Fill(Minv, ptw);

					hCone[1]->Fill(cone, ptw);
					hMinv[1]->Fill(Minv, ptw);
					hEtaPair[1]->Fill(eta_pair, ptw);
					hpTPair[1]->Fill(pT_pair, ptw);
					hPhiRY[1]->Fill(PhiRY, ptw);
					hPhiRB[1]->Fill(PhiRB, ptw);
				}
				if (isJP2)
				{
					hxsecMJP2->Fill(Minv, ptw);

					hCone[2]->Fill(cone, ptw);
					hMinv[2]->Fill(Minv, ptw);
					hEtaPair[2]->Fill(eta_pair, ptw);
					hpTPair[2]->Fill(pT_pair, ptw);
					hPhiRY[2]->Fill(PhiRY, ptw);
					hPhiRB[2]->Fill(PhiRB, ptw);
				}

			} // track 2 loop
		}	  // track loop
		nEventsProcessed++;

	} // event loop
	cout << "Events Passed = " << nentries << "\n";
	cout << "Events Processed = " << nEventsProcessed << endl;
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

	fChain->SetBranchAddress("trigJP0", &trigJP0, &b_trigJP0);
	fChain->SetBranchAddress("trigJP1", &trigJP1, &b_trigJP1);
	fChain->SetBranchAddress("trigJP2", &trigJP2, &b_trigJP2);

	fChain->SetBranchAddress("fVZ", &fVZ, &b_fVZ);
	// fChain->SetBranchAddress("fevTime", &fevTime, &b_fevTime);
	fChain->SetBranchAddress("fverRank", &fverRank, &b_fverRank);
	fChain->SetBranchAddress("fpT", fpT, &b_fpT);
	fChain->SetBranchAddress("fp", fp, &b_fp);
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

	fChain->SetBranchAddress("idTruth", idTruth, &b_idTruth);
	// fChain->SetBranchAddress("mcId", mcId, &b_mcId);
	fChain->SetBranchAddress("fVZ_mc", &fVZ_mc, &b_fVZ_mc);
	fChain->SetBranchAddress("partonicPtBin", &partonicPtBin, &b_partonicPtBin);
	fChain->SetBranchAddress("fpT_pyth", fpT_pyth, &b_fpT_pyth);
	fChain->SetBranchAddress("fp_pyth", fp_pyth, &b_fp_pyth);
	fChain->SetBranchAddress("feta_pyth", feta_pyth, &b_feta_pyth);
	fChain->SetBranchAddress("fphi_pyth", fphi_pyth, &b_fphi_pyth);
	fChain->SetBranchAddress("fpId_pyth", fpId_pyth, &b_fpId_pyth);

	fChain->SetBranchAddress("bPartId", &bPartId, &b_bPartId);
	fChain->SetBranchAddress("bPartX1", &bPartX1, &b_bPartX1);
	fChain->SetBranchAddress("bPartE", &bPartE, &b_bPartE);
	fChain->SetBranchAddress("bPartPhi", &bPartPhi, &b_bPartPhi);
	fChain->SetBranchAddress("bPartEta", &bPartEta, &b_bPartEta);
	fChain->SetBranchAddress("bPartPt", &bPartPt, &b_bPartPt);

	fChain->SetBranchAddress("yPartId", &yPartId, &b_yPartId);
	fChain->SetBranchAddress("yPartX2", &yPartX2, &b_yPartX2);
	fChain->SetBranchAddress("yPartE", &yPartE, &b_yPartE);
	fChain->SetBranchAddress("yPartPhi", &yPartPhi, &b_yPartPhi);
	fChain->SetBranchAddress("yPartEta", &yPartEta, &b_yPartEta);
	fChain->SetBranchAddress("yPartPt", &yPartPt, &b_yPartPt);

	fChain->SetBranchAddress("fPart1Id", &fPart1Id, &b_fPart1Id);
	fChain->SetBranchAddress("fPart1X", &fPart1X, &b_fPart1X);
	fChain->SetBranchAddress("fPart1E", &fPart1E, &b_fPart1E);
	fChain->SetBranchAddress("fPart1Phi", &fPart1Phi, &b_fPart1Phi);
	fChain->SetBranchAddress("fPart1Eta", &fPart1Eta, &b_fPart1Eta);
	fChain->SetBranchAddress("fPart1Pt", &fPart1Pt, &b_fPart1Pt);

	fChain->SetBranchAddress("fPart2Id", &fPart2Id, &b_fPart2Id);
	fChain->SetBranchAddress("fPart2X", &fPart2X, &b_fPart2X);
	fChain->SetBranchAddress("fPart2E", &fPart2E, &b_fPart2E);
	fChain->SetBranchAddress("fPart2Phi", &fPart2Phi, &b_fPart2Phi);
	fChain->SetBranchAddress("fPart2Eta", &fPart2Eta, &b_fPart2Eta);
	fChain->SetBranchAddress("fPart2Pt", &fPart2Pt, &b_fPart2Pt);

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
	hpartonicPt->SetDirectory(fout);
	hpartonicPt->Write();

	hxsecMJP0->SetDirectory(fout);
	hxsecMJP0->Write();
	hxsecMJP1->SetDirectory(fout);
	hxsecMJP1->Write();
	hxsecMJP2->SetDirectory(fout);
	hxsecMJP2->Write();

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
		hPhiRB[i]->SetDirectory(fout);
		hPhiRB[i]->Write();
		hPhiRY[i]->SetDirectory(fout);
		hPhiRY[i]->Write();

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
