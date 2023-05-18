#include "Run12pp200Ana.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include <vector>
#include <algorithm>

using namespace std;
ofstream Output;

void Run12pp200Ana::Loop()
{

	gROOT->Reset();
	TH1::SetDefaultSumw2();

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

	// initialize histograms
	Int_t nTrg = 3;
	Int_t nCh = 2;
	const char *tTrig[3] = {"JP0", "JP1", "JP2"};
	const char *tCh[2] = {"Pos", "Neg"};

	//  updated x-section bins
	const int xnBins = 13;
	Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

	hMinvJP0 = new TH1D("hMinvJP0", "", 200, 0.2, 4.05);
	hMinvJP1 = new TH1D("hMinvJP1", "", 200, 0.2, 4.05);
	hMinvJP2 = new TH1D("hMinvJP2", "", 200, 0.2, 4.05);

	// histograms for trigger efficiency
	hBiasedMinvJP0 = new TH1D("hBiasedMinvJP0", "", xnBins, xnBinsEdges);
	hBiasedMinvJP1 = new TH1D("hBiasedMinvJP1", "", xnBins, xnBinsEdges);
	hBiasedMinvJP2 = new TH1D("hBiasedMinvJP2", "", xnBins, xnBinsEdges);

	// histograms for trigger bias
	// all triggers included
	hQuarksAll = new TH1D("hQuarksAll", "", xnBins, xnBinsEdges);
	hGluonsAll = new TH1D("hGluonsAll", "", xnBins, xnBinsEdges);
	hPartonsAll = new TH1D("hPartonsAll", "", xnBins, xnBinsEdges);
	// jp0 trigger
	hQuarksJP0 = new TH1D("hQuarksJP0", "", xnBins, xnBinsEdges);
	hGluonsJP0 = new TH1D("hGluonsJP0", "", xnBins, xnBinsEdges);
	hPartonsJP0 = new TH1D("hPartonsJP0", "", xnBins, xnBinsEdges);
	// jp1 trigger
	hQuarksJP1 = new TH1D("hQuarksJP1", "", xnBins, xnBinsEdges);
	hGluonsJP1 = new TH1D("hGluonsJP1", "", xnBins, xnBinsEdges);
	hPartonsJP1 = new TH1D("hPartonsJP1", "", xnBins, xnBinsEdges);
	// jp2 trigger
	hQuarksJP2 = new TH1D("hQuarksJP2", "", xnBins, xnBinsEdges);
	hGluonsJP2 = new TH1D("hGluonsJP2", "", xnBins, xnBinsEdges);
	hPartonsJP2 = new TH1D("hPartonsJP2", "", xnBins, xnBinsEdges);
	// jp trigger
	hQuarksJP = new TH1D("hQuarksJP", "", xnBins, xnBinsEdges);
	hGluonsJP = new TH1D("hGluonsJP", "", xnBins, xnBinsEdges);
	hPartonsJP = new TH1D("hPartonsJP", "", xnBins, xnBinsEdges);

	// reconstructed variables
	Double_t p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1;
	Double_t p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
	Double_t cosPhiRB, sinPhiRB, cosPhiRY, sinPhiRY;
	Double_t PhiR_cosB, PhiR_sinB, PhiRB, PhiR_cosY, PhiR_sinY, PhiRY, phiRB, phiDiffB, phi_cos, phi_sin, phi_pair;
	Double_t msqr, msqr1, msqr2;																			  //
	Double_t R1_p1, R2_p1, R1_p2, R2_p2, Rs_tr1, Rs_tr2, sourceParton, sourcePartEta, sourcePartPhi, Rs_pair; // opening angles between parton and tracks and associated parton eta, phi

	// pythia varibales
	Double_t p1x_mc, p2x_mc, p1y_mc, p2y_mc, p1z_mc, p2z_mc, psx_mc, psy_mc, psz_mc, ps_mc, cone_mc, R_mc, Rs_mc, Rx_mc, Ry_mc, Rz_mc, R1_mc, Rx1_mc, Ry1_mc, Rz1_mc;
	Double_t p1_mc, p2_mc, E1_mc, E2_mc, Minv_mc, pT_pair_mc, eta_pair_mc, phi_pair_mc;
	Double_t sourceParton_mc, partonE_mc, E_pair_mc;

	Double_t pi = 3.14159265359;
	Double_t m_pion = 0.1396;	  // GeV
	Bool_t pairMatch = kFALSE;	  // Rs<0.5, Minv<4, and pT_pair<15, pair-source parton matching condition
	Bool_t pairMatch_mc = kFALSE; // Rs<0.5, Minv<4, and pT_pair<15, pair-source parton matching condition
	Double_t pairSep;
	Bool_t geLimit = kFALSE;
	Bool_t pyLimit = kFALSE;

	if (fChain == 0)
		return;

	// calculate weight for each partonic pT bins
	const Int_t nptBin = 13;
	Double_t partPtBin[nptBin + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., 80}; // for trigger eff
	Double_t partPtRange[nptBin + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., -1};
	Double_t fudge[nptBin] = {1.228, 1.051, 1.014, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	Double_t crossSection[nptBin] = {9.0012e+00, 1.46253e+00, 3.54466e-01, 1.51622e-01, 2.49062e-02, 5.84527e-03, 2.30158e-03, 3.42755e-04, 4.57002e-05, 9.72535e-06, 4.69889e-07, 2.69202e-08, 1.43453e-09};
	Double_t numEvents[nptBin] = {3318626, 3301413, 3291662, 3280010, 3282543, 3275693, 3276437, 3276795, 3272804, 2179660, 2183230, 1091927, 1090857}; // number of events in file after processing MuDsts
	Double_t binLumi[nptBin] = {0}, binWt[nptBin] = {0};
	for (int mBin = 0; mBin < nptBin; mBin++)
	{
		binLumi[mBin] = (numEvents[mBin] * crossSection[12]) / (crossSection[mBin] * numEvents[12]); // bin luminosity
		binWt[mBin] = 1. / (binLumi[mBin] * fudge[mBin]);											 // Bin weight for partonic pT weighting.
	}

	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	cout << "Event Loop started. Number of events = " << nentries << endl;
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	// for (Long64_t jentry = 0; jentry < 100; jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		Double_t ptw = 0;
		Double_t extra_weight = 0;

		// soft reweighting from Dmitry
		Double_t p[4] = {1.21733561, -0.32633577, 0.16913723, 0.82134143};
		extra_weight = 1. / (1. + (p[0] + p[1] * (partonicPtBin - 2.) + p[2] * (partonicPtBin - 2.) * (partonicPtBin - 2.)) * exp(-p[3] * (partonicPtBin - 2.)));

		for (int j = 0; j < 13; j++)
		{
			if (j < 12 && partonicPtBin >= partPtRange[j] && partonicPtBin < partPtRange[j + 1])
			{
				ptw = binWt[j] * extra_weight;
			}
			if (j == 12 && partonicPtBin >= partPtRange[j])
			{
				ptw = binWt[j] * extra_weight;
			}
		}

		Bool_t vertexCut = kFALSE;
		Bool_t isJP = kFALSE;
		Bool_t isJP0 = kFALSE;
		Bool_t isJP1 = kFALSE;
		Bool_t isJP2 = kFALSE;

		if (fverRank > 1e6 && fabs(fVZ) < 60)
			vertexCut = kTRUE; // common event vertex selection cut
		if (!vertexCut)
			continue;

		for (unsigned int i = 0; i < ftrigger->size(); i++)
		{
			if (ftrigger->at(i) == 370601)
				isJP0 = true;
			if (ftrigger->at(i) == 370611)
				isJP1 = true;
			if (ftrigger->at(i) == 370621)
				isJP2 = true;
		}

		for (int tr = 0; tr < fmaxpar; tr++)
		{
			if (idTruth[tr] <= 0 || idTruth[tr] != fId_mc[tr])
				continue;
			if (fpT_mc[tr] < 0.5 || fpT_mc[tr] > 15. || fabs(feta_mc[tr]) > 1. || !(fpId_mc[tr] == 8 || fpId_mc[tr] == 9))
				continue;
			//     Track 2 loop
			for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
			{

				if (idTruth[tr2] <= 0 || idTruth[tr2] != fId_mc[tr2])
					continue;

				if (fpT_mc[tr2] < 0.5 || fpT_mc[tr2] > 15. || fabs(feta_mc[tr2]) > 1. || !(fpId_mc[tr2] == 8 || fpId_mc[tr2] == 9))
					continue;

				if (fpId_mc[tr] == fpId_mc[tr2])
					continue;

				// cout << "Pair IDs = " << fpId_mc[tr] << "\t" << fpId_mc[tr2] << endl;

				Double_t phiDiff = fphi_mc[tr] - fphi_mc[tr2];
				if (phiDiff > pi)
					phiDiff -= (2 * pi);
				if (phiDiff < ((-1) * pi))
					phiDiff += (2 * pi);
				if (phiDiff > pi || phiDiff < -pi)
					continue;

				cone = sqrt(pow(feta_mc[tr] - feta_mc[tr2], 2) + pow(phiDiff, 2));
				if (cone >= Cone_Max7 || cone < Cone_Min)
					continue;

				if (fpId_mc[tr] == 8)
				{
					p1x = fpT_mc[tr] * cos(fphi_mc[tr]);
					p2x = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p1y = fpT_mc[tr] * sin(fphi_mc[tr]);
					p2y = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p1z = fpT_mc[tr] * sinh(feta_mc[tr]);
					p2z = fpT_mc[tr2] * sinh(feta_mc[tr2]);
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
					// cout << "Track 1 +ve, " << p1x << ", " << p1y << " , " << p1z << endl;
				}
				if (fpId_mc[tr] == 9)
				{
					p1x = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p2x = fpT_mc[tr] * cos(fphi_mc[tr]);
					p1y = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p2y = fpT_mc[tr] * sin(fphi_mc[tr]);
					p1z = fpT_mc[tr2] * sinh(feta_mc[tr2]);
					p2z = fpT_mc[tr] * sinh(feta_mc[tr]);
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

				if (Minv > Minv_Max || Minv < Minv_Min || fabs(eta_pair) > Eta_Cut || pT_pair > ptPair_Max || pT_pair < ptPair_Min)
					continue;

				if (Minv > 1.5 && partonicPtBin < 5)
					continue; //bad event from lower pt bin having very high pt weight contributing higher mass bin(1 or 2 events).

				// match the dipion to the hard scattered initial state parton(status ==3). Initial state particles are stored in the first few entries along with the beam remanents and the underying event(mother id = 0,1,2)
				double minAngle = 999;
				int matchPartonId = 999;
				for (int i = 0; i < fmaxpar1; i++)
				{
					if (fStatusCode_pyth[i] != 3)
						continue; // selects only the hard scatterd partons
					double angle = sqrt(pow(feta_pyth[i] - eta_pair, 2) + pow(fphi_pyth[i] - phi_pair, 2));
					minAngle = (angle < minAngle) ? angle : minAngle;
					matchPartonId = (angle < minAngle) ? fpId_pyth[i] : matchPartonId;
				}

				Bool_t isParton = fabs(matchPartonId) <= 8 || matchPartonId == 21;
				Bool_t isQuark = fabs(matchPartonId) <= 8;
				Bool_t isGluon = matchPartonId == 21;
				Bool_t isJPtrig = (isJP0 || isJP1 || isJP2);

				if (isJP0)
					hMinvJP0->Fill(Minv, ptw);
				if (isJP1)
					hMinvJP1->Fill(Minv, ptw);
				if (isJP2)
					hMinvJP2->Fill(Minv, ptw);

				for (int i = 0; i < xnBins; i++)
				{
					if (Minv >= xnBinsEdges[i] && Minv < xnBinsEdges[i + 1])
					{

						if (isJP0)
							hBiasedMinvJP0->Fill(Minv, ptw);
						if (isJP1)
							hBiasedMinvJP1->Fill(Minv, ptw);
						if (isJP2)
							hBiasedMinvJP2->Fill(Minv, ptw);

						if (isQuark)
						{

							hQuarksAll->Fill(Minv, ptw);
							if (isJP0)
								hQuarksJP0->Fill(Minv, ptw);
							if (isJP1)
								hQuarksJP1->Fill(Minv, ptw);
							if (isJP2)
								hQuarksJP2->Fill(Minv, ptw);
							if (isJPtrig)
								hQuarksJP->Fill(Minv, ptw);
						}
						if (isGluon)
						{

							hGluonsAll->Fill(Minv, ptw);
							if (isJP0)
								hGluonsJP0->Fill(Minv, ptw);
							if (isJP1)
								hGluonsJP1->Fill(Minv, ptw);
							if (isJP2)
								hGluonsJP2->Fill(Minv, ptw);
							if (isJPtrig)
								hGluonsJP->Fill(Minv, ptw);
						}
						if (isParton)
						{

							hPartonsAll->Fill(Minv, ptw);
							if (isJP0)
								hPartonsJP0->Fill(Minv, ptw);
							if (isJP1)
								hPartonsJP1->Fill(Minv, ptw);
							if (isJP2)
								hPartonsJP2->Fill(Minv, ptw);
							if (isJPtrig)
								hPartonsJP->Fill(Minv, ptw);
						}
					}
				}
			} // 2nd track loop
		}	  // 1st track loop
	}		  // event loop
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
	fChain->SetBranchAddress("fmaxpar1", &fmaxpar1, &b_fmaxpar1);
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
	fChain->SetBranchAddress("idTruth", idTruth, &b_idTruth);
	fChain->SetBranchAddress("ffitPts", ffitPts, &b_ffitPts);
	fChain->SetBranchAddress("ffitPtsPoss", ffitPtsPoss, &b_ffitPtsPoss);
	fChain->SetBranchAddress("fhitsdedx", fhitsdedx, &b_fhitsdedx);
	fChain->SetBranchAddress("fBetaToF", fBetaToF, &b_fBetaToF);
	fChain->SetBranchAddress("fvpdVz", &fvpdVz, &b_fvpdVz);

	fChain->SetBranchAddress("fVZ_mc", &fVZ_mc, &b_fVZ_mc);
	fChain->SetBranchAddress("fpT_mc", fpT_mc, &b_fpT_mc);
	fChain->SetBranchAddress("fp_mc", fp_mc, &b_fp_mc);
	fChain->SetBranchAddress("feta_mc", feta_mc, &b_feta_mc);
	fChain->SetBranchAddress("fphi_mc", fphi_mc, &b_fphi_mc);
	fChain->SetBranchAddress("fpId_mc", fpId_mc, &b_fpId_mc);
	fChain->SetBranchAddress("fId_mc", fId_mc, &b_fId_mc);

	fChain->SetBranchAddress("partonicPtBin", &partonicPtBin, &b_partonicPtBin);
	fChain->SetBranchAddress("fStatusCode_pyth", fStatusCode_pyth, &b_fStatusCode_pyth);
	fChain->SetBranchAddress("fVZ_pyth", &fVZ_pyth, &b_fVZ_pyth);
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
	fout->cd();
	hMinvJP0->Write();
	hMinvJP1->Write();
	hMinvJP2->Write();

	hBiasedMinvJP0->Write();
	hBiasedMinvJP1->Write();
	hBiasedMinvJP2->Write();

	hQuarksAll->Write();
	hGluonsAll->Write();
	hPartonsAll->Write();

	hGluonsJP->Write();
	hPartonsJP->Write();
	hQuarksJP->Write();

	hQuarksJP0->Write();
	hGluonsJP0->Write();
	hPartonsJP0->Write();

	hQuarksJP1->Write();
	hGluonsJP1->Write();
	hPartonsJP1->Write();

	hQuarksJP2->Write();
	hGluonsJP2->Write();
	hPartonsJP2->Write();

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
