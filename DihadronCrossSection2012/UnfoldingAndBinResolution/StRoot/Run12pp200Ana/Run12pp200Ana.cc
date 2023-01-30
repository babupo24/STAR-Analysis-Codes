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
#include <algorithm>

using namespace std;
ofstream Output;

void Run12pp200Ana::Loop()
{

	gROOT->Reset();
	TH1::SetDefaultSumw2();
	//  initialize histograms
	Int_t nTrg = 3;
	Int_t nCh = 2;
	const char *tTrig[3] = {"JP0", "JP1", "JP2"};
	const char *tCh[2] = {"Pos", "Neg"};

	hpartonicPt = new TH1D("hpartonicPt", "", 100, 0, 80);

	// histograms for resolution studies: 12 invariant mass bins with the Minv wifth 0.3 GeV/c^2
	const int nBins = 12;
	Double_t nBinsEdges[nBins + 1] = {0.26, 0.56, 0.86, 1.16, 1.46, 1.76, 2.06, 2.36, 2.66, 2.96, 3.26, 3.56, 4.0};

	const int xnBins = 13;
	Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

	// Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};

	hVZJP0 = new TH1D("hVZ_JP0", "", 100, -65., 65.);
	hVZJP1 = new TH1D("hVZ_JP1", "", 100, -65., 65.);
	hVZJP2 = new TH1D("hVZ_JP2", "", 100, -65., 65.);

	hxsecMinvGen = new TH1D("hxsecMinvGen", "", xnBins, xnBinsEdges);

	hMDiffInt = new TH1D("hMDiffInt", "", 100, -2.0, 2.0);
	for (int i = 0; i < 12; i++)
	{
		hMGenRec[i] = new TH1D(Form("hMGenRec%i", i), "", 100, -2.0, 2.0);
		hMGenRecU[i] = new TH1D(Form("hMGenRecU%i", i), "", 100, -2.0, 2.0);
	}
	// migration matrix with underflow filled. For a reconstructed signal, gen mass > 0. For a reconstructed background, gen mass < 0 (goes to underflow)
	hMGenVsRecJP0 = new TH2D("hMGenVsRecJP0", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP1 = new TH2D("hMGenVsRecJP1", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP2 = new TH2D("hMGenVsRecJP2", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	// migration matrix without underflow/overflow. Reconstructed background will be filled separately and subtracted while unfolding
	hMGenVsRecJP0NoUF = new TH2D("hMGenVsRecJP0NoUF", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP1NoUF = new TH2D("hMGenVsRecJP1NoUF", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP2NoUF = new TH2D("hMGenVsRecJP2NoUF", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	// migration matrix in pseudorapidity bins
	hMGenVsRecJP0Gt = new TH2D("hMGenVsRecJP0Gt", "eta > 0", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP1Gt = new TH2D("hMGenVsRecJP1Gt", "eta > 0", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP2Gt = new TH2D("hMGenVsRecJP2Gt", "eta > 0", 160, 0.27, 4.0, xnBins, xnBinsEdges);

	hMGenVsRecJP0Lt = new TH2D("hMGenVsRecJP0Lt", "eta < 0", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP1Lt = new TH2D("hMGenVsRecJP1Lt", "eta < 0", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMGenVsRecJP2Lt = new TH2D("hMGenVsRecJP2Lt", "eta < 0", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	// histograms for background
	// bkg 1: only one track matched in a pair
	hMRecJP2bk1 = new TH1D("hMRecJP2bk1", "", 160, 0.27, 4.0);
	hMRecJP1bk1 = new TH1D("hMRecJP1bk1", "", 160, 0.27, 4.0);
	hMRecJP0bk1 = new TH1D("hMRecJP0bk1", "", 160, 0.27, 4.0);
	// bkg 2 : both tracks unmatchecd
	hMRecJP2bk2 = new TH1D("hMRecJP2bk2", "", 160, 0.27, 4.0);
	hMRecJP1bk2 = new TH1D("hMRecJP1bk2", "", 160, 0.27, 4.0);
	hMRecJP0bk2 = new TH1D("hMRecJP0bk2", "", 160, 0.27, 4.0);
	// combined background in pseudorapidity bins
	hMRecJP2bkGt = new TH1D("hMRecJP2bkGt", "eta > 0", 160, 0.27, 4.0);
	hMRecJP1bkGt = new TH1D("hMRecJP1bkGt", "eta > 0", 160, 0.27, 4.0);
	hMRecJP0bkGt = new TH1D("hMRecJP0bkGt", "eta > 0", 160, 0.27, 4.0);

	hMRecJP2bkLt = new TH1D("hMRecJP2bkLt", "eta < 0", 160, 0.27, 4.0);
	hMRecJP1bkLt = new TH1D("hMRecJP1bkLt", "eta < 0", 160, 0.27, 4.0);
	hMRecJP0bkLt = new TH1D("hMRecJP0bkLt", "eta < 0", 160, 0.27, 4.0);
	// combined background in x-sec bins
	hMRecJP2bkx = new TH1D("hMRecJP2bkx", "bkg in x-sec bins", xnBins, xnBinsEdges);
	hMRecJP1bkx = new TH1D("hMRecJP1bkx", "", xnBins, xnBinsEdges);
	hMRecJP0bkx = new TH1D("hMRecJP0bkx", "", xnBins, xnBinsEdges);
	// control plots: generated mass distribution
	hMGenJP2 = new TH1D("hMGenJP2", "", xnBins, xnBinsEdges);
	hMGenJP1 = new TH1D("hMGenJP1", "", xnBins, xnBinsEdges);
	hMGenJP0 = new TH1D("hMGenJP0", "", xnBins, xnBinsEdges);

	hMGenJP2f = new TH1D("hMGenJP2f", "", 160, 0.27, 4.0); // generated in fine bins
	hMGenJP1f = new TH1D("hMGenJP1f", "", 160, 0.27, 4.0);
	hMGenJP0f = new TH1D("hMGenJP0f", "", 160, 0.27, 4.0);

	//  reconstructed mass distribution associated with the generated pairs
	hMRecJP2 = new TH1D("hMRecJP2", "fine binning", 160, 0.27, 4.0);
	hMRecJP1 = new TH1D("hMRecJP1", "fine binning", 160, 0.27, 4.0);
	hMRecJP0 = new TH1D("hMRecJP0", "fine binning", 160, 0.27, 4.0);

	hMRecJP2x = new TH1D("hMRecJP2x", "x-sec binning", xnBins, xnBinsEdges);
	hMRecJP1x = new TH1D("hMRecJP1x", "x-sec binning", xnBins, xnBinsEdges);
	hMRecJP0x = new TH1D("hMRecJP0x", "x-sec binning", xnBins, xnBinsEdges);
	// reconstructed mass distributions including background. A fraction of detector level events associated with the generator level events is considered as a matching ratio, which will be multiplied to the real data spectrum before unfolding. This is similar to the background correction.
	hMRecJP2All = new TH1D("hMRecJP2All", "fine binning", 160, 0.27, 4.0);
	hMRecJP1All = new TH1D("hMRecJP1All", "fine binning", 160, 0.27, 4.0);
	hMRecJP0All = new TH1D("hMRecJP0All", "fine binning", 160, 0.27, 4.0);

	hMRecJP2AllLt = new TH1D("hMRecJP2AllLt", "fine binning", 160, 0.27, 4.0);
	hMRecJP1AllLt = new TH1D("hMRecJP1AllLt", "fine binning", 160, 0.27, 4.0);
	hMRecJP0AllLt = new TH1D("hMRecJP0AllLt", "fine binning", 160, 0.27, 4.0);

	hMRecJP2AllGt = new TH1D("hMRecJP2AllGt", "fine binning", 160, 0.27, 4.0);
	hMRecJP1AllGt = new TH1D("hMRecJP1AllGt", "fine binning", 160, 0.27, 4.0);
	hMRecJP0AllGt = new TH1D("hMRecJP0AllGt", "fine binning", 160, 0.27, 4.0);

	hMRecJP2AllXbin = new TH1D("hMRecJP2Allx", "x-sec binning", xnBins, xnBinsEdges);
	hMRecJP1AllXbin = new TH1D("hMRecJP1Allx", "x-sec binning", xnBins, xnBinsEdges);
	hMRecJP0AllXbin = new TH1D("hMRecJP0Allx", "x-sec binning", xnBins, xnBinsEdges);

	// histograms for unfolding test. Half of the reconstructed events will be treated as real data and the other half will be treated as embedding. The first half, the real data will be the input for the unfolding and migration matrix will be construced from the embedding part.
	hMGenVsRecTestM = new TH2D("hMGenVsRecTestM", "", 160, 0.27, 4.0, xnBins, xnBinsEdges);
	hMRecTestInput = new TH1D("hMRecTestInput", "", 160, 0.27, 4.0);
	hMGenTestInput = new TH1D("hMGenTestInput", "", 160, 0.27, 4.0);
	hMGenTestControlPlot1st = new TH1D("hMGenTestControlPlot1st", "", xnBins, xnBinsEdges);
	hMRecTestControlPlot1st = new TH1D("hMRecTestControlPlot1st", "", xnBins, xnBinsEdges);
	hMRecTestControlPlot2nd = new TH1D("hMRecTestControlPlot2nd", "", xnBins, xnBinsEdges);
	hMGenTestControlPlot2nd = new TH1D("hMGenTestControlPlot2nd", "", xnBins, xnBinsEdges);

	// reconstructed variables
	Double_t p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1;
	Double_t p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
	Double_t cosPhiRB, sinPhiRB, cosPhiRY, sinPhiRY;
	Double_t PhiR_cosB, PhiR_sinB, PhiRB, PhiR_cosY, PhiR_sinY, PhiRY, phiRB, phiDiffB, phi_cos, phi_sin, phi_pair;
	Double_t msqr, msqr1, msqr2;																//
	Double_t R1_p1, R2_p1, R1_p2, R2_p2, Rs_tr1, Rs_tr2, sourcePartEta, sourcePartPhi, Rs_pair; // opening angles between parton and tracks and associated parton eta, phi

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
	Long64_t bkgCounts;

	if (fChain == 0)
		return;

	// calculate weight for each partonic pT bins
	const Int_t nptBin = 13;
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
		// parton pT
		hpartonicPt->Fill(partonicPtBin, ptw);

		if (fverRank < 1e6)
			continue;
		if (fabs(fVZ) > 60.)
			continue;

		// if (fabs(fVZ - fVZ_mc) > 2.)
		//	continue;
		// flags to divide the events in two halves for unfolding test
		Bool_t is1st_Half = kFALSE;
		Bool_t is2nd_Half = kFALSE;

		is1st_Half = jentry <= (nentries / 2);
		is2nd_Half = jentry > (nentries / 2);

		// cout << jentry << " " << is1st_Half << " " << is2nd_Half << endl;

		Bool_t isJP = kFALSE;
		Bool_t isJP0 = kFALSE;
		Bool_t isJP1 = kFALSE;
		Bool_t isJP2 = kFALSE;

		if (trigJP0 == 1)
			isJP0 = kTRUE;
		if (trigJP1 == 1)
			isJP1 = kTRUE;
		if (trigJP2 == 1)
			isJP2 = kTRUE;

		if (trigJP0 == 1 || trigJP1 == 1 || trigJP2 == 1)
			isJP = kTRUE;

		if (isJP0)
			hVZJP0->Fill(fVZ, ptw);
		if (isJP1)
			hVZJP1->Fill(fVZ, ptw);
		if (isJP2)
			hVZJP2->Fill(fVZ, ptw);

		if (fPart1Phi > pi)
			fPart1Phi = fPart1Phi - 2. * pi;
		if (fPart2Phi > pi)
			fPart2Phi = fPart2Phi - 2. * pi;
		bkgCounts = 0;
		Bool_t ismc_tr1 = kFALSE;
		Bool_t isrec_tr1 = kFALSE;

		for (int tr = 0; tr < fmaxpar; tr++)
		{
			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];

			if (fpT[tr] > 0.5 && fpT[tr] < 15. && ffitPts[tr] > 15 && fabs(feta[tr]) < 1. && fitPtsRatio_tr > .51 && dcaCut(fpT[tr], fdca[tr]) && fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2.0) // track 1 quality cuts
				isrec_tr1 = kTRUE;
			if (idTruth[tr] > 0 && fpT_mc[tr] > 0.5 && fpT_mc[tr] < 15 && fabs(feta_mc[tr]) < 1. && (fId_mc[tr] == idTruth[tr])) // 2,3 = e+-, 5,6 = mu+-, 8,9=pi+-, 11,12= K+-, 14,15=ppbar
				ismc_tr1 = kTRUE;

			Bool_t isrec_tr2 = kFALSE;
			Bool_t ismc_tr2 = kFALSE;

			//   Track 2 loop
			for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
			{
				Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];

				if (fpT[tr2] > 0.5 && fpT[tr2] < 15. && ffitPts[tr2] > 15 && fabs(feta[tr2]) < 1. && fitPtsRatio_tr2 > .51 && dcaCut(fpT[tr2], fdca[tr2]) && fnSigmaPion[tr2] > -1.0 && fnSigmaPion[tr2] < 2.0)
					isrec_tr2 = kTRUE;
				if (idTruth[tr2] > 0 && fpT_mc[tr2] > 0.5 && fpT_mc[tr2] < 15 && fabs(feta_mc[tr2]) < 1. && (fId_mc[tr2] == idTruth[tr2])) // 2,3 = e+-, 5,6 = mu+-, 8,9=pi+-, 11,12= K+-,
					ismc_tr2 = kTRUE;

				Double_t phiDiff = fphi[tr] - fphi[tr2];
				if (phiDiff > pi)
					phiDiff -= (2 * pi);
				if (phiDiff < ((-1) * pi))
					phiDiff += (2 * pi);

				if (fcharge[tr] == fcharge[tr2])
					continue;

				cone = sqrt(pow(feta[tr] - feta[tr2], 2) + pow(phiDiff, 2));

				if (cone >= 0.7 || !isrec_tr1 || !isrec_tr2)
					continue;
				if (fcharge[tr] > 0)
				{
					// rec. pion
					p1x = fpT[tr] * cos(fphi[tr]);
					p2x = fpT[tr2] * cos(fphi[tr2]);
					p1y = fpT[tr] * sin(fphi[tr]);
					p2y = fpT[tr2] * sin(fphi[tr2]);
					p1z = fpT[tr] * sinh(feta[tr]);
					p2z = fpT[tr2] * sinh(feta[tr2]);
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
				}
				if (fcharge[tr] < 0)
				{
					// rec. pion
					p1x = fpT[tr2] * cos(fphi[tr2]);
					p2x = fpT[tr] * cos(fphi[tr]);
					p1y = fpT[tr2] * sin(fphi[tr2]);
					p2y = fpT[tr] * sin(fphi[tr]);
					p1z = fpT[tr2] * sinh(feta[tr2]);
					p2z = fpT[tr] * sinh(feta[tr]);
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
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
				// phi_cos = acos((pt1ge * cos(phi1ge) + pt2ge * cos(phi2ge)) / pT_pair);
				// phi_sin = asin((pt1ge * sin(phi1ge) + pt2ge * sin(phi2ge)) / pT_pair);
				phi_cos = acos((p1x + p2x) / pT_pair);
				phi_sin = asin((p1y + p2y) / pT_pair);

				if (phi_sin > 0)
					phi_pair = phi_cos;
				if (phi_sin < 0)
					phi_pair = -1. * phi_cos;

				Rs_pair = sqrt(pow(eta_pair - sourcePartEta, 2) + pow(phi_pair - sourcePartPhi, 2));

				if (Minv > 4.0 || fabs(eta_pair) > 1.0 || pT_pair > 15.0)
					continue;

				// histograms for matching ratios: includes background
				if (trigJP0 == 1)
				{
					hMRecJP0All->Fill(Minv, ptw);
					if (eta_pair > 0)
						hMRecJP0AllGt->Fill(Minv, ptw);
					if (eta_pair < 0)
						hMRecJP0AllLt->Fill(Minv, ptw);

					hMRecJP0AllXbin->Fill(Minv, ptw);
				}
				if (trigJP1 == 1)
				{
					hMRecJP1All->Fill(Minv, ptw);
					if (eta_pair > 0)
						hMRecJP1AllGt->Fill(Minv, ptw);
					if (eta_pair < 0)
						hMRecJP1AllLt->Fill(Minv, ptw);

					hMRecJP1AllXbin->Fill(Minv, ptw);
				}
				if (trigJP2 == 1)
				{
					hMRecJP2All->Fill(Minv, ptw);
					if (eta_pair > 0)
						hMRecJP2AllGt->Fill(Minv, ptw);
					if (eta_pair < 0)
						hMRecJP2AllLt->Fill(Minv, ptw);

					hMRecJP2AllXbin->Fill(Minv, ptw);
				}

				// generated pion
				if (ismc_tr1 && ismc_tr2)
				{
					Double_t phiDiff_mc = fphi_mc[tr] - fphi_mc[tr2];
					if (phiDiff_mc > pi)
						phiDiff_mc -= (2 * pi);
					if (phiDiff_mc < ((-1) * pi))
						phiDiff_mc += (2 * pi);

					cone_mc = sqrt(pow(feta_mc[tr] - feta_mc[tr2], 2) + pow(phiDiff_mc, 2));

					if (fcharge[tr] > 0)
					{
						p1x_mc = fpT_mc[tr] * cos(fphi_mc[tr]);
						p2x_mc = fpT_mc[tr2] * cos(fphi_mc[tr2]);
						p1y_mc = fpT_mc[tr] * sin(fphi_mc[tr]);
						p2y_mc = fpT_mc[tr2] * sin(fphi_mc[tr2]);
						p1z_mc = fpT_mc[tr] * sinh(feta_mc[tr]);
						p2z_mc = fpT_mc[tr2] * sinh(feta_mc[tr2]);

						p1_mc = sqrt(p1x_mc * p1x_mc + p1y_mc * p1y_mc + p1z_mc * p1z_mc);
						p2_mc = sqrt(p2x_mc * p2x_mc + p2y_mc * p2y_mc + p2z_mc * p2z_mc);
					}

					if (fcharge[tr] < 0) //  pi-
					{
						p1x_mc = fpT_mc[tr2] * cos(fphi_mc[tr2]);
						p2x_mc = fpT_mc[tr] * cos(fphi_mc[tr]);
						p1y_mc = fpT_mc[tr2] * sin(fphi_mc[tr2]);
						p2y_mc = fpT_mc[tr] * sin(fphi_mc[tr]);
						p1z_mc = fpT_mc[tr2] * sinh(feta_mc[tr2]);
						p2z_mc = fpT_mc[tr] * sinh(feta_mc[tr]);

						p1_mc = sqrt(p1x_mc * p1x_mc + p1y_mc * p1y_mc + p1z_mc * p1z_mc);
						p2_mc = sqrt(p2x_mc * p2x_mc + p2y_mc * p2y_mc + p2z_mc * p2z_mc);
					}
					// generated pair
					psx_mc = p1x_mc + p2x_mc;
					psy_mc = p1y_mc + p2y_mc;
					psz_mc = p1z_mc + p2z_mc;
					ps_mc = sqrt(psx_mc * psx_mc + psy_mc * psy_mc + psz_mc * psz_mc);
					Rx_mc = (p1x_mc - p2x_mc);
					Ry_mc = (p1y_mc - p2y_mc);
					Rz_mc = (p1z_mc - p2z_mc);
					R_mc = sqrt(Rx_mc * Rx_mc + Ry_mc * Ry_mc + Rz_mc * Rz_mc); // R and R1 are same
					E1_mc = sqrt(m_pion * m_pion + p1_mc * p1_mc);
					E2_mc = sqrt(m_pion * m_pion + p2_mc * p2_mc);
					Minv_mc = sqrt(2 * m_pion * m_pion + 2 * (E1_mc * E2_mc - p1x_mc * p2x_mc - p1y_mc * p2y_mc - p1z_mc * p2z_mc));
					eta_pair_mc = TMath::ASinH((p1z_mc + p2z_mc) / pT_pair_mc);
					pT_pair_mc = sqrt((p1x_mc + p2x_mc) * (p1x_mc + p2x_mc) + (p1y_mc + p2y_mc) * (p1y_mc + p2y_mc));

					// if (Minv > 4 || fabs(eta_pair) > 1 || pT_pair > 15)
					//	continue;

					// histograms for bin resolution
					//  for resolution fill histogram for a pair coming from matched tracks
					hMDiffInt->Fill(Minv_mc - Minv);

					for (int i = 0; i < nBins; i++)
					{
						if (Minv >= nBinsEdges[i] && Minv < nBinsEdges[i + 1])
						{
							hMGenRec[i]->Fill(Minv_mc - Minv);
						}
					}
					if (Minv_mc < 4.0 && pT_pair_mc < 15.0 && fabs(eta_pair_mc) < 1.0)
					{
						// unfolding signal histograms: both tracks matched
						if (trigJP0 == 1)
						{
							hMRecJP0->Fill(Minv, ptw);					 // rec control plot
							hMRecJP0x->Fill(Minv, ptw);					 // rec control plot
							hMGenJP0->Fill(Minv_mc, ptw);				 // gen control plot
							hMGenJP0f->Fill(Minv_mc, ptw);				 // gen control plot
							hMGenVsRecJP0NoUF->Fill(Minv, Minv_mc, ptw); // migration matrix with no under/over flow
							hMGenVsRecJP0->Fill(Minv, Minv_mc, ptw);	 // migration matrix with over/under flow
							if (eta_pair > 0)
								hMGenVsRecJP0Gt->Fill(Minv, Minv_mc, ptw);
							if (eta_pair < 0)
								hMGenVsRecJP0Lt->Fill(Minv, Minv_mc, ptw);
						}
						if (trigJP1 == 1)
						{
							hMRecJP1->Fill(Minv, ptw);					 // rec control plot
							hMRecJP1x->Fill(Minv, ptw);					 // rec control plot
							hMGenJP1->Fill(Minv_mc, ptw);				 // gen control plot
							hMGenJP1f->Fill(Minv_mc, ptw);				 // gen control plot
							hMGenVsRecJP1NoUF->Fill(Minv, Minv_mc, ptw); // migration matrix
							hMGenVsRecJP1->Fill(Minv, Minv_mc, ptw);
							if (eta_pair > 0)
								hMGenVsRecJP1Gt->Fill(Minv, Minv_mc, ptw);
							if (eta_pair < 0)
								hMGenVsRecJP1Lt->Fill(Minv, Minv_mc, ptw);
						}
						if (trigJP2 == 1)
						{
							hMRecJP2->Fill(Minv, ptw);					 // rec control plot
							hMRecJP2x->Fill(Minv, ptw);					 // rec control plot
							hMGenJP2->Fill(Minv_mc, ptw);				 // gen control plot
							hMGenJP2f->Fill(Minv_mc, ptw);				 // gen control plot
							hMGenVsRecJP2NoUF->Fill(Minv, Minv_mc, ptw); // migration matrix
							hMGenVsRecJP2->Fill(Minv, Minv_mc, ptw);
							if (eta_pair > 0)
								hMGenVsRecJP2Gt->Fill(Minv, Minv_mc, ptw);
							if (eta_pair < 0)
								hMGenVsRecJP2Lt->Fill(Minv, Minv_mc, ptw);
						}
					}
				}
				if ((ismc_tr1 && !ismc_tr2) || (!ismc_tr1 && ismc_tr2))
				{
					// background type 1: one track matched in a pair
					if (trigJP0 == 1)
					{
						hMRecJP0bk1->Fill(Minv, ptw);
						hMRecJP0bkx->Fill(Minv, ptw);
						hMGenVsRecJP0->Fill(Minv, 0.0, ptw);
						if (eta_pair > 0)
							hMRecJP0bkGt->Fill(Minv, ptw);
						if (eta_pair < 0)
							hMRecJP0bkLt->Fill(Minv, ptw);
					}
					if (trigJP1 == 1)
					{
						hMRecJP1bk1->Fill(Minv, ptw);
						hMRecJP1bkx->Fill(Minv, ptw);
						hMGenVsRecJP1->Fill(Minv, 0.0, ptw);
						if (eta_pair > 0)
							hMRecJP1bkGt->Fill(Minv, ptw);
						if (eta_pair < 0)
							hMRecJP1bkLt->Fill(Minv, ptw);
					}
					if (trigJP2 == 1)
					{
						hMRecJP2bk1->Fill(Minv, ptw);
						hMRecJP2bkx->Fill(Minv, ptw);
						hMGenVsRecJP2->Fill(Minv, 0.0, ptw);
						if (eta_pair > 0)
							hMRecJP2bkGt->Fill(Minv, ptw);
						if (eta_pair < 0)
							hMRecJP2bkLt->Fill(Minv, ptw);
					}
				}
				if (!ismc_tr1 && !ismc_tr2)
				{ // background type 2: both tracks not matched
					if (trigJP0 == 1)
					{
						hMRecJP0bk2->Fill(Minv, ptw);
						hMRecJP0bkx->Fill(Minv, ptw);
						hMGenVsRecJP0->Fill(Minv, 0.0, ptw);
						// same histogram filled again to combine both backgrounds in a single histogram
						if (eta_pair > 0)
							hMRecJP0bkGt->Fill(Minv, ptw);
						if (eta_pair < 0)
							hMRecJP0bkLt->Fill(Minv, ptw);
					}

					if (trigJP1 == 1)
					{
						hMRecJP1bk1->Fill(Minv, ptw);
						hMRecJP1bkx->Fill(Minv, ptw);
						hMGenVsRecJP1->Fill(Minv, 0.0, ptw);
						if (eta_pair > 0)
							hMRecJP1bkGt->Fill(Minv, ptw);
						if (eta_pair < 0)
							hMRecJP1bkLt->Fill(Minv, ptw);
					}
					if (trigJP2 == 1)
					{
						hMRecJP2bk2->Fill(Minv, ptw);
						hMRecJP2bkx->Fill(Minv, ptw);
						hMGenVsRecJP2->Fill(Minv, 0.0, ptw);
						if (eta_pair > 0)
							hMRecJP2bkGt->Fill(Minv, ptw);
						if (eta_pair < 0)
							hMRecJP2bkLt->Fill(Minv, ptw);
					}
				}

				// fill histograms for unfolding test
				// fill the frst half events: treat this as data input
				if (is1st_Half && ismc_tr1 && ismc_tr2)
				{
					hMRecTestInput->Fill(Minv);
					hMGenTestInput->Fill(Minv_mc);
					hMRecTestControlPlot1st->Fill(Minv);
					hMGenTestControlPlot1st->Fill(Minv_mc);
				}
				// fill second half sample as embedidng
				if (is2nd_Half && ismc_tr1 && ismc_tr2)
				{
					hMGenVsRecTestM->Fill(Minv, Minv_mc);
					hMRecTestControlPlot2nd->Fill(Minv);
					hMGenTestControlPlot2nd->Fill(Minv_mc);
				}
			} // track 2 loop
		}	  // track loop
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
	// fChain->SetBranchAddress("isPrim_pyth", isPrim_pyth, &b_isPrim_pyth);
	// fChain->SetBranchAddress("statusCode_pyth", statusCode_pyth, &b_statusCode_pyth);
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

// routine match of all detector pions to the pythia pions......
void Run12pp200Ana::matchFinder(Double_t eta_ge, Double_t phi_ge, Double_t pt_ge, Double_t ch_ge, Int_t pyCounts, Double_t pid_mc[], Double_t eta_mc[], Double_t phi_mc[], Double_t pt_mc[], Int_t statusCode_mc[], Double_t &etage, Double_t &phige, Double_t &pTge, Double_t &chge, Double_t &mAngle, Double_t &etapy, Double_t &phipy, Double_t &pTpy, Double_t &chpy)
{
	Double_t ang;
	vector<Double_t> vAng, vEta_ge, vPhi_ge, vpT_ge, vCh_ge, vEta_mc, vPhi_mc, vpT_mc, vCh_mc;
	Double_t match = -1;
	int index;
	Double_t charge_mc = 0;
	for (int i = 0; i < pyCounts; i++)
	{
		if (statusCode_mc[i] != 1 && abs(pid_mc[i]) != 211)
			continue; // pythia pions
		if (pid_mc[i] == 211)
			charge_mc = 1;
		if (pid_mc[i] == -211)
			charge_mc = -1;
		ang = sqrt(pow(eta_ge - eta_mc[i], 2) + pow(phi_ge - phi_mc[i], 2));
		if (ang > 0.5)
			continue;
		vAng.push_back(ang);
		vpT_ge.push_back(pt_ge);
		vPhi_ge.push_back(phi_ge);
		vEta_ge.push_back(eta_ge);
		vCh_ge.push_back(ch_ge);
		vpT_mc.push_back(pt_mc[i]);
		vPhi_mc.push_back(phi_mc[i]);
		vEta_mc.push_back(eta_mc[i]);
		vCh_mc.push_back(charge_mc);
	}
	if (!vAng.empty())
	{
		match = *min_element(vAng.begin(), vAng.end());
		index = min_element(vAng.begin(), vAng.end()) - vAng.begin();

		pTge = vpT_ge[index];
		phige = vPhi_ge[index];
		etage = vEta_ge[index];
		chge = vCh_ge[index];
		mAngle = match;
		pTpy = vpT_mc[index];
		phipy = vPhi_mc[index];
		etapy = vEta_mc[index];
		chpy = vCh_mc[index];
	}
	vAng.clear();
	vEta_ge.clear();
	vPhi_ge.clear();
	vpT_ge.clear();
	vCh_ge.clear();
	vEta_mc.clear();
	vPhi_mc.clear();
	vpT_mc.clear();
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
	hpartonicPt->Write();

	hVZJP0->Write();
	hVZJP1->Write();
	hVZJP2->Write();

	hMDiffInt->Write();

	for (int i = 0; i < 12; i++)
	{
		hMGenRec[i]->Write();
		hMGenRecU[i]->Write();
	}
	hxsecMinvGen->Write();

	hMGenVsRecJP0->Write();
	hMGenVsRecJP1->Write();
	hMGenVsRecJP2->Write();

	hMGenVsRecJP0NoUF->Write();
	hMGenVsRecJP1NoUF->Write();
	hMGenVsRecJP2NoUF->Write();

	hMGenVsRecJP0Gt->Write();
	hMGenVsRecJP1Gt->Write();
	hMGenVsRecJP2Gt->Write();

	hMGenVsRecJP0Lt->Write();
	hMGenVsRecJP1Lt->Write();
	hMGenVsRecJP2Lt->Write();

	hMGenJP0->Write();
	hMGenJP1->Write();
	hMGenJP2->Write();

	hMGenJP0f->Write();
	hMGenJP1f->Write();
	hMGenJP2f->Write();

	hMRecJP0->Write();
	hMRecJP1->Write();
	hMRecJP2->Write();

	hMRecJP0x->Write();
	hMRecJP1x->Write();
	hMRecJP2x->Write();

	hMRecJP0All->Write();
	hMRecJP1All->Write();
	hMRecJP2All->Write();

	hMRecJP0AllLt->Write();
	hMRecJP1AllLt->Write();
	hMRecJP2AllLt->Write();

	hMRecJP0AllGt->Write();
	hMRecJP1AllGt->Write();
	hMRecJP2AllGt->Write();

	hMRecJP0AllXbin->Write();
	hMRecJP1AllXbin->Write();
	hMRecJP2AllXbin->Write();

	hMRecJP0bkGt->Write();
	hMRecJP1bkGt->Write();
	hMRecJP2bkGt->Write();

	hMRecJP0bkLt->Write();
	hMRecJP1bkLt->Write();
	hMRecJP2bkLt->Write();

	hMRecJP0bk1->Write();
	hMRecJP1bk1->Write();
	hMRecJP2bk1->Write();

	hMRecJP0bk2->Write();
	hMRecJP1bk2->Write();
	hMRecJP2bk2->Write();

	hMRecJP0bkx->Write();
	hMRecJP1bkx->Write();
	hMRecJP2bkx->Write();

	hMRecTestInput->Write();
	hMGenTestInput->Write();
	hMRecTestControlPlot1st->Write();
	hMGenTestControlPlot1st->Write();
	hMRecTestControlPlot2nd->Write();
	hMGenTestControlPlot2nd->Write();
	hMGenVsRecTestM->Write();

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
