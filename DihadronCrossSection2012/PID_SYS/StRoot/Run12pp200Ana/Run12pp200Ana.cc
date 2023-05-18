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
#include <iomanip>

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

	int ww = 20;
	Int_t nTrg = 3;
	Int_t nCh = 2;
	const char *tTrig[3] = {"JP0", "JP1", "JP2"};
	const char *tCh[2] = {"Pos", "Neg"};

	// histograms for resolution studies: 12 invariant mass bins with the Minv wifth 0.3 GeV/c^2
	const int nBins = 12;
	Double_t nBinsEdges[nBins + 1] = {0.26, 0.56, 0.86, 1.16, 1.46, 1.76, 2.06, 2.36, 2.66, 2.96, 3.26, 3.56, 4.0};

	const int xnBins = 13;
	Double_t xnBinsEdges[xnBins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};

	// true pair fraction
	htruePionPairJP0 = new TH1D("htruePionPairJP0", "", xnBins, xnBinsEdges);
	htruePionPairJP1 = new TH1D("htruePionPairJP1", "", xnBins, xnBinsEdges);
	htruePionPairJP2 = new TH1D("htruePionPairJP2", "", xnBins, xnBinsEdges);

	hrecPionPairJP0 = new TH1D("hrecPionPairJP0", "", xnBins, xnBinsEdges);
	hrecPionPairJP1 = new TH1D("hrecPionPairJP1", "", xnBins, xnBinsEdges);
	hrecPionPairJP2 = new TH1D("hrecPionPairJP2", "", xnBins, xnBinsEdges);

	htruePionProtonPairJP0 = new TH1D("htruePionProtonPairJP0", "", xnBins, xnBinsEdges);
	htruePionKaonPairJP0 = new TH1D("htruePionKaonPairJP0", "", xnBins, xnBinsEdges);
	htruePionElectronPairJP0 = new TH1D("htruePionElectronPairJP0", "", xnBins, xnBinsEdges);
	htrueProtonElectronPairJP0 = new TH1D("htrueProtonElectronPairJP0", "", xnBins, xnBinsEdges);
	htrueProtonPairJP0 = new TH1D("htrueProtonPairJP0", "", xnBins, xnBinsEdges);
	htrueProtonKaonPairJP0 = new TH1D("htrueProtonKaonPairJP0", "", xnBins, xnBinsEdges);
	htrueKaonPairJP0 = new TH1D("htrueKaonPairJP0", "", xnBins, xnBinsEdges);
	htrueKaonElectronPairJP0 = new TH1D("htrueKaonElectronPairJP0", "", xnBins, xnBinsEdges);
	htrueElectronPairJP0 = new TH1D("htrueElectronPairJP0", "", xnBins, xnBinsEdges);

	htruePionProtonPairJP1 = new TH1D("htruePionProtonPairJP1", "", xnBins, xnBinsEdges);
	htruePionKaonPairJP1 = new TH1D("htruePionKaonPairJP1", "", xnBins, xnBinsEdges);
	htruePionElectronPairJP1 = new TH1D("htruePionElectronPairJP1", "", xnBins, xnBinsEdges);
	htrueProtonElectronPairJP1 = new TH1D("htrueProtonElectronPairJP1", "", xnBins, xnBinsEdges);
	htrueProtonPairJP1 = new TH1D("htrueProtonPairJP1", "", xnBins, xnBinsEdges);
	htrueProtonKaonPairJP1 = new TH1D("htrueProtonKaonPairJP1", "", xnBins, xnBinsEdges);
	htrueKaonPairJP1 = new TH1D("htrueKaonPairJP1", "", xnBins, xnBinsEdges);
	htrueKaonElectronPairJP1 = new TH1D("htrueKaonElectronPairJP1", "", xnBins, xnBinsEdges);
	htrueElectronPairJP1 = new TH1D("htrueElectronPairJP1", "", xnBins, xnBinsEdges);

	htruePionProtonPairJP2 = new TH1D("htruePionProtonPairJP2", "", xnBins, xnBinsEdges);
	htruePionKaonPairJP2 = new TH1D("htruePionKaonPairJP2", "", xnBins, xnBinsEdges);
	htruePionElectronPairJP2 = new TH1D("htruePionElectronPairJP2", "", xnBins, xnBinsEdges);
	htrueProtonElectronPairJP2 = new TH1D("htrueProtonElectronPairJP2", "", xnBins, xnBinsEdges);
	htrueProtonPairJP2 = new TH1D("htrueProtonPairJP2", "", xnBins, xnBinsEdges);
	htrueProtonKaonPairJP2 = new TH1D("htrueProtonKaonPairJP2", "", xnBins, xnBinsEdges);
	htrueKaonPairJP2 = new TH1D("htrueKaonPairJP2", "", xnBins, xnBinsEdges);
	htrueKaonElectronPairJP2 = new TH1D("htrueKaonElectronPairJP2", "", xnBins, xnBinsEdges);
	htrueElectronPairJP2 = new TH1D("htrueElectronPairJP2", "", xnBins, xnBinsEdges);

	// all true opposite charge hadrons pair
	htrueHadPairJP0 = new TH1D("htrueHadPairJP0", "No nSigma cut", xnBins, xnBinsEdges);
	htrueHadPairJP1 = new TH1D("htrueHadPairJP1", "No nSigma cut", xnBins, xnBinsEdges);
	htrueHadPairJP2 = new TH1D("htrueHadPairJP2", "No nSigma cut", xnBins, xnBinsEdges);

	htruePiPairJP0 = new TH1D("htruePiPairJP0", "true π", xnBins, xnBinsEdges);
	htruePiPairJP1 = new TH1D("htruePiPairJP1", "true π", xnBins, xnBinsEdges);
	htruePiPairJP2 = new TH1D("htruePiPairJP2", "true π", xnBins, xnBinsEdges);

	htrueDetPiPairJP0 = new TH1D("htrueDetPiPairJP0", "true π and det π matched", xnBins, xnBinsEdges);
	htrueDetPiPairJP1 = new TH1D("htrueDetPiPairJP1", "true π and det π matched", xnBins, xnBinsEdges);
	htrueDetPiPairJP2 = new TH1D("htrueDetPiPairJP2", "true π and det π matched", xnBins, xnBinsEdges);

	htrueDetSysCutJP0 = new TH1D("htrueDetSysCutJP0", "true π and det π beyond |nsigmapi|>4", xnBins, xnBinsEdges);
	htrueDetSysCutJP1 = new TH1D("htrueDetSysCutJP1", "true π and det π beyond |nsigmapi|>4", xnBins, xnBinsEdges);
	htrueDetSysCutJP2 = new TH1D("htrueDetSysCutJP2", "true π and det π beyond |nsigmapi|>4", xnBins, xnBinsEdges);

	htrueMissedPiPairJP0 = new TH1D("htrueMissedPiPairJP0", "true π but misidentified by det", xnBins, xnBinsEdges);
	htrueMissedPiPairJP1 = new TH1D("htrueMissedPiPairJP1", "true π but misidentified by det", xnBins, xnBinsEdges);
	htrueMissedPiPairJP2 = new TH1D("htrueMissedPiPairJP2", "true π but misidentified by det", xnBins, xnBinsEdges);

	// pair multiplicities
	hpairMult = new TH1D("hpairMulti", "", 100, 0, 60);
	// nsigmapion in mass bins
	for (int i = 0; i < xnBins; i++)
	{
		hnsigmaPionPiPi[i] = new TH2D(Form("hnsigmaPionPiPi_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionPiP[i] = new TH2D(Form("hnsigmaPionPiP_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionPiK[i] = new TH2D(Form("hnsigmaPionPiK_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionPiE[i] = new TH2D(Form("hnsigmaPionPiE_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionPP[i] = new TH2D(Form("hnsigmaPionPP_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionPK[i] = new TH2D(Form("hnsigmaPionPK_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionPE[i] = new TH2D(Form("hnsigmaPionPE_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionKK[i] = new TH2D(Form("hnsigmaPionKK_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
		hnsigmaPionKE[i] = new TH2D(Form("hnsigmaPionKE_bin%i", i), "Pos vs Neg", 100, -10, 10, 100, -10, 10);
	}

	// hist for combinatorial background systematic
	htrueSysJP0_0 = new TH1D("htrueSysJP0_0", "", xnBins, xnBinsEdges);
	hrecSysJP0_0 = new TH1D("hrecSysJP0_0", "", xnBins, xnBinsEdges);
	hrecSysJP0_1 = new TH1D("hrecSysJP0_1", "", xnBins, xnBinsEdges);
	hrecSysJP0_2 = new TH1D("hrecSysJP0_2", "", xnBins, xnBinsEdges);
	hrecSysJP0_3 = new TH1D("hrecSysJP0_3", "", xnBins, xnBinsEdges);
	hrecSysJP0_4 = new TH1D("hrecSysJP0_4", "", xnBins, xnBinsEdges);

	htrueSysJP1_0 = new TH1D("htrueSysJP1_0", "", xnBins, xnBinsEdges);
	hrecSysJP1_0 = new TH1D("hrecSysJP1_0", "", xnBins, xnBinsEdges);
	hrecSysJP1_1 = new TH1D("hrecSysJP1_1", "", xnBins, xnBinsEdges);
	hrecSysJP1_2 = new TH1D("hrecSysJP1_2", "", xnBins, xnBinsEdges);
	hrecSysJP1_3 = new TH1D("hrecSysJP1_3", "", xnBins, xnBinsEdges);
	hrecSysJP1_4 = new TH1D("hrecSysJP1_4", "", xnBins, xnBinsEdges);

	htrueSysJP2_0 = new TH1D("htrueSysJP2_0", "", xnBins, xnBinsEdges);
	hrecSysJP2_0 = new TH1D("hrecSysJP2_0", "", xnBins, xnBinsEdges);
	hrecSysJP2_1 = new TH1D("hrecSysJP2_1", "", xnBins, xnBinsEdges);
	hrecSysJP2_2 = new TH1D("hrecSysJP2_2", "", xnBins, xnBinsEdges);
	hrecSysJP2_3 = new TH1D("hrecSysJP2_3", "", xnBins, xnBinsEdges);
	hrecSysJP2_4 = new TH1D("hrecSysJP2_4", "", xnBins, xnBinsEdges);

	// reconstructed variables
	Double_t p1x = -999, p2x = -999, p1y = -999, p2y = -999, p1z = -999, p2z = -999, psx = -999, psy = -999, psz = -999, ps = -999, cone = -999, R = -999, Rx = -999, Ry = -999, Rz = -999, R1 = -999;
	Double_t p1 = -999, p2 = -999, E1 = -999, E2 = -999, Minv = -999, pT_pair = -999, pT_min_pair = -999, eta_pair = -999, fitPts_min_pair = -999, phiDiff = -999;
	Double_t cosPhiRB = -999, sinPhiRB = -999, cosPhiRY = -999, sinPhiRY;
	Double_t PhiR_cosB = -999, PhiR_sinB = -999, PhiRB = -999, PhiR_cosY = -999, PhiR_sinY = -999, PhiRY = -999, phiRB = -999, phiDiffB = -999, phi_cos = -999, phi_sin = -999, phi_pair = -999;
	Double_t R1_p1 = -999, R2_p1 = -999, R1_p2 = -999, R2_p2 = -999, Rs_tr1 = -999, Rs_tr2 = -999, sourcePartEta = -999, sourcePartPhi = -999, Rs_pair = -999; // opening angles between parton and tracks and associated parton eta=-999, phi

	// pythia varibales
	Double_t p1x_mc = -999, p2x_mc = -999, p1y_mc = -999, p2y_mc = -999, p1z_mc = -999, p2z_mc = -999, psx_mc = -999, psy_mc = -999, psz_mc = -999, ps_mc = -999, cone_mc = -999, R_mc = -999, Rs_mc = -999, phiDiff_mc = -999, Rx_mc = -999, Ry_mc = -999, Rz_mc = -999, R1_mc = -999, Rx1_mc = -999, Ry1_mc = -999, Rz1_mc = -999;
	Double_t p1_mc = -999, p2_mc = -999, E1_mc = -999, E2_mc = -999, Minv_mc = -999, pT_pair_mc = -999, eta_pair_mc = -999, phi_pair_mc = -999;
	Double_t sourceParton_mc = -999, partonE_mc = -999, E_pair_mc = -999;

	Double_t pi = 3.14159265359;
	Double_t m_pion = 0.1396;	  // GeV
	Bool_t pairMatch = kFALSE;	  // Rs<0.5, Minv<4, and pT_pair<15, pair-source parton matching condition
	Bool_t pairMatch_mc = kFALSE; // Rs<0.5, Minv<4, and pT_pair<15, pair-source parton matching condition
	Double_t pairSep;
	Bool_t geLimit = kFALSE;
	Bool_t pyLimit = kFALSE;
	int pairCounts;

	if (fChain == 0)
		return;

	// calculate weight for each partonic pT bins
	const Int_t nptBin = 13;
	Double_t partPtRange[nptBin + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., -1};
	/*
	Double_t fudge[nptBin] = {1.228, 1.051, 1.014, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	Double_t crossSection[nptBin] = {9.0012e+00, 1.46253e+00, 3.54466e-01, 1.51622e-01, 2.49062e-02, 5.84527e-03, 2.30158e-03, 3.42755e-04, 4.57002e-05, 9.72535e-06, 4.69889e-07, 2.69202e-08, 1.43453e-09};
	Double_t numEvents[nptBin] = {3318626, 3301413, 3291662, 3280010, 3282543, 3275693, 3276437, 3276795, 3272804, 2179660, 2183230, 1091927, 1090857}; // number of events in file after processing MuDsts
	Double_t binLumi[nptBin] = {0}, binWt[nptBin] = {0};
	for (int mBin = 0; mBin < nptBin; mBin++)
	{
		binLumi[mBin] = (numEvents[mBin] * crossSection[12]) / (crossSection[mBin] * numEvents[12]); // bin luminosity
		binWt[mBin] = 1. / (binLumi[mBin] * fudge[mBin]);											 // Bin weight for partonic pT weighting.
	}
	*/
	// calculated bin weights weighted by the highest pt bin weight
	Double_t binWt[13] = {1.67958e+09, 3.20524e+08, 8.07797e+07, 3.51516e+07, 5.76973e+06, 1.35694e+06, 534174, 79541.3, 10618.3, 3392.93, 163.664, 18.7475, 1.0};

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

		Double_t vz_weight = 0.9890 + (fVZ * 0.00110841) + (fVZ * fVZ * 1.79596e-05) + (fVZ * fVZ * fVZ * 5.71566e-07);

		for (int j = 0; j < 13; j++)
		{
			if (j < 12 && partonicPtBin >= partPtRange[j] && partonicPtBin < partPtRange[j + 1])
			{
				ptw = binWt[j] * extra_weight * vz_weight;
			}
			if (j == 12 && partonicPtBin >= partPtRange[j])
			{
				ptw = binWt[j] * extra_weight * vz_weight;
			}
		}

		if (fverRank < 1e6)
			continue;
		if (fabs(fVZ) > 60.)
			continue;

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

		if (fPart1Phi > pi)
			fPart1Phi = fPart1Phi - 2. * pi;
		if (fPart2Phi > pi)
			fPart2Phi = fPart2Phi - 2. * pi;
		pairCounts = 0;

		if (!isJP)
			continue;

		// cout << "tr1 nSigma" << setw(ww) << "tr2 nSigma" << setw(ww) << "case 0" << setw(ww) << "case 1" << setw(ww) << "case 2" << setw(ww) << "case 3" << endl;
		for (int tr = 0; tr < fmaxpar; tr++)
		{

			Bool_t ismc_tr1 = kFALSE;
			Bool_t isrec_tr1 = kFALSE;
			Bool_t ishad_tr1 = kFALSE;

			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];

			if (idTruth[tr] > 0 && fpT_mc[tr] > 0.5 && fpT_mc[tr] < 15 && fabs(feta_mc[tr]) < 1. && (fId_mc[tr] == idTruth[tr])) // 2,3 = e+-, 5,6 = mu+-, 8,9=pi+-, 11,12= K+-, 14,15=ppbar
				ismc_tr1 = kTRUE;

			if (!ismc_tr1)
				continue;

			if (fpT[tr] > 0.5 && fpT[tr] < 15 && fabs(feta[tr]) < 1 && ffitPts[tr] > 15 && fitPtsRatio_tr > 0.51 && dcaCut(fpT[tr], fdca[tr]))
				isrec_tr1 = true;

			if (!isrec_tr1)
				continue;

			ishad_tr1 = (fpId_mc[tr] == 2 || fpId_mc[tr] == 3 || fpId_mc[tr] == 8 || fpId_mc[tr] == 9 || fpId_mc[tr] == 11 || fpId_mc[tr] == 12 || fpId_mc[tr] == 14 || fpId_mc[tr] == 15);

			if (!ishad_tr1)
				continue; // select all hadrons/electrons

			//    Track 2 loop
			for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
			{
				Bool_t isrec_tr2 = kFALSE;
				Bool_t ismc_tr2 = kFALSE;
				Bool_t ishad_tr2 = kFALSE;
				Bool_t isOppPairTrue = kFALSE;
				Bool_t isOppPairDet = kFALSE;

				double nsigmapi_pos = -999, nsigmapi_neg = -999;

				Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];

				if (idTruth[tr2] > 0 && fpT_mc[tr2] > 0.5 && fpT_mc[tr2] < 15 && fabs(feta_mc[tr2]) < 1. && (fId_mc[tr2] == idTruth[tr2])) // 2,3 = e+-, 5,6 = mu+-, 8,9=pi+-, 11,12= K+-,
					ismc_tr2 = kTRUE;

				if (!ismc_tr2)
					continue;

				if (fpT[tr2] > 0.5 && fpT[tr2] < 15 && abs(feta[tr2]) < 1 && ffitPts[tr2] > 15 && fitPtsRatio_tr2 > 0.51 && dcaCut(fpT[tr2], fdca[tr2]))
					isrec_tr2 = kTRUE;

				if (!isrec_tr2)
					continue;

				ishad_tr2 = (fpId_mc[tr2] == 2 || fpId_mc[tr2] == 3 || fpId_mc[tr2] == 8 || fpId_mc[tr2] == 9 || fpId_mc[tr2] == 11 || fpId_mc[tr2] == 12 || fpId_mc[tr2] == 14 || fpId_mc[tr2] == 15);

				if (!ishad_tr2)
					continue; // select all hadrons/electrons

				phiDiff = fphi_mc[tr] - fphi_mc[tr2];
				if (phiDiff > pi)
					phiDiff -= (2 * pi);
				if (phiDiff < ((-1) * pi))
					phiDiff += (2 * pi);

				bool is_pos_tr1 = fpId_mc[tr] == 2 || fpId_mc[tr] == 8 || fpId_mc[tr] == 11 || fpId_mc[tr] == 14;
				bool is_neg_tr1 = fpId_mc[tr] == 3 || fpId_mc[tr] == 9 || fpId_mc[tr] == 12 || fpId_mc[tr] == 15;
				bool is_pos_tr2 = fpId_mc[tr2] == 2 || fpId_mc[tr2] == 8 || fpId_mc[tr2] == 11 || fpId_mc[tr2] == 14;
				bool is_neg_tr2 = fpId_mc[tr2] == 3 || fpId_mc[tr2] == 9 || fpId_mc[tr2] == 12 || fpId_mc[tr2] == 15;

				if ((is_pos_tr1 && is_neg_tr2) || (is_neg_tr1 && is_pos_tr2))
					isOppPairTrue = kTRUE; // reject if same charge pair

				if (fcharge[tr] != fcharge[tr2])
					isOppPairDet = kTRUE;

				phiDiff_mc = fphi_mc[tr] - fphi_mc[tr2];
				if (phiDiff_mc > pi)
					phiDiff_mc -= (2 * pi);
				if (phiDiff_mc < ((-1) * pi))
					phiDiff_mc += (2 * pi);

				cone_mc = sqrt(pow(feta_mc[tr] - feta_mc[tr2], 2) + pow(phiDiff_mc, 2));

				if (cone_mc < 0.02 || cone_mc > 0.7)
					continue;

				bool is_pp_pair = (fpId_mc[tr] == 14 && fpId_mc[tr2] == 15) || (fpId_mc[tr] == 15 && fpId_mc[tr2] == 14);
				bool is_pipi_pair = (fpId_mc[tr] == 8 && fpId_mc[tr2] == 9) || (fpId_mc[tr] == 9 && fpId_mc[tr2] == 8);

				bool is_pip_pair = (fpId_mc[tr] == 8 && fpId_mc[tr2] == 15) || (fpId_mc[tr] == 9 && fpId_mc[tr2] == 14) || (fpId_mc[tr] == 14 && fpId_mc[tr2] == 9) || (fpId_mc[tr] == 15 && fpId_mc[tr2] == 8);

				bool is_pie_pair = (fpId_mc[tr] == 8 && fpId_mc[tr2] == 3) || (fpId_mc[tr] == 9 && fpId_mc[tr2] == 2) || (fpId_mc[tr] == 2 && fpId_mc[tr2] == 9) || (fpId_mc[tr] == 3 && fpId_mc[tr2] == 8);

				bool is_pik_pair = (fpId_mc[tr] == 8 && fpId_mc[tr2] == 12) || (fpId_mc[tr] == 9 && fpId_mc[tr2] == 11) || (fpId_mc[tr] == 11 && fpId_mc[tr2] == 9) || (fpId_mc[tr] == 12 && fpId_mc[tr2] == 8);

				bool is_kk_pair = (fpId_mc[tr] == 11 && fpId_mc[tr2] == 12) || (fpId_mc[tr] == 12 && fpId_mc[tr2] == 11);

				bool is_pk_pair = (fpId_mc[tr] == 14 && fpId_mc[tr2] == 12) || (fpId_mc[tr] == 15 && fpId_mc[tr2] == 11) || (fpId_mc[tr] == 11 && fpId_mc[tr2] == 15) || (fpId_mc[tr] == 12 && fpId_mc[tr2] == 14);

				bool is_pe_pair = (fpId_mc[tr] == 14 && fpId_mc[tr2] == 3) || (fpId_mc[tr] == 15 && fpId_mc[tr2] == 2) || (fpId_mc[tr] == 3 && fpId_mc[tr2] == 14) || (fpId_mc[tr] == 2 && fpId_mc[tr2] == 15);

				bool is_ke_pair = (fpId_mc[tr] == 11 && fpId_mc[tr2] == 3) || (fpId_mc[tr] == 12 && fpId_mc[tr2] == 2) || (fpId_mc[tr] == 3 && fpId_mc[tr2] == 11) || (fpId_mc[tr] == 2 && fpId_mc[tr2] == 12);

				bool is_ee_pair = (fpId_mc[tr] == 2 && fpId_mc[tr2] == 3) || (fpId_mc[tr] == 3 && fpId_mc[tr2] == 2);

				// if (!is_pipi_pair)
				//	continue;

				if (is_pos_tr1)
				{

					p1x_mc = fpT_mc[tr] * cos(fphi_mc[tr]);
					p2x_mc = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p1y_mc = fpT_mc[tr] * sin(fphi_mc[tr]);
					p2y_mc = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p1z_mc = fpT_mc[tr] * sinh(feta_mc[tr]);
					p2z_mc = fpT_mc[tr2] * sinh(feta_mc[tr2]);

					p1_mc = sqrt(p1x_mc * p1x_mc + p1y_mc * p1y_mc + p1z_mc * p1z_mc);
					p2_mc = sqrt(p2x_mc * p2x_mc + p2y_mc * p2y_mc + p2z_mc * p2z_mc);

					nsigmapi_pos = fnSigmaPion[tr];
					nsigmapi_neg = fnSigmaPion[tr2];
				}

				if (is_neg_tr1)
				{
					p1x_mc = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p2x_mc = fpT_mc[tr] * cos(fphi_mc[tr]);
					p1y_mc = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p2y_mc = fpT_mc[tr] * sin(fphi_mc[tr]);
					p1z_mc = fpT_mc[tr2] * sinh(feta_mc[tr2]);
					p2z_mc = fpT_mc[tr] * sinh(feta_mc[tr]);

					p1_mc = sqrt(p1x_mc * p1x_mc + p1y_mc * p1y_mc + p1z_mc * p1z_mc);
					p2_mc = sqrt(p2x_mc * p2x_mc + p2y_mc * p2y_mc + p2z_mc * p2z_mc);

					nsigmapi_pos = fnSigmaPion[tr2];
					nsigmapi_neg = fnSigmaPion[tr];
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

				if (Minv_mc > 1.5 && partonicPtBin < 5.) // filter to "bad" event
					continue;

				// for combinatorial background correction
				// check if the rec pair is from the true pions
				bool true_pair = (fpId_mc[tr] == 8 && fpId_mc[tr2] == 9) || (fpId_mc[tr] == 9 && fpId_mc[tr2] == 8);
				bool det_pair = (fnSigmaPion[tr] > -1 && fnSigmaPion[tr] < 2) && (fnSigmaPion[tr2] > -1 && fnSigmaPion[tr2] < 2);

				// for combinatorial background correction systematic
				bool sys_cut_default = fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2 && fnSigmaPion[tr2] > -1.0 && fnSigmaPion[tr2] < 2;
				bool sys_cut0 = fnSigmaPion[tr] > -1.1 && fnSigmaPion[tr] < 2.1 && fnSigmaPion[tr2] > -1.1 && fnSigmaPion[tr2] < 2.1;
				bool sys_cut1 = fnSigmaPion[tr] > -1.2 && fnSigmaPion[tr] < 2.2 && fnSigmaPion[tr2] > -1.2 && fnSigmaPion[tr2] < 2.2;
				bool sys_cut2 = fnSigmaPion[tr] > -1.3 && fnSigmaPion[tr] < 2.3 && fnSigmaPion[tr2] > -1.3 && fnSigmaPion[tr2] < 2.3;
				bool sys_cut3 = fnSigmaPion[tr] > -1.4 && fnSigmaPion[tr] < 2.4 && fnSigmaPion[tr2] > -1.4 && fnSigmaPion[tr2] < 2.4;
				bool sys_cut4 = fnSigmaPion[tr] > -1.5 && fnSigmaPion[tr] < 2.5 && fnSigmaPion[tr2] > -1.5 && fnSigmaPion[tr2] < 2.5;
				bool sys_cut = fabs(fnSigmaPion[tr]) > 2.0 || fabs(fnSigmaPion[tr2]) > 2.0;

				// cout << fnSigmaPion[tr] << setw(ww) << fnSigmaPion[tr2] << setw(ww) << sys_cut0 << setw(ww) << sys_cut1 << setw(ww) << sys_cut2 << setw(ww) << sys_cut3 << endl;
				//  all true hadron pairs

				if (isJP0)
					htrueHadPairJP0->Fill(Minv_mc, ptw);
				if (isJP1)
					htrueHadPairJP1->Fill(Minv_mc, ptw);
				if (isJP2)
					htrueHadPairJP2->Fill(Minv_mc, ptw);

				for (int i = 0; i < 13; i++)
				{
					if (Minv_mc >= xnBinsEdges[i] && Minv_mc < xnBinsEdges[i + 1])
					{
						if (is_pipi_pair)
							hnsigmaPionPiPi[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_pip_pair)
							hnsigmaPionPiP[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_pik_pair)
							hnsigmaPionPiK[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_pie_pair)
							hnsigmaPionPiE[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_pp_pair)
							hnsigmaPionPP[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_pk_pair)
							hnsigmaPionPK[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_pe_pair)
							hnsigmaPionPE[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_kk_pair)
							hnsigmaPionKK[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
						if (is_ke_pair)
							hnsigmaPionKE[i]->Fill(nsigmapi_neg, nsigmapi_pos, ptw);
					}
				}

				if (!is_pipi_pair)
					continue;

				// if (abs(fnSigmaPion[tr2]) > 4 || abs(fnSigmaPion[tr]) > 4)
				//	cout << "Found one---" << endl;

				// now only pipi pair
				if (isJP0)
				{
					htruePiPairJP0->Fill(Minv_mc, ptw);
					if (sys_cut)
						htrueDetSysCutJP0->Fill(Minv_mc, ptw);
					if (sys_cut_default)
						htrueDetPiPairJP0->Fill(Minv_mc, ptw);
					else
						htrueMissedPiPairJP0->Fill(Minv_mc, ptw);
				}
				if (isJP1)
				{
					htruePiPairJP1->Fill(Minv_mc, ptw);
					if (sys_cut)
						htrueDetSysCutJP1->Fill(Minv_mc, ptw);
					if (sys_cut_default)
						htrueDetPiPairJP1->Fill(Minv_mc, ptw);
					else
						htrueMissedPiPairJP1->Fill(Minv_mc, ptw);
				}
				if (isJP2)
				{
					htruePiPairJP2->Fill(Minv_mc, ptw);
					if (sys_cut)
						htrueDetSysCutJP2->Fill(Minv_mc, ptw);
					if (sys_cut_default)
						htrueDetPiPairJP2->Fill(Minv_mc, ptw);
					else
						htrueMissedPiPairJP2->Fill(Minv_mc, ptw);
				}

				// for systematic studies
				if (trigJP0 == 1)
				{
					htrueSysJP0_0->Fill(Minv_mc, ptw);
					if (sys_cut0)
						hrecSysJP0_0->Fill(Minv_mc, ptw);
					if (sys_cut1)
						hrecSysJP0_1->Fill(Minv_mc, ptw);
					if (sys_cut2)
						hrecSysJP0_2->Fill(Minv_mc, ptw);
					if (sys_cut3)
						hrecSysJP0_3->Fill(Minv_mc, ptw);
					if (sys_cut4)
						hrecSysJP0_4->Fill(Minv_mc, ptw);
				}
				if (trigJP1 == 1)
				{
					htrueSysJP1_0->Fill(Minv_mc, ptw);
					if (sys_cut0)
						hrecSysJP1_0->Fill(Minv_mc, ptw);
					if (sys_cut1)
						hrecSysJP1_1->Fill(Minv_mc, ptw);
					if (sys_cut2)
						hrecSysJP1_2->Fill(Minv_mc, ptw);
					if (sys_cut3)
						hrecSysJP1_3->Fill(Minv_mc, ptw);
					if (sys_cut4)
						hrecSysJP1_4->Fill(Minv_mc, ptw);
				}
				if (trigJP2 == 1)
				{
					htrueSysJP2_0->Fill(Minv_mc, ptw);
					if (sys_cut0)
						hrecSysJP2_0->Fill(Minv_mc, ptw);
					if (sys_cut1)
						hrecSysJP2_1->Fill(Minv_mc, ptw);
					if (sys_cut2)
						hrecSysJP2_2->Fill(Minv_mc, ptw);
					if (sys_cut3)
						hrecSysJP2_3->Fill(Minv_mc, ptw);
					if (sys_cut4)
						hrecSysJP2_4->Fill(Minv_mc, ptw);
				}

				pairCounts++;
			} // track 2 loop
		}	  // track loop
		hpairMult->Fill(pairCounts);
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

	htruePiPairJP0->Write();
	htrueDetPiPairJP0->Write();
	htrueMissedPiPairJP0->Write();

	htruePiPairJP1->Write();
	htrueDetPiPairJP1->Write();
	htrueMissedPiPairJP1->Write();

	htruePiPairJP2->Write();
	htrueDetPiPairJP2->Write();
	htrueMissedPiPairJP2->Write();

	htrueSysJP0_0->Write();
	hrecSysJP0_0->Write();
	hrecSysJP0_1->Write();
	hrecSysJP0_2->Write();
	hrecSysJP0_3->Write();
	hrecSysJP0_4->Write();

	htrueDetSysCutJP0->Write();
	htrueDetSysCutJP1->Write();
	htrueDetSysCutJP2->Write();

	htrueSysJP1_0->Write();
	hrecSysJP1_0->Write();
	hrecSysJP1_1->Write();
	hrecSysJP1_2->Write();
	hrecSysJP1_3->Write();
	hrecSysJP1_4->Write();

	htrueSysJP2_0->Write();
	hrecSysJP2_0->Write();
	hrecSysJP2_1->Write();
	hrecSysJP2_2->Write();
	hrecSysJP2_3->Write();
	hrecSysJP2_4->Write();

	for (int i = 0; i < 13; i++)
	{
		hnsigmaPionPiPi[i]->Write();
		hnsigmaPionPiP[i]->Write();
		hnsigmaPionPiK[i]->Write();
		hnsigmaPionPiE[i]->Write();
		hnsigmaPionPP[i]->Write();
		hnsigmaPionPK[i]->Write();
		hnsigmaPionPE[i]->Write();
		hnsigmaPionKK[i]->Write();
		hnsigmaPionKE[i]->Write();
	}

	hpairMult->Write();

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
