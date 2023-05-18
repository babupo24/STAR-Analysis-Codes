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

	const int nBin = 5;
	const int nBinpTi = 9;
	// eficiency dependence on mass bins
	Double_t mRange[6] = {0.27, 0.45, 0.75, 1.35, 2.20, 4.0};
	// eficiency dependence on pt bins as a function of eta
	Double_t ptRange[6] = {0.2, 0.5, 1.0, 2.0, 4.0, 15.0};
	Double_t ptRangei[nBinpTi + 1] = {0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 3.0, 5.0, 8.0, 15.0};
	// eficiency dependence on Cone bins
	Double_t coneRange[6] = {0., 0.1, 0.2, 0.35, 0.5, 0.7};

	hpairGenJP0C = new TH1D("hGenJP0C", "", xnBins, xnBinsEdges);
	hpairGenRecJP0C = new TH1D("hGenRecJP0C", "", xnBins, xnBinsEdges);
	hpairGenJP1C = new TH1D("hGenJP1C", "", xnBins, xnBinsEdges);
	hpairGenRecJP1C = new TH1D("hGenRecJP1C", "", xnBins, xnBinsEdges);
	hpairGenJP2C = new TH1D("hGenJP2C", "", xnBins, xnBinsEdges);
	hpairGenRecJP2C = new TH1D("hGenRecJP2C", "", xnBins, xnBinsEdges);

	// in each of the mass bins
	for (int ntrg = 0; ntrg < nTrg; ntrg++)
	{
		for (int nch = 0; nch < nCh; nch++)
		{
			for (int nbin = 0; nbin < xnBins; nbin++)
			{
				// x-section level cuts
				hGenRecEta_Mx[ntrg][nch][nbin] = new TH1D(Form("hGenRecEta_Mx_%s%s%i", tTrig[ntrg], tCh[nch], nbin), "numerator", 60, -1.5, 1.5);
				hGenEta_Mx[ntrg][nch][nbin] = new TH1D(Form("hGenEta_Mx_%s%s%i", tTrig[ntrg], tCh[nch], nbin), "denominator", 60, -1.5, 1.5);
				hGenRecpT_Mx[ntrg][nch][nbin] = new TH1D(Form("hGenRecpT_Mx_%s%s%i", tTrig[ntrg], tCh[nch], nbin), "numerator", 60, 0.2, 15.5);
				hGenpT_Mx[ntrg][nch][nbin] = new TH1D(Form("hGenpT_Mx_%s%s%i", tTrig[ntrg], tCh[nch], nbin), "denominator", 60, 0.2, 15.5);

				// histograms for efficiency includeing acceptance effect
				hGenRecpT[ntrg][nch][nbin] = new TH1D(Form("hGenRecpT_%s%s%i", tTrig[ntrg], tCh[nch], nbin), " 4 accept. effect ", 60, 0.2, 15.5);
				hGenpT[ntrg][nch][nbin] = new TH1D(Form("hGenpT_%s%s%i", tTrig[ntrg], tCh[nch], nbin), "fpr accept. effect", 60, 0.2, 15.5);
			}
		}
	}

	for (int nbin = 0; nbin < xnBins; nbin++)
	{
		hpairGenJP0X[nbin] = new TH1D(Form("hpairGenJP0X_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
		hpairRecJP0X[nbin] = new TH1D(Form("hpairRecJP0X_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);

		hpairGenJP1X[nbin] = new TH1D(Form("hpairGenJP1X_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
		hpairRecJP1X[nbin] = new TH1D(Form("hpairRecJP1X_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);

		hpairGenJP2X[nbin] = new TH1D(Form("hpairGenJP2X_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
		hpairRecJP2X[nbin] = new TH1D(Form("hpairRecJP2X_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);

		hpairGenJP0[nbin] = new TH1D(Form("hpairGenJP0_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
		hpairRecJP0[nbin] = new TH1D(Form("hpairRecJP0_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);

		hpairGenJP1[nbin] = new TH1D(Form("hpairGenJP1_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
		hpairRecJP1[nbin] = new TH1D(Form("hpairRecJP1_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);

		hpairGenJP2[nbin] = new TH1D(Form("hpairGenJP2_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
		hpairRecJP2[nbin] = new TH1D(Form("hpairRecJP2_%i", nbin), "", 20, xnBinsEdges[nbin], xnBinsEdges[nbin + 1]);
	}
	// test histograms

	// reconstructed variables
	Double_t p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1;
	Double_t p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
	Double_t phi_cos, phi_sin, phi_pair;

	// monte carlo varibales
	Double_t p1x_mc, p2x_mc, p1y_mc, p2y_mc, p1z_mc, p2z_mc, psx_mc, psy_mc, psz_mc, ps_mc, cone_mc, R_mc, Rs_mc, Rx_mc, Ry_mc, Rz_mc, R1_mc, Rx1_mc, Ry1_mc, Rz1_mc;
	Double_t p1_mc, p2_mc, E1_mc, E2_mc, Minv_mc, pT_pair_mc, eta_pair_mc, phi_pair_mc;

	double pt_pos = -999;
	double pt_neg = -999;
	int pion_type1 = 0;
	int pion_type2 = 0;

	Double_t pi = 3.14159265359;
	Double_t m_pion = 0.1396; // GeV

	if (fChain == 0)
		return;

	// calculate weight for each partonic pT bins
	const Int_t nptBin = 13;
	Double_t partPtBin[nptBin + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., 80}; // for trigger eff
	Double_t partPtRange[nptBin + 1] = {2., 3., 4., 5., 7., 9., 11., 15., 20., 25., 35., 45., 55., -1};
	Double_t fudge[nptBin] = {1.228, 1.051, 1.014, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	/*
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

		Bool_t vertexCut = kFALSE;
		Bool_t isJP = kFALSE;
		Bool_t isJP0 = kFALSE;
		Bool_t isJP1 = kFALSE;
		Bool_t isJP2 = kFALSE;

		if (fverRank > 1e6 && fabs(fVZ) < 60)
			vertexCut = kTRUE; // common event vertex selection cut
		if (!vertexCut)
			continue;

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

			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];

			if (idTruth[tr] <= 0 || fId_mc[tr] != idTruth[tr])
				continue;

			Bool_t tr1_gen = kFALSE;
			Bool_t tr1_genX = kFALSE;

			if (fabs(feta_mc[tr]) < 2.5 && fpT_mc[tr] > 0.2 && (fpId_mc[tr] == 8 || fpId_mc[tr] == 9))
				tr1_gen = kTRUE;
			if (!tr1_gen)
				continue;
			if (fabs(feta_mc[tr]) < 1. && fpT_mc[tr] > 0.5 && fpT_mc[tr] < 15.0 && (fpId_mc[tr] == 8 || fpId_mc[tr] == 9))
				tr1_genX = kTRUE;

			//    Track 2 loop
			for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
			{
				if (idTruth[tr2] <= 0 || fId_mc[tr2] == -999)
					continue;

				Bool_t tr2_gen = kFALSE;
				Bool_t tr2_genX = kFALSE;
				Bool_t pipos_gen_rec = kFALSE;
				Bool_t pipos_gen_recX = kFALSE;
				Bool_t pineg_gen_rec = kFALSE;
				Bool_t pineg_gen_recX = kFALSE;

				Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];
				// generated tracks flag.
				if (fabs(feta_mc[tr2]) < 2.5 && fpT_mc[tr2] > 0.2 && (fpId_mc[tr2] == 8 || fpId_mc[tr2] == 9))
					tr2_gen = kTRUE;
				if (!tr2_gen)
					continue;

				if (fabs(feta_mc[tr2]) < 1. && fpT_mc[tr2] > 0.5 && fpT_mc[tr2] < 15.0 && (fpId_mc[tr2] == 8 || fpId_mc[tr2] == 9))
					tr2_genX = kTRUE;

				if (fpId_mc[tr] == fpId_mc[tr2])
					continue;

				Double_t phiDiff = fphi_mc[tr] - fphi_mc[tr2];
				if (phiDiff > pi)
					phiDiff -= (2 * pi);
				if (phiDiff < ((-1) * pi))
					phiDiff += (2 * pi);
				if (phiDiff == 0 || phiDiff > pi || phiDiff < -pi)
					continue;

				cone = sqrt(pow(feta_mc[tr] - feta_mc[tr2], 2) + pow(phiDiff, 2));

				if (cone >= Cone_Max7)
					continue;

				if (fpId_mc[tr] == 8) // pi+
				{
					// reco pair
					p1x = fpT_mc[tr] * cos(fphi_mc[tr]);
					p2x = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p1y = fpT_mc[tr] * sin(fphi_mc[tr]);
					p2y = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p1z = fpT_mc[tr] * sinh(feta_mc[tr]);
					p2z = fpT_mc[tr2] * sinh(feta[tr2]);
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
					// mc pair
					p1x_mc = fpT_mc[tr] * cos(fphi_mc[tr]);
					p2x_mc = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p1y_mc = fpT_mc[tr] * sin(fphi_mc[tr]);
					p2y_mc = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p1z_mc = fpT_mc[tr] * sinh(feta_mc[tr]);
					p2z_mc = fpT_mc[tr2] * sinh(feta_mc[tr2]);
					p1_mc = sqrt(p1x_mc * p1x_mc + p1y_mc * p1y_mc + p1z_mc * p1z_mc);
					p2_mc = sqrt(p2x_mc * p2x_mc + p2y_mc * p2y_mc + p2z_mc * p2z_mc);

					if (tr1_gen && fpT[tr] > 0.5 && fpT[tr] < 15.0 && fabs(feta[tr]) < 1. && ffitPts[tr] > 15 && dcaCut(fpT[tr], fdca[tr]) && fitPtsRatio_tr > 0.51)
						pipos_gen_rec = true;

					if (tr1_genX && fpT[tr] > 0.5 && fpT[tr] < 15.0 && fabs(feta[tr]) < 1. && ffitPts[tr] > 15 && dcaCut(fpT[tr], fdca[tr]) && fitPtsRatio_tr > 0.51)
						pipos_gen_recX = kTRUE;

					if (tr2_gen && fpT[tr2] > 0.5 && fpT[tr2] < 15.0 && fabs(feta[tr2]) < 1. && ffitPts[tr2] > 15 && dcaCut(fpT[tr2], fdca[tr2]) && fitPtsRatio_tr2 > 0.51)
						pineg_gen_rec = kTRUE;

					if (tr2_genX && fpT[tr2] > 0.5 && fpT[tr2] < 15.0 && fabs(feta[tr2]) < 1. && ffitPts[tr2] > 15 && dcaCut(fpT[tr2], fdca[tr2]) && fitPtsRatio_tr2 > 0.51)
						pineg_gen_recX = kTRUE;

					pt_pos = fpT_mc[tr];
					pt_neg = fpT_mc[tr2];

					// cout << "Track 1 +ve, " << p1x << ", " << p1y << " , " << p1z << endl;
				}
				if (fpId_mc[tr] == 9) // pi-
				{
					// reco pair
					p1x = fpT[tr2] * cos(fphi[tr2]);
					p2x = fpT[tr] * cos(fphi[tr]);
					p1y = fpT[tr2] * sin(fphi[tr2]);
					p2y = fpT[tr] * sin(fphi[tr]);
					p1z = fpT[tr2] * sinh(feta[tr2]);
					p2z = fpT[tr] * sinh(feta[tr]);
					p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
					p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
					// mc pair
					p1x_mc = fpT_mc[tr2] * cos(fphi_mc[tr2]);
					p2x_mc = fpT_mc[tr] * cos(fphi_mc[tr]);
					p1y_mc = fpT_mc[tr2] * sin(fphi_mc[tr2]);
					p2y_mc = fpT_mc[tr] * sin(fphi_mc[tr]);
					p1z_mc = fpT_mc[tr2] * sinh(feta_mc[tr2]);
					p2z_mc = fpT_mc[tr] * sinh(feta_mc[tr]);
					p1_mc = sqrt(p1x_mc * p1x_mc + p1y_mc * p1y_mc + p1z_mc * p1z_mc);
					p2_mc = sqrt(p2x_mc * p2x_mc + p2y_mc * p2y_mc + p2z_mc * p2z_mc);

					if (tr1_gen && fpT[tr] > 0.5 && fpT[tr] < 15.0 && fabs(feta[tr]) < 1. && ffitPts[tr] > 15 && dcaCut(fpT[tr], fdca[tr]) && fitPtsRatio_tr > 0.51)
						pineg_gen_rec = true;

					if (tr1_genX && fpT[tr] > 0.5 && fpT[tr] < 15.0 && fabs(feta[tr]) < 1. && ffitPts[tr] > 15 && dcaCut(fpT[tr], fdca[tr]) && fitPtsRatio_tr > 0.51)
						pineg_gen_recX = true;

					if (tr2_gen && fpT[tr2] > 0.5 && fpT[tr2] < 15.0 && fabs(feta[tr2]) < 1. && ffitPts[tr2] > 15 && dcaCut(fpT[tr2], fdca[tr2]) && fitPtsRatio_tr2 > 0.51)
						pipos_gen_rec = true;

					if (tr2_genX && fpT[tr2] > 0.5 && fpT[tr2] < 15.0 && fabs(feta[tr2]) < 1. && ffitPts[tr2] > 15 && dcaCut(fpT[tr2], fdca[tr2]) && fitPtsRatio_tr2 > 0.51)
						pipos_gen_recX = true;

					pt_pos = fpT_mc[tr2];
					pt_neg = fpT_mc[tr];
					// cout << "Track 1 -ve, " << p1x << ", " << p1y << " , " << p1z << endl;
				}
				// Components of sum vector
				// reco
				psx = p1x + p2x;
				psy = p1y + p2y;
				psz = p1z + p2z;
				// mc
				psx_mc = p1x_mc + p2x_mc;
				psy_mc = p1y_mc + p2y_mc;
				psz_mc = p1z_mc + p2z_mc;
				// Sum vector
				// reco
				ps = sqrt(psx * psx + psy * psy + psz * psz);
				// mc
				ps_mc = sqrt(psx_mc * psx_mc + psy_mc * psy_mc + psz_mc * psz_mc);
				// Relative momentum
				// reco
				Rx = (p1x - p2x);
				Ry = (p1y - p2y);
				Rz = (p1z - p2z);
				// mc
				Rx_mc = (p1x_mc - p2x_mc);
				Ry_mc = (p1y_mc - p2y_mc);
				Rz_mc = (p1z_mc - p2z_mc);

				// R-mag
				// rec
				R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
				// rec
				R_mc = sqrt(Rx_mc * Rx_mc + Ry_mc * Ry_mc + Rz_mc * Rz_mc);

				// calculate M,pt,eta
				// pion energy
				// reco
				E1 = sqrt(m_pion * m_pion + p1 * p1);
				E2 = sqrt(m_pion * m_pion + p2 * p2);
				// mc
				E1_mc = sqrt(m_pion * m_pion + p1_mc * p1_mc);
				E2_mc = sqrt(m_pion * m_pion + p2_mc * p2_mc);
				// invariant mass
				// reco
				Minv = sqrt(2 * m_pion * m_pion + 2 * (E1 * E2 - p1x * p2x - p1y * p2y - p1z * p2z));
				// mc
				Minv_mc = sqrt(2 * m_pion * m_pion + 2 * (E1_mc * E2_mc - p1x_mc * p2x_mc - p1y_mc * p2y_mc - p1z_mc * p2z_mc));
				// cout << "tr1 = " << fpT_mc[tr] << " tr2 = " << fpT_mc[tr2] << " Mass Gen = " << Minv_mc << " rec = " << Minv << endl;

				// p_T of pair
				// reco
				pT_pair = sqrt((p1x + p2x) * (p1x + p2x) + (p1y + p2y) * (p1y + p2y));
				// mc
				pT_pair_mc = sqrt((p1x_mc + p2x_mc) * (p1x_mc + p2x_mc) + (p1y_mc + p2y_mc) * (p1y_mc + p2y_mc));
				// eta of pair
				// reco
				eta_pair = TMath::ASinH((p1z + p2z) / pT_pair);
				// mc
				eta_pair_mc = TMath::ASinH((p1z_mc + p2z_mc) / pT_pair_mc);
				// phi_cos = acos((fpT[tr] * cos(fphi[tr]) + fpT[tr2] * cos(fphi[tr2])) / pT_pair);
				// phi_sin = asin((fpT[tr] * sin(fphi[tr]) + fpT[tr2] * sin(fphi[tr2])) / pT_pair);
				phi_cos = acos((p1x + p2x) / pT_pair);
				phi_sin = asin((p1y + p2y) / pT_pair);
				// cout << "Cosine: " << (Double_t)acos((p1x + p2x) / pT_pair) << "  " << phi_cos << "\t Sin  ";
				// cout << (Double_t)asin((p1y + p2y) / pT_pair) << "  " << phi_sin << endl;

				if (Minv_mc < Minv_Min || Minv_mc > Minv_Max || fabs(eta_pair_mc) > Eta_Cut || pT_pair_mc > ptPair_Max || pT_pair_mc < ptPair_Min)
					continue;

				if (Minv_mc > 1.50 && partonicPtBin < 5.)
					continue;

				// cout << pipos_gen_rec << "  " << pineg_gen_rec << " " << fpId_mc[tr] << "  " << fpId_mc[tr2] << "  " << fpT_mc[tr] << "  " << fpT_mc[tr2] << "  " << pt_pos << "  " << pt_neg << endl;

				if (tr1_gen && tr2_gen) // select generated pair at wide fiducial volume
				{
					if (isJP0)
					{
						// in x-sec bins
						for (int i = 0; i < xnBins; i++)
						{
							if (Minv_mc >= xnBinsEdges[i] && Minv_mc < xnBinsEdges[i + 1])
							{
								// pair level efficieicy with min cone cut
								if (cone > Cone_Min)
								{
									hpairGenJP0C->Fill(Minv_mc, ptw);
									if (pipos_gen_rec && pineg_gen_rec)
										hpairGenRecJP0C->Fill(Minv_mc, ptw);
								}
								// pair level efficieicy
								hpairGenJP0[i]->Fill(Minv_mc, ptw);
								if (pipos_gen_rec && pineg_gen_rec)
									hpairRecJP0[i]->Fill(Minv_mc, ptw);

								// track level efficiency
								//  fill positive pion eta
								hGenpT[0][0][i]->Fill(pt_pos, ptw);
								// fill negative pion eta
								hGenpT[0][1][i]->Fill(pt_neg, ptw);

								// fill generated and reconstructed
								if (pipos_gen_rec)
									hGenRecpT[0][0][i]->Fill(pt_pos, ptw);
								if (pineg_gen_rec)
									hGenRecpT[0][1][i]->Fill(pt_neg, ptw);

								if (tr1_genX && tr2_genX) // selects generated pair with same fiducial cut as rec.
								{
									// pair level efficieicy
									hpairGenJP0X[i]->Fill(Minv_mc, ptw);
									if (pipos_gen_recX && pineg_gen_recX)
										hpairRecJP0X[i]->Fill(Minv_mc, ptw);

									// fill positive pion eta
									hGenpT_Mx[0][0][i]->Fill(pt_pos, ptw);
									// fill negative pion eta
									hGenpT_Mx[0][1][i]->Fill(pt_neg, ptw);
									// fill generated and reconstructed
									if (pipos_gen_recX)
										hGenRecpT_Mx[0][0][i]->Fill(pt_pos, ptw);
									if (pineg_gen_recX)
										hGenRecpT_Mx[0][1][i]->Fill(pt_neg, ptw);
								}
							}
						}
					} // jp0 selection ends
					if (isJP1)
					{
						// in x-sec bins
						for (int i = 0; i < xnBins; i++)
						{
							if (Minv_mc >= xnBinsEdges[i] && Minv_mc < xnBinsEdges[i + 1])
							{

								if (cone > Cone_Min)
								{
									hpairGenJP1C->Fill(Minv_mc, ptw);
									if (pipos_gen_rec && pineg_gen_rec)
										hpairGenRecJP1C->Fill(Minv_mc, ptw);
								}
								// pair level efficieicy
								hpairGenJP1[i]->Fill(Minv_mc, ptw);
								if (pipos_gen_rec && pineg_gen_rec)
									hpairRecJP1[i]->Fill(Minv_mc, ptw);

								// track level efficiency
								//  fill positive pion eta
								hGenpT[1][0][i]->Fill(pt_pos, ptw);
								// fill negative pion eta
								hGenpT[1][1][i]->Fill(pt_neg, ptw);

								// fill generated and reconstructed
								if (pipos_gen_rec)
									hGenRecpT[1][0][i]->Fill(pt_pos, ptw);
								if (pineg_gen_rec)
									hGenRecpT[1][1][i]->Fill(pt_neg, ptw);

								if (tr1_genX && tr2_genX) // selects generated pair with same fiducial cut as rec.
								{
									// pair level efficieicy
									hpairGenJP1X[i]->Fill(Minv_mc, ptw);
									if (pipos_gen_recX && pineg_gen_recX)
										hpairRecJP1X[i]->Fill(Minv_mc, ptw);

									// fill positive pion eta
									hGenpT_Mx[1][0][i]->Fill(pt_pos, ptw);
									// fill negative pion eta
									hGenpT_Mx[1][1][i]->Fill(pt_neg, ptw);
									// fill generated and reconstructed
									if (pipos_gen_recX)
										hGenRecpT_Mx[1][0][i]->Fill(pt_pos, ptw);
									if (pineg_gen_recX)
										hGenRecpT_Mx[1][1][i]->Fill(pt_neg, ptw);
								}
							}
						}
					}
					if (isJP2)
					{
						// in x-sec bins
						for (int i = 0; i < xnBins; i++)
						{
							if (Minv_mc >= xnBinsEdges[i] && Minv_mc < xnBinsEdges[i + 1])
							{
								if (cone > Cone_Min)
								{
									hpairGenJP2C->Fill(Minv_mc, ptw);
									if (pipos_gen_rec && pineg_gen_rec)
										hpairGenRecJP2C->Fill(Minv_mc, ptw);
								}

								// pair level efficieicy
								hpairGenJP2[i]->Fill(Minv_mc, ptw);
								if (pipos_gen_rec && pineg_gen_rec)
									hpairRecJP2[i]->Fill(Minv_mc, ptw);

								// track level efficiency
								//  fill positive pion eta
								hGenpT[2][0][i]->Fill(pt_pos, ptw);
								// fill negative pion eta
								hGenpT[2][1][i]->Fill(pt_neg, ptw);

								// fill generated and reconstructed
								if (pipos_gen_rec)
									hGenRecpT[2][0][i]->Fill(pt_pos, ptw);
								if (pineg_gen_rec)
									hGenRecpT[2][1][i]->Fill(pt_neg, ptw);

								if (tr1_genX && tr2_genX) // selects generated pair with same fiducial cut as rec.
								{
									// pair level efficieicy
									hpairGenJP2X[i]->Fill(Minv_mc, ptw);
									if (pipos_gen_recX && pineg_gen_recX)
										hpairRecJP2X[i]->Fill(Minv_mc, ptw);

									// fill positive pion eta
									hGenpT_Mx[2][0][i]->Fill(pt_pos, ptw);
									// fill negative pion eta
									hGenpT_Mx[2][1][i]->Fill(pt_neg, ptw);
									// fill generated and reconstructed
									if (pipos_gen_recX)
										hGenRecpT_Mx[2][0][i]->Fill(pt_pos, ptw);
									if (pineg_gen_recX)
										hGenRecpT_Mx[2][1][i]->Fill(pt_neg, ptw);
								}
							}
						}
					}
				} // ends filling histograms with the same fiducial cuts
			}	  // track 2 loop
		}		  // track 1 loop
	}			  // event loop
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

	for (int ntrg = 0; ntrg < 3; ntrg++)
	{
		for (int nch = 0; nch < 2; nch++)
		{
			for (int nbin = 0; nbin < 13; nbin++)
			{
				hGenRecpT_Mx[ntrg][nch][nbin]->Write();
				hGenpT_Mx[ntrg][nch][nbin]->Write();

				hGenRecpT[ntrg][nch][nbin]->Write();
				hGenpT[ntrg][nch][nbin]->Write();
			}
		}
	}
	for (int nbin = 0; nbin < 13; nbin++)
	{
		hpairGenJP0[nbin]->Write();
		hpairRecJP0[nbin]->Write();
		hpairGenJP1[nbin]->Write();
		hpairRecJP1[nbin]->Write();
		hpairGenJP2[nbin]->Write();
		hpairRecJP2[nbin]->Write();

		hpairGenJP0X[nbin]->Write();
		hpairRecJP0X[nbin]->Write();
		hpairGenJP1X[nbin]->Write();
		hpairRecJP1X[nbin]->Write();
		hpairGenJP2X[nbin]->Write();
		hpairRecJP2X[nbin]->Write();
	}
	hpairGenJP0C->Write();
	hpairGenRecJP0C->Write();
	hpairGenJP1C->Write();
	hpairGenRecJP1C->Write();
	hpairGenJP2C->Write();
	hpairGenRecJP2C->Write();

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
