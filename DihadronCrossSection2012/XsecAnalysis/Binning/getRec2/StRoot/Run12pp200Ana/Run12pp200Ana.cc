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
// function to match the reconstructed track to the generated

void Run12pp200Ana::Loop()
{

	gROOT->Reset();
	TH1::SetDefaultSumw2();
	// initialize histograms
	Int_t nTrg = 3;
	Int_t nCh = 2;
	const char *tTrig[3] = {"JP0", "JP1", "JP2"};
	const char *tCh[2] = {"Pos", "Neg"};

	hpartonicPt = new TH1D("hpartonicPt", "", 100, 0, 80);

	// histograms for resolution studies: 12 invariant mass bins with the Minv wifth 0.3 GeV/c^2
	const int nBins = 12;
	Double_t nBinsEdges[nBins + 1] = {0.26, 0.56, 0.86, 1.16, 1.46, 1.76, 2.06, 2.36, 2.66, 2.96, 3.26, 3.56, 4.0};
	const int xnBins = 16;
	Double_t xnBinsEdges[xnBins + 1] = {0.28, 0.34, 0.40, 0.46, 0.52, 0.59, 0.67, 0.76, 0.86, 0.97, 1.10, 1.30, 1.6, 2.0, 2.5, 3.2, 4.0};

	hxsecMinvGen = new TH1D("hxsecMinvGen", "", xnBins, xnBinsEdges);

	hxsecMinvRecJP0 = new TH1D("hxsecMinvRecJP0", "", xnBins, xnBinsEdges);
	hxsecMinvRecJP1 = new TH1D("hxsecMinvRecJP1", "", xnBins, xnBinsEdges);
	hxsecMinvRecJP2 = new TH1D("hxsecMinvRecJP2", "", xnBins, xnBinsEdges);

	for (int i = 0; i < 12; i++)
	{
		hMGenRec[i] = new TH1D(Form("hMGenRec%i", i), "", 100, -2.0, 2.0);
		hMGenRecU[i] = new TH1D(Form("hMGenRecU%i", i), "", 100, -2.0, 2.0);
	}

	// reconstructed variables
	Double_t p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1;
	Double_t p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
	Double_t cosPhiRB, sinPhiRB, cosPhiRY, sinPhiRY;
	Double_t PhiR_cosB, PhiR_sinB, PhiRB, PhiR_cosY, PhiR_sinY, PhiRY, phiRB, phiDiffB, phi_cos, phi_sin, phi_pair;
	Double_t msqr, msqr1, msqr2;																//
	Double_t R1_p1, R2_p1, R1_p2, R2_p2, Rs_tr1, Rs_tr2, sourcePartEta, sourcePartPhi, Rs_pair; // opening angles between parton and tracks and associated parton eta, phi

	// pythia varibales
	Double_t p1x_py, p2x_py, p1y_py, p2y_py, p1z_py, p2z_py, psx_py, psy_py, psz_py, ps_py, cone_py, R_py, Rs_py, Rx_py, Ry_py, Rz_py, R1_py, Rx1_py, Ry1_py, Rz1_py;
	Double_t p1_py, p2_py, E1_py, E2_py, Minv_py, pT_pair_py, eta_pair_py, phi_pair_py;
	Double_t sourceParton_py, partonE_py, E_pair_py;

	Double_t pi = 3.14159265359;
	Double_t m_pion = 0.1396;	  // GeV
	Bool_t pairMatch = kFALSE;	  // Rs<0.5, Minv<4, and pT_pair<15, pair-source parton matching condition
	Bool_t pairMatch_py = kFALSE; // Rs<0.5, Minv<4, and pT_pair<15, pair-source parton matching condition
	Double_t pairSep;
	Bool_t geLimit = kFALSE;
	Bool_t pyLimit = kFALSE;

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
	// vectors to find unique pairs
	vector<Int_t> v_track1;
	vector<Int_t> v_track2;

	vector<Double_t> v_pairMinv;
	vector<Double_t> v_pairCone;
	vector<Double_t> v_pairEta;
	vector<Double_t> v_pairPt;

	vector<Double_t> v_pairMinv_py;
	vector<Double_t> v_pairCone_py;
	vector<Double_t> v_pairEta_py;
	vector<Double_t> v_pairPt_py;

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
		if (fabs(fVZ - fVZ_pyth) > 2.)
			continue;

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

		if (!isJP)
			continue; // only JP0, JP1, and JP2 triggers

		if (fPart1Phi > pi)
			fPart1Phi = fPart1Phi - 2. * pi;
		if (fPart2Phi > pi)
			fPart2Phi = fPart2Phi - 2. * pi;

		// angle between rec. track and gen. track.
		Double_t matchangle1;
		Double_t matchangle2;
		// pt, eta and phi of matched tracks (ge = rec., py = pythia)
		Double_t pt1ge, eta1ge, ch1ge, phi1ge;
		Double_t pt2ge, eta2ge, ch2ge, phi2ge;

		Double_t pt1py, eta1py, ch1py, phi1py;
		Double_t pt2py, eta2py, ch2py, phi2py;

		Int_t pid1py, pid2py;

		v_track1.clear();
		v_track2.clear();

		v_pairMinv.clear();
		v_pairCone.clear();
		v_pairEta.clear();
		v_pairPt.clear();

		v_pairMinv_py.clear();
		v_pairCone_py.clear();
		v_pairEta_py.clear();
		v_pairPt_py.clear();
		// Reconstructed Track Loop
		for (int tr = 0; tr < fmaxpar; tr++)
		{
			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];

			if (fpT[tr] > 0.5 && fpT[tr] < 15. && ffitPts[tr] > 15 && fabs(feta[tr]) < 1. && fitPtsRatio_tr > .51 && fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2.0) // track 1 quality cuts
			{
				// find source parton for track 1
				R1_p1 = sqrt(pow(feta[tr] - fPart1Eta, 2) + pow(fphi[tr] - fPart1Phi, 2));
				R1_p2 = sqrt(pow(feta[tr] - fPart2Eta, 2) + pow(fphi[tr] - fPart2Phi, 2));
				if (R1_p1 >= R1_p2)
				{
					Rs_tr1 = R1_p2;
					sourcePartEta = fPart2Eta;
					sourcePartPhi = fPart2Phi;
				}
				else if (R1_p1 < R1_p2)
				{
					Rs_tr1 = R1_p1;
					sourcePartEta = fPart1Eta;
					sourcePartPhi = fPart1Phi;
				}
				// if (Rs_tr1 > 0.5)
				//	continue;

				if (!dcaCut(fpT[tr], fdca[tr])) // pt dependent dca cut
					continue;

				// find MC track that is closest to the Rec. track.
				matchFinder(feta[tr], fphi[tr], fpT[tr], fcharge[tr], fmaxpar1, statusCode_pyth, fpId_pyth, feta_pyth, fphi_pyth, fpT_pyth, eta1ge, phi1ge, pt1ge, ch1ge, matchangle1, eta1py, phi1py, pt1py, ch1py, pid1py);

				//  Track 2 loop
				for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
				{
					Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];
					if (fpT[tr2] > 0.5 && fpT[tr2] < 15. && ffitPts[tr2] > 15 && fabs(feta[tr2]) < 1. && fnSigmaPion[tr2] > -1.0 && fnSigmaPion[tr2] < 2.0 && fitPtsRatio_tr2 > .51)
					{
						// find source parton for track 2
						R2_p1 = sqrt(pow(feta[tr2] - fPart1Eta, 2) + pow(fphi[tr2] - fPart1Phi, 2));
						R2_p2 = sqrt(pow(feta[tr2] - fPart2Eta, 2) + pow(fphi[tr2] - fPart2Phi, 2));
						if (R2_p1 >= R2_p2)
						{
							Rs_tr2 = R2_p2;
							// source parton should be the same for both tracks
						}
						else if (R2_p1 < R2_p2)
						{
							Rs_tr2 = R2_p1;
						}

						// if (Rs_tr2 > 0.5)
						//	continue;

						if (!dcaCut(fpT[tr2], fdca[tr2]))
							continue;

						Double_t phiDiff = fphi[tr] - fphi[tr2];
						if (phiDiff > pi)
							phiDiff -= (2 * pi);
						if (phiDiff < ((-1) * pi))
							phiDiff += (2 * pi);
						Double_t phiDiff_py = phi1py - phi2py;
						if (phiDiff_py > pi)
							phiDiff_py -= (2 * pi);
						if (phiDiff_py < ((-1) * pi))
							phiDiff_py += (2 * pi);

						// find match for track 2
						matchFinder(feta[tr2], fphi[tr2], fpT[tr2], fcharge[tr2], fmaxpar1, statusCode_pyth, fpId_pyth, feta_pyth, fphi_pyth, fpT_pyth, eta2ge, phi2ge, pt2ge, ch2ge, matchangle2, eta2py, phi2py, pt2py, ch2py, pid2py);

						// if ((ch1ge == ch2ge) && (ch1py == ch2py))
						if (ch1ge == ch2ge)
							continue;
						if (ch1py == ch2py)
							continue;
						cone = sqrt(pow(eta1ge - eta2ge, 2) + pow(phiDiff, 2));
						cone_py = sqrt(pow(eta1py - eta2py, 2) + pow(phiDiff, 2));
						// cout << "Cone: " << cone << endl;
						if (cone >= 0.7)
							continue;
						// cout << "tr1 pid = " << pid1py << " , tr2 pid = " << pid2py << endl;
						if (ch1ge > 0)
						{
							// rec. pion
							p1x = pt1ge * cos(phi1ge);
							p2x = pt2ge * cos(phi2ge);
							p1y = pt1ge * sin(phi1ge);
							p2y = pt2ge * sin(phi2ge);
							p1z = pt1ge * sinh(eta1ge);
							p2z = pt2ge * sinh(eta2ge);
							p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
							p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
						}
						// gen pion
						if (ch1py > 0)
						{
							p1x_py = pt1py * cos(phi1py);
							p2x_py = pt2py * cos(phi2py);
							p1y_py = pt1py * sin(phi1py);
							p2y_py = pt2py * sin(phi2py);
							p1z_py = pt1py * sinh(eta1py);
							p2z_py = pt2py * sinh(eta2py);
							p1_py = sqrt(p1x_py * p1x_py + p1y_py * p1y_py + p1z_py * p1z_py);
							p2_py = sqrt(p2x_py * p2x_py + p2y_py * p2y_py + p2z_py * p2z_py);
						}
						if (ch1ge < 0)
						{
							p1x = pt2ge * cos(phi2ge);
							p2x = pt1ge * cos(phi1ge);
							p1y = pt2ge * sin(phi2ge);
							p2y = pt1ge * sin(phi1ge);
							p1z = pt2ge * sinh(eta2ge);
							p2z = pt1ge * sinh(eta1ge);
							p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
							p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
						}
						if (ch1py < 0)
						{
							// gen pion
							p1x_py = pt2py * cos(phi2py);
							p2x_py = pt1py * cos(phi1py);
							p1y_py = pt2py * sin(phi2py);
							p2y_py = pt1py * sin(phi1py);
							p1z_py = pt2py * sinh(eta2py);
							p2z_py = pt1py * sinh(eta1py);
							p1_py = sqrt(p1x_py * p1x_py + p1y_py * p1y_py + p1z_py * p1z_py);
							p2_py = sqrt(p2x_py * p2x_py + p2y_py * p2y_py + p2z_py * p2z_py);
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
						// cout << "Rs_tr1 = " << Rs_tr1 << ", Rs_tr2 =  " << Rs_tr2 << ",  Rs_pair = " << Rs_pair << endl;

						// if (Rs_pair > 0.5) // || Minv > 4.0 || fabs(eta_pair) > 1.0 || pT_pair > 15.0)
						if (Minv > 4.0 || fabs(eta_pair) > 1.0 || pT_pair > 15.0)
							continue;

						// generated pair
						psx_py = p1x_py + p2x_py;
						psy_py = p1y_py + p2y_py;
						psz_py = p1z_py + p2z_py;
						ps_py = sqrt(psx_py * psx_py + psy_py * psy_py + psz_py * psz_py);
						Rx_py = (p1x_py - p2x_py);
						Ry_py = (p1y_py - p2y_py);
						Rz_py = (p1z_py - p2z_py);
						R_py = sqrt(Rx_py * Rx_py + Ry_py * Ry_py + Rz_py * Rz_py); // R and R1 are same
						E1_py = sqrt(m_pion * m_pion + p1_py * p1_py);
						E2_py = sqrt(m_pion * m_pion + p2_py * p2_py);
						Minv_py = sqrt(2 * m_pion * m_pion + 2 * (E1_py * E2_py - p1x_py * p2x_py - p1y_py * p2y_py - p1z_py * p2z_py));
						eta_pair_py = TMath::ASinH((p1z_py + p2z_py) / pT_pair_py);
						pT_pair_py = sqrt((p1x_py + p2x_py) * (p1x_py + p2x_py) + (p1y_py + p2y_py) * (p1y_py + p2y_py));

						// cout << "Gen: " << Minv_py << ", Rec: " << Minv << ", Diff: " << Minv_py - Minv << endl;

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

						// cout << "Pythia - Reco. = " << Minv_py - Minv << endl;

						for (int i = 0; i < nBins; i++)
						{
							if (Minv >= nBinsEdges[i] && Minv < nBinsEdges[i + 1])
							{
								hMGenRec[i]->Fill(Minv_py - Minv);
							}
						}

						v_track1.push_back(tr);
						v_track2.push_back(tr2);

						v_pairCone.push_back(cone);
						v_pairMinv.push_back(Minv);
						v_pairPt.push_back(pT_pair);
						v_pairEta.push_back(eta_pair);

						v_pairCone_py.push_back(cone_py);
						v_pairMinv_py.push_back(Minv_py);
						v_pairPt_py.push_back(pT_pair_py);
						v_pairEta_py.push_back(eta_pair_py);

					} // track 2 cuts
				}	  // track 2 loop
			}		  // track 1 PID and selection cuts
		}			  // track loop

		if (v_pairCone.size() == 0)
			continue;

		/*cout << "Event no : " << jentry << "\n";
		cout << "Old vector: "
			 << "\n";

		for (int k = 0; k < v_track1.size(); k++)
		{
			cout << "tr 1= " << v_track1[k] << ", tr2 =  " << v_track2[k] << " Cone" << v_pairCone[k] << ", PairPt = " << v_pairPt[k] << " pairEta = " << v_pairEta[k] << " pairMass = " << v_pairMinv[k] << endl;
		}
		*/
		// remove repeated track 1 keeping one with the smallest cone
		if (v_pairCone.size() > 1 && adjacent_find(v_track1.begin(), v_track1.end()) != v_track1.end())
		{
			for (int i = 0; i < v_track1.size(); i++)
			{
				for (int j = i + 1; j < v_track1.size(); j++)
				{
					if (v_track1[i] == v_track1[j] && v_pairCone[j] >= v_pairCone[i])
					{

						v_track1.erase(v_track1.begin() + j);
						v_track2.erase(v_track2.begin() + j);

						v_pairCone.erase(v_pairCone.begin() + j);
						v_pairMinv.erase(v_pairMinv.begin() + j);
						v_pairEta.erase(v_pairEta.begin() + j);
						v_pairPt.erase(v_pairPt.begin() + j);

						v_pairCone_py.erase(v_pairCone_py.begin() + j);
						v_pairMinv_py.erase(v_pairMinv_py.begin() + j);
						v_pairEta_py.erase(v_pairEta_py.begin() + j);
						v_pairPt_py.erase(v_pairPt_py.begin() + j);

						i--;
						continue;
					}
					else if (v_track1[i] == v_track1[j] && v_pairCone[j] < v_pairCone[i])
					{
						v_track1.erase(v_track1.begin() + i);
						v_track2.erase(v_track2.begin() + i);

						v_pairCone.erase(v_pairCone.begin() + i);
						v_pairMinv.erase(v_pairMinv.begin() + i);
						v_pairEta.erase(v_pairEta.begin() + i);
						v_pairPt.erase(v_pairPt.begin() + i);

						v_pairCone_py.erase(v_pairCone_py.begin() + j);
						v_pairMinv_py.erase(v_pairMinv_py.begin() + j);
						v_pairEta_py.erase(v_pairEta_py.begin() + j);
						v_pairPt_py.erase(v_pairPt_py.begin() + j);

						i--;
						continue;
					}
				}
			}
		}

		// reomve events having commom tracks from track 2
		vector<Int_t> cv_track2 = v_track2;
		sort(cv_track2.begin(), cv_track2.end());
		// if (v_pairCone.size() > 1 && adjacent_find(v_track2.begin(), v_track2.end()) != v_track2.end())
		if (v_pairCone.size() > 1 && adjacent_find(cv_track2.begin(), cv_track2.end()) != cv_track2.end())
		{
			for (int i = 0; i < v_track2.size(); i++)
			{
				for (int j = i + 1; j < v_track2.size(); j++)
				{
					if (v_track2[i] == v_track2[j] && v_pairCone[j] >= v_pairCone[i])
					{

						v_pairCone.erase(v_pairCone.begin() + j);
						v_track1.erase(v_track1.begin() + j);
						v_track2.erase(v_track2.begin() + j);
						v_pairMinv.erase(v_pairMinv.begin() + j);
						v_pairEta.erase(v_pairEta.begin() + j);
						v_pairPt.erase(v_pairPt.begin() + j);

						v_pairCone_py.erase(v_pairCone_py.begin() + j);
						v_pairMinv_py.erase(v_pairMinv_py.begin() + j);
						v_pairEta_py.erase(v_pairEta_py.begin() + j);
						v_pairPt_py.erase(v_pairPt_py.begin() + j);

						i--;
						continue;
					}
					else if (v_track2[i] == v_track2[j] && v_pairCone[j] < v_pairCone[i])
					{

						v_track1.erase(v_track1.begin() + i);
						v_track2.erase(v_track2.begin() + i);

						v_pairCone.erase(v_pairCone.begin() + i);
						v_pairMinv.erase(v_pairMinv.begin() + i);
						v_pairEta.erase(v_pairEta.begin() + i);
						v_pairPt.erase(v_pairPt.begin() + i);

						v_pairCone_py.erase(v_pairCone_py.begin() + j);
						v_pairMinv_py.erase(v_pairMinv_py.begin() + j);
						v_pairEta_py.erase(v_pairEta_py.begin() + j);
						v_pairPt_py.erase(v_pairPt_py.begin() + j);

						i--;
						continue;
					}
				}
			}
		}

		// check if track 1 and track 2 has common tracks and remove one with larger cone
		// copy vectors for both tr1 and tr2 in to a single vector

		vector<int> tr_ids;
		tr_ids.insert(tr_ids.begin(), v_track1.begin(), v_track1.end());
		tr_ids.insert(tr_ids.end(), v_track2.begin(), v_track2.end());
		sort(tr_ids.begin(), tr_ids.end());
		if (v_pairCone.size() > 1 && adjacent_find(tr_ids.begin(), tr_ids.end()) != tr_ids.end())
		{
			for (int i = 0; i < v_track1.size(); i++)
			{
				for (int j = 0; j < v_track2.size(); j++)
				{
					if (v_track1[i] == v_track2[j] && v_pairCone[j] >= v_pairCone[i])
					{

						v_pairCone.erase(v_pairCone.begin() + j);
						v_track1.erase(v_track1.begin() + j);
						v_track2.erase(v_track2.begin() + j);

						v_pairMinv.erase(v_pairMinv.begin() + j);
						v_pairEta.erase(v_pairEta.begin() + j);
						v_pairPt.erase(v_pairPt.begin() + j);

						i--;
						j--;
						continue;
					}
					else if (v_track1[i] == v_track2[j] && v_pairCone[j] < v_pairCone[i])
					{

						v_track1.erase(v_track1.begin() + i);
						v_track2.erase(v_track2.begin() + i);

						v_pairCone.erase(v_pairCone.begin() + i);
						v_pairMinv.erase(v_pairMinv.begin() + i);
						v_pairEta.erase(v_pairEta.begin() + i);
						v_pairPt.erase(v_pairPt.begin() + i);

						i--;
						j--;
						continue;
					}
				}
			}
		}

		/*cout << "Final unique pair: \n " << endl;
		for (int k = 0; k < v_track1.size(); k++)
		{
			cout << "tr 1= " << v_track1[k] << ", tr2 =  " << v_track2[k] << " Cone" << v_pairCone[k] << ", PairPt= " << v_pairPt[k] << " pairEta = " << v_pairEta[k] << " pairMass = " << v_pairMinv[k] << endl;
		}
		*/
		for (int j = 0; j < v_track1.size(); j++)
		{
			for (int i = 0; i < nBins; i++)
			{
				if (v_pairMinv[j] >= nBinsEdges[i] && v_pairMinv[j] < nBinsEdges[i + 1])
				{
					hMGenRecU[i]->Fill(v_pairMinv_py[j] - v_pairMinv[j]);
				}
			}
		}
		for (int j = 0; j < v_track1.size(); j++)
		{
			hxsecMinvGen->Fill(v_pairMinv_py[j]);
			for (int i = 0; i < xnBins; i++)
			{
				if (v_pairMinv[j] >= xnBinsEdges[i] && v_pairMinv[j] < xnBinsEdges[i + 1])
				{
					if (isJP0)
						hxsecMinvRecJP0->Fill(v_pairMinv[j]);
					if (isJP1)
						hxsecMinvRecJP1->Fill(v_pairMinv[j]);
					if (isJP2)
						hxsecMinvRecJP2->Fill(v_pairMinv[j]);
				}
			}
		}

		tr_ids.clear();
		cv_track2.clear();
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
	fChain->SetBranchAddress("fevTime", &fevTime, &b_fevTime);
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

	fChain->SetBranchAddress("fVZ_pyth", &fVZ_pyth, &b_fVZ_pyth);
	fChain->SetBranchAddress("partonicPtBin", &partonicPtBin, &b_partonicPtBin);
	fChain->SetBranchAddress("statusCode_pyth", statusCode_pyth, &b_statusCode_pyth);
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

// routine match of detector pion to the pythia pions......
void Run12pp200Ana::matchFinder(Double_t eta_ge, Double_t phi_ge, Double_t pt_ge, Double_t ch_ge, int pyCounts, int statusCode[], Double_t pid_py[], Double_t eta_py[], Double_t phi_py[], Double_t pt_py[], Double_t &etage, Double_t &phige, Double_t &pTge, Double_t &chge, Double_t &mAngle, Double_t &etapy, Double_t &phipy, Double_t &pTpy, Double_t &chpy, int &pidpy)
{
	Double_t ang;
	vector<Double_t> vAng, vEta_ge, vPhi_ge, vpT_ge, vCh_ge, vEta_py, vPhi_py, vpT_py, vCh_py;
	vector<Int_t> vPid_py;
	Double_t match = -1;
	int index;
	Double_t charge_py = 0;
	for (int i = 0; i < pyCounts; i++)
	{
		if (abs(pid_py[i]) != 211 && abs(pid_py[i]) != 321 && abs(pid_py[i]) != 2212)
			continue; // pythia pions
		if (pid_py[i] > 0)
			charge_py = 1;
		if (pid_py[i] < 0)
			charge_py = -1;
		ang = sqrt(pow(eta_ge - eta_py[i], 2) + pow(phi_ge - phi_py[i], 2));
		if (ang > 0.5)
			continue;
		vAng.push_back(ang);
		vpT_ge.push_back(pt_ge);
		vPhi_ge.push_back(phi_ge);
		vEta_ge.push_back(eta_ge);
		vCh_ge.push_back(ch_ge);
		vpT_py.push_back(pt_py[i]);
		vPhi_py.push_back(phi_py[i]);
		vEta_py.push_back(eta_py[i]);
		vCh_py.push_back(charge_py);
		vPid_py.push_back(pid_py[i]);
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
		pTpy = vpT_py[index];
		phipy = vPhi_py[index];
		etapy = vEta_py[index];
		chpy = vCh_py[index];
		pidpy = vPid_py[index];
	}
	vAng.clear();
	vEta_ge.clear();
	vPhi_ge.clear();
	vpT_ge.clear();
	vCh_ge.clear();
	vEta_py.clear();
	vPhi_py.clear();
	vpT_py.clear();
	vPid_py.clear();
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

	for (int i = 0; i < 12; i++)
	{
		hMGenRec[i]->SetDirectory(fout);
		hMGenRec[i]->Write();
		hMGenRecU[i]->SetDirectory(fout);
		hMGenRecU[i]->Write();
	}
	hxsecMinvGen->SetDirectory(fout);
	hxsecMinvGen->Write();

	hxsecMinvRecJP0->SetDirectory(fout);
	hxsecMinvRecJP0->Write();
	hxsecMinvRecJP1->SetDirectory(fout);
	hxsecMinvRecJP1->Write();
	hxsecMinvRecJP2->SetDirectory(fout);
	hxsecMinvRecJP2->Write();

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
