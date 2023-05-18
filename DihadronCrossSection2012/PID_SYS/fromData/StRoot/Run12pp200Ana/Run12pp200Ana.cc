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
	Double_t p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1;
	Double_t p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
	Double_t cosPhiRB, sinPhiRB, cosPhiRY, sinPhiRY;
	Double_t PhiR_cosB, PhiR_sinB, PhiRB, PhiR_cosY, PhiR_sinY, PhiRY, phiRB, phiDiffB, phi_cos, phi_sin, phi_pair;
	Double_t msqr, msqr1, msqr2; //
	Double_t pi = 3.14159265359;
	Double_t CONE_CUT_MIN = 0.05;
	Double_t m_pion = 0.1396; // GeV

	// for sub-divided PID in momentum bins within mass bins
	Double_t pBinEdges[6] = {0.5, 0.6, 1.0, 1.5, 2.5, 10};
	const char *pBin[5] = {"p0", "p1", "p2", "p3", "p4"};
	for (int i = 0; i < 5; i++)
	{
		hSigmaPivsK_M0[i] = new TH2D(Form("hSigmaPivsK_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPivsP_M0[i] = new TH2D(Form("hSigmaPivsP_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPivsE_M0[i] = new TH2D(Form("hSigmaPivsE_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPion_M0[i] = new TH1D(Form("hSigmaPion_M0%s", pBin[i]), "", 100, -10, 10);

		hSigmaPivsKtof_M0[i] = new TH2D(Form("hSigmaPivsKtof_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPivsPtof_M0[i] = new TH2D(Form("hSigmaPivsPtof_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPivsEtof_M0[i] = new TH2D(Form("hSigmaPivsEtof_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPiontof_M0[i] = new TH1D(Form("hSigmaPiontof_M0%s", pBin[i]), "", 100, -10, 10);

		hSigmaPivsKboth_M0[i] = new TH2D(Form("hSigmaPivsKboth_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPivsPboth_M0[i] = new TH2D(Form("hSigmaPivsPboth_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPivsEboth_M0[i] = new TH2D(Form("hSigmaPivsEboth_M0%s", pBin[i]), "", 100, -1, 1, 100, -5, 5);
		hSigmaPionboth_M0[i] = new TH1D(Form("hSigmaPionboth_M0%s", pBin[i]), "", 100, -10, 10);
	}
	//----------------------------
	const int nBins = 17;
	const int ncharge = 2;
	const char *charge[ncharge] = {"Pos", "Neg"};
	// updated binning
	Double_t mBinsEdges[nBins + 1] = {0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.58, 0.64, 0.72, 0.82, 0.95, 1.10, 1.40, 1.90, 2.5, 3.1, 4.0};

	for (int j = 0; j < ncharge; j++)
	{
		for (int i = 0; i < nBins; i++)
		{
			hSigmaPivsK[j][i] = new TH2D(Form("hSigmaPivsK_%s_M%i", charge[j], i), "", 100, -1, 1, 100, -5, 5);
			hSigmaPivsP[j][i] = new TH2D(Form("hSigmaPivsP_%s_M%i", charge[j], i), "", 100, -1, 1, 100, -5, 5);
			hSigmaPivsE[j][i] = new TH2D(Form("hSigmaPivsE_%s_M%i", charge[j], i), "", 100, -1, 1, 100, -5, 5);

			hPidEdxVsP[j][i] = new TH2D(Form("hPidEdxVsP_%s_M%i", charge[j], i), "", 100, 0.0, 5.0, 100, 0.0, 10.0);
			hdEdxVsP[j][i] = new TH2D(Form("hdEdxVsP_%s_M%i", charge[j], i), "", 100, 0.0, 5.0, 100, 0.0, 10.0);

			hSigmaPionTPCvsTOF[j][i] = new TH2D(Form("hSigmaPionTPCvsTOF_%s_M%i", charge[j], i), "", 100, -10.0, 10.0, 100, -10.0, 10.0);

			hSigmaPion[j][i] = new TH1D(Form("hSigmaPion_%s_M%i", charge[j], i), "", 100, -10, 10);

			hSigmaPionTOF[j][i] = new TH1D(Form("hSigmaPionTOF_%s_M%i", charge[j], i), "Only TOF tracks", 100, -10, 10);
			hSigmaPionTPC[j][i] = new TH1D(Form("hSigmaPionTPC_%s_M%i", charge[j], i), "Only TOF tracks", 100, -10, 10);

			hp[j][i] = new TH1D(Form("hp_%s_M%i", charge[j], i), "", 100, 0.45, 15.0);
		}
	}

	if (fChain == 0)
		return;

	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	cout << "Event loop sarted, Number of Events: " << nentries << endl;
	// Event Loop

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
		if (!isJP)
			continue; // only JP0, JP1, and JP2 triggers

		// Track Loop
		for (int tr = 0; tr < fmaxpar; tr++)
		{
			Double_t fitPtsRatio_tr = (Double_t)ffitPts[tr] / (Double_t)ffitPtsPoss[tr];
			if (fpT[tr] > 0.5 && ffitPts[tr] > 15 && abs(feta[tr]) < 1. && fitPtsRatio_tr > .51 && fabs(fnSigmaPion[tr]) < 10.0 && fhitsdedx[tr] > 20)
			{

				// if (fBetaToF[tr] != -999 && fabs(fnSigmaPionTof[tr]) > 2.0)
				if (fBetaToF[tr] != -999 && fabs(fnSigmaPionTof[tr]) > 10.0)
					continue; // if track has TOF info
				Double_t inverse_pion_beta1 = sqrt(pow(m_pion, 2) / pow(fp[tr], 2) + 1);
				if (fBetaToF[tr] != -999)
					msqr1 = fp[tr] * fp[tr] * (1 / (fBetaToF[tr] * fBetaToF[tr]) - 1);

				// if (!(fnSigmaPion[tr] > -1.0 && fnSigmaPion[tr] < 2.0))
				//	continue;

				// Track 2 loop
				for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
				{
					Double_t fitPtsRatio_tr2 = (Double_t)ffitPts[tr2] / (Double_t)ffitPtsPoss[tr2];
					if (fpT[tr2] > 0.5 && ffitPts[tr2] > 15 && fabs(feta[tr2]) < 1. && fabs(fnSigmaPion[tr2]) < 10.0 && fitPtsRatio_tr2 > .51 && fhitsdedx[tr2] > 20)
					{
						// if (fBetaToF[tr2] != -999 && fabs(fnSigmaPionTof[tr2]) > 2.0)
						if (fBetaToF[tr2] != -999 && fabs(fnSigmaPionTof[tr2]) > 10.0)
							continue; // if track has TOF info

						Double_t inverse_pion_beta2 = sqrt(pow(m_pion, 2) / pow(fp[tr2], 2) + 1);

						if (fBetaToF[tr2] != -999)
							msqr2 = fp[tr2] * fp[tr2] * (1 / (fBetaToF[tr2] * fBetaToF[tr2]) - 1);

						Double_t phiDiff = fphi[tr] - fphi[tr2];
						if (phiDiff > pi)
							phiDiff -= (2 * pi);
						if (phiDiff < ((-1) * pi))
							phiDiff += (2 * pi);
						if (phiDiff > pi || phiDiff < -pi)
							continue;

						if (fcharge[tr] == fcharge[tr2])
							continue;

						cone = sqrt(pow(feta[tr] - feta[tr2], 2) + pow(phiDiff, 2)); // cone cut < 0.3

						if (cone >= 0.7)
							continue;

						if (fcharge[tr] > 0)
						{
							p1x = fpX[tr];
							p2x = fpX[tr2];
							p1y = fpY[tr];
							p2y = fpY[tr2];
							p1z = fpZ[tr];
							p2z = fpZ[tr2];
							p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
							p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
						}
						if (fcharge[tr] < 0)
						{
							p1x = fpX[tr2];
							p2x = fpX[tr];
							p1y = fpY[tr2];
							p2y = fpY[tr];
							p1z = fpZ[tr2];
							p2z = fpZ[tr];
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
						double sigmaPionP = -999, sigmaProtonP = -999, sigmaElectronP = -999, sigmaKaonP = -999, sigmaPionPtof = -999, sigmaKaonPtof = -999, sigmaProtonPtof = -999, sigmaElectronPtof = -999, sigmaPionNtof = -999, sigmaKaonNtof = -999, sigmaProtonNtof = -999, sigmaElectronNtof = -999, pP = -999, ptP = -999, sigmaPionN = -999, sigmaProtonN = -999, sigmaElectronN = -999, sigmaKaonN = -999, pN = -999, ptN = -999, dEdxN = -999, dEdxP = -999, tofBetaP = -999, tofBetaN = -999;

						if (fcharge[tr] > 0)
						{
							pP = fp[tr];
							ptP = fpT[tr];
							sigmaPionP = fnSigmaPion[tr];
							sigmaProtonP = fnSigmaProton[tr];
							sigmaElectronP = fnSigmaElectron[tr];
							sigmaKaonP = fnSigmaKaon[tr];
							sigmaPionPtof = fnSigmaPionTof[tr];
							sigmaProtonPtof = fnSigmaProtonTof[tr];
							sigmaElectronPtof = fnSigmaElectronTof[tr];
							sigmaKaonPtof = fnSigmaKaonTof[tr];
							dEdxP = fdEdx[tr];
							tofBetaP = fBetaToF[tr];

							pN = fp[tr2];
							ptN = fpT[tr2];
							sigmaPionN = fnSigmaPion[tr2];
							sigmaProtonN = fnSigmaProton[tr2];
							sigmaElectronN = fnSigmaElectron[tr2];
							sigmaKaonN = fnSigmaKaon[tr2];
							sigmaPionNtof = fnSigmaPionTof[tr2];
							sigmaProtonNtof = fnSigmaProtonTof[tr2];
							sigmaElectronNtof = fnSigmaElectronTof[tr2];
							sigmaKaonNtof = fnSigmaKaonTof[tr2];
							dEdxN = fdEdx[tr2];
							tofBetaN = fBetaToF[tr2];
						}
						else if (fcharge[tr] < 0)
						{
							pP = fp[tr2];
							ptP = fpT[tr2];
							sigmaPionP = fnSigmaPion[tr2];
							sigmaProtonP = fnSigmaProton[tr2];
							sigmaElectronP = fnSigmaElectron[tr2];
							sigmaKaonP = fnSigmaKaon[tr2];
							sigmaPionPtof = fnSigmaPionTof[tr2];
							sigmaProtonPtof = fnSigmaProtonTof[tr2];
							sigmaElectronPtof = fnSigmaElectronTof[tr2];
							sigmaKaonPtof = fnSigmaKaonTof[tr2];
							dEdxP = fdEdx[tr2];
							tofBetaP = fBetaToF[tr2];

							pN = fp[tr];
							ptN = fpT[tr];
							sigmaPionN = fnSigmaPion[tr];
							sigmaProtonN = fnSigmaProton[tr];
							sigmaElectronN = fnSigmaElectron[tr];
							sigmaKaonN = fnSigmaKaon[tr];
							sigmaPionNtof = fnSigmaPionTof[tr];
							sigmaProtonNtof = fnSigmaProtonTof[tr];
							sigmaElectronNtof = fnSigmaElectronTof[tr];
							sigmaKaonNtof = fnSigmaKaonTof[tr];
							dEdxN = fdEdx[tr];
							tofBetaN = fBetaToF[tr];
						}
						for (int nbin = 0; nbin < nBins; nbin++)
						{
							if (Minv >= mBinsEdges[nbin] && Minv < mBinsEdges[nbin + 1])
							{
								hSigmaPion[0][nbin]->Fill(sigmaPionP);
								hSigmaPion[1][nbin]->Fill(sigmaPionN);

								if (tofBetaP != -999)
								{
									hSigmaPionTPCvsTOF[0][nbin]->Fill(sigmaPionPtof, sigmaPionP);
									hSigmaPionTOF[0][nbin]->Fill(sigmaPionPtof);
									hSigmaPionTPC[0][nbin]->Fill(sigmaPionP);
								}
								if (tofBetaN != -999)
								{
									hSigmaPionTPCvsTOF[1][nbin]->Fill(sigmaPionNtof, sigmaPionN);
									hSigmaPionTOF[1][nbin]->Fill(sigmaPionNtof);
									hSigmaPionTPC[1][nbin]->Fill(sigmaPionN);
								}

								hSigmaPion[1][nbin]->Fill(sigmaPionN);

								hSigmaPivsP[0][nbin]->Fill(sigmaProtonP, sigmaPionP);
								hSigmaPivsP[1][nbin]->Fill(sigmaProtonN, sigmaPionN);

								hSigmaPivsK[0][nbin]->Fill(sigmaKaonP, sigmaPionP);
								hSigmaPivsK[1][nbin]->Fill(sigmaKaonN, sigmaPionN);

								hSigmaPivsE[0][nbin]->Fill(sigmaElectronP, sigmaPionP);
								hSigmaPivsE[1][nbin]->Fill(sigmaElectronN, sigmaPionN);

								hdEdxVsP[0][nbin]->Fill(pP, dEdxP);
								hdEdxVsP[1][nbin]->Fill(pN, dEdxN);

								if (sigmaPionP > -1 && sigmaPionP < 2)
									hPidEdxVsP[0][nbin]->Fill(pP, dEdxP);
								if (sigmaPionN > -1 && sigmaPionN < 2)
									hPidEdxVsP[1][nbin]->Fill(pN, dEdxN);

								hp[0][nbin]->Fill(pP);
								hp[1][nbin]->Fill(pN);

								// only in first Minv bin and for +ve tracks
								if (nbin == 0)
								{
									for (int pbin = 0; pbin < 5; pbin++)
									{
										if (pP >= pBinEdges[pbin] && pP < pBinEdges[pbin + 1])
										{
											hSigmaPivsP_M0[pbin]->Fill(sigmaProtonP, sigmaPionP);
											hSigmaPivsK_M0[pbin]->Fill(sigmaKaonP, sigmaPionP);
											hSigmaPivsE_M0[pbin]->Fill(sigmaElectronP, sigmaPionP);
											hSigmaPion_M0[pbin]->Fill(sigmaPionP);
											if (tofBetaP != -999 && abs(sigmaPionPtof) < 2.0)
											{ // tpc and tof cut all tof tracks
												hSigmaPivsPtof_M0[pbin]->Fill(sigmaProtonP, sigmaPionP);
												hSigmaPivsKtof_M0[pbin]->Fill(sigmaKaonP, sigmaPionP);
												hSigmaPivsEtof_M0[pbin]->Fill(sigmaElectronP, sigmaPionP);
												hSigmaPiontof_M0[pbin]->Fill(sigmaPionP);

												hSigmaPivsPboth_M0[pbin]->Fill(sigmaProtonP, sigmaPionP);
												hSigmaPivsKboth_M0[pbin]->Fill(sigmaKaonP, sigmaPionP);
												hSigmaPivsEboth_M0[pbin]->Fill(sigmaElectronP, sigmaPionP);
												hSigmaPionboth_M0[pbin]->Fill(sigmaPionP);
											}
											else if(tofBetaP == 999)
											{
												hSigmaPivsPboth_M0[pbin]->Fill(sigmaProtonP, sigmaPionP);
												hSigmaPivsKboth_M0[pbin]->Fill(sigmaKaonP, sigmaPionP);
												hSigmaPivsEboth_M0[pbin]->Fill(sigmaElectronP, sigmaPionP);
												hSigmaPionboth_M0[pbin]->Fill(sigmaPionP);
											}
										}
									}
								}
							}
						}

					} // track 2 cuts
				}	  // track 2 loop
			}		  // track 1 PID and selection cuts
		}			  // track loop
	}				  // event loop
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

	fChain->SetBranchAddress("fnSigmaPionTof", fnSigmaPionTof, &b_fnSigmaPionTof);
	fChain->SetBranchAddress("fnSigmaKaonTof", fnSigmaKaonTof, &b_fnSigmaKaonTof);
	fChain->SetBranchAddress("fnSigmaProtonTof", fnSigmaProtonTof, &b_fnSigmaProtonTof);
	fChain->SetBranchAddress("fnSigmaElectronTof", fnSigmaElectronTof, &b_fnSigmaElectronTof);

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
	// The Notify() functionTof is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}
void Run12pp200Ana::Finish(char *ofile)
{
	TFile *fout = new TFile(ofile, "recreate");
	for (int nch = 0; nch < 2; nch++)
	{
		for (int mbin = 0; mbin < 17; mbin++)
		{
			hSigmaPion[nch][mbin]->SetDirectory(fout);
			hSigmaPion[nch][mbin]->Write();
			hSigmaPionTOF[nch][mbin]->SetDirectory(fout);
			hSigmaPionTOF[nch][mbin]->Write();
			hSigmaPionTPC[nch][mbin]->SetDirectory(fout);
			hSigmaPionTPC[nch][mbin]->Write();

			hSigmaPionTPCvsTOF[nch][mbin]->SetDirectory(fout);
			hSigmaPionTPCvsTOF[nch][mbin]->Write();

			hSigmaPivsP[nch][mbin]->SetDirectory(fout);
			hSigmaPivsP[nch][mbin]->Write();

			hSigmaPivsK[nch][mbin]->SetDirectory(fout);
			hSigmaPivsK[nch][mbin]->Write();

			hSigmaPivsE[nch][mbin]->SetDirectory(fout);
			hSigmaPivsE[nch][mbin]->Write();

			hPidEdxVsP[nch][mbin]->SetDirectory(fout);
			hPidEdxVsP[nch][mbin]->Write();

			hdEdxVsP[nch][mbin]->SetDirectory(fout);
			hdEdxVsP[nch][mbin]->Write();

			hp[nch][mbin]->SetDirectory(fout);
			hp[nch][mbin]->Write();
		}
	}

	for (int pbin = 0; pbin < 5; pbin++)
	{
		hSigmaPivsP_M0[pbin]->SetDirectory(fout);
		hSigmaPivsP_M0[pbin]->Write();
		hSigmaPivsK_M0[pbin]->SetDirectory(fout);
		hSigmaPivsK_M0[pbin]->Write();
		hSigmaPivsE_M0[pbin]->SetDirectory(fout);
		hSigmaPivsE_M0[pbin]->Write();
		hSigmaPion_M0[pbin]->SetDirectory(fout);
		hSigmaPion_M0[pbin]->Write();

		hSigmaPivsPtof_M0[pbin]->SetDirectory(fout);
		hSigmaPivsPtof_M0[pbin]->Write();
		hSigmaPivsKtof_M0[pbin]->SetDirectory(fout);
		hSigmaPivsKtof_M0[pbin]->Write();
		hSigmaPivsEtof_M0[pbin]->SetDirectory(fout);
		hSigmaPivsEtof_M0[pbin]->Write();
		hSigmaPiontof_M0[pbin]->SetDirectory(fout);
		hSigmaPiontof_M0[pbin]->Write();

		hSigmaPivsPboth_M0[pbin]->SetDirectory(fout);
		hSigmaPivsPboth_M0[pbin]->Write();
		hSigmaPivsKboth_M0[pbin]->SetDirectory(fout);
		hSigmaPivsKboth_M0[pbin]->Write();
		hSigmaPivsEboth_M0[pbin]->SetDirectory(fout);
		hSigmaPivsEboth_M0[pbin]->Write();
		hSigmaPionboth_M0[pbin]->SetDirectory(fout);
		hSigmaPionboth_M0[pbin]->Write();
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
