

/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 22 17:04:53 2010 by ROOT version 5.22/00
// from TTree ftree/LRC
// found on file: 00E5DE07942C320D8D35F5F8AD40C261_7total.root
//////////////////////////////////////////////////////////

#ifndef Run12pp200Ana_h
#define Run12pp200Ana_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include <fstream>

using namespace std;

class Run12pp200Ana
{
public:
        TChain *fChain; // pointer to the analyzed TTree or TChain
        Int_t fCurrent; // current Tree number in a TChain

        // Declaration of leaf types
        Int_t frefmult;
        Int_t fmaxpar;
        Int_t ffillNum;
        Int_t frunNum;
        Int_t fspinconfig;
        vector<int> *ftrigger;

        Int_t trigJP0;
        Int_t trigJP1;
        Int_t trigJP2;
        Int_t mcId[2422];
        Int_t idTruth[2422];

        Double_t fVZ;
        Double_t fVZ_mc;
        Double_t fvpdVz;
        unsigned int fevTime;
        Double_t fverRank;
        Double_t fpT[2422];
        Double_t fp[2422];
        Double_t fpX[2422];
        Double_t fpY[2422];
        Double_t fpZ[2422];
        Double_t feta[2422];
        Double_t fphi[2422];
        Short_t fcharge[2422];
        Double_t fnSigmaPion[2422];
        Double_t fnSigmaKaon[2422];
        Double_t fnSigmaProton[2422];
        Double_t fnSigmaElectron[2422];
        Double_t fdEdx[2422];
        Double_t fdca[2422];
        UShort_t ffitPts[2422];
        UShort_t ffitPtsPoss[2422];
        UShort_t fhitsdedx[2422];
        Double_t fBetaToF[2422];

        // Part level info
        Double_t partonicPtBin;
        Int_t bPartId;
        Double_t bPartX1;
        Double_t bPartE;
        Double_t bPartEta;
        Double_t bPartPt;
        Double_t bPartPhi;

        Int_t yPartId;
        Double_t yPartX2;
        Double_t yPartE;
        Double_t yPartPhi;
        Double_t yPartEta;
        Double_t yPartPt;

        Int_t fPart1Id;
        Double_t fPart1Pt;
        Double_t fPart1Eta;
        Double_t fPart1Phi;
        Double_t fPart1E;
        Double_t fPart1X;

        Int_t fPart2Id;
        Double_t fPart2Pt;
        Double_t fPart2Eta;
        Double_t fPart2Phi;
        Double_t fPart2E;
        Double_t fPart2X;
        Double_t ParticPtBin;
        Double_t fVZ_pyth;
        Double_t fpT_pyth[2422];
        Double_t fp_pyth[2422];
        Double_t feta_pyth[2422];
        Double_t fphi_pyth[2422];
        Double_t fpId_pyth[2422];
        // List of branches
        TBranch *b_frefmult;
        TBranch *b_fmaxpar;
        TBranch *b_ffillNum;
        TBranch *b_frunNum;
        TBranch *b_fspinconfig;
        TBranch *b_ftrigger;
        TBranch *b_trigJP0;
        TBranch *b_trigJP1;
        TBranch *b_trigJP2;
        TBranch *b_fVZ;
        TBranch *b_fvpdVz;
        TBranch *b_fevTime;
        TBranch *b_fverRank;
        TBranch *b_fpT;
        TBranch *b_fp;

        TBranch *b_feta;
        TBranch *b_fphi;
        TBranch *b_fcharge;
        TBranch *b_fnSigmaPion;
        TBranch *b_fnSigmaKaon;
        TBranch *b_fnSigmaProton;
        TBranch *b_fnSigmaElectron;
        TBranch *b_fdEdx;
        TBranch *b_fdca;
        TBranch *b_ffitPts;
        TBranch *b_ffitPtsPoss;
        TBranch *b_fhitsdedx;
        TBranch *b_fBetaToF;

        TBranch *b_fVZ_pyth;
        TBranch *b_fVZ_mc;
        TBranch *b_bPartId;
        TBranch *b_bPartEta;
        TBranch *b_bPartPt;
        TBranch *b_bPartPhi;
        TBranch *b_bPartX1;
        TBranch *b_bPartE;

        TBranch *b_yPartId;
        TBranch *b_yPartX2;
        TBranch *b_yPartPhi;
        TBranch *b_yPartEta;

        TBranch *b_yPartE;
        TBranch *b_yPartPt;

        TBranch *b_fPart1Id;
        TBranch *b_fPart1Eta;
        TBranch *b_fPart1Phi;
        TBranch *b_fPart1E;
        TBranch *b_fPart1X;
        TBranch *b_fPart1Pt;

        TBranch *b_fPart2X;
        TBranch *b_fPart2Id;
        TBranch *b_fPart2Eta;
        TBranch *b_fPart2Phi;
        TBranch *b_fPart2E;
        TBranch *b_fPart2Pt;
        TBranch *b_partonicPtBin;
        TBranch *b_fpT_pyth;
        TBranch *b_fp_pyth;
        TBranch *b_feta_pyth;
        TBranch *b_fphi_pyth;
        TBranch *b_fpId_pyth;
        TBranch *b_idTruth;
        TBranch *b_mcId;

        Run12pp200Ana(char *ifile);
        virtual ~Run12pp200Ana();
        virtual Int_t Cut(Long64_t entry);
        virtual Int_t GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void Init();
        virtual void Loop();
        virtual Bool_t Notify();
        virtual void Finish(char *ofile);
        virtual void Show(Long64_t entry = -1);
        virtual Bool_t dcaCut(Double_t, Double_t);
        inline void setHome(TDirectory *pHome)
        {
                cout << "setting home.." << endl;
                home = pHome;
        }

        static const int n = 23;
        int trigger;
        TH1D *hpartonicPt;
        TH1D *hxsecMJP0;
        TH1D *hxsecMJP1;
        TH1D *hxsecMJP2;
        // All tracks
        TH1D *hVZ[3];

        TH1D *htrackpT[3];
        TH1D *htrackEta[3]; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *htrackPhi[3];
        TH1D *htracknSigmaPion[3];

        TH1D *hpionpT[3];
        TH1D *hpionEta[3]; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *hpionPhi[3];
        // for unique pair
        TH1D *hpionpTPosU[3];
        TH1D *hpionEtaPosU[3]; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *hpionPhiPosU[3];
        TH1D *hpionpTNegU[3];
        TH1D *hpionEtaNegU[3]; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *hpionPhiNegU[3];
        // pair quantities
        TH1D *hCone[3];
        TH1D *hMinv[3];
        TH1D *hpTPair[3];
        TH1D *hEtaPair[3];
        TH1D *hPhiRY[3];
        TH1D *hPhiRB[3];
        // for unique pair
        TH1D *hConeU[3];
        TH1D *hMinvU[3];
        TH1D *hpTPairU[3];
        TH1D *hEtaPairU[3];

protected:
        TDirectory *home;

private:
        ClassDef(Run12pp200Ana, 1);
};

#endif
